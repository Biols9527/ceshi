import logging
import os
import uuid
import shutil
from collections import defaultdict
# Bio.SeqIO is imported where needed to avoid issues if not always available during module load
# Bio.Blast.NCBIXML is imported inside select_best_by_alignment_quality

logger = logging.getLogger(__name__) # cds_processor.duplicate_handler

def process_fasta_batch_for_duplicates(
    batch_records,
    species_seqs_dict,
    species_ids_sequences_dict,
    species_counts_dict,
    duplicate_strategy,
    fasta_file_path,
    preprocess_function,
    blast_params=None,
    config=None
):
    from .fasta_parser import extract_species_name # Assuming this is safe or will be moved

    for record in batch_records:
        seq_id = record.id
        species = extract_species_name(record.description)
        raw_sequence = str(record.seq)

        processed_sequence, _ = preprocess_function(
            raw_sequence, seq_id, fasta_file_path, blast_params, config
        )

        if duplicate_strategy == "alignment_quality":
            if species not in species_ids_sequences_dict:
                species_ids_sequences_dict[species] = {}
            species_ids_sequences_dict[species][seq_id] = processed_sequence
            species_counts_dict[species] += 1
        else:
            if species in species_seqs_dict:
                species_counts_dict[species] += 1
                current_len = len(species_seqs_dict[species])
                new_len = len(processed_sequence)

                if duplicate_strategy == "longest":
                    if new_len > current_len:
                        species_seqs_dict[species] = processed_sequence
                elif duplicate_strategy == "first":
                    pass # Keep existing
                elif duplicate_strategy == "rename":
                    new_species_name = f"{species}_dup{species_counts_dict[species]}"
                    species_seqs_dict[new_species_name] = processed_sequence
                # Implicitly, if not longer or other strategy, current is kept or warning issued by caller
            else:
                species_seqs_dict[species] = processed_sequence
                species_counts_dict[species] = 1

def select_best_by_alignment_quality(
    species_ids_sequences_dict,
    blast_params,
    utils_module, # For things like hash_file_content, not used in this version
    blast_handler_module,
    config=None
):
    logger.info("Attempting to select best sequences by alignment quality.")
    final_sequences = {}

    try:
        from Bio.Blast import NCBIXML as BiopythonNCBIXMLParser
        from Bio import SeqIO as BiopythonSeqIO
    except ImportError:
        logger.critical("Biopython NCBIXML/SeqIO not found. Cannot perform BLAST-based duplicate selection. Falling back to 'longest'.")
        for species, id_to_seq_map in species_ids_sequences_dict.items():
            if not id_to_seq_map: continue
            final_sequences[species] = max(id_to_seq_map.values(), key=len)
        return final_sequences

    use_blast_for_dups = blast_params and blast_params.get("use_blast_for_duplicates", False)
    if not use_blast_for_dups:
        logger.warning("Alignment quality strategy: BLAST for duplicates not enabled. Falling back to 'longest'.")
        for species, id_to_seq_map in species_ids_sequences_dict.items():
            if not id_to_seq_map: continue
            final_sequences[species] = max(id_to_seq_map.values(), key=len)
        return final_sequences

    blastn_exe = blast_params.get("blastn_exe")
    makeblastdb_exe = blast_params.get("makeblastdb_exe")
    base_temp_dir = blast_params.get("temp_dir", "./temp_processing_data") # Default temp dir

    selection_temp_root = os.path.join(base_temp_dir, "dup_sel_blast_" + uuid.uuid4().hex)
    try:
        os.makedirs(selection_temp_root, exist_ok=True)
    except OSError as e:
        logger.error(f"Cannot create root temp dir {selection_temp_root} for duplicate selection: {e}. Fallback to 'longest'.")
        for species, id_to_seq_map in species_ids_sequences_dict.items():
            if not id_to_seq_map: continue
            final_sequences[species] = max(id_to_seq_map.values(), key=len)
        return final_sequences

    if not all([blastn_exe, makeblastdb_exe]):
        logger.error("Missing BLAST executables for duplicate selection. Fallback to 'longest'.")
        # Code for fallback
        for species, id_to_seq_map in species_ids_sequences_dict.items():
            if not id_to_seq_map: continue
            final_sequences[species] = max(id_to_seq_map.values(), key=len)
        shutil.rmtree(selection_temp_root, ignore_errors=True)
        return final_sequences

    for species, id_to_seq_map in species_ids_sequences_dict.items():
        if not id_to_seq_map or len(id_to_seq_map) <= 1:
            if id_to_seq_map: # Only one sequence
                final_sequences[species] = next(iter(id_to_seq_map.values()))
            continue

        logger.info(f"Processing {species} ({len(id_to_seq_map)} duplicates) with all-vs-all BLAST.")

        # Create a dedicated temp dir for this species' BLAST operations
        species_temp_dir = os.path.join(selection_temp_root, f"{species}_{uuid.uuid4().hex}")
        os.makedirs(species_temp_dir, exist_ok=True)

        group_fasta_path = os.path.join(species_temp_dir, "group.fasta")

        temp_records = []
        for seq_id, sequence_str in id_to_seq_map.items():
            # Create a minimal SeqRecord-like object for BiopythonSeqIO.write
            # This requires Bio.Seq and Bio.SeqRecord, or manual FASTA string writing
            # For simplicity with BiopythonSeqIO.write, actual SeqRecord might be better if available
            # Manual writing:
             with open(group_fasta_path, "a") as f_group_append: # Append to group fasta
                 f_group_append.write(f">{seq_id}\n{sequence_str}\n")

        db_prefix = None
        try:
            # DB is created inside species_temp_dir to keep it contained
            db_prefix = blast_handler_module.create_temp_blast_db(
                group_fasta_path, species_temp_dir, makeblastdb_exe
            )
        except Exception as e_db_create:
            logger.error(f"Exception creating BLAST DB for {species}: {e_db_create}")

        if not db_prefix:
            logger.error(f"Failed to create BLAST DB for {species}. Fallback to 'longest'.")
            final_sequences[species] = max(id_to_seq_map.values(), key=len)
        else:
            sequence_scores = {seq_id: 0.0 for seq_id in id_to_seq_map.keys()}

            for query_id, query_sequence in id_to_seq_map.items():
                query_single_fasta_path = os.path.join(species_temp_dir, f"{query_id}_query.fasta")
                with open(query_single_fasta_path, "w") as f_query:
                    f_query.write(f">{query_id}\n{query_sequence}\n")

                # BLAST output also goes into the species-specific temp dir
                blast_output_xml = blast_handler_module.run_blastn(
                    query_single_fasta_path, db_prefix, blastn_exe, species_temp_dir
                )

                if blast_output_xml and os.path.exists(blast_output_xml):
                    try:
                        with open(blast_output_xml) as result_handle:
                            for record in BiopythonNCBIXMLParser.parse(result_handle):
                                for alignment in record.alignments:
                                    if alignment.hit_id == query_id: continue # Ignore self-hits
                                    for hsp in alignment.hsps:
                                        sequence_scores[query_id] += hsp.score # Summing bit scores
                    except Exception as e_parse:
                        logger.error(f"Error parsing BLAST XML for {query_id} (species {species}): {e_parse}")
                    finally:
                        if os.path.exists(blast_output_xml): os.remove(blast_output_xml)
                if os.path.exists(query_single_fasta_path): os.remove(query_single_fasta_path)

            if not sequence_scores or all(s == 0 for s in sequence_scores.values()):
                logger.warning(f"No BLAST scores determined for {species}. Fallback to longest.")
                best_seq_id = max(id_to_seq_map, key=lambda k: len(id_to_seq_map[k]))
            else:
                best_seq_id = max(sequence_scores, key=sequence_scores.get)

            final_sequences[species] = id_to_seq_map[best_seq_id]
            logger.info(f"For {species}, selected {best_seq_id} by BLAST (score: {sequence_scores.get(best_seq_id, 'N/A')}).")

        # Clean up species-specific temp dir after processing
        shutil.rmtree(species_temp_dir, ignore_errors=True)

    # Clean up the root temp dir for duplicate selection
    shutil.rmtree(selection_temp_root, ignore_errors=True)
    return final_sequences

def convert_duplicates_to_single_sequences(species_to_seqs_map_of_maps, strategy='longest'):
    logger.debug(f"Converting nested duplicate dict with strategy: {strategy}")
    single_seqs = {}
    for species, id_to_seq_map in species_to_seqs_map_of_maps.items():
        if not id_to_seq_map: continue

        if len(id_to_seq_map) == 1:
            single_seqs[species] = next(iter(id_to_seq_map.values()))
        else:
            if strategy == 'longest':
                best_id = max(id_to_seq_map, key=lambda k: len(id_to_seq_map[k]))
                single_seqs[species] = id_to_seq_map[best_id]
            elif strategy == 'first':
                first_id = next(iter(id_to_seq_map.keys()))
                single_seqs[species] = id_to_seq_map[first_id]
            else:
                logger.warning(f"Unsupported conversion strategy '{strategy}' for {species}. Defaulting to 'longest'.")
                best_id = max(id_to_seq_map, key=lambda k: len(id_to_seq_map[k]))
                single_seqs[species] = id_to_seq_map[best_id]
    return single_seqs
