import logging
import os
import subprocess
import uuid
import shutil # For cleanup of non-cached DBs if needed

# from .utils import hash_file_content # For DB caching
# from .file_manager import cleanup_temp_files # For cleaning up temp query/output files
# from Bio.Blast import NCBIXML # For parsing BLAST output

logger = logging.getLogger(__name__) # cds_processor.blast_handler

# Global cache for BLAST database paths to avoid recreating them for the same file content
# Key: file_hash, Value: db_prefix_path
BLAST_DB_CACHE = {}
BLAST_DB_CACHE_SIZE = 5 # Max number of DBs to cache

def create_temp_blast_db(target_fasta_path, temp_dir, makeblastdb_executable_path, config=None):
    """Creates a temporary BLAST database from a FASTA file."""
    logger.info(f"Creating temporary BLAST database for: {target_fasta_path} in {temp_dir}")

    # Ensure makeblastdb_executable_path is provided and valid
    if not makeblastdb_executable_path or not os.path.isfile(makeblastdb_executable_path):
        logger.error(f"makeblastdb executable not found or invalid path: {makeblastdb_executable_path}")
        return None

    try:
        os.makedirs(temp_dir, exist_ok=True)

        # Use a unique name for the database to avoid collisions if temp_dir is shared or not cleaned properly
        db_name = f"temp_blastdb_{uuid.uuid4().hex}"
        db_prefix = os.path.join(temp_dir, db_name)

        # Command for makeblastdb
        cmd = [
            makeblastdb_executable_path,
            "-in", target_fasta_path,
            "-dbtype", "nucl",
            "-out", db_prefix,
            "-title", db_name
        ]
        logger.debug(f"makeblastdb command: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True, check=False) # check=False to handle errors manually

        if result.returncode != 0:
            logger.error(f"makeblastdb failed for {target_fasta_path}. Error: {result.stderr}")
            # Attempt to clean up partial DB files if creation failed
            for ext in [".nhr", ".nin", ".nsq", ".fasta", ".nal", ".log"]:
                if os.path.exists(db_prefix + ext):
                    try:
                        os.remove(db_prefix + ext)
                    except OSError:
                        pass # Ignore if removal fails, main error is already logged
            return None

        logger.info(f"Successfully created BLAST database: {db_prefix}")
        return db_prefix # Return the prefix (path without extension)

    except Exception as e:
        logger.error(f"Exception during BLAST DB creation for {target_fasta_path}: {str(e)}", exc_info=True)
        return None

def create_or_get_blast_db(target_fasta_path, temp_dir, makeblastdb_executable_path, utils_module, config=None):
    """Creates a BLAST DB or retrieves it from cache if content matches."""
    # utils_module is passed to access hash_file_content
    global BLAST_DB_CACHE, BLAST_DB_CACHE_SIZE

    file_hash = utils_module.hash_file_content(target_fasta_path)
    if not file_hash:
        logger.warning(f"Could not hash file {target_fasta_path}. Proceeding to create DB without caching.")
        # Fallback to creating DB without caching for this specific call if hashing fails
        return create_temp_blast_db(target_fasta_path, temp_dir, makeblastdb_executable_path, config)

    if file_hash in BLAST_DB_CACHE:
        db_prefix = BLAST_DB_CACHE[file_hash]
        # Verify that the database files still exist (e.g., not cleaned up by another process or OS)
        if os.path.exists(db_prefix + ".nsq"): # .nsq is a common BLAST db file extension
            logger.info(f"Using cached BLAST database for {target_fasta_path}: {db_prefix}")
            return db_prefix
        else:
            logger.warning(f"Cached BLAST DB for {target_fasta_path} (hash {file_hash}) not found at {db_prefix}. Recreating.")
            del BLAST_DB_CACHE[file_hash] # Remove stale entry

    # If not in cache or files are missing, create it
    db_prefix = create_temp_blast_db(target_fasta_path, temp_dir, makeblastdb_executable_path, config)
    if db_prefix:
        # Manage cache size: if full, remove the oldest entry (FIFO for simplicity here)
        if len(BLAST_DB_CACHE) >= BLAST_DB_CACHE_SIZE:
            oldest_hash_key = next(iter(BLAST_DB_CACHE)) # Get the first key (oldest in Python 3.7+ dicts)
            oldest_db_prefix_to_remove = BLAST_DB_CACHE.pop(oldest_hash_key)
            logger.info(f"BLAST DB cache full. Removing oldest entry: {oldest_db_prefix_to_remove} (hash {oldest_hash_key})")
            # Actively clean up files of the removed database to save space
            # This assumes the db_prefix is just the path and name, extensions are standard
            for ext in [".nhr", ".nin", ".nsq", ".fasta", ".nal", ".log"]: # Add other extensions if used by makeblastdb
                db_file_to_remove = oldest_db_prefix_to_remove + ext
                if os.path.exists(db_file_to_remove):
                    try:
                        os.remove(db_file_to_remove)
                        logger.debug(f"Removed cached DB file: {db_file_to_remove}")
                    except OSError as e_remove:
                        logger.warning(f"Could not remove cached DB file {db_file_to_remove}: {e_remove}")

        BLAST_DB_CACHE[file_hash] = db_prefix
        logger.info(f"Cached new BLAST database for {target_fasta_path} (hash {file_hash}): {db_prefix}")
    return db_prefix

def run_blastn(query_sequence_path, db_prefix_path, blastn_executable_path, output_dir, config=None):
    """Runs BLASTN and returns the path to the XML output file."""
    logger.info(f"Running BLASTN for query {query_sequence_path} against DB {db_prefix_path}")

    if not blastn_executable_path or not os.path.isfile(blastn_executable_path):
        logger.error(f"blastn executable not found or invalid path: {blastn_executable_path}")
        return None
    if not os.path.exists(query_sequence_path):
        logger.error(f"Query sequence file not found: {query_sequence_path}")
        return None
    if not os.path.exists(db_prefix_path + ".nsq"): # Check if DB seems valid
        logger.error(f"BLAST database not found or incomplete: {db_prefix_path}")
        return None

    try:
        os.makedirs(output_dir, exist_ok=True)
        # Unique output file name
        blast_output_filename = f"blast_out_{uuid.uuid4().hex}.xml"
        blast_output_path = os.path.join(output_dir, blast_output_filename)

        # Standard BLASTN command parameters
        cmd = [
            blastn_executable_path,
            "-query", query_sequence_path,
            "-db", db_prefix_path,
            "-out", blast_output_path,
            "-outfmt", "5",  # XML output format
            "-task", "blastn", # As per legacy script
            "-word_size", "7", # As per legacy script
            "-evalue", "0.01" # As per legacy script
        ]
        # Add other params from config if needed (e.g., num_threads)
        # if config and config.get("blast_threads"):
        # cmd.extend(["-num_threads", str(config.get("blast_threads"))])

        logger.debug(f"BLASTN command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)

        if result.returncode != 0:
            logger.error(f"BLASTN run failed for query {query_sequence_path}. Error: {result.stderr}")
            # Cleanup output file if it was created but empty/error
            if os.path.exists(blast_output_path):
                try: os.remove(blast_output_path)
                except OSError: pass
            return None

        logger.info(f"BLASTN completed. Output at: {blast_output_path}")
        return blast_output_path

    except Exception as e:
        logger.error(f"Exception during BLASTN run for {query_sequence_path}: {str(e)}", exc_info=True)
        return None

# This is the function that was complex in the legacy script.
# It needs access to sequence_validator.handle_non_triplet_sequence for fallback.
# It also needs Bio.Blast.NCBIXML for parsing.
def handle_non_triplet_with_blast(
    sequence, seq_id, original_file_path,
    blastn_exe, makeblastdb_exe, temp_dir,
    utils_module, # For hash_file_content
    sequence_validator_module, # For handle_non_triplet_sequence (fallback)
    biopython_blast_parser_module, # For NCBIXML
    config=None):
    """
    Uses BLAST to determine padding for non-triplet sequences.
    This is a placeholder. The full logic needs careful migration.
    """
    logger.info(f"Attempting BLAST-based non-triplet handling for seq ID: {seq_id}, length: {len(sequence)}")

    if len(sequence) % 3 == 0:
        return sequence # Already a triplet

    # Ensure necessary modules and paths are provided
    if not all([blastn_exe, makeblastdb_exe, temp_dir, original_file_path,
                utils_module, sequence_validator_module, biopython_blast_parser_module]):
        logger.warning(f"Missing parameters for BLAST non-triplet handling of {seq_id}. Using default padding.")
        return sequence_validator_module.handle_non_triplet_sequence(sequence, seq_id, position="smart", context=None)

    # 1. Create or get BLAST DB for the original_file_path
    db_prefix = create_or_get_blast_db(original_file_path, os.path.join(temp_dir, "blast_db_cache"), makeblastdb_exe, utils_module, config)
    if not db_prefix:
        logger.warning(f"Failed to create/get BLAST DB for {original_file_path}. Using default padding for {seq_id}.")
        return sequence_validator_module.handle_non_triplet_sequence(sequence, seq_id, position="smart", context=None)

    # 2. Create temporary FASTA file for the current sequence
    query_fasta_path = os.path.join(temp_dir, f"query_{seq_id}_{uuid.uuid4().hex}.fasta")
    with open(query_fasta_path, "w") as f_query:
        f_query.write(f">{seq_id}\n{sequence}\n")

    # 3. Run BLASTN
    blast_xml_output = run_blastn(query_fasta_path, db_prefix, blastn_exe, os.path.join(temp_dir, "blast_results"), config)

    padded_sequence = None
    if blast_xml_output and os.path.exists(blast_xml_output):
        try:
            with open(blast_xml_output) as result_handle:
                blast_records = list(biopython_blast_parser_module.parse(result_handle))

            if blast_records and blast_records[0].alignments:
                logger.info(f"BLAST hit found for {seq_id}. (Detailed parsing logic to be ported)")
                # --- Start of refined logic block ---
                best_alignment = blast_records[0].alignments[0]
                hsp = best_alignment.hsps[0]  # Best HSP

                # Extract alignment details
                query_start = hsp.query_start
                query_end = hsp.query_end
                # Subject start/end and frame are crucial for frame consistency logic
                sbjct_start = hsp.sbjct_start
                # sbjct_end = hsp.sbjct_end # Not directly used in legacy padding decision, but good for context

                remainder = len(sequence) % 3
                padding_needed = 3 - remainder
                padding_Ns = "N" * padding_needed

                # Decision logic from legacy script:
                # Based on query alignment coverage and frame consistency with subject

                # Case 1: Significant unaligned region at the 5-prime end of the query
                # The legacy used "> 3" which is a bit small. Let's use padding_needed as a more dynamic threshold.
                if query_start > padding_needed + 1: # query_start is 1-based
                    padded_sequence = padding_Ns + sequence
                    logger.info(f"BLAST for {seq_id}: Query aligns significantly downstream ({query_start=}). Padding 5-prime end with {padding_needed} Ns.")

                # Case 2: Significant unaligned region at the 3-prime end of the query
                elif (len(sequence) - query_end) > padding_needed:
                    padded_sequence = sequence + padding_Ns
                    logger.info(f"BLAST for {seq_id}: Query aligns significantly upstream ({query_end=}/{len(sequence)}). Padding 3-prime end with {padding_needed} Ns.")

                # Case 3: Alignment covers most of the sequence, check frame consistency
                else:
                    # Determine frame of subject and query based on their start positions in the alignment
                    # Frame is (start_position - 1) % 3
                    subject_frame_offset = (sbjct_start - 1) % 3
                    query_frame_offset_in_alignment = (query_start - 1) % 3 # query's frame relative to its own start in alignment

                    # We want the original sequence's actual starting frame to match the subject's frame.
                    # The part of the query *before* query_start is (query_start - 1) bases long.
                    # If we pad at the start, we shift the original sequence.
                    # The goal is that after padding, the *new* start of the sequence, when it aligns,
                    # implies a frame consistent with subject_frame_offset.

                    if subject_frame_offset != query_frame_offset_in_alignment:
                        # Frames are inconsistent. Attempt to fix by padding the 5-prime end.
                        # This part of the legacy logic was:
                        # if current_offset < frame_offset: padding = 'N' * (frame_offset - current_offset)
                        # else: padding = 'N' * (3 - (current_offset - frame_offset))
                        # This logic seems to adjust the query's frame to match the subject's.
                        # Let's re-verify. If query is frame 0 (starts at base 1 of alignment) and subject is frame 1 (starts at base 2 of its sequence),
                        # we need to make query behave as if it's frame 1.
                        # query_frame_offset_in_alignment = (query_start - 1) % 3.
                        # This is the frame of the *aligned part* of the query.
                        # The actual sequence `sequence` starts at frame 0 of itself.
                        # The legacy `current_offset` was `(query_start - 1) % 3`.
                        # The legacy `frame_offset` was `(sbjct_start - 1) % 3`.

                        # Let current_actual_frame_start_of_sequence = 0 (implicit for the raw sequence)
                        # We want new_frame_start_of_sequence such that
                        # (new_frame_start_of_sequence + query_start -1) % 3 == subject_frame_offset % 3
                        # (padding_length + query_start -1) % 3 == subject_frame_offset % 3

                        # The legacy logic seems more direct: adjust the existing query_frame_offset_in_alignment
                        # to match subject_frame_offset by adding Ns at the beginning of `sequence`.
                        required_padding_for_frame_match = 0
                        if query_frame_offset_in_alignment < subject_frame_offset:
                            required_padding_for_frame_match = subject_frame_offset - query_frame_offset_in_alignment
                        else: # query_frame_offset_in_alignment > subject_frame_offset
                            required_padding_for_frame_match = 3 - (query_frame_offset_in_alignment - subject_frame_offset)

                        # This padding is for frame alignment. It must also satisfy the triplet requirement.
                        # This is where it gets tricky. The original padding_Ns (for triplet) vs this.
                        # The legacy script's frame adjustment padding:
                        #   `padding = 'N' * (frame_offset - current_offset)` or `padding = 'N' * (3 - (current_offset - frame_offset))`
                        #   This `padding` was then prepended. This could make the sequence *not* a multiple of 3.
                        # This interpretation needs to be careful. The primary goal is to make len(sequence) % 3 == 0.
                        # The frame alignment is secondary or informs *where* the padding_Ns (of length padding_needed) go.

                        # Re-evaluating legacy: it seems it prepends Ns to align frames, then if the result is still not triplet,
                        # it might fall through to the 'smart' handler or default padding.
                        # This part of the legacy logic is the most complex and potentially ambiguous.

                        # Let's stick to the primary goal: use padding_Ns (of length padding_needed).
                        # Where should it go if frames are inconsistent?
                        # If subject_frame_offset = 1, query_frame_offset_in_alignment = 0.
                        # Query is "too early" in frame. Prepending one N would make query_frame_offset_in_alignment = 1.
                        # If padding_needed = 1, then this is perfect: pad 5-prime.
                        # If padding_needed = 2, and we add one N for frame, we still need one N. Add to 3-prime?

                        # A simpler interpretation from legacy: if frames misalign, prefer 5-prime padding using the `padding_Ns`
                        # determined by `padding_needed` for triplet length.
                        padded_sequence = padding_Ns + sequence
                        logger.info(f"BLAST for {seq_id}: Query frame ({query_frame_offset_in_alignment}) and subject frame ({subject_frame_offset}) mismatch. Padding 5-prime with {padding_needed} Ns (heuristic).")

                    else: # Frames are consistent
                        # If frames match, and not covered by Case 1 or 2, the alignment is good.
                        # Decide padding based on start/stop codons in the original sequence (passed via context to smart handler)
                        # or default to 3-prime. The legacy script had a block for this:
                        # has_start = sequence[:3].upper() == 'ATG'
                        # has_stop = sequence[-3:].upper() in ['TAA', 'TAG', 'TGA']
                        # if has_start and not has_stop: padded_seq = sequence + padding
                        # elif not has_start and has_stop: padded_seq = padding + sequence
                        # else: padded_seq = sequence + padding (defaulting to 3-prime)
                        # This means we should call the sequence_validator's smart handler.
                        # The `context` for smart handler is not directly available here.
                        # The legacy called `handle_non_triplet_sequence(sequence, seq_id, position='smart')`
                        # which internally re-checks start/stop.
                        logger.info(f"BLAST for {seq_id}: Query and subject frames match. Alignment good. Using smart fallback for padding.")
                        # This call needs the `sequence_validator_module`
                        padded_sequence = sequence_validator_module.handle_non_triplet_sequence(sequence, seq_id, position="smart", context=None)
                # --- End of refined logic block ---
            else:
                logger.warning(f"No BLAST alignments found for {seq_id}. Using default padding.")
                padded_sequence = sequence_validator_module.handle_non_triplet_sequence(sequence, seq_id, position="smart", context=None)
        except Exception as e_parse:
            logger.error(f"Error parsing BLAST XML for {seq_id}: {e_parse}", exc_info=True)
            padded_sequence = sequence_validator_module.handle_non_triplet_sequence(sequence, seq_id, position="smart", context=None)
        finally:
            if os.path.exists(blast_xml_output): os.remove(blast_xml_output)
    else:
        logger.warning(f"BLAST run failed or produced no output for {seq_id}. Using default padding.")
        padded_sequence = sequence_validator_module.handle_non_triplet_sequence(sequence, seq_id, position="smart", context=None)

    if os.path.exists(query_fasta_path): os.remove(query_fasta_path)

    return padded_sequence if padded_sequence else sequence_validator_module.handle_non_triplet_sequence(sequence, seq_id, position="smart", context=None)
