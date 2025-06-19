import logging
from Bio import SeqIO
from collections import defaultdict
# utils will be needed for logger, GENETIC_CODE etc.
# blast_handler might be needed if parse_fasta calls preprocess_cds_sequence which uses BLAST
# sequence_validator for check_coding_sequence
# orf_analyzer for analyze_codon_structure

logger = logging.getLogger(__name__) # cds_processor.fasta_parser

def extract_species_name(header):
    """Extracts species name from FASTA header."""
    logger.debug(f"Extracting species name from header: {header}")
    if header.startswith(">"):
        header = header[1:]
    parts = header.split(" ", 1)
    if len(parts) > 1:
        species_name = parts[1].strip()
        logger.debug(f"Extracted species name: {species_name}")
        return species_name
    else:
        # If no space, use the whole identifier as species name (or decide on an error/default)
        species_name = header.strip()
        logger.warning(f"No space in header '{header}', using full ID as species name: {species_name}")
        return species_name

# Placeholder for the main parsing function. Will be complex.
def parse_fasta_file(file_path, duplicate_strategy="longest", blast_params=None, cds_processor_config=None):
    """
    Parses a FASTA file, processes sequences, and handles duplicates.
    This is a simplified placeholder. The full logic from the legacy script
    (parse_fasta, process_fasta_batch, parse_fasta_with_duplicates, process_duplicate_batch,
    preprocess_cds_sequence) will need to be integrated here or in helper functions.
    """
    logger.info(f"Starting FASTA parsing for {file_path} with strategy: {duplicate_strategy}")

    # This will eventually call preprocess_cds_sequence, which in turn might call
    # sequence_validator.check_coding_sequence, orf_analyzer.analyze_codon_structure,
    # and blast_handler.handle_non_triplet_with_blast or sequence_validator.handle_non_triplet_sequence

    species_seqs = {}
    species_counts = defaultdict(int)
    # original_ids = {} # If needed for strategies like rename

    try:
        for record in SeqIO.parse(file_path, "fasta"):
            species = extract_species_name(record.description) # Use .description for full header
            sequence = str(record.seq)
            seq_id = record.id

            # Simplified preprocessing call (actual one is more complex)
            # processed_sequence, _ = preprocess_cds_sequence(sequence, seq_id, file_path, blast_params, cds_processor_config)
            # For now, just use the raw sequence for basic migration
            processed_sequence = sequence

            # Simplified duplicate handling (actual one is in process_fasta_batch)
            if species in species_seqs:
                species_counts[species] += 1
                if duplicate_strategy == "longest":
                    if len(processed_sequence) > len(species_seqs[species]):
                        logger.debug(f"Replacing sequence for {species} (new length {len(processed_sequence)} > old length {len(species_seqs[species])})")
                        species_seqs[species] = processed_sequence
                elif duplicate_strategy == "first":
                    logger.debug(f"Keeping first sequence for {species}, ignoring duplicate ID {seq_id}")
                    pass # Keep the first one already there
                # Add other strategies like "rename", "alignment_quality" later
                else:
                    logger.warning(f"Duplicate species {species} (ID: {seq_id}) found, strategy '{duplicate_strategy}' leads to keeping current.")
                    # By default, or if strategy is unclear, might keep first or error
                    pass # Default to keeping first or implement error
            else:
                species_seqs[species] = processed_sequence
                species_counts[species] = 1

        logger.info(f"Successfully parsed {file_path}. Found {len(species_seqs)} unique species.")
        return species_seqs

    except FileNotFoundError:
        logger.error(f"FASTA file not found: {file_path}")
        return {}
    except Exception as e:
        logger.error(f"Error parsing FASTA file {file_path}: {str(e)}", exc_info=True)
        return {}
