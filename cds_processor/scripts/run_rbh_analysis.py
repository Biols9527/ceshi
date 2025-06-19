import argparse
import logging
import os
import sys
import shutil # For cleanup
import uuid # For unique temp dir names

# Adjust path to import from the cds_processor package
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR) # Should be "cds_processor" main directory
PACKAGE_ROOT = os.path.join(PROJECT_ROOT, "cds_processor") # "cds_processor/cds_processor"

if PACKAGE_ROOT not in sys.path:
    sys.path.insert(0, PACKAGE_ROOT)
if PROJECT_ROOT not in sys.path: # To allow finding cds_processor.cds_processor
    sys.path.insert(0, PROJECT_ROOT)

# Import modules from the cds_processor package
try:
    from cds_processor import utils
    from cds_processor import rbh_handler
    from cds_processor import blast_handler # Needed by rbh_handler.perform_blast_for_rbh
    from cds_processor import fasta_parser # For extract_species_name
    from Bio import SeqIO
except ImportError as e:
    print(f"Error importing cds_processor modules or Biopython: {e}", file=sys.stderr)
    sys.exit(1)

logger = None # Will be configured in main

def extract_sequences_for_species(original_fasta_file, species_name_target, output_temp_fasta_path, fasta_parser_module):
    """
    Extracts sequences for a specific species from a FASTA file.
    If species_name_target is None, all sequences are copied.
    Writes them to output_temp_fasta_path.
    Returns a dictionary of {seq_id: length} for the extracted sequences.
    """
    selected_seq_lengths = {}
    records_to_write = []
    logger.info(f"Extracting sequences from {original_fasta_file} for species '{species_name_target if species_name_target else 'all'}'.")

    try:
        for record in SeqIO.parse(original_fasta_file, "fasta"):
            if species_name_target:
                current_species = fasta_parser_module.extract_species_name(record.description)
                if current_species == species_name_target:
                    records_to_write.append(record)
                    selected_seq_lengths[str(record.id)] = len(record.seq)
            else:
                # No specific species, take all
                records_to_write.append(record)
                selected_seq_lengths[str(record.id)] = len(record.seq)

        if not records_to_write:
            logger.warning(f"No sequences found for species '{species_name_target}' in {original_fasta_file}.")
            # Create an empty fasta file to prevent downstream errors if one input is empty
            with open(output_temp_fasta_path, "w") as f_out:
                pass # Write empty file
            return selected_seq_lengths, False # Indicate no sequences found

        SeqIO.write(records_to_write, output_temp_fasta_path, "fasta")
        logger.info(f"Wrote {len(records_to_write)} sequences to {output_temp_fasta_path}.")
        return selected_seq_lengths, True

    except FileNotFoundError:
        logger.error(f"FASTA file not found: {original_fasta_file}")
        with open(output_temp_fasta_path, "w") as f_out: pass # Write empty file
        return {}, False
    except Exception as e:
        logger.error(f"Error processing FASTA file {original_fasta_file} for species extraction: {e}", exc_info=True)
        with open(output_temp_fasta_path, "w") as f_out: pass # Write empty file
        return {}, False


def main():
    global logger
    parser = argparse.ArgumentParser(description="Reciprocal Best Hit (RBH) Analysis for CDS sequences.")

    parser.add_argument("--input_fasta_A", required=True, help="Path to the first input FASTA file (Set A).")
    parser.add_argument("--input_fasta_B", required=True, help="Path to the second input FASTA file (Set B).")
    parser.add_argument("--species_A", default=None, help="Specific species name from input_fasta_A to use (optional).")
    parser.add_argument("--species_B", default=None, help="Specific species name from input_fasta_B to use (optional).")
    parser.add_argument("--output_rbh_file", required=True, help="Path to the output TSV file for RBH pairs.")

    parser.add_argument("--blastn_path", required=True, help="Path to blastn executable.")
    parser.add_argument("--makeblastdb_path", required=True, help="Path to makeblastdb executable.")

    parser.add_argument("--temp_dir", default=None, help="Directory for temporary files. Defaults to a new dir in output_rbh_file's directory.")
    parser.add_argument("--log_level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], default='INFO', help="Logging level.")
    parser.add_argument("--clean_temp", action="store_true", default=False, help="Clean up temporary directory after execution.")
    parser.add_argument("--blast_evalue", type=float, default=0.01, help="E-value threshold for BLAST searches.") # Example for future config

    args = parser.parse_args()

    # Configure logger
    log_file_name = os.path.splitext(args.output_rbh_file)[0] + "_rbh_analysis.log"
    logger = utils.setup_logging(log_level_str=args.log_level, log_file=log_file_name)
    logger.info(f"Starting RBH Analysis with arguments: {args}")

    # Setup temp directory
    if args.temp_dir:
        temp_dir_base = args.temp_dir
    else:
        temp_dir_base = os.path.join(os.path.dirname(os.path.abspath(args.output_rbh_file)), "rbh_temp_" + uuid.uuid4().hex)

    try:
        os.makedirs(temp_dir_base, exist_ok=True)
        logger.info(f"Temporary files will be stored in: {temp_dir_base}")
    except OSError as e:
        logger.error(f"Failed to create temporary directory {temp_dir_base}: {e}")
        sys.exit(1)

    # ---- Main Logic ----
    # 1. Prepare input sequences for A and B
    temp_fasta_A_path = os.path.join(temp_dir_base, "input_A_selected.fasta")
    temp_fasta_B_path = os.path.join(temp_dir_base, "input_B_selected.fasta")

    lengths_A, found_A = extract_sequences_for_species(args.input_fasta_A, args.species_A, temp_fasta_A_path, fasta_parser)
    lengths_B, found_B = extract_sequences_for_species(args.input_fasta_B, args.species_B, temp_fasta_B_path, fasta_parser)

    if not found_A or not lengths_A:
        logger.error(f"No sequences selected for Set A from {args.input_fasta_A} (species: {args.species_A}). Cannot proceed.")
        if args.clean_temp: shutil.rmtree(temp_dir_base, ignore_errors=True)
        sys.exit(1)
    if not found_B or not lengths_B:
        logger.error(f"No sequences selected for Set B from {args.input_fasta_B} (species: {args.species_B}). Cannot proceed.")
        if args.clean_temp: shutil.rmtree(temp_dir_base, ignore_errors=True)
        sys.exit(1)

    # 2. Run reciprocal BLAST
    #    (blast_handler_module, utils_module are already imported as blast_handler, utils)
    logger.info("Performing BLAST: Set A vs Set B...")
    xml_A_vs_B_path = rbh_handler.perform_blast_for_rbh(
        query_fasta_path=temp_fasta_A_path,
        target_fasta_path=temp_fasta_B_path,
        blastn_exe=args.blastn_path,
        makeblastdb_exe=args.makeblastdb_path,
        blast_handler_module=blast_handler,
        temp_dir_base=os.path.join(temp_dir_base, "A_vs_B_blast"), # Specific sub-temp for this BLAST run
        utils_module=utils
    )
    if not xml_A_vs_B_path:
        logger.error("BLAST A vs B failed. See logs for details. Exiting.")
        if args.clean_temp: shutil.rmtree(temp_dir_base, ignore_errors=True)
        sys.exit(1)

    logger.info("Performing BLAST: Set B vs Set A...")
    xml_B_vs_A_path = rbh_handler.perform_blast_for_rbh(
        query_fasta_path=temp_fasta_B_path,
        target_fasta_path=temp_fasta_A_path,
        blastn_exe=args.blastn_path,
        makeblastdb_exe=args.makeblastdb_path,
        blast_handler_module=blast_handler,
        temp_dir_base=os.path.join(temp_dir_base, "B_vs_A_blast"), # Specific sub-temp
        utils_module=utils
    )
    if not xml_B_vs_A_path:
        logger.error("BLAST B vs A failed. See logs for details. Exiting.")
        if args.clean_temp: shutil.rmtree(temp_dir_base, ignore_errors=True)
        sys.exit(1)

    # 3. Parse BLAST results
    logger.info("Parsing BLAST results for A vs B...")
    hits_A_vs_B = rbh_handler.parse_blast_xml(xml_A_vs_B_path, lengths_A, lengths_B)
    logger.info("Parsing BLAST results for B vs A...")
    hits_B_vs_A = rbh_handler.parse_blast_xml(xml_B_vs_A_path, lengths_B, lengths_A)

    # 4. Find RBH pairs
    logger.info("Finding RBH pairs...")
    rbh_pairs = rbh_handler.find_rbh_pairs(hits_A_vs_B, hits_B_vs_A)

    # 5. Write output
    logger.info(f"Writing {len(rbh_pairs)} RBH pairs to {args.output_rbh_file}...")
    try:
        with open(args.output_rbh_file, "w") as f_out:
            f_out.write("SeqA_ID\tSeqB_ID\tEvalue_A_vs_B\tBitscore_A_vs_B\tPIdent_A_vs_B\tAlignLength_A_vs_B\t"
                        "Evalue_B_vs_A\tBitscore_B_vs_A\tPIdent_B_vs_A\tAlignLength_B_vs_A\t"
                        "SeqA_Length\tSeqB_Length\tQueryA_Start_End\tSubjectB_Start_End\t"
                        "QueryB_Start_End\tSubjectA_Start_End\n")
            for pair in rbh_pairs:
                id_A, id_B, hit_A_to_B, hit_B_to_A = pair
                f_out.write(f"{id_A}\t{id_B}\t"
                            f"{hit_A_to_B.evalue:.2e}\t{hit_A_to_B.bitscore:.1f}\t{hit_A_to_B.pident:.2f}\t{hit_A_to_B.length}\t"
                            f"{hit_B_to_A.evalue:.2e}\t{hit_B_to_A.bitscore:.1f}\t{hit_B_to_A.pident:.2f}\t{hit_B_to_A.length}\t"
                            f"{hit_A_to_B.query_len}\t{hit_A_to_B.subject_len}\t" # Note: subject_len for A_to_B is B's length
                            f"{hit_A_to_B.query_start}-{hit_A_to_B.query_end}\t{hit_A_to_B.subject_start}-{hit_A_to_B.subject_end}\t"
                            f"{hit_B_to_A.query_start}-{hit_B_to_A.query_end}\t{hit_B_to_A.subject_start}-{hit_B_to_A.subject_end}\n")
        logger.info("RBH analysis complete.")
    except IOError as e:
        logger.error(f"Failed to write RBH output file {args.output_rbh_file}: {e}")
        if args.clean_temp: shutil.rmtree(temp_dir_base, ignore_errors=True)
        sys.exit(1)

    # 6. Cleanup
    if args.clean_temp:
        logger.info(f"Cleaning up temporary directory: {temp_dir_base}")
        shutil.rmtree(temp_dir_base, ignore_errors=True)
    else:
        logger.info(f"Temporary files retained at {temp_dir_base}")

if __name__ == "__main__":
    main()
