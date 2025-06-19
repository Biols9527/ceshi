import argparse
import logging
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

# Adjust path to import from the parent directory (cds_processor)
# This assumes process_cds.py is in cds_processor/scripts/
# and the package is cds_processor/cds_processor/
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR) # This should be the root "cds_processor" dir
PACKAGE_ROOT = os.path.join(PROJECT_ROOT, "cds_processor") # This is "cds_processor/cds_processor"

if PACKAGE_ROOT not in sys.path:
    sys.path.insert(0, PACKAGE_ROOT)
if PROJECT_ROOT not in sys.path: # To allow finding cds_processor.cds_processor
    sys.path.insert(0, PROJECT_ROOT)

# Import modules from the cds_processor package
from cds_processor import utils
from cds_processor import file_manager
from cds_processor import fasta_parser # Will contain the main parsing and preprocessing logic
from cds_processor import duplicate_handler # For alignment_quality final selection
# Other modules like sequence_validator, orf_analyzer, blast_handler will be used by fasta_parser internally

# Initialize a logger for the script itself
# The setup_logging from utils will configure the root logger and cds_processor package logger
logger = utils.setup_logging(log_level_str="INFO", log_file="cds_main_processing.log")

# Placeholder for a global config object if needed later
CDS_PROCESSOR_CONFIG = {
    "min_repeat_report": utils.DEFAULT_MIN_REPEAT_COUNT_TO_REPORT,
    "blast_db_cache_size": 5, # From legacy
    # Paths for BLAST executables will come from args
}

def process_single_file_wrapper(args_tuple):
    """
    Wrapper for process_single_cds_file to be used with ThreadPoolExecutor.
    Unpacks arguments for the actual processing function.
    """
    input_file_path, output_file_path, duplicate_strategy, base_blast_params, config_params, main_preprocess_func = args_tuple

    # Create a specific blast_params dict for this file if needed (e.g., if temp_dir is file-specific)
    # For now, base_blast_params might be sufficient if temp_dir is shared or managed inside blast_handler
    current_blast_params = base_blast_params.copy() if base_blast_params else {}

    try:
        logger.info(f"Processing file: {input_file_path}")

        # This is where the main parsing and processing logic for a single file resides.
        # It will call functions from fasta_parser, which in turn use other modules.
        # The legacy `parse_fasta` and its helpers had this logic.
        # We need a function in `fasta_parser` that encapsulates this.
        # Let's assume `fasta_parser.parse_and_process_file` exists or will be created.

        # For now, using the placeholder from fasta_parser.py
        # This will need the full preprocess_cds_sequence logic integrated within it.
        all_species_data = fasta_parser.parse_fasta_file(
            input_file_path,
            duplicate_strategy=duplicate_strategy,
            blast_params=current_blast_params, # Contains paths to BLAST tools, temp_dir, etc.
            cds_processor_config=config_params # Global settings like repeat count
        )

        # If duplicate_strategy was "alignment_quality", all_species_data would be
        # {species: {id: sequence}}. It needs further processing by duplicate_handler.
        if duplicate_strategy == "alignment_quality":
            if not current_blast_params or not current_blast_params.get("use_blast_for_duplicates"):
                logger.warning(f"Alignment quality strategy for {input_file_path} but BLAST for duplicates not enabled. Defaulting to longest among duplicates.")
            # This function will use BLAST (if enabled) or fallback to longest
            processed_species_data = duplicate_handler.select_best_by_alignment_quality(
                all_species_data, # This should be the {species: {id:seq}} map
                blast_params=current_blast_params,
                config=config_params
            )
        else:
            # For other strategies, parse_fasta_file should return {species: sequence}
            processed_species_data = all_species_data

        if not processed_species_data:
            logger.error(f"No sequences processed from {input_file_path}")
            return None

        success = file_manager.write_output(processed_species_data, output_file_path)
        if not success:
            logger.error(f"Failed to write output for {input_file_path} to {output_file_path}")
            return None

        logger.info(f"Successfully processed {input_file_path} -> {output_file_path} ({len(processed_species_data)} species)")
        return {
            "input_file": input_file_path,
            "output_file": output_file_path,
            "species_count": len(processed_species_data),
            # "sequences": processed_species_data # Potentially large, maybe omit from summary
        }
    except Exception as e:
        logger.error(f"Error processing single file {input_file_path}: {str(e)}", exc_info=True)
        return None

def main():
    parser = argparse.ArgumentParser(description="CDS Sequence Processor - Refactored")
    parser.add_argument("--input_dir", required=True, help="Directory containing input .cds files.")
    parser.add_argument("--output_dir", required=True, help="Directory to save processed files.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing.")
    parser.add_argument("--log_level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], default='INFO', help="Logging level.")

    # Duplicate handling strategy
    parser.add_argument("--duplicate_strategy",
                        choices=['longest', 'first', 'rename', 'alignment_quality'],
                        default='longest',
                        help="Strategy for handling duplicate species in a file.")

    # BLAST related arguments
    parser.add_argument("--use_blast_for_non_triplets", action="store_true", help="Use BLAST to aid fixing non-triplet sequences.")
    parser.add_argument("--use_blast_for_duplicates", action="store_true", help="Use BLAST for 'alignment_quality' duplicate strategy (Not fully implemented).")
    parser.add_argument("--blastn_path", help="Path to blastn executable.")
    parser.add_argument("--makeblastdb_path", help="Path to makeblastdb executable.")
    # parser.add_argument("--blast_db_cache_size", type=int, default=5, help="Max BLAST DBs to cache in memory.") # Handled by blast_handler global

    # Other config from legacy
    parser.add_argument("--min_repeat_report", type=int, default=utils.DEFAULT_MIN_REPEAT_COUNT_TO_REPORT,
                        help="Min repeat count to report as a quality issue.")
    parser.add_argument("--clean_temp_on_finish", action="store_true", default=True, help="Clean up temp directory after processing.")

    args = parser.parse_args()

    # Update logger level based on args
    logger = utils.setup_logging(log_level_str=args.log_level, log_file=os.path.join(args.output_dir, "cds_processing_run.log"))
    logger.info(f"Starting CDS processing with arguments: {args}")
    start_time = time.time()

    # Update global config based on args
    CDS_PROCESSOR_CONFIG["min_repeat_report"] = args.min_repeat_report
    # CDS_PROCESSOR_CONFIG["blast_db_cache_size"] = args.blast_db_cache_size # Global in blast_handler

    # Setup directories
    processed_output_dir = os.path.join(args.output_dir, "processed_cds")
    temp_processing_dir = os.path.join(args.output_dir, "temp_processing_data")
    try:
        os.makedirs(processed_output_dir, exist_ok=True)
        os.makedirs(temp_processing_dir, exist_ok=True)
        logger.info(f"Output will be in: {processed_output_dir}")
        logger.info(f"Temporary files in: {temp_processing_dir}")
    except OSError as e:
        logger.error(f"Failed to create output/temp directories: {e}")
        sys.exit(1)

    # Prepare BLAST parameters dictionary to pass around
    blast_params = None
    if args.use_blast_for_non_triplets or args.use_blast_for_duplicates:
        if not args.blastn_path or not args.makeblastdb_path:
            logger.error("BLAST usage requested, but --blastn_path or --makeblastdb_path not provided. Exiting.")
            sys.exit(1)
        blast_params = {
            "use_blast_for_non_triplets": args.use_blast_for_non_triplets,
            "use_blast_for_duplicates": args.use_blast_for_duplicates,
            "blastn_exe": os.path.abspath(args.blastn_path),
            "makeblastdb_exe": os.path.abspath(args.makeblastdb_path),
            "temp_dir": temp_processing_dir, # Main temp dir for BLAST operations
            # "db_cache_path": os.path.join(temp_processing_dir, "blast_db_cache") # If disk cache for DBs needed
        }
        logger.info(f"BLAST operations enabled with params: {blast_params}")

    input_cds_files = file_manager.find_cds_files(args.input_dir)
    if not input_cds_files:
        logger.warning(f"No .cds files found in {args.input_dir}. Exiting.")
        sys.exit(0)

    logger.info(f"Found {len(input_cds_files)} CDS files to process.")

    # Prepare arguments for parallel processing
    # The `main_preprocess_func` is the core logic that would have been in `preprocess_cds_sequence`
    # and then `parse_fasta` in the legacy script. This needs to be available in `fasta_parser`.
    # For now, this is a conceptual link.
    # Assume `fasta_parser.full_sequence_processing_pipeline` is this function.
    # This function would take (raw_sequence, seq_id, file_path_context, blast_params, config)
    # and return (processed_sequence_string, analysis_dict).
    # This is a GAPING HOLE in the current migration - the main `preprocess_cds_sequence` from legacy
    # is not yet fully refactored into `fasta_parser` or a similar orchestrator module.
    # The current `fasta_parser.parse_fasta_file` is too simplistic.
    # For this step, we will assume it exists and is passed conceptually.

    # A placeholder for what `fasta_parser.preprocess_cds_sequence` would be:
    def placeholder_preprocess_function(sequence, seq_id, file_path, blast_params_for_call, config_for_call):
        # This function needs to call:
        # 1. orf_analyzer.analyze_codon_structure
        # 2. Adjust frame: sequence = sequence[structure['best_frame']:]
        # 3. Handle non-triplet:
        #    if blast_params_for_call and blast_params_for_call.get('use_blast_for_non_triplets'):
        #        sequence = blast_handler.handle_non_triplet_with_blast(...)
        #    else:
        #        sequence = sequence_validator.handle_non_triplet_sequence(...)
        # 4. sequence_validator.assess_sequence_quality (optional, for logging)
        # 5. sequence_validator.check_coding_sequence (optional, for logging)
        logger.debug(f"Placeholder preprocessing for {seq_id} from {file_path}")
        # Simulate basic processing: ensure it's a string
        # Actual processing is much more involved.
        if len(sequence) % 3 != 0:
             # Basic padding if not multiple of 3, mimicking non-BLAST fallback
             from cds_processor import sequence_validator # Temp import
             sequence = sequence_validator.handle_non_triplet_sequence(sequence, seq_id, position="end")
        return str(sequence), {"best_frame": 0} # Return string and dummy dict

    tasks_args = []
    for cds_file_name in input_cds_files:
        input_file = os.path.join(args.input_dir, cds_file_name)
        output_file = os.path.join(processed_output_dir, f"{os.path.splitext(cds_file_name)[0]}_processed.fasta")
        tasks_args.append((
            input_file,
            output_file,
            args.duplicate_strategy,
            blast_params,
            CDS_PROCESSOR_CONFIG,
            placeholder_preprocess_function # This is the critical function to pass
        ))

    results_summary = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        # Using map to preserve order for results might be good if needed, but future returns out of order
        futures = [executor.submit(process_single_file_wrapper, task_arg_set) for task_arg_set in tasks_args]

        for i, future in enumerate(futures):
            try:
                result = future.result() # Blocks until this future is done
                if result:
                    results_summary.append(result)
                logger.info(f"Progress: {i+1}/{len(input_cds_files)} files processed.")
            except Exception as e_future:
                logger.error(f"Exception processing file for task {tasks_args[i][0]}: {e_future}", exc_info=True)

    successful_files = len(results_summary)
    logger.info(f"Processing complete. {successful_files}/{len(input_cds_files)} files processed successfully.")

    if args.clean_temp_on_finish:
        logger.info(f"Cleaning up temporary directory: {temp_processing_dir}")
        file_manager.cleanup_temp_directory(temp_processing_dir)
    else:
        logger.info(f"Temporary files retained at: {temp_processing_dir}")

    end_time = time.time()
    logger.info(f"Total execution time: {end_time - start_time:.2f} seconds.")
    logger.info(f"Processed files are in: {processed_output_dir}")

if __name__ == "__main__":
    # Ensure Biopython is available if actual parsing is done
    try:
        from Bio import SeqIO
    except ImportError:
        print("Biopython (Bio module) is not installed. Please install it: pip install biopython", file=sys.stderr)
        # logger might not be configured yet if this fails early
        # logger.critical("Biopython (Bio module) not found. Please install it.")
        sys.exit(1)

    # For BLAST parsing
    try:
        from Bio.Blast import NCBIXML
    except ImportError:
        print("Biopython (Bio.Blast.NCBIXML) is not available. Please ensure full Biopython install.", file=sys.stderr)
        sys.exit(1)

    main()
