import os
import shutil
import logging

logger = logging.getLogger(__name__)

def find_cds_files(directory="."):
    """Finds all .cds files in the given directory."""
    logger.info(f"Searching for .cds files in directory: {directory}")
    try:
        if not os.path.isdir(directory):
            logger.error(f"Provided path is not a directory: {directory}")
            return []
        cds_files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and f.endswith(".cds")]
        logger.info(f"Found {len(cds_files)} .cds files in {directory}: {cds_files}")
        return sorted(cds_files)
    except Exception as e:
        logger.error(f"Error finding CDS files in {directory}: {str(e)}", exc_info=True)
        return []

def write_output(species_seqs_data, output_file_path):
    """Writes processed sequences to a FASTA file."""
    logger.info(f"Preparing to write sequences to {output_file_path}")
    if not species_seqs_data:
        logger.warning(f"No sequences provided to write to {output_file_path}")
        return False

    try:
        output_dir = os.path.dirname(os.path.abspath(output_file_path))
        if output_dir and not os.path.exists(output_dir):
            logger.info(f"Creating output directory: {output_dir}")
            os.makedirs(output_dir, exist_ok=True)

        sequences_to_write = species_seqs_data

        with open(output_file_path, "w") as f_out:
            for species_name, sequence_str in sorted(sequences_to_write.items()):
                if not isinstance(sequence_str, str):
                    logger.error(f"Sequence for species {species_name} is not a string ({type(sequence_str)}). Skipping.")
                    continue
                f_out.write(f">{species_name}\n")
                for i in range(0, len(sequence_str), 60):
                    f_out.write(sequence_str[i:i+60] + "\n")

        num_written = len(sequences_to_write)
        logger.info(f"Successfully wrote {num_written} sequences to {output_file_path}")
        return True
    except Exception as e:
        logger.error(f"Error writing output to {output_file_path}: {str(e)}", exc_info=True)
        return False

def cleanup_temp_directory(temp_dir_to_remove):
    """Removes the specified temporary directory and all its contents."""
    logger.info(f"Attempting to clean up temporary directory: {temp_dir_to_remove}")
    if not temp_dir_to_remove or not isinstance(temp_dir_to_remove, str):
        logger.warning(f"Invalid temporary directory path provided: {temp_dir_to_remove}. Skipping cleanup.")
        return

    if not os.path.exists(temp_dir_to_remove):
        logger.info(f"Temporary directory {temp_dir_to_remove} does not exist. Nothing to clean.")
        return

    if not os.path.isdir(temp_dir_to_remove):
        logger.warning(f"Path {temp_dir_to_remove} is not a directory. Skipping cleanup.")
        return

    try:
        shutil.rmtree(temp_dir_to_remove)
        logger.info(f"Successfully removed temporary directory: {temp_dir_to_remove}")
    except Exception as e:
        logger.error(f"Error cleaning up temporary directory {temp_dir_to_remove}: {str(e)}", exc_info=True)
