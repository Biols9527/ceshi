import logging
import hashlib
import os
import sys # Added sys for logger stream handler
from functools import lru_cache

# Define 4-fold degenerate codons
FOURFOLD_CODONS = {
    "GC": True, "CG": True, "GG": True, "CT": True,
    "CC": True, "TC": True, "AC": True, "GT": True
}

# Define genetic code table
GENETIC_CODE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M", "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K", "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L", "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q", "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V", "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E", "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S", "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*", "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
}

DEFAULT_MIN_REPEAT_COUNT_TO_REPORT = 70

logger = logging.getLogger(__name__)

def setup_logging(log_level_str="INFO", log_file="cds_processing.log"):
    log_level = getattr(logging, log_level_str.upper(), logging.INFO)

    root_logger_instance = logging.getLogger()
    package_logger_instance = logging.getLogger("cds_processor")

    for lg in [root_logger_instance, package_logger_instance]:
        if lg.hasHandlers():
            lg.handlers.clear()

    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] [%(name)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

    package_logger_instance.setLevel(log_level)

    logger.info(f"Logging setup complete. Level: {log_level_str}, File: {log_file}")
    return package_logger_instance

def hash_file_content(file_path):
    logger.debug(f"Calculating MD5 hash for file: {file_path}")
    try:
        hasher = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                hasher.update(chunk)
        hex_digest = hasher.hexdigest()
        logger.debug(f"MD5 hash for {file_path}: {hex_digest}")
        return hex_digest
    except FileNotFoundError:
        logger.error(f"File not found for hashing: {file_path}")
        return None
    except Exception as e:
        logger.error(f"Error calculating file hash for {file_path}: {str(e)}")
        return None
