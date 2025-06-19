import logging
from collections import defaultdict
# from .utils import GENETIC_CODE, DEFAULT_MIN_REPEAT_COUNT_TO_REPORT # Will be needed

logger = logging.getLogger(__name__) # cds_processor.sequence_validator

# Moved from legacy script, MIN_REPEAT_COUNT_TO_REPORT will come from config or utils
# MIN_REPEAT_COUNT_TO_REPORT = 70 # Placeholder, should be configurable

def check_coding_sequence(sequence, genetic_code_map=None):
    """Checks if a sequence appears to be a valid coding sequence."""
    # GENETIC_CODE will be imported from utils or passed as arg
    if genetic_code_map is None:
        # This is a placeholder, proper GENETIC_CODE should be used
        genetic_code_map = {"TAA": "*", "TAG": "*", "TGA": "*"}

    logger.debug(f"Checking coding sequence of length {len(sequence)}")
    if not sequence:
        return False, "Empty sequence"
    if len(sequence) % 3 != 0:
        return False, f"Sequence length {len(sequence)} is not a multiple of 3"

    # Check for internal stop codons
    for i in range(0, len(sequence) - 3, 3): # Exclude the last codon
        codon = sequence[i:i+3].upper()
        if codon in genetic_code_map and genetic_code_map[codon] == "*":
            return False, f"Internal stop codon {codon} found at position {i}"

    # Check for invalid characters (more comprehensive set)
    valid_chars = set("ACGTRYKMSWBDHVN-") # DNA + IUPAC ambiguities + gap
    # Convert sequence to set of unique characters for efficient checking
    sequence_chars = set(sequence.upper())
    invalid_chars = sequence_chars - valid_chars
    if invalid_chars:
        return False, f"Sequence contains invalid characters: {invalid_chars}"

    return True, "Valid coding sequence"

def assess_sequence_quality(sequence, seq_id="Unknown", min_repeat_report=70, config=None):
    """Assesses sequence quality for potential issues."""
    # config might pass MIN_REPEAT_COUNT_TO_REPORT
    # min_repeat_report can also be directly from utils.DEFAULT_MIN_REPEAT_COUNT_TO_REPORT

    logger.debug(f"Assessing quality for sequence ID: {seq_id}, length: {len(sequence)}")
    issues = []
    clean_seq = sequence.upper().replace("-", "")
    if not clean_seq: # Handle empty sequence after cleaning
        issues.append("Sequence is empty after cleaning gaps.")
        return issues

    if len(clean_seq) < 30: # Arbitrary short length threshold (e.g., < 10 codons)
        issues.append(f"Sequence is very short (length: {len(clean_seq)})")

    # GC content (only if sequence is not empty)
    gc_count = clean_seq.count("G") + clean_seq.count("C")
    gc_content = gc_count / len(clean_seq) if len(clean_seq) > 0 else 0
    if len(clean_seq) > 0 and (gc_content < 0.25 or gc_content > 0.75):
        issues.append(f"Atypical GC content: {gc_content:.2f}")

    # N ratio (unknown bases)
    n_count = clean_seq.count("N")
    n_ratio = n_count / len(clean_seq) if len(clean_seq) > 0 else 0
    if n_ratio > 0.1: # More than 10% Ns
        issues.append(f"High proportion of unknown bases (N): {n_ratio:.2f}")

    # Simplified repeat check (from legacy, can be expanded)
    # Only check if sequence is long enough
    if len(clean_seq) > 100:
        # Check for non-triplet repeats (e.g., 7-mers, 10-mers) that might indicate frameshifts
        for repeat_len in [7, 10]:
            repeat_counts = defaultdict(int)
            for i in range(len(clean_seq) - repeat_len + 1):
                fragment = clean_seq[i:i+repeat_len]
                if "N" not in fragment: # Ignore fragments with Ns
                    repeat_counts[fragment] += 1

            # Count fragments that repeat more than once
            num_repeating_fragments = sum(1 for count in repeat_counts.values() if count > 1)
            if num_repeating_fragments >= min_repeat_report: # Use the passed or default threshold
                issues.append(f"High number of {repeat_len}-mer repeats ({num_repeating_fragments}), suggesting potential quality issues or frameshifts.")

    if issues:
        logger.debug(f"Quality issues for {seq_id}: {issues}")
    else:
        logger.debug(f"No major quality issues detected for {seq_id}")
    return issues

def handle_non_triplet_sequence(sequence, seq_id, position="end", context=None):
    """Handles sequences not of length multiple of 3 by padding with Ns."""
    logger.debug(f"Handling non-triplet sequence: {seq_id}, length: {len(sequence)}, position: {position}")
    remainder = len(sequence) % 3
    if remainder == 0:
        return sequence

    padding_needed = 3 - remainder
    padding = "N" * padding_needed

    # Smart padding based on context (start/stop codons)
    if position == "smart" and context:
        has_start = context.get("has_start", False) or sequence[:3].upper() == "ATG"
        has_stop = context.get("has_stop", False) or sequence[-3:].upper() in ["TAA", "TAG", "TGA"]

        if has_start and not has_stop:
            logger.info(f"Smart padding {seq_id} at 3-prime end (has start, no stop). Adding {padding_needed} Ns.")
            return sequence + padding
        elif not has_start and has_stop:
            logger.info(f"Smart padding {seq_id} at 5-prime end (no start, has stop). Adding {padding_needed} Ns.")
            return padding + sequence
        else:
            logger.info(f"Smart padding {seq_id} at 3-prime end (default for ambiguous context). Adding {padding_needed} Ns.")
            return sequence + padding # Default to end if context is ambiguous
    elif position == "start":
        logger.info(f"Padding {seq_id} at 5-prime end. Adding {padding_needed} Ns.")
        return padding + sequence
    else: # Default to "end"
        logger.info(f"Padding {seq_id} at 3-prime end. Adding {padding_needed} Ns.")
        return sequence + padding
