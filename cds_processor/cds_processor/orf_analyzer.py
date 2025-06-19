import logging
# from .utils import GENETIC_CODE # Will be needed
from functools import lru_cache

logger = logging.getLogger(__name__) # cds_processor.orf_analyzer

# The lru_cache for safe_translate_codon was in the legacy script global scope.
# If safe_translate_codon is defined here, the cache can be applied directly.
@lru_cache(maxsize=1024)
def safe_translate_codon(codon, genetic_code_map=None):
    """Translates a codon to an amino acid, handling ambiguities with X."""
    # GENETIC_CODE will be imported from utils or passed as arg
    if genetic_code_map is None:
        # This is a placeholder, proper GENETIC_CODE should be used
        # For now, provide a minimal map to avoid NameError if not passed
        genetic_code_map = {
            "ATG": "M", "TAA": "*", "TAG": "*", "TGA": "*", # Minimal for basic ORF finding
            # Add more if needed for full translation, or ensure GENETIC_CODE from utils is used
        }

    codon_upper = codon.upper().strip()
    if codon_upper in genetic_code_map:
        return genetic_code_map[codon_upper]
    if "N" in codon_upper or any(base not in "ACGT" for base in codon_upper):
        return "X" # Standard ambiguity character
    logger.warning(f"Unknown codon encountered: {codon_upper}. Returning X.")
    return "X"

def analyze_codon_structure(sequence, genetic_code_map=None):
    """Analyzes codon structure, identifies ORFs, and determines the best reading frame."""
    # GENETIC_CODE will be imported from utils or passed as arg
    logger.debug(f"Analyzing codon structure for sequence of length {len(sequence)}")
    clean_seq = sequence.upper().replace("-", "").replace("N", "") # Remove gaps and Ns for ORF analysis
    if not clean_seq or len(clean_seq) < 3:
        logger.warning("Sequence too short or empty for ORF analysis after cleaning.")
        return {"best_frame": 0, "longest_orf": None, "has_start": False, "has_stop": False, "orf_length": 0, "start_positions": []}

    best_frame = 0
    longest_orf_info = {"codons": [], "start_codon_pos_in_orf": -1, "stop_codon_pos_in_orf": -1}
    all_start_positions = [] # Store all ATG positions found across frames

    for frame_offset in range(3):
        codons_in_frame = []
        current_orf_codons = []
        in_orf = False
        orf_start_pos_global = -1 # Start position of current ORF in the original cleaned sequence

        for i in range(frame_offset, len(clean_seq) - 2, 3):
            codon = clean_seq[i:i+3]
            codons_in_frame.append(codon)
            aa = safe_translate_codon(codon, genetic_code_map) # Use the cached version

            if codon == "ATG": # Potential start codon
                all_start_positions.append({"frame": frame_offset, "position": i, "codon_idx_in_frame": len(codons_in_frame)-1})
                if not in_orf:
                    in_orf = True
                    current_orf_codons = [codon]
                    orf_start_pos_global = i
                else: # Already in ORF, ATG is just another codon
                    current_orf_codons.append(codon)
            elif aa == "*" and in_orf: # Stop codon and currently in an ORF
                # ORF ended here (exclusive of stop codon for length, but include for context if needed)
                if len(current_orf_codons) > len(longest_orf_info["codons"]):
                    longest_orf_info = {
                        "codons": list(current_orf_codons),
                        "start_codon_pos_in_orf": 0, # Assuming ORF started with the first codon captured
                        "stop_codon_pos_in_orf": len(current_orf_codons) * 3, # Position after last codon
                        "frame": frame_offset,
                        "start_pos_global": orf_start_pos_global
                    }
                    best_frame = frame_offset
                in_orf = False
                current_orf_codons = []
                orf_start_pos_global = -1
            elif in_orf: # Codon is part of the current ORF
                current_orf_codons.append(codon)

        # After iterating through codons in a frame, if still in_orf (no stop codon found)
        if in_orf and len(current_orf_codons) > len(longest_orf_info["codons"]):
            longest_orf_info = {
                "codons": list(current_orf_codons),
                "start_codon_pos_in_orf": 0,
                "stop_codon_pos_in_orf": -1, # No stop codon found
                "frame": frame_offset,
                "start_pos_global": orf_start_pos_global
            }
            best_frame = frame_offset

    has_start_codon = False
    if longest_orf_info["codons"] and longest_orf_info["codons"][0] == "ATG":
        has_start_codon = True

    # Check if the longest ORF ends with a stop codon (implicitly, if stop_codon_pos_in_orf is not -1)
    # The original script had has_stop as False initially, this is a bit more direct.
    has_stop_codon_defined = (longest_orf_info.get("stop_codon_pos_in_orf", -1) != -1 and longest_orf_info["codons"])

    result = {
        "best_frame": best_frame,
        "longest_orf": longest_orf_info["codons"], # Just the list of codons
        "orf_start_pos_global": longest_orf_info.get("start_pos_global", -1),
        "orf_length": len(longest_orf_info["codons"]) * 3,
        "has_start": has_start_codon, # Whether the selected longest ORF starts with ATG
        "has_stop": has_stop_codon_defined, # Whether the selected longest ORF had a defined stop codon
        "start_positions": all_start_positions # All ATGs found
    }
    logger.debug(f"Codon structure analysis result: {result}")
    return result
