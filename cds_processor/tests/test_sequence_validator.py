import unittest
import os
import sys

# Add project root to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from cds_processor.cds_processor import sequence_validator # CORRECTED IMPORT
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class TestSequenceValidator(unittest.TestCase):

    def test_is_valid_dna_sequence_valid(self):
        # Valid sequence (length multiple of 3, only ATCGN)
        seq_rec = SeqRecord(Seq("ATGCGT"), id="valid_seq") # CORRECTED LENGTH (6)
        self.assertTrue(sequence_validator.is_valid_dna_sequence(seq_rec))

    def test_is_valid_dna_sequence_invalid_chars(self):
        seq_rec = SeqRecord(Seq("ATGCGXN"), id="invalid_chars") # X is not valid
        self.assertFalse(sequence_validator.is_valid_dna_sequence(seq_rec))

    def test_is_valid_dna_sequence_invalid_length(self):
        seq_rec = SeqRecord(Seq("ATGCG"), id="invalid_length") # Length 5, not multiple of 3
        self.assertFalse(sequence_validator.is_valid_dna_sequence(seq_rec))

    def test_is_valid_dna_sequence_lowercase(self):
        # Should convert to uppercase and validate
        seq_rec = SeqRecord(Seq("atgcgt"), id="lowercase_seq") # CORRECTED LENGTH (6)
        self.assertTrue(sequence_validator.is_valid_dna_sequence(seq_rec))

    def test_is_valid_dna_sequence_empty(self):
        # Empty sequence, length 0 (multiple of 3) but might be considered invalid by some criteria
        # Current implementation allows it if it passes regex (it does, as ^[ATCGN]+$ matches empty if not careful)
        # Let's refine the regex in `is_valid_dna_sequence` to `^[ATCGN]+$` which means at least one char.
        # An empty sequence will fail `^[ATCGN]+$`.
        # If an empty string should be valid, regex should be `^[ATCGN]*$`.
        # Given these are CDS, empty is probably not valid.
        seq_rec = SeqRecord(Seq(""), id="empty_seq")
        self.assertFalse(sequence_validator.is_valid_dna_sequence(seq_rec)) # Fails due to non-matching regex
                                                                          # and length % 3 == 0 is true.

    def test_validate_sequences_mixed(self):
        sequences = {
            "s1": SeqRecord(Seq("ATGCGTAGC"), id="s1"),      # Valid
            "s2": SeqRecord(Seq("ATGCGTX"), id="s2"),       # Invalid char
            "s3": SeqRecord(Seq("ATGCGTA"), id="s3"),       # Invalid length (7)
            "s4": SeqRecord(Seq("NNNCGNNNN"), id="s4"),    # Valid
            "s5": SeqRecord(Seq(""), id="s5")               # Invalid (empty)
        }
        valid_sequences = sequence_validator.validate_sequences(sequences)
        self.assertIn("s1", valid_sequences)
        self.assertNotIn("s2", valid_sequences)
        self.assertNotIn("s3", valid_sequences)
        self.assertIn("s4", valid_sequences)
        self.assertNotIn("s5", valid_sequences)
        self.assertEqual(len(valid_sequences), 2)

    def test_validate_sequences_all_valid(self):
        sequences = {
            "s1": SeqRecord(Seq("ATGCGTAGC"), id="s1"),
            "s2": SeqRecord(Seq("NNNCGNNNN"), id="s2")
        }
        valid_sequences = sequence_validator.validate_sequences(sequences)
        self.assertEqual(len(valid_sequences), 2)

    def test_validate_sequences_all_invalid(self):
        sequences = {
            "s1": SeqRecord(Seq("ATGCGX"), id="s1"),
            "s2": SeqRecord(Seq("ATGCG"), id="s2")
        }
        valid_sequences = sequence_validator.validate_sequences(sequences)
        self.assertEqual(len(valid_sequences), 0)

    def test_validate_sequences_empty_input(self):
        sequences = {}
        valid_sequences = sequence_validator.validate_sequences(sequences)
        self.assertEqual(len(valid_sequences), 0)

if __name__ == "__main__":
    unittest.main()
