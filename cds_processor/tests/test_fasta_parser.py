import unittest
import os
import sys

# Add project root to sys.path to allow importing cds_processor modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from cds_processor.cds_processor import fasta_parser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class TestFastaParser(unittest.TestCase):

    def setUp(self):
        # Create a dummy FASTA file for testing
        self.test_fasta_content = ">seq1 description1\nATGCGT\n>seq2 description2\nGTACNA\n"
        self.test_fasta_filepath = "test_input.fasta"
        with open(self.test_fasta_filepath, "w") as f:
            f.write(self.test_fasta_content)

        self.empty_fasta_filepath = "empty.fasta"
        with open(self.empty_fasta_filepath, "w") as f:
            f.write("")

        self.malformed_fasta_filepath = "malformed.fasta"
        with open(self.malformed_fasta_filepath, "w") as f:
            f.write(">seq1\nATGC\n>seq2\nGTAC\nSEQ3WITHOUTHEADER\n")


    def tearDown(self):
        # Clean up the dummy FASTA file
        if os.path.exists(self.test_fasta_filepath):
            os.remove(self.test_fasta_filepath)
        if os.path.exists(self.empty_fasta_filepath):
            os.remove(self.empty_fasta_filepath)
        if os.path.exists(self.malformed_fasta_filepath):
            os.remove(self.malformed_fasta_filepath)

    def test_parse_valid_fasta_file(self):
        sequences = fasta_parser.parse_fasta_file(self.test_fasta_filepath)
        self.assertEqual(len(sequences), 2)
        self.assertIn("seq1", sequences)
        self.assertIn("seq2", sequences)
        self.assertEqual(str(sequences["seq1"].seq), "ATGCGT")
        self.assertEqual(sequences["seq1"].description, "seq1 description1") # SeqIO takes the whole header line
        self.assertEqual(str(sequences["seq2"].seq), "GTACNA")
        self.assertEqual(sequences["seq2"].id, "seq2")

    def test_parse_non_existent_file(self):
        sequences = fasta_parser.parse_fasta_file("non_existent.fasta")
        self.assertEqual(len(sequences), 0) # Expect an empty dict and an error message (printed by func)

    def test_parse_empty_fasta_file(self):
        sequences = fasta_parser.parse_fasta_file(self.empty_fasta_filepath)
        self.assertEqual(len(sequences), 0)

    def test_parse_malformed_fasta_file(self):
        # SeqIO can be quite robust, let's see how it handles this.
        # It might parse the valid parts or raise an error depending on strictness (not set here).
        # The current fasta_parser catches general exceptions.
        sequences = fasta_parser.parse_fasta_file(self.malformed_fasta_filepath)
        # Depending on Biopython's SeqIO behavior, it might parse the first two.
        # If it errors out, the function returns {}. Let's assume it parses what it can.
        # Based on typical SeqIO behavior, it should parse seq1 and seq2.
        self.assertIn("seq1", sequences)
        self.assertIn("seq2", sequences)
        self.assertNotIn("SEQ3WITHOUTHEADER", sequences) # This part is malformed.
        # If strict parsing is desired, the underlying SeqIO.parse call would need adjustment,
        # or more specific error handling. The current wrapper prints an error and returns {}.

if __name__ == "__main__":
    unittest.main()
