import unittest
import os
import sys

# Add project root to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from cds_processor.cds_processor import duplicate_handler
from Bio.Seq import Seq

class TestDuplicateHandler(unittest.TestCase):

    def setUp(self):
        self.orf_results_simple_duplicates = {
            "seq1": [Seq("ATGCATGC"), Seq("GTGTGT"), Seq("ATGCATGC")], # ORF 0 and 2 are duplicates
            "seq2": [Seq("GTGTGT"), Seq("AGAGAG")]
        }
        # seq1_orf_0: ATGCATGC
        # seq1_orf_1: GTGTGT
        # seq1_orf_2: ATGCATGC (duplicate of seq1_orf_0)
        # seq2_orf_0: GTGTGT (duplicate of seq1_orf_1)
        # seq2_orf_1: AGAGAG

        self.expected_unique_orf_results = {
            "seq1": [Seq("ATGCATGC"), Seq("GTGTGT")],
            "seq2": [Seq("AGAGAG")]
        }
        self.expected_duplicate_info = {
            "seq1_orf_2": "seq1_orf_0", # ATGCATGC is dup of first ATGCATGC
            "seq2_orf_0": "seq1_orf_1"  # GTGTGT is dup of first GTGTGT (which was in seq1)
        }

    def test_find_duplicates_by_sequence(self):
        unique_orfs, dup_info = duplicate_handler.find_duplicates_by_sequence(self.orf_results_simple_duplicates)

        # Check unique ORFs
        self.assertIn("seq1", unique_orfs)
        self.assertEqual(len(unique_orfs["seq1"]), 2) # Should have 2 unique ORFs
        self.assertIn(Seq("ATGCATGC"), unique_orfs["seq1"])
        self.assertIn(Seq("GTGTGT"), unique_orfs["seq1"])

        self.assertIn("seq2", unique_orfs)
        self.assertEqual(len(unique_orfs["seq2"]), 1) # Should have 1 unique ORF
        self.assertIn(Seq("AGAGAG"), unique_orfs["seq2"])

        # Check duplicate info
        self.assertEqual(len(dup_info), 2)
        self.assertEqual(dup_info["seq1_orf_2"], "seq1_orf_0")
        self.assertEqual(dup_info["seq2_orf_0"], "seq1_orf_1")


    def test_find_duplicates_by_sequence_no_duplicates(self):
        no_dup_orf_results = {
            "seqA": [Seq("ATGCG"), Seq("CGTAC")],
            "seqB": [Seq("TTTTT")]
        }
        unique_orfs, dup_info = duplicate_handler.find_duplicates_by_sequence(no_dup_orf_results)
        self.assertEqual(len(unique_orfs["seqA"]), 2)
        self.assertEqual(len(unique_orfs["seqB"]), 1)
        self.assertEqual(len(dup_info), 0) # No duplicates

    def test_find_duplicates_by_sequence_all_duplicates_of_one(self):
        all_dup_orf_results = {
            "s1": [Seq("AAAAA"), Seq("AAAAA")],
            "s2": [Seq("AAAAA")]
        }
        # Expected: s1 has [Seq("AAAAA")], s2 is empty or not present if all its ORFs were duplicates
        # The current implementation keeps the entry for s2 if it initially had ORFs, even if all become dups.
        # Let's check the output: s1: [AAAAA], s2: [] (if it had ORFs initially)
        # The current code structure: if unique_orfs_for_seq: unique_orfs_dict[original_seq_id] = unique_orfs_for_seq
        # So if all ORFs in seq2 are duplicates of ORFs in seq1, seq2 might be missing from unique_orfs.

        unique_orfs, dup_info = duplicate_handler.find_duplicates_by_sequence(all_dup_orf_results)

        self.assertIn("s1", unique_orfs)
        self.assertEqual(len(unique_orfs["s1"]), 1)
        self.assertEqual(str(unique_orfs["s1"][0]), "AAAAA")

        self.assertNotIn("s2", unique_orfs) # Because all its ORFs were duplicates of s1's ORF

        self.assertEqual(dup_info["s1_orf_1"], "s1_orf_0")
        self.assertEqual(dup_info["s2_orf_0"], "s1_orf_0")


    def test_find_duplicates_by_blast_placeholder(self):
        # Placeholder test
        blast_hits_placeholder = {"orf1": [], "orf2": []}
        result = duplicate_handler.find_duplicates_by_blast(blast_hits_placeholder)
        self.assertTrue(isinstance(result, dict)) # Placeholder returns a dict

    def test_filter_duplicates_strategy_sequence(self):
        filtered_orfs = duplicate_handler.filter_duplicates(self.orf_results_simple_duplicates, strategy="sequence")
        # This should be the same as unique_orfs from the specific test.
        self.assertEqual(len(filtered_orfs["seq1"]), 2)
        self.assertEqual(len(filtered_orfs["seq2"]), 1)

    def test_filter_duplicates_strategy_blast_placeholder(self):
        # Placeholder test, currently returns original if strategy is 'blast'
        filtered_orfs = duplicate_handler.filter_duplicates(self.orf_results_simple_duplicates, strategy="blast", blast_results={})
        self.assertEqual(filtered_orfs, self.orf_results_simple_duplicates) # As it's a placeholder

    def test_filter_duplicates_unknown_strategy(self):
        # Should return original and print a message
        filtered_orfs = duplicate_handler.filter_duplicates(self.orf_results_simple_duplicates, strategy="unknown")
        self.assertEqual(filtered_orfs, self.orf_results_simple_duplicates)

if __name__ == "__main__":
    unittest.main()
