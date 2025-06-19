import unittest
import os
import sys

# Add project root to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from cds_processor.cds_processor import orf_analyzer
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class TestOrfAnalyzer(unittest.TestCase):

    def test_find_orfs_simple_case_one_orf(self):
        # ATG...(stop)
        seq_rec = SeqRecord(Seq("CGATGCGGCCGTGACTACG"), id="simple_orf") # ORF: ATGCGGCCGTGA
        orfs = orf_analyzer.find_orfs(seq_rec)
        self.assertEqual(len(orfs), 1)
        self.assertEqual(str(orfs[0]), "ATGCGGCCGTGA")

    def test_find_orfs_multiple_orfs_same_frame(self):
        # ATG...(stop)...ATG...(stop)
        seq_rec = SeqRecord(Seq("ATGCGGTAGCCATGCCCTAA"), id="multi_orf_same_frame")
        # ORF1: ATGCGGTAG
        # ORF2: ATGCCCTAA
        orfs = orf_analyzer.find_orfs(seq_rec)
        self.assertEqual(len(orfs), 2)
        self.assertIn(Seq("ATGCGGTAG"), orfs)
        self.assertIn(Seq("ATGCCCTAA"), orfs)

    def test_find_orfs_different_frames(self):
        # Frame 0: ATG...TAG
        # Frame 1: .CATG...TAA
        seq_str = "ATGCCCTAG" + "C" + "CATGCGGTAA" # ORF1: ATGCCCTAG, ORF2: ATGCGGTAA (from 2nd ATG)
        seq_rec = SeqRecord(Seq(seq_str), id="multi_frame_orf")
        orfs = orf_analyzer.find_orfs(seq_rec)

        expected_orf1 = Seq("ATGCCCTAG")
        expected_orf2 = Seq("ATGCGGTAA") # This is from the C ATGCGGTAA part, frame 1 of original

        self.assertEqual(len(orfs), 2)
        self.assertIn(expected_orf1, orfs)
        self.assertIn(expected_orf2, orfs)


    def test_find_orfs_no_start_codon(self):
        seq_rec = SeqRecord(Seq("CGTACGTACGTAGCTA"), id="no_start")
        orfs = orf_analyzer.find_orfs(seq_rec)
        self.assertEqual(len(orfs), 0)

    def test_find_orfs_no_stop_codon(self):
        # ATG followed by sequence without stop codon
        seq_rec = SeqRecord(Seq("ATGCGCCGCGCG"), id="no_stop")
        orfs = orf_analyzer.find_orfs(seq_rec)
        self.assertEqual(len(orfs), 0) # Current implementation requires a stop codon

    def test_find_orfs_internal_atg(self):
        # ATG...ATG...(stop) - should find the outermost ORF.
        # The simplified model might also find internal ORFs if they are in frame.
        # Current model: finds first ATG, then looks for stop. Then continues search after that ATG.
        # So, ATG (ATG TAG) -> ATATGTAG is found.
        # And after that, the second ATG is found, resulting in ATGTAG.
        seq_rec = SeqRecord(Seq("CGATGATGTAGCC"), id="internal_atg") # ORF1: ATGATGTAG, ORF2: ATGTAG
        orfs = orf_analyzer.find_orfs(seq_rec)
        self.assertEqual(len(orfs), 2)
        self.assertIn(Seq("ATGATGTAG"), orfs)
        self.assertIn(Seq("ATGTAG"), orfs)


    def test_find_orfs_overlapping_orfs_different_frames(self):
        # Example: AGATGCxxxxxxTAGxxxx (Frame 1: ATGCxxxxxxTAG)
        #          A GATGCGTGTAG (Frame 1)
        #           G ATGCGTGTA G (Frame 2, no start)
        #            A TGCGTGTAG (Frame 0 of subseq, no start)
        # Complex cases like this test frame handling thoroughly.
        # seq_str = "AGATGCGTTAG" # ORF in frame 1: ATGCGTTAG (length 9)
        # Frame 0: AGA TGC GTT AG (no orf)
        # Frame 1: G ATG CGT TAG (ATG CGT TAG)
        # Frame 2: A TGC GTT AG (no orf)
        seq_rec = SeqRecord(Seq("AGATGCGTTAG"), id="overlap_orf")
        orfs = orf_analyzer.find_orfs(seq_rec)
        self.assertEqual(len(orfs), 1)
        self.assertEqual(str(orfs[0]), "ATGCGTTAG")

    def test_find_orfs_minimal_orf(self):
        # ATGTAG (start, stop, 6nt)
        seq_rec = SeqRecord(Seq("ATGTAG"), id="min_orf")
        orfs = orf_analyzer.find_orfs(seq_rec)
        self.assertEqual(len(orfs), 1)
        self.assertEqual(str(orfs[0]), "ATGTAG")

    def test_analyze_orfs_multiple_sequences(self):
        sequences = {
            "s1": SeqRecord(Seq("ATGCGGCCGTGA"), id="s1"), # ORF: ATGCGGCCGTGA
            "s2": SeqRecord(Seq("NOATGCRAPTAA"), id="s2"), # No ORF
            "s3": SeqRecord(Seq("ATGTAACATGTGA"), id="s3") # ORF1: ATGTAA, ORF2: ATGTGA
        }
        orf_results = orf_analyzer.analyze_orfs(sequences)
        self.assertIn("s1", orf_results)
        self.assertEqual(len(orf_results["s1"]), 1)
        self.assertEqual(str(orf_results["s1"][0]), "ATGCGGCCGTGA")

        self.assertNotIn("s2", orf_results) # No ORFs, so key s2 should not be in result

        self.assertIn("s3", orf_results)
        self.assertEqual(len(orf_results["s3"]), 2)
        self.assertIn(Seq("ATGTAA"), orf_results["s3"])
        self.assertIn(Seq("ATGTGA"), orf_results["s3"])

    def test_analyze_orfs_no_sequences(self):
        sequences = {}
        orf_results = orf_analyzer.analyze_orfs(sequences)
        self.assertEqual(len(orf_results), 0)

    def test_analyze_orfs_no_orfs_found(self):
        sequences = {"s1": SeqRecord(Seq("ACGTACGT"), id="s1")}
        orf_results = orf_analyzer.analyze_orfs(sequences)
        self.assertEqual(len(orf_results), 0)


if __name__ == "__main__":
    unittest.main()
