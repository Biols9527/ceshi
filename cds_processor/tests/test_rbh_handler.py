import unittest
import os
import tempfile

# Adjust path to import from cds_processor package
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR) # This should be "cds_processor"
PACKAGE_ROOT = os.path.join(PROJECT_ROOT, "cds_processor")

if PACKAGE_ROOT not in sys.path:
    sys.path.insert(0, PACKAGE_ROOT)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from cds_processor import rbh_handler

class TestRBHHit(unittest.TestCase):
    def test_rbh_hit_creation(self):
        hit = rbh_handler.RBHHit(
            query_id="q1", subject_id="s1", pident=90.0, length=100,
            evalue=1e-10, bitscore=200.0, query_len=110, subject_len=120,
            query_start=1, query_end=100, subject_start=10, subject_end=109
        )
        self.assertEqual(hit.query_id, "q1")
        self.assertEqual(hit.subject_id, "s1")
        self.assertEqual(hit.pident, 90.0)
        self.assertEqual(hit.length, 100)
        self.assertEqual(hit.evalue, 1e-10)
        self.assertEqual(hit.bitscore, 200.0)
        self.assertEqual(hit.query_len, 110)
        self.assertEqual(hit.subject_len, 120)
        self.assertEqual(hit.query_start, 1)
        self.assertEqual(hit.query_end, 100)
        self.assertEqual(hit.subject_start, 10)
        self.assertEqual(hit.subject_end, 109)

class TestRBHHitGroup(unittest.TestCase):
    def setUp(self):
        self.group = rbh_handler.RBHHitGroup("queryA")
        # score (desc), evalue (asc)
        self.hit1 = rbh_handler.RBHHit("queryA", "s1", 90, 100, 1e-20, 300)
        self.hit2 = rbh_handler.RBHHit("queryA", "s2", 95, 110, 1e-30, 400) # Best
        self.hit3 = rbh_handler.RBHHit("queryA", "s3", 80, 90, 1e-10, 200)
        self.hit4_tie_score_better_evalue = rbh_handler.RBHHit("queryA", "s4", 98, 120, 1e-40, 400) # Tied best with hit2
        self.hit5_tie_score_worse_evalue = rbh_handler.RBHHit("queryA", "s5", 92, 115, 1e-25, 400) # Tied score, but worse e-val than hit4_tie_score_better_evalue

    def test_add_hit(self):
        self.group.add_hit(self.hit1)
        self.assertEqual(len(self.group.hits), 1)
        with self.assertRaises(ValueError):
            self.group.add_hit(rbh_handler.RBHHit("wrongQuery", "s", 0,0,0,0))

    def test_best_hits_empty(self):
        self.assertEqual(self.group.best_hits, [])

    def test_best_hits_single_best(self):
        self.group.add_hit(self.hit1)
        self.group.add_hit(self.hit2) # hit2 is clearly best
        self.group.add_hit(self.hit3)
        best = self.group.best_hits
        self.assertEqual(len(best), 1)
        self.assertIs(best[0], self.hit2)

    def test_best_hits_tie_in_score_and_evalue(self):
        # hit2 and hit4_tie_score_better_evalue have same score (400)
        # Let's make hit4_tie_score_better_evalue have same evalue as hit2 for a true tie
        self.hit4_tie_score_better_evalue.evalue = self.hit2.evalue # 1e-30

        self.group.add_hit(self.hit1)
        self.group.add_hit(self.hit2)
        self.group.add_hit(self.hit3)
        self.group.add_hit(self.hit4_tie_score_better_evalue)

        best = self.group.best_hits
        self.assertEqual(len(best), 2)
        self.assertIn(self.hit2, best)
        self.assertIn(self.hit4_tie_score_better_evalue, best)

    def test_best_hits_tie_in_score_different_evalue(self):
        # hit2 (400, 1e-30), hit4 (400, 1e-40), hit5 (400, 1e-25)
        # hit4 should be the sole best because of better e-value
        self.group.add_hit(self.hit2)
        self.group.add_hit(self.hit4_tie_score_better_evalue) # Best: score 400, evalue 1e-40
        self.group.add_hit(self.hit5_tie_score_worse_evalue) # Score 400, evalue 1e-25 (worse than hit4)

        best = self.group.best_hits
        self.assertEqual(len(best), 1)
        self.assertIs(best[0], self.hit4_tie_score_better_evalue)

    def test_sorting_logic(self):
        # hit2 (400, 1e-30)
        # hit1 (300, 1e-20)
        # hit3 (200, 1e-10)
        # hit4 (400, 1e-40) -> best overall
        # hit5 (400, 1e-25) -> second best group (tied score with hit2, but worse evalue)

        self.group.add_hit(self.hit1) # 300
        self.group.add_hit(self.hit2) # 400, 1e-30
        self.group.add_hit(self.hit3) # 200
        self.group.add_hit(self.hit4_tie_score_better_evalue) # 400, 1e-40 (BEST)
        self.group.add_hit(self.hit5_tie_score_worse_evalue) # 400, 1e-25

        # Access best_hits to trigger sorting
        _ = self.group.best_hits

        # Check internal order of hits list after sorting
        # Expected order: hit4, hit2, hit5, hit1, hit3 (based on __lt__ in RBHHit)
        self.assertEqual(self.group.hits[0], self.hit4_tie_score_better_evalue) # Best score, best e-value
        self.assertEqual(self.group.hits[1], self.hit2) # Same score as hit4, but worse e-value
        self.assertEqual(self.group.hits[2], self.hit5_tie_score_worse_evalue) # Same score as hit2, but worse e-value
        self.assertEqual(self.group.hits[3], self.hit1) # Lower score
        self.assertEqual(self.group.hits[4], self.hit3) # Lowest score

class TestGetSequenceLengthsFromFasta(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.fasta_file_valid = os.path.join(self.temp_dir, "valid.fasta")
        with open(self.fasta_file_valid, "w") as f:
            f.write(">seq1 description1
AAACCCGGG
") # len 9
            f.write(">seq2
TTTT
") # len 4
            f.write(">seq3nospace
AAAAACCCCC
") # len 10

        self.fasta_file_empty = os.path.join(self.temp_dir, "empty.fasta")
        with open(self.fasta_file_empty, "w") as f:
            pass # Empty file

        self.fasta_file_malformed = os.path.join(self.temp_dir, "malformed.fasta")
        with open(self.fasta_file_malformed, "w") as f:
            f.write("This is not fasta
")


    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_valid_fasta(self):
        lengths = rbh_handler.get_sequence_lengths_from_fasta(self.fasta_file_valid)
        self.assertEqual(len(lengths), 3)
        self.assertEqual(lengths["seq1"], 9)
        self.assertEqual(lengths["seq2"], 4)
        self.assertEqual(lengths["seq3nospace"], 10)

    def test_empty_fasta(self):
        lengths = rbh_handler.get_sequence_lengths_from_fasta(self.fasta_file_empty)
        self.assertEqual(len(lengths), 0)

    def test_non_existent_fasta(self):
        lengths = rbh_handler.get_sequence_lengths_from_fasta("non_existent.fasta")
        self.assertEqual(len(lengths), 0) # Should log error and return empty

    def test_malformed_fasta(self):
        # Biopython's SeqIO.parse might raise an error or return nothing depending on malformation
        # Expecting it to handle gracefully by returning no sequences or raising an error caught by the function
        lengths = rbh_handler.get_sequence_lengths_from_fasta(self.fasta_file_malformed)
        self.assertEqual(len(lengths), 0) # Should log error and return empty

# (Keep existing TestRBHHit, TestRBHHitGroup, TestGetSequenceLengthsFromFasta)

# --- Sample BLAST XML Data for testing parse_blast_xml ---
SAMPLE_BLAST_XML_MULTIPLE_HITS = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.10.1+</BlastOutput_version>
  <BlastOutput_reference>Stephen F. Altschul, et al.</BlastOutput_reference>
  <BlastOutput_db>test_db</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>queryA_seq1</BlastOutput_query-def>
  <BlastOutput_query-len>100</BlastOutput_query-len>
  <BlastOutput_param><Parameters> <Parameters_expect>0.01</Parameters_expect> </Parameters></BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>queryA_seq1</Iteration_query-ID>
      <Iteration_query-def>queryA_seq1</Iteration_query-def>
      <Iteration_query-len>100</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>subjectB_seq1</Hit_id>
          <Hit_def>subjectB_seq1</Hit_def>
          <Hit_accession>subjectB_seq1</Hit_accession>
          <Hit_len>120</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>250.0</Hsp_bit-score>
              <Hsp_score>130</Hsp_score>
              <Hsp_evalue>1e-25</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>90</Hsp_query-to>
              <Hsp_hit-from>10</Hsp_hit-from>
              <Hsp_hit-to>99</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>85</Hsp_identity>
              <Hsp_positive>85</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>90</Hsp_align-len>
              <Hsp_qseq>...</Hsp_qseq>
              <Hsp_hseq>...</Hsp_hseq>
              <Hsp_midline>...</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>subjectB_seq2</Hit_id>
          <Hit_def>subjectB_seq2</Hit_def>
          <Hit_accession>subjectB_seq2</Hit_accession>
          <Hit_len>110</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>150.0</Hsp_bit-score>
              <Hsp_score>70</Hsp_score>
              <Hsp_evalue>1e-15</Hsp_evalue>
              <Hsp_query-from>5</Hsp_query-from>
              <Hsp_query-to>80</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>76</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>70</Hsp_identity>
              <Hsp_positive>70</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>76</Hsp_align-len>
              <Hsp_qseq>...</Hsp_qseq>
              <Hsp_hseq>...</Hsp_hseq>
              <Hsp_midline>...</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics> <Statistics_db-num>10</Statistics_db-num> <Statistics_db-len>5000</Statistics_db-len> <Statistics_hsp-len>20</Statistics_hsp-len> <Statistics_eff-space>100000</Statistics_eff-space> <Statistics_kappa>0.41</Statistics_kappa> <Statistics_lambda>0.625</Statistics_lambda> <Statistics_entropy>0.78</Statistics_entropy> </Statistics>
      </Iteration_stat>
    </Iteration>
    <Iteration>
      <Iteration_iter-num>2</Iteration_iter-num>
      <Iteration_query-ID>queryA_seq2</Iteration_query-ID>
      <Iteration_query-def>queryA_seq2</Iteration_query-def>
      <Iteration_query-len>80</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>subjectB_seq1</Hit_id>
          <Hit_def>subjectB_seq1</Hit_def>
          <Hit_accession>subjectB_seq1</Hit_accession>
          <Hit_len>120</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>180.0</Hsp_bit-score>
              <Hsp_score>90</Hsp_score>
              <Hsp_evalue>1e-18</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>70</Hsp_query-to>
              <Hsp_hit-from>20</Hsp_hit-from>
              <Hsp_hit-to>89</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>65</Hsp_identity>
              <Hsp_positive>65</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>70</Hsp_align-len>
              <Hsp_qseq>...</Hsp_qseq>
              <Hsp_hseq>...</Hsp_hseq>
              <Hsp_midline>...</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat><Statistics/></Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""

SAMPLE_BLAST_XML_NO_HITS = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>query_no_hit</BlastOutput_query-def>
  <BlastOutput_query-len>50</BlastOutput_query-len>
  <BlastOutput_param><Parameters><Parameters_expect>10</Parameters_expect></Parameters></BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>query_no_hit</Iteration_query-ID>
      <Iteration_query-def>query_no_hit</Iteration_query-def>
      <Iteration_query-len>50</Iteration_query-len>
      <Iteration_hits></Iteration_hits>
      <Iteration_stat><Statistics/></Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""

class TestParseBlastXML(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.xml_file_multi = os.path.join(self.temp_dir, "multi.xml")
        with open(self.xml_file_multi, "w") as f:
            f.write(SAMPLE_BLAST_XML_MULTIPLE_HITS)

        self.xml_file_none = os.path.join(self.temp_dir, "none.xml")
        with open(self.xml_file_none, "w") as f:
            f.write(SAMPLE_BLAST_XML_NO_HITS)

        self.xml_file_empty_content = os.path.join(self.temp_dir, "empty_content.xml")
        with open(self.xml_file_empty_content, "w") as f:
            f.write("") # Empty content

        self.query_lengths = {"queryA_seq1": 100, "queryA_seq2": 80, "query_no_hit": 50}
        self.subject_lengths = {"subjectB_seq1": 120, "subjectB_seq2": 110}

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_parse_valid_xml_multiple_hits(self):
        hit_groups = rbh_handler.parse_blast_xml(self.xml_file_multi, self.query_lengths, self.subject_lengths)
        self.assertEqual(len(hit_groups), 2) # queryA_seq1, queryA_seq2

        # Check queryA_seq1
        self.assertIn("queryA_seq1", hit_groups)
        group1 = hit_groups["queryA_seq1"]
        self.assertEqual(group1.query_id, "queryA_seq1")
        self.assertEqual(len(group1.hits), 2)

        # First hit of queryA_seq1 (to subjectB_seq1)
        hit1_1 = next(h for h in group1.hits if h.subject_id == "subjectB_seq1")
        self.assertEqual(hit1_1.query_id, "queryA_seq1")
        self.assertEqual(hit1_1.subject_id, "subjectB_seq1")
        self.assertAlmostEqual(hit1_1.pident, (85/90)*100)
        self.assertEqual(hit1_1.length, 90)
        self.assertEqual(hit1_1.evalue, 1e-25)
        self.assertEqual(hit1_1.bitscore, 250.0)
        self.assertEqual(hit1_1.query_len, 100)
        self.assertEqual(hit1_1.subject_len, 120)
        self.assertEqual(hit1_1.query_start, 1)
        self.assertEqual(hit1_1.query_end, 90)
        self.assertEqual(hit1_1.subject_start, 10)
        self.assertEqual(hit1_1.subject_end, 99)

        # Second hit of queryA_seq1 (to subjectB_seq2)
        hit1_2 = next(h for h in group1.hits if h.subject_id == "subjectB_seq2")
        self.assertEqual(hit1_2.subject_id, "subjectB_seq2")
        self.assertEqual(hit1_2.bitscore, 150.0)
        self.assertEqual(hit1_2.evalue, 1e-15)
        self.assertEqual(hit1_2.query_len, 100)
        self.assertEqual(hit1_2.subject_len, 110)


        # Check queryA_seq2
        self.assertIn("queryA_seq2", hit_groups)
        group2 = hit_groups["queryA_seq2"]
        self.assertEqual(group2.query_id, "queryA_seq2")
        self.assertEqual(len(group2.hits), 1)
        hit2_1 = group2.hits[0]
        self.assertEqual(hit2_1.subject_id, "subjectB_seq1")
        self.assertEqual(hit2_1.bitscore, 180.0)
        self.assertEqual(hit2_1.query_len, 80)
        self.assertEqual(hit2_1.subject_len, 120)


    def test_parse_xml_no_hits(self):
        hit_groups = rbh_handler.parse_blast_xml(self.xml_file_none, self.query_lengths, self.subject_lengths)
        self.assertIn("query_no_hit", hit_groups) # Group is created for the query
        self.assertEqual(len(hit_groups["query_no_hit"].hits), 0) # But it has no hits

    def test_parse_xml_empty_file_content(self):
        # Biopython's NCBIXML.parse on an empty file will likely raise an exception
        # The function should catch this and return an empty dict
        hit_groups = rbh_handler.parse_blast_xml(self.xml_file_empty_content, self.query_lengths, self.subject_lengths)
        self.assertEqual(len(hit_groups), 0)

    def test_parse_non_existent_xml(self):
        hit_groups = rbh_handler.parse_blast_xml("non_existent.xml", self.query_lengths, self.subject_lengths)
        self.assertEqual(len(hit_groups), 0)


class TestFindRBHPairs(unittest.TestCase):
    def _create_hit(self, q_id, s_id, score, evalue, q_len=100, s_len=100):
        return rbh_handler.RBHHit(q_id, s_id, 100, 100, evalue, score, q_len, s_len)

    def _create_group(self, q_id, hits_list):
        group = rbh_handler.RBHHitGroup(q_id)
        for hit in hits_list:
            group.add_hit(hit)
        return group

    def test_simple_one_to_one_rbh(self):
        # A1's best hit is B1; B1's best hit is A1
        hits_A_vs_B = {
            "A1": self._create_group("A1", [self._create_hit("A1", "B1", 200, 1e-20)])
        }
        hits_B_vs_A = {
            "B1": self._create_group("B1", [self._create_hit("B1", "A1", 200, 1e-20)])
        }
        rbh_pairs = rbh_handler.find_rbh_pairs(hits_A_vs_B, hits_B_vs_A)
        self.assertEqual(len(rbh_pairs), 1)
        self.assertEqual(rbh_pairs[0][0], "A1")
        self.assertEqual(rbh_pairs[0][1], "B1")

    def test_no_rbh_different_reciprocal_best(self):
        # A1's best hit is B1; B1's best hit is A2
        hits_A_vs_B = {
            "A1": self._create_group("A1", [self._create_hit("A1", "B1", 200, 1e-20)])
        }
        hits_B_vs_A = {
            "B1": self._create_group("B1", [self._create_hit("B1", "A2", 200, 1e-20)]) # B1 hits A2 better
        }
        rbh_pairs = rbh_handler.find_rbh_pairs(hits_A_vs_B, hits_B_vs_A)
        self.assertEqual(len(rbh_pairs), 0)

    def test_one_to_many_rbh_tie(self):
        # A1's best hits are B1, B2 (tie). B1's best hit is A1. B2's best hit is A1.
        hit_A1_B1 = self._create_hit("A1", "B1", 200, 1e-20)
        hit_A1_B2 = self._create_hit("A1", "B2", 200, 1e-20) # Same score & evalue
        hits_A_vs_B = {
            "A1": self._create_group("A1", [hit_A1_B1, hit_A1_B2])
        }

        hit_B1_A1 = self._create_hit("B1", "A1", 190, 1e-19)
        hit_B2_A1 = self._create_hit("B2", "A1", 210, 1e-21)
        hits_B_vs_A = {
            "B1": self._create_group("B1", [hit_B1_A1]),
            "B2": self._create_group("B2", [hit_B2_A1])
        }
        rbh_pairs = rbh_handler.find_rbh_pairs(hits_A_vs_B, hits_B_vs_A)
        self.assertEqual(len(rbh_pairs), 2)
        pair_ids = sorted([(p[0], p[1]) for p in rbh_pairs])
        self.assertIn(("A1", "B1"), pair_ids)
        self.assertIn(("A1", "B2"), pair_ids)

    def test_one_sided_hits(self):
        # A1 hits B1, but B1 has no hits back to A
        hits_A_vs_B = {
            "A1": self._create_group("A1", [self._create_hit("A1", "B1", 200, 1e-20)])
        }
        hits_B_vs_A = {
            "B1": self._create_group("B1", []) # B1 has no hits
        }
        rbh_pairs = rbh_handler.find_rbh_pairs(hits_A_vs_B, hits_B_vs_A)
        self.assertEqual(len(rbh_pairs), 0)

    def test_empty_inputs(self):
        rbh_pairs_empty_A = rbh_handler.find_rbh_pairs({}, {"B1": self._create_group("B1", [])})
        self.assertEqual(len(rbh_pairs_empty_A), 0)
        rbh_pairs_empty_B = rbh_handler.find_rbh_pairs({"A1": self._create_group("A1", [])}, {})
        self.assertEqual(len(rbh_pairs_empty_B), 0)
        rbh_pairs_both_empty = rbh_handler.find_rbh_pairs({}, {})
        self.assertEqual(len(rbh_pairs_both_empty), 0)

    def test_complex_ties_mixed_reciprocity(self):
        # A1 best hits are B1, B2 (tie)
        # B1 best hit is A1
        # B2 best hit is A2 (not A1)
        # Expected: (A1,B1) is RBH, (A1,B2) is not.
        hit_A1_B1 = self._create_hit("A1", "B1", 200, 1e-20)
        hit_A1_B2 = self._create_hit("A1", "B2", 200, 1e-20) # Tie with B1
        hits_A_vs_B = {
            "A1": self._create_group("A1", [hit_A1_B1, hit_A1_B2])
        }

        hit_B1_A1 = self._create_hit("B1", "A1", 190, 1e-19) # B1's best hit to A1
        hit_B2_A2 = self._create_hit("B2", "A2", 210, 1e-21) # B2's best hit to A2
        hits_B_vs_A = {
            "B1": self._create_group("B1", [hit_B1_A1]),
            "B2": self._create_group("B2", [hit_B2_A2])
        }
        rbh_pairs = rbh_handler.find_rbh_pairs(hits_A_vs_B, hits_B_vs_A)
        self.assertEqual(len(rbh_pairs), 1)
        self.assertEqual(rbh_pairs[0][0], "A1")
        self.assertEqual(rbh_pairs[0][1], "B1")

    def test_reciprocal_but_not_best_overall(self):
        # A1 hits B1 (best), B2 (worse)
        # B1 hits A1 (best)
        # Should still be an RBH (A1,B1)
        hit_A1_B1_best = self._create_hit("A1", "B1", 200, 1e-20)
        hit_A1_B2_worse = self._create_hit("A1", "B2", 100, 1e-10)
        hits_A_vs_B = {
            "A1": self._create_group("A1", [hit_A1_B1_best, hit_A1_B2_worse])
        }

        hit_B1_A1_best = self._create_hit("B1", "A1", 190, 1e-19)
        hits_B_vs_A = {
            "B1": self._create_group("B1", [hit_B1_A1_best])
        }
        rbh_pairs = rbh_handler.find_rbh_pairs(hits_A_vs_B, hits_B_vs_A)
        self.assertEqual(len(rbh_pairs), 1)
        self.assertEqual(rbh_pairs[0][0], "A1")
        self.assertEqual(rbh_pairs[0][1], "B1")

    def test_subject_not_in_return_hits(self):
        # A1 hits B1 (best)
        # B1 is not a query in hits_B_vs_A (e.g. B1 was filtered out or had no hits)
        hits_A_vs_B = {
            "A1": self._create_group("A1", [self._create_hit("A1", "B1", 200, 1e-20)])
        }
        hits_B_vs_A = { # B1 is missing as a query key
            "B2": self._create_group("B2", [self._create_hit("B2", "A1", 200, 1e-20)])
        }
        rbh_pairs = rbh_handler.find_rbh_pairs(hits_A_vs_B, hits_B_vs_A)
        self.assertEqual(len(rbh_pairs), 0)

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
