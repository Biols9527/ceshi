import unittest
from unittest import mock
import os
import sys
import tempfile
import shutil
import argparse

# Adjust path to import from cds_processor package
SCRIPT_DIR_TEST = os.path.dirname(os.path.abspath(__file__)) # cds_processor/tests
PROJECT_ROOT_TEST = os.path.dirname(SCRIPT_DIR_TEST) # cds_processor
SCRIPT_MODULE_DIR = os.path.join(PROJECT_ROOT_TEST, "scripts") # cds_processor/scripts
PACKAGE_ROOT_TEST = os.path.join(PROJECT_ROOT_TEST, "cds_processor") # cds_processor/cds_processor

# Ensure scripts directory and package root are in path for imports by the script itself
if SCRIPT_MODULE_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_MODULE_DIR)
if PACKAGE_ROOT_TEST not in sys.path: # For cds_processor.utils etc.
    sys.path.insert(0, PACKAGE_ROOT_TEST)
if PROJECT_ROOT_TEST not in sys.path: # For cds_processor.cds_processor
     sys.path.insert(0, PROJECT_ROOT_TEST)


# Now import the script we want to test
# This is tricky because it's a script, not a module designed for import.
# We can use importlib if it's structured well, or run it as a subprocess.
# For now, let's try to import its functions if possible, or mock heavily.
# Assuming run_rbh_analysis can be imported or its main function parts can be tested.
from cds_processor.scripts import run_rbh_analysis
from cds_processor import fasta_parser # For extract_species_name, used by run_rbh_analysis.extract_sequences_for_species


class TestExtractSequencesForSpecies(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.fasta_file_A = os.path.join(self.temp_dir, "species_test.fasta")
        with open(self.fasta_file_A, "w") as f:
            f.write(">seq1 species_X gene1
AAACC
")
            f.write(">seq2 species_Y gene2
GGGTT
")
            f.write(">seq3 species_X gene3
TTTAAA
")

        self.output_fasta = os.path.join(self.temp_dir, "extracted.fasta")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    @mock.patch('cds_processor.scripts.run_rbh_analysis.logger') # Mock logger inside the script
    def test_extract_specific_species(self, mock_logger):
        lengths, found = run_rbh_analysis.extract_sequences_for_species(
            self.fasta_file_A, "species_X", self.output_fasta, fasta_parser
        )
        self.assertTrue(found)
        self.assertEqual(len(lengths), 2)
        self.assertIn("seq1", lengths)
        self.assertEqual(lengths["seq1"], 5)
        self.assertIn("seq3", lengths)
        self.assertEqual(lengths["seq3"], 6)

        # Verify content of output file
        with open(self.output_fasta, "r") as f:
            content = f.read()
            self.assertIn(">seq1", content)
            self.assertIn(">seq3", content)
            self.assertNotIn(">seq2", content)

    @mock.patch('cds_processor.scripts.run_rbh_analysis.logger')
    def test_extract_all_species(self, mock_logger):
        lengths, found = run_rbh_analysis.extract_sequences_for_species(
            self.fasta_file_A, None, self.output_fasta, fasta_parser
        )
        self.assertTrue(found)
        self.assertEqual(len(lengths), 3)
        self.assertIn("seq1", lengths)
        self.assertIn("seq2", lengths)
        self.assertIn("seq3", lengths)

    @mock.patch('cds_processor.scripts.run_rbh_analysis.logger')
    def test_extract_species_not_found(self, mock_logger):
        lengths, found = run_rbh_analysis.extract_sequences_for_species(
            self.fasta_file_A, "species_Z", self.output_fasta, fasta_parser
        )
        self.assertFalse(found) # Found should be False
        self.assertEqual(len(lengths), 0)
        mock_logger.warning.assert_called_with(f"No sequences found for species 'species_Z' in {self.fasta_file_A}.")


class TestRunRBHAnalysisMain(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.input_A = os.path.join(self.temp_dir, "inputA.fasta")
        self.input_B = os.path.join(self.temp_dir, "inputB.fasta")
        self.output_rbh = os.path.join(self.temp_dir, "output.rbh.tsv")

        with open(self.input_A, "w") as f:
            f.write(">A1 species_A
AAA
")
        with open(self.input_B, "w") as f:
            f.write(">B1 species_B
CCC
")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    @mock.patch('cds_processor.scripts.run_rbh_analysis.utils.setup_logging')
    @mock.patch('cds_processor.scripts.run_rbh_analysis.rbh_handler.perform_blast_for_rbh')
    @mock.patch('cds_processor.scripts.run_rbh_analysis.rbh_handler.parse_blast_xml')
    @mock.patch('cds_processor.scripts.run_rbh_analysis.rbh_handler.find_rbh_pairs')
    @mock.patch('builtins.open', new_callable=mock.mock_open) # Mocks file writing for output
    def test_main_flow_successful_run(self, mock_file_open, mock_find_rbh, mock_parse_xml, mock_perform_blast, mock_setup_logging):
        # Configure mocks
        mock_setup_logging.return_value = mock.MagicMock() # Mock logger object
        mock_perform_blast.side_effect = ["/mock/path/A_vs_B.xml", "/mock/path/B_vs_A.xml"] # Return different paths for two calls

        mock_hits_A_vs_B = {"A1": mock.MagicMock()}
        mock_hits_B_vs_A = {"B1": mock.MagicMock()}
        mock_parse_xml.side_effect = [mock_hits_A_vs_B, mock_hits_B_vs_A]

        mock_rbh_pair = ("A1", "B1", mock.MagicMock(), mock.MagicMock())
        mock_rbh_pair[2].evalue=1e-10; mock_rbh_pair[2].bitscore=100; mock_rbh_pair[2].pident=90; mock_rbh_pair[2].length=3
        mock_rbh_pair[2].query_len=3; mock_rbh_pair[2].subject_len=3
        mock_rbh_pair[2].query_start=1; mock_rbh_pair[2].query_end=3; mock_rbh_pair[2].subject_start=1; mock_rbh_pair[2].subject_end=3

        mock_rbh_pair[3].evalue=1e-10; mock_rbh_pair[3].bitscore=100; mock_rbh_pair[3].pident=90; mock_rbh_pair[3].length=3
        mock_rbh_pair[3].query_start=1; mock_rbh_pair[3].query_end=3; mock_rbh_pair[3].subject_start=1; mock_rbh_pair[3].subject_end=3


        mock_find_rbh.return_value = [mock_rbh_pair]

        # Prepare args for the script's main function
        args = [
            "--input_fasta_A", self.input_A,
            "--input_fasta_B", self.input_B,
            "--output_rbh_file", self.output_rbh,
            "--blastn_path", "/fake/blastn",
            "--makeblastdb_path", "/fake/makeblastdb",
            "--temp_dir", os.path.join(self.temp_dir, "run_temp"),
            "--log_level", "INFO",
            # "--clean_temp" # Not testing cleanup here
        ]

        with mock.patch.object(sys, 'argv', ['run_rbh_analysis.py'] + args):
            run_rbh_analysis.main()

        # Assertions
        mock_setup_logging.assert_called()
        self.assertEqual(mock_perform_blast.call_count, 2)
        # Check calls to perform_blast_for_rbh
        # Call 1 (A vs B)
        call_A_vs_B_args = mock_perform_blast.call_args_list[0][1] # Get kwargs of first call
        self.assertTrue(call_A_vs_B_args['query_fasta_path'].endswith("input_A_selected.fasta"))
        self.assertTrue(call_A_vs_B_args['target_fasta_path'].endswith("input_B_selected.fasta"))

        # Call 2 (B vs A)
        call_B_vs_A_args = mock_perform_blast.call_args_list[1][1] # Get kwargs of second call
        self.assertTrue(call_B_vs_A_args['query_fasta_path'].endswith("input_B_selected.fasta"))
        self.assertTrue(call_B_vs_A_args['target_fasta_path'].endswith("input_A_selected.fasta"))

        self.assertEqual(mock_parse_xml.call_count, 2)
        mock_parse_xml.assert_any_call("/mock/path/A_vs_B.xml", {"A1": 3}, {"B1": 3})
        mock_parse_xml.assert_any_call("/mock/path/B_vs_A.xml", {"B1": 3}, {"A1": 3})

        mock_find_rbh.assert_called_once_with(mock_hits_A_vs_B, mock_hits_B_vs_A)

        mock_file_open.assert_called_once_with(self.output_rbh, "w")
        # Check if header and data line were written (simplified check)
        handle = mock_file_open()
        self.assertIn("SeqA_ID\tSeqB_ID", handle.write.call_args_list[0][0][0]) # Header
        self.assertIn("A1\tB1", handle.write.call_args_list[1][0][0]) # Data

    @mock.patch('cds_processor.scripts.run_rbh_analysis.utils.setup_logging')
    @mock.patch('cds_processor.scripts.run_rbh_analysis.extract_sequences_for_species')
    @mock.patch('sys.exit') # Mock sys.exit to prevent test runner from exiting
    def test_main_flow_no_sequences_A(self, mock_sys_exit, mock_extract_seqs, mock_setup_logging):
        mock_setup_logging.return_value = mock.MagicMock()
        # Simulate extract_sequences_for_species returning no sequences for A
        mock_extract_seqs.side_effect = [({}, False), ({"B1":3}, True)]


        args = [
            "--input_fasta_A", self.input_A,
            "--input_fasta_B", self.input_B,
            "--output_rbh_file", self.output_rbh,
            "--blastn_path", "/fake/blastn",
            "--makeblastdb_path", "/fake/makeblastdb",
        ]
        with mock.patch.object(sys, 'argv', ['run_rbh_analysis.py'] + args):
            run_rbh_analysis.main()

        mock_sys_exit.assert_called_once_with(1)


if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
