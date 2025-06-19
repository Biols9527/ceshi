import unittest
import os
import sys
import shutil # For cleanup in tests

# Add project root to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from cds_processor.cds_processor import file_manager
# from Bio.Seq import Seq # Not used by new file_manager directly
# from Bio.SeqRecord import SeqRecord # Not used by new file_manager directly
# from Bio import SeqIO # Not used by new file_manager directly for writing

class TestFileManagerNew(unittest.TestCase):

    def setUp(self):
        self.test_dir = "temp_fm_test_dir"
        os.makedirs(self.test_dir, exist_ok=True)

        # Create some dummy .cds files for find_cds_files
        with open(os.path.join(self.test_dir, "file1.cds"), "w") as f: f.write("content1")
        with open(os.path.join(self.test_dir, "file2.cds"), "w") as f: f.write("content2")
        with open(os.path.join(self.test_dir, "file.txt"), "w") as f: f.write("textfile")

        self.output_file_path = os.path.join(self.test_dir, "output.fasta")
        self.temp_cleanup_dir = os.path.join(self.test_dir, "to_be_cleaned")
        os.makedirs(self.temp_cleanup_dir, exist_ok=True)
        with open(os.path.join(self.temp_cleanup_dir, "dummy.txt"), "w") as f: f.write("dummy")


    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def test_find_cds_files(self):
        cds_files = file_manager.find_cds_files(self.test_dir)
        self.assertEqual(len(cds_files), 2)
        self.assertIn("file1.cds", cds_files)
        self.assertIn("file2.cds", cds_files)
        self.assertNotIn("file.txt", cds_files)

    def test_find_cds_files_non_existent_dir(self):
        cds_files = file_manager.find_cds_files("non_existent_dir_for_fm")
        self.assertEqual(len(cds_files), 0)

    def test_find_cds_files_empty_dir(self):
        empty_dir = os.path.join(self.test_dir, "empty_subdir")
        os.makedirs(empty_dir, exist_ok=True)
        cds_files = file_manager.find_cds_files(empty_dir)
        self.assertEqual(len(cds_files), 0)

    def test_write_output_valid(self):
        species_data = {
            "speciesA": "ATGCATGC",
            "speciesB": "CGTACGTA"
        }
        result = file_manager.write_output(species_data, self.output_file_path)
        self.assertTrue(result)
        self.assertTrue(os.path.exists(self.output_file_path))
        with open(self.output_file_path, "r") as f:
            content = f.read()
        self.assertIn(">speciesA\nATGCATGC\n", content)
        self.assertIn(">speciesB\nCGTACGTA\n", content)

    def test_write_output_empty_data(self):
        result = file_manager.write_output({}, self.output_file_path)
        self.assertFalse(result) # No sequences provided
        self.assertFalse(os.path.exists(self.output_file_path)) # Should not create file for empty data

    def test_write_output_invalid_sequence_type(self):
        species_data = {
            "speciesC": 12345 # Not a string
        }
        # This should log an error and skip the sequence.
        # If other valid sequences exist, it might still return True.
        # For this test, let's assume it's the only sequence.
        # Depending on strictness, it could return False or True but not write this seq.
        # The current file_manager.py skips it and continues.
        result = file_manager.write_output(species_data, self.output_file_path)
        self.assertTrue(result) # File is written, but speciesC is skipped.
        if os.path.exists(self.output_file_path): # Check content if file was written
            with open(self.output_file_path, "r") as f:
                content = f.read()
            self.assertNotIn(">speciesC", content)


    def test_cleanup_temp_directory_valid(self):
        self.assertTrue(os.path.exists(self.temp_cleanup_dir))
        file_manager.cleanup_temp_directory(self.temp_cleanup_dir)
        self.assertFalse(os.path.exists(self.temp_cleanup_dir))

    def test_cleanup_temp_directory_non_existent(self):
        non_existent_dir = os.path.join(self.test_dir, "does_not_exist")
        # Should not raise error
        file_manager.cleanup_temp_directory(non_existent_dir)

    def test_cleanup_temp_directory_is_file(self):
        file_path_as_dir = os.path.join(self.test_dir, "a_file.txt")
        with open(file_path_as_dir, "w") as f: f.write("i am a file")
        # Should not raise error, should log warning
        file_manager.cleanup_temp_directory(file_path_as_dir)
        self.assertTrue(os.path.exists(file_path_as_dir)) # File should still be there


if __name__ == "__main__":
    unittest.main()
