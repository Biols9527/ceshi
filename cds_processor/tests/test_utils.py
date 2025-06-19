import unittest
import os
import sys
import logging
import hashlib

# Add project root to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from cds_processor.cds_processor import utils

class TestUtils(unittest.TestCase):

    def setUp(self):
        self.test_log_dir = "temp_test_logs" # Renamed to avoid conflict with actual logs
        self.test_log_file = "test_app.log"
        self.full_log_path = os.path.join(self.test_log_dir, self.test_log_file)

        # Create the log directory if it doesn't exist
        os.makedirs(self.test_log_dir, exist_ok=True)

        self.test_hash_file = os.path.join(self.test_log_dir, "sample_for_hash.txt")
        with open(self.test_hash_file, "w") as f:
            f.write("This is a test file for hashing.")

    def tearDown(self):
        # Clean up log files and directory
        # Close all handlers to release the log file
        # Get the specific logger instance returned by setup_logging if possible, or root
        logger_instance = logging.getLogger("cds_processor") # or logging.getLogger()
        for handler in logger_instance.handlers[:]:
            handler.close()
            logger_instance.removeHandler(handler)
        for handler in logging.getLogger().handlers[:]: # Clear root handlers too
            handler.close()
            logging.getLogger().removeHandler(handler)

        if os.path.exists(self.full_log_path):
            os.remove(self.full_log_path)
        if os.path.exists(self.test_hash_file):
            os.remove(self.test_hash_file)
        if os.path.exists(self.test_log_dir) and not os.listdir(self.test_log_dir):
            try:
                os.rmdir(self.test_log_dir)
            except OSError as e: # Might fail if logs from other tests ended up here
                print(f"Could not remove {self.test_log_dir}: {e}")


    def test_setup_logging_new(self):
        # Use a unique log file for this test to avoid interference
        specific_log_file = os.path.join(self.test_log_dir, "specific_test_run.log")
        if os.path.exists(specific_log_file):
            os.remove(specific_log_file)

        logger = utils.setup_logging(log_level_str="DEBUG", log_file=specific_log_file)
        self.assertTrue(os.path.exists(specific_log_file))
        self.assertEqual(logger.name, "cds_processor")
        self.assertEqual(logger.level, logging.DEBUG)

        test_message = "This is a debug test message for the new logging setup."
        logger.debug(test_message)

        with open(specific_log_file, "r") as f:
            log_content = f.read()
        self.assertIn(test_message, log_content)
        self.assertIn("[DEBUG]", log_content)
        self.assertIn("[cds_processor]", log_content) # Check package logger name in format

        # Test if root logger also got configured (though we mostly care about the package logger)
        root_logger_test_message = "Root logger test message."
        logging.info(root_logger_test_message) # Use root logger directly

        # Explicitly flush handlers before reading the log file
        for handler in logging.getLogger().handlers:
            handler.flush()
        for handler in logging.getLogger("cds_processor").handlers:
            handler.flush()

        with open(specific_log_file, "r") as f:
            log_content = f.read()
        self.assertIn(test_message, log_content) # Re-check this too, just in case
        self.assertIn(root_logger_test_message, log_content)


    def test_hash_file_content_valid_file(self):
        expected_hash = hashlib.md5("This is a test file for hashing.".encode()).hexdigest()
        actual_hash = utils.hash_file_content(self.test_hash_file)
        self.assertEqual(actual_hash, expected_hash)

    def test_hash_file_content_not_found(self):
        actual_hash = utils.hash_file_content("non_existent_file.txt")
        self.assertIsNone(actual_hash)
        # Check logs for error message (optional, requires log capturing or specific log file)

    def test_hash_file_content_empty_file(self):
        empty_file = os.path.join(self.test_log_dir, "empty_for_hash.txt")
        with open(empty_file, "w") as f:
            pass
        expected_hash = hashlib.md5("".encode()).hexdigest()
        actual_hash = utils.hash_file_content(empty_file)
        self.assertEqual(actual_hash, expected_hash)
        os.remove(empty_file)

if __name__ == "__main__":
    # Ensure logs dir exists for the default logger if tests are run directly
    # Though setup_logging in tests now uses a specific file.
    if not os.path.exists("logs"): # Default log dir used by setup_logging if no file is passed
        os.makedirs("logs")
    unittest.main()
