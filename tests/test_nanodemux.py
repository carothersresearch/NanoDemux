#!/usr/bin/env python3
"""
Tests for the comprehensive nanodemux workflow script.
"""

import unittest
import tempfile
import os
import shutil
import subprocess
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TestNanodemuxWorkflow(unittest.TestCase):
    """Test the comprehensive nanodemux workflow script."""
    
    # Test barcode sequences
    TEST_ROW_PRIMER = "AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC"
    TEST_COL_PRIMER = "CAAGCAGAAGACGGCATACGAGATATTACTCGGTCTCGTGGGCTCGG"
    TEST_INSERT = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 5
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.output_dir = os.path.join(self.temp_dir, "output")
        
        # Create test barcode CSV
        self.barcode_csv = os.path.join(self.temp_dir, "test_barcodes.csv")
        self._create_barcode_csv()
        
        # Create test FASTQ file
        self.fastq_file = os.path.join(self.temp_dir, "test.fastq")
        self._create_test_fastq()
        
        # Path to nanodemux script
        self.nanodemux_script = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "nanodemux"
        )
    
    def tearDown(self):
        """Clean up temporary files."""
        shutil.rmtree(self.temp_dir)
    
    def _create_barcode_csv(self):
        """Create a test barcode CSV file."""
        barcodes = [
            ("A1", "R oDA373.D501", self.TEST_ROW_PRIMER),
            ("A1", "F oDA361.D701", self.TEST_COL_PRIMER),
            ("A2", "R oDA373.D501", self.TEST_ROW_PRIMER),
            ("A2", "F oDA362.D702", "CAAGCAGAAGACGGCATACGAGATTCCGGAGAGTCTCGTGGGCTCGG"),
            ("B1", "R oDA374.D502", "AATGATACGGCGACCACCGAGATCTACACATAGAGGCTCGTCGGCAGCGTC"),
            ("B1", "F oDA361.D701", self.TEST_COL_PRIMER),
        ]
        
        df = pd.DataFrame(barcodes, columns=["Well Position", "Sequence Name", "Sequence"])
        df.to_csv(self.barcode_csv, index=False)
    
    def _create_test_fastq(self):
        """Create a test FASTQ file with realistic sequences."""
        records = []
        
        test_cases = [
            (self.TEST_ROW_PRIMER + self.TEST_INSERT + self.TEST_COL_PRIMER, "read_both_primers"),
            (self.TEST_ROW_PRIMER + self.TEST_INSERT, "read_row_only"),
            (self.TEST_COL_PRIMER + self.TEST_INSERT, "read_col_only"),
            (self.TEST_INSERT, "read_no_primers"),
            (self.TEST_ROW_PRIMER + self.TEST_INSERT + self.TEST_COL_PRIMER, "read_both_primers_2"),
        ]
        
        for seq_str, read_id in test_cases:
            rec = SeqRecord(Seq(seq_str), id=read_id)
            rec.letter_annotations["phred_quality"] = [40] * len(seq_str)
            records.append(rec)
        
        SeqIO.write(records, self.fastq_file, "fastq")
    
    def test_nanodemux_help(self):
        """Test that nanodemux --help works."""
        result = subprocess.run(
            [sys.executable, self.nanodemux_script, "--help"],
            capture_output=True,
            text=True
        )
        
        self.assertEqual(result.returncode, 0)
        self.assertIn("NanoDemux Comprehensive Workflow", result.stdout)
        self.assertIn("Raw Quality Analysis Options", result.stdout)
        self.assertIn("Demultiplexing Options", result.stdout)
        self.assertIn("Alignment Options", result.stdout)
    
    def test_nanodemux_basic_workflow(self):
        """Test basic nanodemux workflow with minimal options."""
        result = subprocess.run(
            [
                sys.executable, self.nanodemux_script,
                self.fastq_file,
                self.barcode_csv,
                "--output", self.output_dir,
                "--cpus", "1"
            ],
            capture_output=True,
            text=True
        )
        
        # Print output for debugging if test fails
        if result.returncode != 0:
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)
        
        self.assertEqual(result.returncode, 0)
        
        # Check that output directory exists
        self.assertTrue(os.path.isdir(self.output_dir))
        
        # Check that all expected subdirectories exist
        self.assertTrue(os.path.isdir(os.path.join(self.output_dir, "1_raw_quality_report")))
        self.assertTrue(os.path.isdir(os.path.join(self.output_dir, "2_demultiplexed")))
        self.assertTrue(os.path.isdir(os.path.join(self.output_dir, "3_demux_quality_report")))
        self.assertTrue(os.path.isdir(os.path.join(self.output_dir, "4_aligned")))
        
        # Check comprehensive report exists
        report_file = os.path.join(self.output_dir, "comprehensive_report.html")
        self.assertTrue(os.path.exists(report_file))
        
        # Verify report content
        with open(report_file, 'r') as f:
            content = f.read()
            self.assertIn("NanoDemux Comprehensive Workflow Report", content)
            self.assertIn("Raw Quality Report", content)
            self.assertIn("Demultiplexed Data", content)
            self.assertIn("Aligned Sequences", content)
    
    def test_nanodemux_with_skip_options(self):
        """Test nanodemux with skip options."""
        result = subprocess.run(
            [
                sys.executable, self.nanodemux_script,
                self.fastq_file,
                self.barcode_csv,
                "--output", self.output_dir,
                "--skip-raw-qc",
                "--skip-alignment",
                "--cpus", "1"
            ],
            capture_output=True,
            text=True
        )
        
        self.assertEqual(result.returncode, 0)
        
        # Check that skipped directories don't exist
        self.assertFalse(os.path.isdir(os.path.join(self.output_dir, "1_raw_quality_report")))
        self.assertFalse(os.path.isdir(os.path.join(self.output_dir, "4_aligned")))
        
        # Check that non-skipped directories exist
        self.assertTrue(os.path.isdir(os.path.join(self.output_dir, "2_demultiplexed")))
        self.assertTrue(os.path.isdir(os.path.join(self.output_dir, "3_demux_quality_report")))
    
    def test_nanodemux_with_custom_parameters(self):
        """Test nanodemux with custom parameters for each step."""
        result = subprocess.run(
            [
                sys.executable, self.nanodemux_script,
                self.fastq_file,
                self.barcode_csv,
                "--output", self.output_dir,
                "--min-length", "40",
                "--max-penalty", "70",
                "--flank", "80",
                "--min-quality", "15",
                "--cpus", "1"
            ],
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)
        
        self.assertEqual(result.returncode, 0)
        
        # Check that comprehensive report exists
        report_file = os.path.join(self.output_dir, "comprehensive_report.html")
        self.assertTrue(os.path.exists(report_file))
    
    def test_nanodemux_invalid_input(self):
        """Test nanodemux with invalid input files."""
        result = subprocess.run(
            [
                sys.executable, self.nanodemux_script,
                "nonexistent.fastq",
                self.barcode_csv,
                "--output", self.output_dir
            ],
            capture_output=True,
            text=True
        )
        
        # Should fail with non-zero exit code
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not found", result.stderr)
    
    def test_nanodemux_creates_barcode_stats(self):
        """Test that demultiplexing creates barcode statistics."""
        result = subprocess.run(
            [
                sys.executable, self.nanodemux_script,
                self.fastq_file,
                self.barcode_csv,
                "--output", self.output_dir,
                "--cpus", "1"
            ],
            capture_output=True,
            text=True
        )
        
        self.assertEqual(result.returncode, 0)
        
        # Check barcode stats file exists
        stats_file = os.path.join(self.output_dir, "2_demultiplexed", "barcode_stats.csv")
        self.assertTrue(os.path.exists(stats_file))
        
        # Verify it's a valid CSV that can be read
        df = pd.read_csv(stats_file)
        # The CSV is a matrix format with rows (A-H) and columns (1-12, X)
        self.assertGreater(len(df), 0, "Barcode stats CSV should not be empty")


if __name__ == '__main__':
    unittest.main()
