#!/usr/bin/env python3
"""
Integration tests for the full demultiplexing pipeline.
"""

import unittest
import tempfile
import os
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from demux_barcodes import main, load_barcodes, process_chunk


class TestDemultiplexingPipeline(unittest.TestCase):
    """Integration tests for the full demultiplexing pipeline."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.output_dir = os.path.join(self.temp_dir, "output")
        
        # Create test barcode CSV with real-looking sequences
        self.barcode_csv = os.path.join(self.temp_dir, "test_barcodes.csv")
        self._create_barcode_csv()
        
        # Create test FASTQ file
        self.fastq_file = os.path.join(self.temp_dir, "test.fastq")
        self._create_test_fastq()

    def tearDown(self):
        """Clean up temporary files."""
        shutil.rmtree(self.temp_dir)

    def _create_barcode_csv(self):
        """Create a test barcode CSV file."""
        barcodes = [
            ("A1", "R oDA373.D501", "AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC"),
            ("A1", "F oDA361.D701", "CAAGCAGAAGACGGCATACGAGATATTACTCGGTCTCGTGGGCTCGG"),
            ("A2", "R oDA373.D501", "AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC"),
            ("A2", "F oDA362.D702", "CAAGCAGAAGACGGCATACGAGATTCCGGAGAGTCTCGTGGGCTCGG"),
            ("B1", "R oDA374.D502", "AATGATACGGCGACCACCGAGATCTACACATAGAGGCTCGTCGGCAGCGTC"),
            ("B1", "F oDA361.D701", "CAAGCAGAAGACGGCATACGAGATATTACTCGGTCTCGTGGGCTCGG"),
        ]
        
        df = pd.DataFrame(barcodes, columns=["Well Position", "Sequence Name", "Sequence"])
        df.to_csv(self.barcode_csv, index=False)

    def _create_test_fastq(self):
        """Create a test FASTQ file with realistic sequences."""
        records = []
        
        # Create a read with both row and column barcodes
        row_primer = "AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC"
        col_primer = "CAAGCAGAAGACGGCATACGAGATATTACTCGGTCTCGTGGGCTCGG"
        insert = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 5
        
        test_cases = [
            (row_primer + insert + col_primer, "read_both_primers"),
            (row_primer + insert, "read_row_only"),
            (col_primer + insert, "read_col_only"),
            ("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 5, "read_no_primers"),
            ("ATCG" * 5, "read_short"),
        ]
        
        for seq_str, read_id in test_cases:
            rec = SeqRecord(Seq(seq_str), id=read_id)
            rec.letter_annotations["phred_quality"] = [40] * len(seq_str)
            records.append(rec)
        
        SeqIO.write(records, self.fastq_file, "fastq")

    def test_pipeline_creates_output_dir(self):
        """Test that pipeline creates output directory."""
        main(
            self.fastq_file,
            self.barcode_csv,
            self.output_dir,
            min_length=50,
            max_penalty=60,
            cpus=1,
            flank=100
        )
        
        self.assertTrue(os.path.isdir(self.output_dir), "Output directory should be created")

    def test_pipeline_creates_stats_file(self):
        """Test that pipeline creates barcode_stats.csv."""
        main(
            self.fastq_file,
            self.barcode_csv,
            self.output_dir,
            min_length=50,
            max_penalty=60,
            cpus=1,
            flank=100
        )
        
        stats_file = os.path.join(self.output_dir, "barcode_stats.csv")
        self.assertTrue(os.path.isfile(stats_file), "Stats file should be created")

    def test_pipeline_stats_contents(self):
        """Test that stats file has expected structure."""
        main(
            self.fastq_file,
            self.barcode_csv,
            self.output_dir,
            min_length=50,
            max_penalty=60,
            cpus=1,
            flank=100
        )
        
        stats_file = os.path.join(self.output_dir, "barcode_stats.csv")
        df = pd.read_csv(stats_file, index_col=0)
        
        # Check for expected rows (A1-A2, B1, X for row_only, Stats)
        self.assertIn("A", df.index, "Should have row A")
        self.assertIn("Stats", df.index, "Should have Stats row")

    def test_pipeline_creates_fastq_outputs(self):
        """Test that pipeline creates FASTQ files for matched reads."""
        main(
            self.fastq_file,
            self.barcode_csv,
            self.output_dir,
            min_length=50,
            max_penalty=60,
            cpus=1,
            flank=100
        )
        
        # Check if any well-specific FASTQ files were created
        fastq_files = [f for f in os.listdir(self.output_dir) if f.endswith("_reads.fastq")]
        self.assertGreater(len(fastq_files), 0, "Should create at least one FASTQ file")

    def test_barcode_loading_consistency(self):
        """Test that loaded barcodes are consistent."""
        row_map, col_map = load_barcodes(self.barcode_csv)
        
        # Row primers should be unique
        self.assertEqual(len(set(row_map.values())), 2, "Should have 2 unique rows (A, B)")
        
        # Column primers should be unique
        self.assertEqual(len(set(col_map.values())), 2, "Should have 2 unique columns")


class TestProcessChunk(unittest.TestCase):
    """Test the process_chunk function in isolation."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.barcode_csv = os.path.join(self.temp_dir, "test_barcodes.csv")
        self._create_barcode_csv()
        self.row_map, self.col_map = load_barcodes(self.barcode_csv)

    def tearDown(self):
        """Clean up temporary files."""
        shutil.rmtree(self.temp_dir)

    def _create_barcode_csv(self):
        """Create a test barcode CSV file."""
        barcodes = [
            ("A1", "R test", "AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC"),
            ("A1", "F test", "CAAGCAGAAGACGGCATACGAGATATTACTCGGTCTCGTGGGCTCGG"),
        ]
        df = pd.DataFrame(barcodes, columns=["Well Position", "Sequence Name", "Sequence"])
        df.to_csv(self.barcode_csv, index=False)

    def test_process_chunk_returns_tuple(self):
        """Test that process_chunk returns stats and reads."""
        rec = SeqRecord(Seq("ATCG" * 50), id="test")
        rec.letter_annotations["phred_quality"] = [40] * 200
        records = [rec]
        
        args = (records, self.row_map, self.col_map, 50, 60, 100, [])
        stats, reads = process_chunk(args)
        
        self.assertIsInstance(stats, dict, "Stats should be a dict")
        self.assertIsInstance(reads, dict, "Reads should be a dict")

    def test_process_chunk_counts_total(self):
        """Test that process_chunk counts total reads."""
        records = []
        for i in range(2):
            rec = SeqRecord(Seq("ATCG" * 50), id=f"test{i+1}")
            rec.letter_annotations["phred_quality"] = [40] * 200
            records.append(rec)
        
        args = (records, self.row_map, self.col_map, 50, 60, 100, [])
        stats, reads = process_chunk(args)
        
        self.assertEqual(stats['GLOBAL']['total'], 2, "Should count 2 reads")


class TestMultiFileProcessing(unittest.TestCase):
    """Test processing multiple FASTQ files from a directory."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.input_dir = os.path.join(self.temp_dir, "input_data")
        os.makedirs(self.input_dir)
        
        # Create test barcode CSV
        self.barcode_csv = os.path.join(self.temp_dir, "test_barcodes.csv")
        self._create_barcode_csv()
        
        # Create multiple test FASTQ files
        self._create_test_fastq_files()

    def tearDown(self):
        """Clean up temporary files."""
        shutil.rmtree(self.temp_dir)

    def _create_barcode_csv(self):
        """Create a test barcode CSV file."""
        barcodes = [
            ("A1", "R oDA373.D501", "AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC"),
            ("A1", "F oDA361.D701", "CAAGCAGAAGACGGCATACGAGATATTACTCGGTCTCGTGGGCTCGG"),
        ]
        df = pd.DataFrame(barcodes, columns=["Well Position", "Sequence Name", "Sequence"])
        df.to_csv(self.barcode_csv, index=False)

    def _create_test_fastq_files(self):
        """Create multiple test FASTQ files."""
        for file_num in [1, 2]:
            fastq_file = os.path.join(self.input_dir, f"test_{file_num}.fastq")
            records = []
            row_primer = "AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC"
            col_primer = "CAAGCAGAAGACGGCATACGAGATATTACTCGGTCTCGTGGGCTCGG"
            insert = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 5
            
            seq_str = row_primer + insert + col_primer
            rec = SeqRecord(Seq(seq_str), id=f"read_{file_num}")
            rec.letter_annotations["phred_quality"] = [40] * len(seq_str)
            records.append(rec)
            
            SeqIO.write(records, fastq_file, "fastq")

    def test_directory_processing(self):
        """Test that directory processing creates correct output structure."""
        main(
            self.input_dir,
            self.barcode_csv,
            "demuxed",
            min_length=50,
            max_penalty=60,
            cpus=1,
            flank=100
        )
        
        # Check that output directory was created with correct structure
        expected_base = os.path.join("demplex_data", "input_data")
        self.assertTrue(os.path.isdir(expected_base), "Base output directory should exist")
        
        # Check that subdirectories for each file were created
        for file_num in [1, 2]:
            file_dir = os.path.join(expected_base, f"test_{file_num}")
            self.assertTrue(os.path.isdir(file_dir), f"Output directory for test_{file_num} should exist")
            
            # Check that stats file was created
            stats_file = os.path.join(file_dir, "barcode_stats.csv")
            self.assertTrue(os.path.isfile(stats_file), f"Stats file should exist for test_{file_num}")
        
        # Cleanup
        shutil.rmtree("demplex_data")

    def test_single_file_default_output(self):
        """Test that single file processing uses demplex_data/<filename> structure."""
        fastq_file = os.path.join(self.input_dir, "test_1.fastq")
        
        main(
            fastq_file,
            self.barcode_csv,
            "demuxed",  # Default value
            min_length=50,
            max_penalty=60,
            cpus=1,
            flank=100
        )
        
        # Check that output directory was created with correct structure
        expected_dir = os.path.join("demplex_data", "test_1")
        self.assertTrue(os.path.isdir(expected_dir), "Output directory should exist in demplex_data")
        
        # Check that stats file was created
        stats_file = os.path.join(expected_dir, "barcode_stats.csv")
        self.assertTrue(os.path.isfile(stats_file), "Stats file should exist")
        
        # Cleanup
        shutil.rmtree("demplex_data")

    def test_single_file_custom_output(self):
        """Test that single file processing respects custom outdir."""
        fastq_file = os.path.join(self.input_dir, "test_1.fastq")
        custom_outdir = os.path.join(self.temp_dir, "custom_output")
        
        main(
            fastq_file,
            self.barcode_csv,
            custom_outdir,
            min_length=50,
            max_penalty=60,
            cpus=1,
            flank=100
        )
        
        # Check that custom output directory was used
        self.assertTrue(os.path.isdir(custom_outdir), "Custom output directory should exist")
        
        # Check that stats file was created
        stats_file = os.path.join(custom_outdir, "barcode_stats.csv")
        self.assertTrue(os.path.isfile(stats_file), "Stats file should exist in custom directory")


if __name__ == "__main__":
    unittest.main()
