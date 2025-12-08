#!/usr/bin/env python3
"""
Unit tests for well alignment and consensus generation.
"""

import unittest
import tempfile
import os
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import pandas as pd

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from align_wells import (
    calculate_quality_weighted_consensus,
    process_well_fastq,
    write_aligned_fastq,
    process_demux_directory,
    write_consensus_csv
)


class TestConsensusCalculation(unittest.TestCase):
    """Test consensus sequence calculation."""
    
    def test_single_sequence(self):
        """Test consensus with a single sequence."""
        sequences = ["ATCGATCG"]
        qualities = [[30, 30, 30, 30, 30, 30, 30, 30]]
        
        consensus_seq, consensus_qual, depth = calculate_quality_weighted_consensus(
            sequences, qualities
        )
        
        self.assertEqual(consensus_seq, "ATCGATCG")
        self.assertEqual(len(consensus_qual), 8)
        self.assertEqual(depth, [1] * 8)
    
    def test_identical_sequences(self):
        """Test consensus with identical sequences."""
        sequences = ["ATCG", "ATCG", "ATCG"]
        qualities = [[30] * 4, [35] * 4, [40] * 4]
        
        consensus_seq, consensus_qual, depth = calculate_quality_weighted_consensus(
            sequences, qualities
        )
        
        self.assertEqual(consensus_seq, "ATCG")
        self.assertEqual(depth, [3] * 4)
    
    def test_quality_weighted_consensus(self):
        """Test that higher quality bases are preferred."""
        # First sequence has high quality A at position 0
        # Second sequence has low quality T at position 0
        sequences = ["ATCG", "TTCG"]
        qualities = [[40, 30, 30, 30], [20, 30, 30, 30]]
        
        consensus_seq, consensus_qual, depth = calculate_quality_weighted_consensus(
            sequences, qualities
        )
        
        # Should prefer A due to higher quality
        self.assertEqual(consensus_seq[0], "A")
        self.assertEqual(depth, [2] * 4)
    
    def test_low_quality_bases(self):
        """Test that low quality bases are ignored."""
        sequences = ["ATCG", "ATCG"]
        qualities = [[30, 30, 30, 30], [10, 10, 10, 10]]  # Second seq has low quality
        
        consensus_seq, consensus_qual, depth = calculate_quality_weighted_consensus(
            sequences, qualities, min_quality=20
        )
        
        # Only first sequence should contribute
        self.assertEqual(consensus_seq, "ATCG")
        # Depth should be 1 at all positions (only first sequence counts)
        self.assertEqual(depth, [1] * 4)
    
    def test_no_coverage(self):
        """Test position with no coverage produces N."""
        sequences = ["AT"]
        qualities = [[5, 5]]  # All low quality
        
        consensus_seq, consensus_qual, depth = calculate_quality_weighted_consensus(
            sequences, qualities, min_quality=20
        )
        
        # Should produce Ns due to low quality
        self.assertEqual(consensus_seq, "NN")
        self.assertEqual(depth, [0, 0])


class TestWellProcessing(unittest.TestCase):
    """Test processing of well FASTQ files."""
    
    def setUp(self):
        """Create temporary directory for test files."""
        self.test_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.test_dir)
    
    def test_process_single_read(self):
        """Test processing a well with one read."""
        # Create test FASTQ file
        fastq_path = os.path.join(self.test_dir, "A1_reads.fastq")
        record = SeqRecord(
            Seq("ATCGATCG"),
            id="test_read",
            description="",
            letter_annotations={'phred_quality': [30] * 8}
        )
        SeqIO.write([record], fastq_path, "fastq")
        
        # Process the well
        data = process_well_fastq(fastq_path, align_sequences=False)
        
        self.assertEqual(data['num_reads'], 1)
        self.assertEqual(data['consensus_seq'], "ATCGATCG")
        self.assertEqual(len(data['consensus_qual']), 8)
    
    def test_process_multiple_reads(self):
        """Test processing a well with multiple reads."""
        # Create test FASTQ file with multiple identical reads
        fastq_path = os.path.join(self.test_dir, "A1_reads.fastq")
        records = [
            SeqRecord(
                Seq("ATCGATCG"),
                id=f"read_{i}",
                description="",
                letter_annotations={'phred_quality': [30] * 8}
            )
            for i in range(3)
        ]
        SeqIO.write(records, fastq_path, "fastq")
        
        # Process the well
        data = process_well_fastq(fastq_path, align_sequences=False)
        
        self.assertEqual(data['num_reads'], 3)
        self.assertEqual(data['consensus_seq'], "ATCGATCG")
    
    def test_max_reads_limit(self):
        """Test that max_reads parameter limits processing."""
        # Create test FASTQ file with many reads
        fastq_path = os.path.join(self.test_dir, "A1_reads.fastq")
        records = [
            SeqRecord(
                Seq("ATCGATCG"),
                id=f"read_{i}",
                description="",
                letter_annotations={'phred_quality': [30] * 8}
            )
            for i in range(10)
        ]
        SeqIO.write(records, fastq_path, "fastq")
        
        # Process with limit
        data = process_well_fastq(fastq_path, max_reads=5, align_sequences=False)
        
        self.assertEqual(data['num_reads'], 5)
    
    def test_empty_fastq(self):
        """Test processing an empty FASTQ file."""
        fastq_path = os.path.join(self.test_dir, "A1_reads.fastq")
        with open(fastq_path, 'w') as f:
            f.write("")
        
        data = process_well_fastq(fastq_path, align_sequences=False)
        
        self.assertEqual(data['num_reads'], 0)
        self.assertEqual(data['consensus_seq'], "")


class TestOutputGeneration(unittest.TestCase):
    """Test output file generation."""
    
    def setUp(self):
        """Create temporary directory for test files."""
        self.test_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.test_dir)
    
    def test_write_aligned_fastq(self):
        """Test writing aligned FASTQ file."""
        # Create well data
        records = [
            SeqRecord(
                Seq("ATCGATCG"),
                id="read1",
                description="",
                letter_annotations={'phred_quality': [30] * 8}
            )
        ]
        well_data = {
            'sequences': ["ATCGATCG"],
            'qualities': [[30] * 8],
            'records': records,
            'aligned_sequences': ["ATCGATCG"],
            'aligned_qualities': [[30] * 8],
            'consensus_seq': "ATCGATCG",
            'consensus_qual': [30] * 8,
            'alignment_depth': [1] * 8,
            'num_reads': 1
        }
        
        output_path = os.path.join(self.test_dir, "A1_aligned.fastq")
        write_aligned_fastq(well_data, output_path, "A1")
        
        # Check output file exists
        self.assertTrue(os.path.exists(output_path))
        
        # Read and verify output
        output_records = list(SeqIO.parse(output_path, "fastq"))
        self.assertEqual(len(output_records), 2)  # Consensus + original
        self.assertEqual(output_records[0].id, "A1_consensus")
        self.assertEqual(str(output_records[0].seq), "ATCGATCG")
    
    def test_write_consensus_csv(self):
        """Test writing consensus CSV file."""
        well_data = {
            'A1': {
                'num_reads': 3,
                'consensus_seq': "ATCGATCG",
                'consensus_qual': [30] * 8,
                'alignment_depth': [3] * 8
            },
            'A2': {
                'num_reads': 5,
                'consensus_seq': "GCTAGCTA",
                'consensus_qual': [35] * 8,
                'alignment_depth': [5] * 8
            }
        }
        
        csv_path = os.path.join(self.test_dir, "consensus.csv")
        write_consensus_csv(well_data, csv_path)
        
        # Check output file exists
        self.assertTrue(os.path.exists(csv_path))
        
        # Read and verify CSV
        df = pd.read_csv(csv_path)
        self.assertEqual(len(df), 2)
        self.assertIn('Well', df.columns)
        self.assertIn('Consensus_Sequence', df.columns)
        self.assertIn('Num_Reads', df.columns)
        
        # Check specific values
        a1_row = df[df['Well'] == 'A1'].iloc[0]
        self.assertEqual(a1_row['Num_Reads'], 3)
        self.assertEqual(a1_row['Consensus_Sequence'], "ATCGATCG")


class TestDirectoryProcessing(unittest.TestCase):
    """Test processing entire demux directories."""
    
    def setUp(self):
        """Create temporary directories for test files."""
        self.test_dir = tempfile.mkdtemp()
        self.output_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up temporary directories."""
        shutil.rmtree(self.test_dir)
        shutil.rmtree(self.output_dir)
    
    def test_process_empty_directory(self):
        """Test processing directory with no FASTQ files."""
        well_data = process_demux_directory(
            self.test_dir,
            self.output_dir,
            align=False
        )
        
        self.assertEqual(len(well_data), 0)
    
    def test_process_directory_with_wells(self):
        """Test processing directory with multiple wells."""
        # Create test FASTQ files for different wells
        for well in ['A1', 'A2', 'B1']:
            fastq_path = os.path.join(self.test_dir, f"{well}_reads.fastq")
            record = SeqRecord(
                Seq("ATCGATCG"),
                id=f"{well}_read",
                description="",
                letter_annotations={'phred_quality': [30] * 8}
            )
            SeqIO.write([record], fastq_path, "fastq")
        
        # Process directory
        well_data = process_demux_directory(
            self.test_dir,
            self.output_dir,
            align=False
        )
        
        self.assertEqual(len(well_data), 3)
        self.assertIn('A1', well_data)
        self.assertIn('A2', well_data)
        self.assertIn('B1', well_data)
        
        # Check output files were created
        for well in ['A1', 'A2', 'B1']:
            output_file = os.path.join(self.output_dir, f"{well}_aligned.fastq")
            self.assertTrue(os.path.exists(output_file))


if __name__ == '__main__':
    unittest.main()
