#!/usr/bin/env python3
"""
Tests for the quality report generation functionality.
"""

import unittest
import tempfile
import os
import shutil
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import generate_quality_report

# Test barcode sequences (extracted for clarity)
TEST_ROW_BARCODE = "AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC"
TEST_COL_BARCODE = "CAAGCAGAAGACGGCATACGAGATATTACTCGGTCTCGTGGGCTCGG"
TEST_READ_SPACER = "N" * 100


class TestQualityReport(unittest.TestCase):
    """Test quality report generation."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.demux_dir = os.path.join(self.temp_dir, "demux_output")
        os.makedirs(self.demux_dir)
        
        # Create test barcode CSV
        self.barcode_csv = os.path.join(self.temp_dir, "test_barcodes.csv")
        self._create_barcode_csv()
        
        # Create test demuxed FASTQ files
        self._create_test_fastq_files()
        
        # Create test barcode_stats.csv
        self._create_stats_csv()

    def tearDown(self):
        """Clean up temporary files."""
        shutil.rmtree(self.temp_dir)

    def _create_barcode_csv(self):
        """Create a test barcode CSV file."""
        content = """Well Position,Sequence Name,Sequence
A1,R oDA373.D501,AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC
A1,F oDA361.D701,CAAGCAGAAGACGGCATACGAGATATTACTCGGTCTCGTGGGCTCGG
B1,R oDA374.D502,AATGATACGGCGACCACCGAGATCTACACATAGAGGCTCGTCGGCAGCGTC
B1,F oDA361.D701,CAAGCAGAAGACGGCATACGAGATATTACTCGGTCTCGTGGGCTCGG
"""
        with open(self.barcode_csv, 'w') as f:
            f.write(content)

    def _create_test_fastq_files(self):
        """Create test FASTQ files for a few wells."""
        wells = ['A1', 'B1']
        for well in wells:
            records = []
            for i in range(5):
                # Construct read with test barcodes and spacer
                seq = TEST_ROW_BARCODE + TEST_READ_SPACER + TEST_COL_BARCODE
                qual = "I" * len(seq)
                record = SeqRecord(
                    Seq(seq),
                    id=f"{well}_read_{i}",
                    description="",
                    letter_annotations={'phred_quality': [40] * len(seq)}
                )
                records.append(record)
            
            fastq_path = os.path.join(self.demux_dir, f"{well}_reads.fastq")
            SeqIO.write(records, fastq_path, "fastq")

    def _create_stats_csv(self):
        """Create a test barcode_stats.csv file."""
        content = """,1,2,3,4,5,6,7,8,9,10,11,12,X
A,5,0,0,0,0,0,0,0,0,0,0,0,0
B,5,0,0,0,0,0,0,0,0,0,0,0,0
C,0,0,0,0,0,0,0,0,0,0,0,0,0
D,0,0,0,0,0,0,0,0,0,0,0,0,0
E,0,0,0,0,0,0,0,0,0,0,0,0,0
F,0,0,0,0,0,0,0,0,0,0,0,0,0
G,0,0,0,0,0,0,0,0,0,0,0,0,0
H,0,0,0,0,0,0,0,0,0,0,0,0,0
X,0,0,0,0,0,0,0,0,0,0,0,0,0
Stats,total,length_ok,mapped,single,adapter_only,no_match,too_short,ambiguous_multiple_cols,ambiguous_multiple_rows,ambiguous_multiple_rows_and_cols,,,
,100,95,10,80,0,5,5,0,0,0,,,
"""
        stats_path = os.path.join(self.demux_dir, "barcode_stats.csv")
        with open(stats_path, 'w') as f:
            f.write(content)

    def test_load_barcode_map(self):
        """Test loading barcode map from CSV."""
        barcodes = generate_quality_report.load_barcode_map(self.barcode_csv)
        self.assertGreater(len(barcodes), 0)
        # Check that we have both row and column barcodes
        types = set(info['type'] for info in barcodes.values())
        self.assertIn('row', types)
        self.assertIn('col', types)

    def test_find_barcode_positions(self):
        """Test finding barcode positions in sequences."""
        barcodes = generate_quality_report.load_barcode_map(self.barcode_csv)
        seq = "AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC" + "N" * 100
        positions = generate_quality_report.find_barcode_positions(seq, barcodes)
        # Should find at least one barcode
        self.assertGreater(len(positions), 0)

    def test_plot_read_length_distribution(self):
        """Test read length distribution plotting."""
        output_file = os.path.join(self.temp_dir, "test_length_dist.png")
        stats = generate_quality_report.plot_read_length_distribution(
            self.demux_dir, output_file
        )
        self.assertTrue(os.path.exists(output_file))
        self.assertIn('total_reads', stats)
        self.assertEqual(stats['total_reads'], 10)  # 5 reads in 2 wells

    def test_plot_msa_layout(self):
        """Test MSA layout visualization."""
        barcodes = generate_quality_report.load_barcode_map(self.barcode_csv)
        output_file = os.path.join(self.temp_dir, "test_msa.png")
        generate_quality_report.plot_msa_layout(
            self.demux_dir, barcodes, output_file, max_reads=10
        )
        self.assertTrue(os.path.exists(output_file))

    def test_plot_barcode_position_heatmap(self):
        """Test barcode position heatmap."""
        barcodes = generate_quality_report.load_barcode_map(self.barcode_csv)
        output_file = os.path.join(self.temp_dir, "test_barcode_pos.png")
        generate_quality_report.plot_barcode_position_heatmap(
            self.demux_dir, barcodes, output_file
        )
        self.assertTrue(os.path.exists(output_file))

    def test_plot_barcode_presence_summary(self):
        """Test barcode presence summary."""
        stats_csv = os.path.join(self.demux_dir, "barcode_stats.csv")
        output_file = os.path.join(self.temp_dir, "test_barcode_summary.png")
        generate_quality_report.plot_barcode_presence_summary(stats_csv, output_file)
        self.assertTrue(os.path.exists(output_file))

    def test_generate_html_report(self):
        """Test HTML report generation."""
        report_dir = os.path.join(self.temp_dir, "report")
        os.makedirs(report_dir)
        stats = {
            'total_reads': 10,
            'mean_length': 200,
            'median_length': 180,
            'min_length': 100,
            'max_length': 300,
            'std_length': 50
        }
        html_file = generate_quality_report.generate_html_report(
            self.demux_dir, report_dir, stats
        )
        self.assertTrue(os.path.exists(html_file))
        
        # Check HTML content
        with open(html_file, 'r') as f:
            content = f.read()
            self.assertIn('NanoDemux Quality Report', content)
            self.assertIn('TOTAL READS', content)
            self.assertIn('10', content)


if __name__ == '__main__':
    unittest.main()
