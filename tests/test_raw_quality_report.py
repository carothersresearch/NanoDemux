#!/usr/bin/env python3
"""
Tests for the raw quality report generation functionality.
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

import generate_raw_quality_report


class TestRawQualityReport(unittest.TestCase):
    """Test raw quality report generation."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.test_fastq = os.path.join(self.temp_dir, "test_raw.fastq")
        self._create_test_fastq()

    def tearDown(self):
        """Clean up temporary files."""
        shutil.rmtree(self.temp_dir)

    def _create_test_fastq(self):
        """Create a test FASTQ file with known properties."""
        records = []
        
        # Create reads with varying lengths and quality
        for i in range(50):
            length = 100 + i * 10  # Varying lengths from 100 to 590
            seq = "ATCG" * (length // 4)  # Simple repeating pattern
            
            # Pad sequence if needed to match exact length
            if len(seq) < length:
                seq += "ATCG"[:length - len(seq)]
            seq = seq[:length]  # Trim to exact length
            
            # Vary quality scores
            qual = [30 + (i % 10)] * len(seq)  # Quality scores from 30 to 39
            
            record = SeqRecord(
                Seq(seq),
                id=f"test_read_{i}",
                description="",
                letter_annotations={'phred_quality': qual}
            )
            records.append(record)
        
        SeqIO.write(records, self.test_fastq, "fastq")

    def test_calculate_n50(self):
        """Test N50 calculation."""
        lengths = [100, 200, 300, 400, 500]
        n50 = generate_raw_quality_report.calculate_n50(lengths)
        # Total length = 1500, half = 750
        # Cumulative: 500, 500+400=900 >= 750, so N50 = 400
        self.assertEqual(n50, 400)

    def test_analyze_fastq(self):
        """Test FASTQ analysis."""
        data = generate_raw_quality_report.analyze_fastq(self.test_fastq)
        
        # Check that data structure is correct
        self.assertIn('stats', data)
        self.assertIn('lengths', data)
        self.assertIn('qualities', data)
        self.assertIn('base_counts', data)
        self.assertIn('quality_by_position', data)
        
        # Check statistics
        stats = data['stats']
        self.assertEqual(stats['total_reads'], 50)
        self.assertGreater(stats['mean_length'], 0)
        self.assertGreater(stats['n50'], 0)
        self.assertGreater(stats['mean_quality'], 0)
        
        # Check that we have base counts
        base_counts = data['base_counts']
        self.assertGreater(base_counts['A'], 0)
        self.assertGreater(base_counts['T'], 0)
        self.assertGreater(base_counts['G'], 0)
        self.assertGreater(base_counts['C'], 0)

    def test_analyze_fastq_with_max_reads(self):
        """Test FASTQ analysis with max_reads limit."""
        data = generate_raw_quality_report.analyze_fastq(self.test_fastq, max_reads=10)
        
        stats = data['stats']
        self.assertEqual(stats['total_reads'], 10)

    def test_plot_read_length_distribution(self):
        """Test read length distribution plotting."""
        data = generate_raw_quality_report.analyze_fastq(self.test_fastq)
        output_file = os.path.join(self.temp_dir, "test_length_dist.png")
        
        generate_raw_quality_report.plot_read_length_distribution(data, output_file)
        self.assertTrue(os.path.exists(output_file))
        self.assertGreater(os.path.getsize(output_file), 0)

    def test_plot_quality_distribution(self):
        """Test quality distribution plotting."""
        data = generate_raw_quality_report.analyze_fastq(self.test_fastq)
        output_file = os.path.join(self.temp_dir, "test_quality_dist.png")
        
        generate_raw_quality_report.plot_quality_distribution(data, output_file)
        self.assertTrue(os.path.exists(output_file))
        self.assertGreater(os.path.getsize(output_file), 0)

    def test_plot_base_composition(self):
        """Test base composition plotting."""
        data = generate_raw_quality_report.analyze_fastq(self.test_fastq)
        output_file = os.path.join(self.temp_dir, "test_base_comp.png")
        
        generate_raw_quality_report.plot_base_composition(data, output_file)
        self.assertTrue(os.path.exists(output_file))
        self.assertGreater(os.path.getsize(output_file), 0)

    def test_generate_html_report(self):
        """Test HTML report generation."""
        data = generate_raw_quality_report.analyze_fastq(self.test_fastq)
        report_dir = os.path.join(self.temp_dir, "report")
        os.makedirs(report_dir)
        
        html_file = generate_raw_quality_report.generate_html_report(
            self.test_fastq, report_dir, data
        )
        
        self.assertTrue(os.path.exists(html_file))
        
        # Check HTML content
        with open(html_file, 'r') as f:
            content = f.read()
            self.assertIn('NanoDemux Raw Data Quality Report', content)
            self.assertIn('TOTAL READS', content)
            self.assertIn('50', content)  # Should have 50 reads
            self.assertIn('GC CONTENT', content)

    def test_empty_fastq(self):
        """Test handling of empty FASTQ file."""
        empty_fastq = os.path.join(self.temp_dir, "empty.fastq")
        with open(empty_fastq, 'w') as f:
            f.write("")
        
        data = generate_raw_quality_report.analyze_fastq(empty_fastq)
        stats = data['stats']
        
        # Should handle empty file gracefully
        self.assertEqual(stats['total_reads'], 0)
        self.assertEqual(stats['total_bases'], 0)


if __name__ == '__main__':
    unittest.main()
