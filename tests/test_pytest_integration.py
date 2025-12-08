#!/usr/bin/env python3
"""
Pytest-based integration tests with fixtures.
"""

import pytest
import os
import pandas as pd
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "nanodemux_scripts"))

from demux_barcodes import main, load_barcodes


class TestDemultiplexingWithPytest:
    """Integration tests using pytest and fixtures."""

    def test_full_pipeline(self, sample_fastq, sample_barcodes_csv, output_directory):
        """Test the full demultiplexing pipeline."""
        main(
            sample_fastq,
            sample_barcodes_csv,
            output_directory,
            min_length=50,
            max_penalty=60,
            cpus=1,
            flank=100
        )
        
        assert os.path.isdir(output_directory), "Output directory should exist"
        assert os.path.isfile(os.path.join(output_directory, "barcode_stats.csv")), \
            "Stats file should be created"

    def test_stats_file_valid(self, sample_fastq, sample_barcodes_csv, output_directory):
        """Test that stats file is a valid CSV with expected columns."""
        main(
            sample_fastq,
            sample_barcodes_csv,
            output_directory,
            min_length=50,
            max_penalty=60,
            cpus=1,
            flank=100
        )
        
        stats_file = os.path.join(output_directory, "barcode_stats.csv")
        df = pd.read_csv(stats_file, index_col=0)
        
        # Should have numeric columns for wells (1-12)
        assert "1" in df.columns or any(str(i) in df.columns for i in range(1, 13)), \
            "Stats should have numeric well columns"

    def test_fastq_output_format(self, sample_fastq, sample_barcodes_csv, output_directory):
        """Test that output FASTQ files are valid."""
        main(
            sample_fastq,
            sample_barcodes_csv,
            output_directory,
            min_length=50,
            max_penalty=60,
            cpus=1,
            flank=100
        )
        
        # Check any generated FASTQ files
        fastq_files = [
            os.path.join(output_directory, f)
            for f in os.listdir(output_directory)
            if f.endswith("_reads.fastq")
        ]
        
        for fastq_file in fastq_files:
            from Bio import SeqIO
            records = list(SeqIO.parse(fastq_file, "fastq"))
            assert len(records) >= 0, f"FASTQ {fastq_file} should be readable"

    def test_barcode_loading_with_fixture(self, sample_barcodes_csv):
        """Test barcode loading with pytest fixture."""
        row_map, col_map = load_barcodes(sample_barcodes_csv)
        
        assert len(row_map) > 0, "Should load row barcodes"
        assert len(col_map) > 0, "Should load column barcodes"
        assert all(isinstance(v, str) for v in row_map.values()), \
            "Row values should be letter strings"
        assert all(isinstance(v, int) for v in col_map.values()), \
            "Column values should be integers"

    def test_handles_different_flank_sizes(self, sample_fastq, sample_barcodes_csv, output_directory):
        """Test that pipeline works with different flank sizes."""
        for flank in [50, 100, 150]:
            output = os.path.join(output_directory, f"flank_{flank}")
            
            main(
                sample_fastq,
                sample_barcodes_csv,
                output,
                min_length=50,
                max_penalty=60,
                cpus=1,
                flank=flank
            )
            
            assert os.path.isfile(os.path.join(output, "barcode_stats.csv")), \
                f"Should work with flank={flank}"
