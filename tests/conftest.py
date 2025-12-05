#!/usr/bin/env python3
"""
Configuration and fixtures for pytest.
"""

import pytest
import tempfile
import shutil
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd


@pytest.fixture
def temp_directory():
    """Provide a temporary directory that's cleaned up after the test."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)


@pytest.fixture
def sample_barcodes_csv(temp_directory):
    """Create a sample barcodes CSV file."""
    barcodes = [
        ("A1", "R oDA373.D501", "AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC"),
        ("A1", "F oDA361.D701", "CAAGCAGAAGACGGCATACGAGATATTACTCGGTCTCGTGGGCTCGG"),
        ("A2", "R oDA373.D501", "AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC"),
        ("A2", "F oDA362.D702", "CAAGCAGAAGACGGCATACGAGATTCCGGAGAGTCTCGTGGGCTCGG"),
        ("B1", "R oDA374.D502", "AATGATACGGCGACCACCGAGATCTACACATAGAGGCTCGTCGGCAGCGTC"),
        ("B1", "F oDA361.D701", "CAAGCAGAAGACGGCATACGAGATATTACTCGGTCTCGTGGGCTCGG"),
    ]
    
    csv_file = os.path.join(temp_directory, "test_barcodes.csv")
    df = pd.DataFrame(barcodes, columns=["Well Position", "Sequence Name", "Sequence"])
    df.to_csv(csv_file, index=False)
    
    return csv_file


@pytest.fixture
def sample_fastq(temp_directory):
    """Create a sample FASTQ file for testing."""
    records = []
    
    # Read with both barcodes
    row_primer = "AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC"
    col_primer = "CAAGCAGAAGACGGCATACGAGATATTACTCGGTCTCGTGGGCTCGG"
    insert = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 5
    
    test_cases = [
        (row_primer + insert + col_primer, "read_both"),
        (row_primer + insert, "read_row"),
        ("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 5, "read_none"),
        ("ATCG" * 5, "read_short"),
    ]
    
    for seq_str, read_id in test_cases:
        rec = SeqRecord(Seq(seq_str), id=read_id)
        rec.letter_annotations["phred_quality"] = [40] * len(seq_str)
        records.append(rec)
    
    fastq_file = os.path.join(temp_directory, "test.fastq")
    SeqIO.write(records, fastq_file, "fastq")
    
    return fastq_file


@pytest.fixture
def output_directory(temp_directory):
    """Provide an output directory."""
    return os.path.join(temp_directory, "output")
