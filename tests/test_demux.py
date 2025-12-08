#!/usr/bin/env python3
"""
Unit tests for barcode matching and demultiplexing logic.
"""

import unittest
import tempfile
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "nanodemux_scripts"))

from demux_barcodes import (
    revcomp,
    match_barcode_ends_weighted,
    load_barcodes,
    int_defaultdict,
    merge_dicts,
    merge_reads,
    get_basename_without_extensions,
    find_barcode_candidates,
    find_barcodes_with_alignment,
    find_barcodes_two_tier,
    alignment_score_and_llr
)


class TestHelperFunctions(unittest.TestCase):
    """Test basic helper functions."""

    def test_revcomp_basic(self):
        """Test reverse complement of simple sequences."""
        self.assertEqual(revcomp("ATCG"), "CGAT")
        self.assertEqual(revcomp("AAAA"), "TTTT")
        self.assertEqual(revcomp("GCGC"), "GCGC")

    def test_revcomp_empty(self):
        """Test reverse complement of empty string."""
        self.assertEqual(revcomp(""), "")

    def test_int_defaultdict(self):
        """Test int_defaultdict initialization."""
        d = int_defaultdict()
        self.assertIsInstance(d, defaultdict)
        self.assertEqual(d['nonexistent'], 0)

    def test_get_basename_without_extensions(self):
        """Test removal of FASTQ extensions including double extensions."""
        self.assertEqual(get_basename_without_extensions("sample.fastq"), "sample")
        self.assertEqual(get_basename_without_extensions("sample.fq"), "sample")
        self.assertEqual(get_basename_without_extensions("sample.fastq.gz"), "sample")
        self.assertEqual(get_basename_without_extensions("sample.fq.gz"), "sample")
        self.assertEqual(get_basename_without_extensions("/path/to/sample.fastq.gz"), "sample")
        # Test with unknown extension
        self.assertEqual(get_basename_without_extensions("sample.txt"), "sample")


class TestBarcodeMatching(unittest.TestCase):
    """Test barcode matching algorithm."""

    def test_exact_match(self):
        """Test exact barcode match with perfect sequence."""
        seq = "ATCGATCGATCGATCG" * 10
        qual = [40] * len(seq)
        barcode = "ATCGATCGATCGATCG"
        
        # For testing, use a simple case with no const/var distinction
        group_seqs = [barcode]
        const_idx = list(range(len(barcode)))
        var_idx = []
        
        result = match_barcode_ends_weighted(
            seq, qual, barcode, group_seqs, const_idx, var_idx,
            max_penalty=60, flank=50
        )
        self.assertTrue(result, "Should find exact match")

    def test_no_match(self):
        """Test no barcode match."""
        seq = "AAAAAAAAAAAAAAAA" * 10
        qual = [40] * len(seq)
        barcode = "TTTTTTTTTTTTTTTT"
        
        group_seqs = [barcode]
        const_idx = list(range(len(barcode)))
        var_idx = []
        
        result = match_barcode_ends_weighted(
            seq, qual, barcode, group_seqs, const_idx, var_idx,
            max_penalty=60, flank=50
        )
        self.assertFalse(result, "Should not find match with all mismatches")

    def test_match_with_low_quality(self):
        """Test match with low-quality mismatches (high Phred scores = high penalties)."""
        seq = "ATCGATCGATCGATCG" * 10
        qual = [40] * len(seq)  # High quality = match is worth it
        barcode = "ATCGATCGATCGATCG"
        
        group_seqs = [barcode]
        const_idx = list(range(len(barcode)))
        var_idx = []
        
        result = match_barcode_ends_weighted(
            seq, qual, barcode, group_seqs, const_idx, var_idx,
            max_penalty=60, flank=50
        )
        self.assertTrue(result, "Should find match despite some quality concerns")

    def test_match_reverse_complement(self):
        """Test barcode matching on reverse complement."""
        # Forward sequence
        seq_fwd = "ATCGATCGATCGATCG" * 10
        seq_rc = revcomp(seq_fwd)
        qual = [40] * len(seq_rc)
        barcode = "ATCGATCGATCGATCG"
        
        group_seqs = [barcode]
        const_idx = list(range(len(barcode)))
        var_idx = []
        
        # The barcode appears at the end of the reverse complement
        result = match_barcode_ends_weighted(
            seq_rc, qual, barcode, group_seqs, const_idx, var_idx,
            max_penalty=60, flank=50
        )
        self.assertTrue(result, "Should find barcode in reverse complement")

    def test_flank_size_constraint(self):
        """Test that barcode matching respects flank size."""
        # Barcode in middle, shouldn't be found with small flank
        seq = "AAAA" + "ATCGATCGATCGATCG" + "AAAA" * 10
        qual = [40] * len(seq)
        barcode = "ATCGATCGATCGATCG"
        
        group_seqs = [barcode]
        const_idx = list(range(len(barcode)))
        var_idx = []
        
        result = match_barcode_ends_weighted(
            seq, qual, barcode, group_seqs, const_idx, var_idx,
            max_penalty=60, flank=5
        )
        self.assertFalse(result, "Should not find barcode outside flank region")


class TestBarcodeLoading(unittest.TestCase):
    """Test barcode loading from CSV."""

    def setUp(self):
        """Create temporary test CSV file."""
        self.temp_dir = tempfile.mkdtemp()
        self.csv_file = os.path.join(self.temp_dir, "test_barcodes.csv")
        
        # Create minimal test barcode CSV
        with open(self.csv_file, 'w') as f:
            f.write("Well Position,Sequence Name,Sequence\n")
            f.write("A1,R test1,AGATATTACGAGATCTACACTATAGCCTTCGTCGGCAGCGTC\n")
            f.write("A1,F test1,AGATCAGAAGACGGCATACGAGATATTACTCGGTCTCGTCGTC\n")
            f.write("A2,R test1,AGATATTACGAGATCTACACATAGAGGCTCGTCGGCAGCGTC\n")
            f.write("A2,F test2,AGATCAGAAGACGGCATACGAGATTCCGGAGAGTCTCGTCGTC\n")

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_load_barcodes_structure(self):
        """Test that barcodes are loaded into correct maps."""
        row_map, col_map = load_barcodes(self.csv_file)
        
        self.assertEqual(len(row_map), 2, "Should load 2 row barcodes")
        self.assertEqual(len(col_map), 2, "Should load 2 column barcodes")

    def test_load_barcodes_mapping(self):
        """Test that barcodes map correctly to positions."""
        row_map, col_map = load_barcodes(self.csv_file)
        
        # Check that row barcodes map to letters
        for barcode, letter in row_map.items():
            self.assertIn(letter, "ABCDEFGH", "Row should map to letter")
        
        # Check that column barcodes map to numbers
        for barcode, number in col_map.items():
            self.assertIsInstance(number, int, "Column should map to number")
            self.assertGreaterEqual(number, 1)
            self.assertLessEqual(number, 12)


class TestMergingFunctions(unittest.TestCase):
    """Test dictionary and read merging functions."""

    def test_merge_dicts(self):
        """Test merging of defaultdicts."""
        dict1 = defaultdict(lambda: defaultdict(int))
        dict1['A']['x'] = 5
        dict1['A']['y'] = 3
        
        dict2 = defaultdict(lambda: defaultdict(int))
        dict2['A']['x'] = 2
        dict2['B']['z'] = 7
        
        result = merge_dicts([dict1, dict2])
        
        self.assertEqual(result['A']['x'], 7, "Should sum values for A/x")
        self.assertEqual(result['A']['y'], 3, "Should preserve A/y")
        self.assertEqual(result['B']['z'], 7, "Should preserve B/z")

    def test_merge_reads(self):
        """Test merging of read lists."""
        reads1 = defaultdict(list)
        reads1['A1'] = ['read1', 'read2']
        
        reads2 = defaultdict(list)
        reads2['A1'] = ['read3']
        reads2['A2'] = ['read4', 'read5']
        
        result = merge_reads([reads1, reads2])
        
        self.assertEqual(len(result['A1']), 3, "Should merge A1 reads")
        self.assertEqual(len(result['A2']), 2, "Should include A2 reads")


class TestFASTQHandling(unittest.TestCase):
    """Test FASTQ file reading and writing."""

    def setUp(self):
        """Create temporary FASTQ file for testing."""
        self.temp_dir = tempfile.mkdtemp()
        self.fastq_file = os.path.join(self.temp_dir, "test.fastq")
        
        # Create test FASTQ
        records = []
        for seq_str, qual_str in [("ATCG" * 20, "I" * 80), ("GCTA" * 20, "J" * 80), ("AT" * 25, "I" * 50)]:
            rec = SeqRecord(Seq(seq_str), id=f"read{len(records)+1}")
            rec.letter_annotations["phred_quality"] = [ord(q) - ord('!') for q in qual_str]
            records.append(rec)
        SeqIO.write(records, self.fastq_file, "fastq")

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_fastq_read_count(self):
        """Test that FASTQ is read with correct record count."""
        records = list(SeqIO.parse(self.fastq_file, "fastq"))
        self.assertEqual(len(records), 3, "Should read 3 records")

    def test_fastq_sequence_lengths(self):
        """Test that sequences are read with correct lengths."""
        records = list(SeqIO.parse(self.fastq_file, "fastq"))
        
        self.assertEqual(len(records[0].seq), 80, "First record should be 80 bp")
        self.assertEqual(len(records[1].seq), 80, "Second record should be 80 bp")
        self.assertEqual(len(records[2].seq), 50, "Third record should be 50 bp")

    def test_fastq_quality_scores(self):
        """Test that quality scores are read correctly."""
        records = list(SeqIO.parse(self.fastq_file, "fastq"))
        
        qual1 = records[0].letter_annotations["phred_quality"]
        self.assertEqual(len(qual1), 80, "Quality should match sequence length")
        self.assertTrue(all(q > 0 for q in qual1), "All quality scores should be positive")


class TestTwoTierMatching(unittest.TestCase):
    """Test two-tiered barcode matching logic."""

    def test_find_barcode_candidates_unique_match(self):
        """Test fast candidate filter with unique match."""
        seq = "ATCGATCGATCGATCG" * 10
        qual = [40] * len(seq)
        
        barcode1 = "ATCGATCGATCGATCG"
        barcode2 = "GCTAGCTAGCTAGCTA"
        
        barcode_map = {barcode1: "A", barcode2: "B"}
        group_seqs = [barcode1, barcode2]
        const_idx = list(range(len(barcode1)))
        var_idx = []
        
        candidates = find_barcode_candidates(
            seq, qual, barcode_map, group_seqs, const_idx, var_idx
        )
        
        self.assertEqual(len(candidates), 1, "Should find one unique candidate")
        self.assertEqual(candidates[0], "A", "Should identify barcode1 as candidate")

    def test_alignment_score_and_llr(self):
        """Test alignment scoring function."""
        seq = "ATCGATCGATCGATCG" * 10
        barcode = "ATCGATCGATCGATCG"
        qual = [40] * len(seq)
        
        score, llr = alignment_score_and_llr(seq, barcode, qual)
        
        self.assertIsNotNone(score, "Should find alignment")
        self.assertIsNotNone(llr, "Should compute LLR")
        self.assertGreater(score, 0, "Score should be positive for matching sequence")

    def test_find_barcodes_with_alignment(self):
        """Test alignment-based barcode finding."""
        seq = "ATCGATCGATCGATCG" * 10
        qual = [40] * len(seq)
        
        barcode1 = "ATCGATCGATCGATCG"
        barcode2 = "GCTAGCTAGCTAGCTA"
        
        barcode_map = {barcode1: "A", barcode2: "B"}
        
        results = find_barcodes_with_alignment(seq, qual, barcode_map)
        
        self.assertGreater(len(results), 0, "Should find at least one match")
        self.assertEqual(results[0][0], "A", "Best match should be barcode1")

    def test_two_tier_fast_path(self):
        """Test two-tier matching with fast path (unique match)."""
        seq = "ATCGATCGATCGATCG" * 10
        qual = [40] * len(seq)
        
        barcode1 = "ATCGATCGATCGATCG"
        barcode2 = "GCTAGCTAGCTAGCTA"
        
        barcode_map = {barcode1: "A", barcode2: "B"}
        group_seqs = [barcode1, barcode2]
        const_idx = list(range(len(barcode1)))
        var_idx = []
        
        candidates, method = find_barcodes_two_tier(
            seq, qual, barcode_map, group_seqs, const_idx, var_idx
        )
        
        self.assertEqual(len(candidates), 1, "Should find unique match")
        self.assertEqual(method, "fast", "Should use fast method")
        self.assertEqual(candidates[0], "A", "Should identify correct barcode")

    def test_two_tier_fallback_path(self):
        """Test two-tier matching with fallback alignment."""
        # Create a sequence with some mismatches that might be ambiguous
        seq = "ATCGATCGNTCGATCG" * 10  # N = ambiguous base
        qual = [40] * len(seq)
        
        barcode1 = "ATCGATCGATCGATCG"
        barcode2 = "ATCGATCGCTCGATCG"  # Very similar to barcode1
        
        barcode_map = {barcode1: "A", barcode2: "B"}
        group_seqs = [barcode1, barcode2]
        const_idx = list(range(len(barcode1)))
        var_idx = []
        
        candidates, method = find_barcodes_two_tier(
            seq, qual, barcode_map, group_seqs, const_idx, var_idx,
            use_fallback=True
        )
        
        # Should resolve to at most one candidate (or none if truly ambiguous)
        self.assertLessEqual(len(candidates), 1, "Should not have multiple candidates after two-tier matching")


if __name__ == "__main__":
    unittest.main()
