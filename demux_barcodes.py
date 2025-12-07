#!/usr/bin/env python3

"""
Demultiplex nanopore FASTQ by row/column barcodes with:
 - barcode analysis (constant vs variable positions)
 - quality-weighted matching on constant positions
 - variable-region disambiguation tolerant of low-quality bases (but requires uniqueness)
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from collections import defaultdict
import os
from multiprocessing import Pool
import glob
import sys


import itertools
import parasail
import numpy as np

# ---------------------------------
# Helpers
# ---------------------------------

def revcomp(seq):
    return str(Seq(seq).reverse_complement())

def int_defaultdict():
    return defaultdict(int)

def hamming(a, b):
    """Hamming distance of equal-length strings."""
    if len(a) != len(b):
        raise ValueError("Hamming requires equal lengths")
    return sum(x != y for x, y in zip(a, b))

def get_basename_without_extensions(filename):
    """
    Remove common FASTQ extensions from filename, handling double extensions.
    Examples:
        sample.fastq -> sample
        sample.fq -> sample
        sample.fastq.gz -> sample
        sample.fq.gz -> sample
    """
    basename = os.path.basename(filename)
    # Remove known FASTQ extensions
    for ext in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']:
        if basename.endswith(ext):
            return basename[:-len(ext)]
    # Fallback to standard splitext if no known extension found
    return os.path.splitext(basename)[0]

# ---------------------------------
# Barcode analysis
# ---------------------------------

def analyze_barcode_set(sequences):
    """
    Given list of equal-length barcode sequences (strings), return metadata:
      - seqs: list of sequences
      - length: k
      - const_idx: indices that are constant across sequences
      - var_idx: indices that vary across sequences
      - pairwise: list of (i, j, hamming)
    Also prints a short summary (min/mean/max Hamming distance and variable positions).
    """
    seqs = list(sequences)
    if not seqs:
        return {'seqs': [], 'length': 0, 'const_idx': [], 'var_idx': [], 'pairwise': []}
    k = len(seqs[0])
    for s in seqs:
        if len(s) != k:
            raise ValueError("All barcodes in a set must be same length for analysis")

    # positions variable vs constant
    const_idx = []
    var_idx = []
    for i in range(k):
        bases = {s[i] for s in seqs}
        if len(bases) == 1:
            const_idx.append(i)
        else:
            var_idx.append(i)

    # pairwise distances
    pairwise = []
    distances = []
    for (i,s1),(j,s2) in itertools.combinations(list(enumerate(seqs)), 2):
        d = hamming(s1, s2)
        pairwise.append((i, j, d))
        distances.append(d)

    if distances:
        mn = min(distances); mx = max(distances); avg = sum(distances)/len(distances)
    else:
        mn = mx = avg = 0

    print(f"[barcode_analysis] {len(seqs)} sequences of length {k}.")
    print(f"  variable positions: {var_idx}")
    print(f"  constant positions: {const_idx}")
    print(f"  pairwise Hamming distances: min={mn}, mean={avg:.2f}, max={mx}")

    return {
        'seqs': seqs,
        'length': k,
        'const_idx': const_idx,
        'var_idx': var_idx,
        'pairwise': pairwise
    }

# ---------------------------------
# Fast candidate filter (Method 1)
# ---------------------------------

def find_barcode_candidates(seq, qual, barcode_map, group_seqs, const_idx, var_idx,
                            max_penalty=60, flank=100, var_q=10):
    """
    Fast candidate filter: find all barcodes that pass the quality-weighted constant check
    and variable disambiguation.
    
    Parameters:
        seq, qual: sequence string and list of Phred ints
        barcode_map: dict mapping barcode sequences to identifiers (e.g., row letters or col numbers)
        group_seqs: list of all barcode sequences in the group
        const_idx, var_idx: constant and variable position indices
        max_penalty: allowed sum of Phred scores of mismatches at constant positions
        flank: bases from each end to search
        var_q: Phred cutoff for confident variable-base match
    
    Returns:
        list of identifiers for barcodes that match (empty if none match)
    """
    candidates = []
    for barcode_seq, identifier in barcode_map.items():
        if match_barcode_ends_weighted(seq, qual, barcode_seq, group_seqs,
                                      const_idx, var_idx, max_penalty, flank, var_q):
            candidates.append(identifier)
    return candidates

# ---------------------------------
# Constants for alignment scoring
# ---------------------------------

# Parasail alignment matrix (match=2, mismatch=-1)
# Cached to avoid recreating on every alignment call
_ALIGNMENT_MATRIX = None

def _get_alignment_matrix():
    """Get or create the cached alignment matrix."""
    global _ALIGNMENT_MATRIX
    if _ALIGNMENT_MATRIX is None:
        _ALIGNMENT_MATRIX = parasail.matrix_create("ACGT", 2, -1)
    return _ALIGNMENT_MATRIX

# LLR calculation constants
RANDOM_BASE_PROBABILITY = 0.25  # Probability of random base match (1/4 bases)
MATCH_SCORE = 2  # Match score used in alignment

# Alignment uniqueness threshold
MIN_SCORE_DIFFERENCE = 5  # Minimum score difference for unique alignment

# ---------------------------------
# Fallback alignment (Method 2)
# ---------------------------------

def alignment_score_and_llr(seq, barcode_seq, qual, flank=100):
    """
    Compute local alignment score using parasail and calculate a log-likelihood ratio.
    
    Parameters:
        seq: read sequence (string)
        barcode_seq: barcode sequence (string)
        qual: quality scores (list of Phred ints)
        flank: bases from each end to search
    
    Returns:
        tuple: (alignment_score, llr) or (None, None) if no alignment found
    """
    # Extract both ends for alignment
    k = len(barcode_seq)
    ends = [seq[:flank], seq[-flank:]]
    
    best_score = -float('inf')
    best_llr = -float('inf')
    
    # Use Smith-Waterman local alignment with cached matrix
    # Match=2, Mismatch=-1, Gap open=2, Gap extend=1 (typical nanopore-friendly scoring)
    matrix = _get_alignment_matrix()
    
    for end_seq in ends:
        if len(end_seq) < k:
            continue
        
        # Perform Smith-Waterman alignment
        result = parasail.sw_trace_scan_16(barcode_seq, end_seq, 2, 1, matrix)
        
        if result.score > best_score:
            best_score = result.score
            
            # Calculate LLR based on alignment quality
            # LLR = log(P(alignment|correct) / P(alignment|incorrect))
            # Simplified: use alignment score normalized by barcode length
            # and penalized by expected random match
            expected_random_score = k * RANDOM_BASE_PROBABILITY * MATCH_SCORE
            best_llr = (result.score - expected_random_score) / k
    
    return (best_score, best_llr) if best_score > 0 else (None, None)

def find_barcodes_with_alignment(seq, qual, barcode_map, min_score=20, min_llr=0.5, flank=100):
    """
    Find barcodes using alignment-based approach for ambiguous cases.
    
    Parameters:
        seq, qual: sequence and quality scores
        barcode_map: dict mapping barcode sequences to identifiers
        min_score: minimum alignment score to accept
        min_llr: minimum log-likelihood ratio to accept
        flank: bases from each end to search
    
    Returns:
        list of tuples: [(identifier, score, llr), ...] sorted by score descending
    """
    results = []
    for barcode_seq, identifier in barcode_map.items():
        score, llr = alignment_score_and_llr(seq, barcode_seq, qual, flank)
        if score is not None and score >= min_score and llr >= min_llr:
            results.append((identifier, score, llr))
    
    # Sort by score (descending), then by llr (descending)
    results.sort(key=lambda x: (x[1], x[2]), reverse=True)
    return results

# ---------------------------------
# Two-tiered matching wrapper
# ---------------------------------

def find_barcodes_two_tier(seq, qual, barcode_map, group_seqs, const_idx, var_idx,
                          max_penalty=60, flank=100, var_q=10,
                          min_alignment_score=20, min_llr=0.5,
                          use_fallback=True):
    """
    Two-tiered barcode matching:
    1. Try fast candidate filter first
    2. If ambiguous (0 matches or >1 matches), optionally use alignment fallback
    
    Parameters:
        seq, qual: sequence and quality scores
        barcode_map: dict mapping barcode sequences to identifiers
        group_seqs: list of all barcode sequences in the group
        const_idx, var_idx: constant and variable position indices
        max_penalty: allowed sum of Phred scores of mismatches at constant positions
        flank: bases from each end to search
        var_q: Phred cutoff for confident variable-base match
        min_alignment_score: minimum alignment score for fallback
        min_llr: minimum log-likelihood ratio for fallback
        use_fallback: whether to use alignment fallback for ambiguous cases
    
    Returns:
        tuple: (list of identifiers, method_used)
               method_used: 'fast' or 'alignment' or 'none'
    """
    # Method 1: Fast candidate filter
    candidates = find_barcode_candidates(seq, qual, barcode_map, group_seqs,
                                        const_idx, var_idx, max_penalty, flank, var_q)
    
    # If unique match found, return immediately (fast path)
    if len(candidates) == 1:
        return candidates, 'fast'
    
    # If no ambiguity and fast filter worked, return
    if len(candidates) == 0 and not use_fallback:
        return [], 'fast'
    
    # Method 2: Fallback alignment for ambiguous cases (0 or >1 candidates)
    if use_fallback and (len(candidates) == 0 or len(candidates) > 1):
        alignment_results = find_barcodes_with_alignment(
            seq, qual, barcode_map, min_alignment_score, min_llr, flank
        )
        
        if len(alignment_results) > 0:
            # Take best alignment (first in sorted list)
            best_match = alignment_results[0]
            # Only return if significantly better than second-best
            if len(alignment_results) == 1:
                return [best_match[0]], 'alignment'
            else:
                # Check if best is significantly better than second
                score_diff = best_match[1] - alignment_results[1][1]
                if score_diff > MIN_SCORE_DIFFERENCE:  # Require significant score difference for uniqueness
                    return [best_match[0]], 'alignment'
        
        # Alignment didn't help disambiguate
        return [], 'alignment'
    
    # Multiple candidates from fast filter, no fallback
    return candidates, 'fast'

# ---------------------------------
# Improved matching function (kept for backwards compatibility)
# ---------------------------------

def match_barcode_ends_weighted(seq, qual, barcode_seq,
                                group_seqs, const_idx, var_idx,
                                max_penalty=60, flank=20, var_q=10):
    """
    Match 'barcode_seq' against either end of 'seq' (searching windows inside 'flank'),
    using:
      - quality-weighted penalties for mismatches **only at constant positions** (sum of Phred)
      - variable positions used to disambiguate candidate barcode vs other barcodes in the same group
        (we count how many variable positions match for each candidate barcode,
         and require the target barcode to be strictly better than any other barcode)

    Parameters:
      seq, qual: sequence string and list of Phred ints (same length)
      barcode_seq: candidate barcode sequence (length k)
      group_seqs: list of all barcode sequences in the same group (including barcode_seq)
      const_idx: indices (0..k-1) that are constant across group_seqs
      var_idx: indices that vary across group_seqs
      max_penalty: allowed sum of Phred scores of mismatches at constant positions
      flank: how many bases from each end to search
      var_q: Phred cutoff to count a variable-base match as "confident"

    Returns:
      True if at least one window matches under these rules, else False.
    """
    k = len(barcode_seq)
    ends = [
        (seq[:flank], qual[:flank]),
        (seq[-flank:], qual[-flank:])
    ]

    # Precompute dict for group sequences -> list of chars (for fast index)
    group_arrays = {g: list(g) for g in group_seqs}

    for region_seq, region_qual in ends:
        L = len(region_seq)
        if L < k:
            continue
        for offset in range(L - k + 1):
            window = region_seq[offset:offset + k]
            window_qual = region_qual[offset:offset + k]

            # constant-position penalty: sum Phred where mismatch at const positions
            penalty_const = 0
            for idx in const_idx:
                if window[idx] != barcode_seq[idx]:
                    penalty_const += window_qual[idx]

            if penalty_const > max_penalty:
                # too many constant-position mismatches; reject this window fast
                continue

            # variable-position support counts for every barcode in group
            var_support = {}
            var_confident_support = {}
            for g in group_seqs:
                # count number of variable positions that equal this barcode base
                sarr = group_arrays[g]
                sup = 0
                csup = 0
                for idx in var_idx:
                    if window[idx] == sarr[idx]:
                        sup += 1
                        if window_qual[idx] >= var_q:
                            csup += 1
                var_support[g] = sup
                var_confident_support[g] = csup

            # Determine if target barcode is unambiguously best on var positions
            target_sup = var_support.get(barcode_seq, 0)
            target_csup = var_confident_support.get(barcode_seq, 0)

            # find best other
            other_sups = [v for k_, v in var_support.items() if k_ != barcode_seq]
            other_csups = [v for k_, v in var_confident_support.items() if k_ != barcode_seq]
            best_other = max(other_sups) if other_sups else -1
            best_other_conf = max(other_csups) if other_csups else -1

            # uniqueness tests:
            # If there are no variable positions, accept if constant positions match
            if not var_idx:
                if penalty_const <= max_penalty:
                    return True
            else:
                # - either variable support strictly greater than best other
                # - or confident support strictly greater than best other confident support
                unique_by_var = (target_sup > best_other)
                unique_by_conf = (target_csup > best_other_conf)

                if (unique_by_var or unique_by_conf) and (penalty_const <= max_penalty):
                    # accepted: constant region ok AND variable region uniquely identifies barcode
                    return True

    return False

# ---------------------------------
# Load adapters from Python file
# ---------------------------------

def load_adapters(adapter_file):
    """
    Load adapter sequences from a Python file that defines ADAPTERS list.
    
    Parameters:
        adapter_file (str): Path to Python file containing ADAPTERS = [...]
    
    Returns:
        list: List of adapter sequences (uppercase)
    """
    if not adapter_file:
        return []
    
    print(f"Loading adapters from {adapter_file}...")
    adapters_module = {}
    with open(adapter_file, 'r') as f:
        exec(f.read(), adapters_module)
    
    adapters = adapters_module.get('ADAPTERS', [])
    print(f"Loaded {len(adapters)} adapter sequences.")
    return adapters

# ---------------------------------
# Load barcodes from CSV
# ---------------------------------

def load_barcodes(csv_file):
    print(f"Loading barcodes from {csv_file}...")
    df = pd.read_csv(csv_file, engine="python")
    row_map = {}   # idx8 -> row_letter
    col_map = {}   # idx8 -> col_number
    for _, r in df.iterrows():
        well = r["Well Position"]  # e.g. A1
        letter, number = well[0], int(well[1:])
        full = r["Sequence"].strip().upper()
        # Removed AGAT validation to support different barcode formats
        idx8 = full
        if r["Sequence Name"].startswith("R"):
            row_map[idx8] = letter
        else:
            col_map[idx8] = number
    print(f"Loaded {len(row_map)} row primers and {len(col_map)} column primers.")
    return row_map, col_map

# ---------------------------------
# Process a chunk of reads
# ---------------------------------

def process_chunk(args):
    """
    args is a tuple:
      (chunk, row_map, col_map, row_meta, col_meta, min_length, max_penalty, flank, adapters, var_q, use_fallback)
    """
    chunk, row_map, col_map, row_meta, col_meta, min_length, max_penalty, flank, adapters, var_q, use_fallback = args
    stats = defaultdict(int_defaultdict)
    grouped_reads = defaultdict(list)

    # convenience lists for group sequences (for match function)
    row_group = row_meta['seqs']
    col_group = col_meta['seqs']

    for rec in chunk:
        seq = str(rec.seq).upper()
        rc_seq = revcomp(seq)
        qual = rec.letter_annotations.get("phred_quality", [])
        if not qual:
            # fallback for rare cases where phred_quality is not set
            qual = [ord(c) - 33 for c in rec.format("fastq").splitlines()[3].strip()]
        rc_qual = qual[::-1]

        stats['GLOBAL']['total'] += 1
        if len(seq) < min_length:
            stats['GLOBAL']['too_short'] += 1
            continue

        stats['GLOBAL']['length_ok'] += 1

        # Two-tiered matching: try forward and reverse complement for both row and column
        # Find row primers using two-tiered approach
        found_rows_fwd, method_row_fwd = find_barcodes_two_tier(
            seq, qual, row_map, row_group,
            row_meta['const_idx'], row_meta['var_idx'],
            max_penalty, flank, var_q, use_fallback=use_fallback
        )
        found_rows_rc, method_row_rc = find_barcodes_two_tier(
            rc_seq, rc_qual, row_map, row_group,
            row_meta['const_idx'], row_meta['var_idx'],
            max_penalty, flank, var_q, use_fallback=use_fallback
        )
        
        # Combine forward and reverse complement results
        found_rows = list(set(found_rows_fwd + found_rows_rc))
        row_method = method_row_fwd if found_rows_fwd else method_row_rc

        # Find column primers using two-tiered approach
        found_cols_fwd, method_col_fwd = find_barcodes_two_tier(
            seq, qual, col_map, col_group,
            col_meta['const_idx'], col_meta['var_idx'],
            max_penalty, flank, var_q, use_fallback=use_fallback
        )
        found_cols_rc, method_col_rc = find_barcodes_two_tier(
            rc_seq, rc_qual, col_map, col_group,
            col_meta['const_idx'], col_meta['var_idx'],
            max_penalty, flank, var_q, use_fallback=use_fallback
        )
        
        # Combine forward and reverse complement results
        found_cols = list(set(found_cols_fwd + found_cols_rc))
        col_method = method_col_fwd if found_cols_fwd else method_col_rc

        # Track which method was used for successful mappings
        if row_method == 'alignment':
            stats['method']['row_alignment'] += 1
        elif row_method == 'fast' and len(found_rows) > 0:
            stats['method']['row_fast'] += 1
            
        if col_method == 'alignment':
            stats['method']['col_alignment'] += 1
        elif col_method == 'fast' and len(found_cols) > 0:
            stats['method']['col_fast'] += 1

        # Assignment logic
        if len(found_rows) == 1 and len(found_cols) == 1:
            well = f"{found_rows[0]}{found_cols[0]}"
            stats[well]['both'] += 1
            stats[well]['total'] += 1
            stats['GLOBAL']['mapped'] += 1
            grouped_reads[well].append(rec)

        # Exactly one row but multiple columns → ambiguous columns
        elif len(found_rows) == 1 and len(found_cols) > 1:
            stats['ambiguous']['multiple_cols'] += 1

        # Exactly one column but multiple rows → ambiguous rows
        elif len(found_cols) == 1 and len(found_rows) > 1:
            stats['ambiguous']['multiple_rows'] += 1

        # Multiple rows and multiple columns → ambiguous both
        elif len(found_rows) > 1 and len(found_cols) > 1:
            stats['ambiguous']['multiple_rows_and_cols'] += 1

        elif len(found_rows) == 1:
            stats['row_only'][found_rows[0]] += 1
            stats['GLOBAL']['single'] += 1

        elif len(found_cols) == 1:
            stats['col_only'][found_cols[0]] += 1
            stats['GLOBAL']['single'] += 1

        else:
            # Look for adapters (if provided)
            if adapters:
                adapter_found = any(
                    match_barcode_ends_weighted(seq, qual, adapter,
                                                adapters, [], [],  # adapters: no const/var distinction
                                                max_penalty, flank, var_q) or
                    match_barcode_ends_weighted(rc_seq, rc_qual, adapter,
                                                adapters, [], [], max_penalty, flank, var_q)
                    for adapter in adapters
                )
                if adapter_found:
                    stats['GLOBAL']['adapters_only'] += 1
                else:
                    stats['GLOBAL']['no_match'] += 1
            else:
                stats['GLOBAL']['no_match'] += 1

    return stats, grouped_reads

# ---------------------------------
# Read chunking and merging
# ---------------------------------

def chunk_reads(fastq_file, chunk_size=5000):
    chunk = []
    for rec in SeqIO.parse(fastq_file, "fastq"):
        chunk.append(rec)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk

def merge_dicts(dicts):
    result = defaultdict(lambda: defaultdict(int))
    for d in dicts:
        for key, sub in d.items():
            for k, v in sub.items():
                result[key][k] += v
    return result

def merge_reads(read_sets):
    result = defaultdict(list)
    for d in read_sets:
        for well, reads in d.items():
            result[well].extend(reads)
    return result

# ---------------------------------
# Write barcode match stats grid
# ---------------------------------

def write_barcode_grid_csv(stats, outpath):
    rows = []
    for letter in "ABCDEFGH":
        row = [stats[f"{letter}{col}"]["both"] for col in range(1, 13)]
        row.append(stats['row_only'][letter])
        rows.append(row)
    row_i = [stats['col_only'][col] for col in range(1, 13)] + [0]
    rows.append(row_i)

    # Add ambiguous counts row
    ambiguous = stats.get('ambiguous', {})

    g = stats['GLOBAL']
    total_row_labels = ['total', 'length_ok', 'mapped', 'single','adapter_only', 'no_match', 'too_short',
                        'ambiguous_multiple_cols', 'ambiguous_multiple_rows', 'ambiguous_multiple_rows_and_cols']
    total_row_values = [
        g['total'], g['length_ok'], g['mapped'], g['single'], g['adapters_only'], g['no_match'], g['too_short'],
        ambiguous.get('multiple_cols', 0),
        ambiguous.get('multiple_rows', 0),
        ambiguous.get('multiple_rows_and_cols', 0)
    ]
    total_row_values += [''] * (13 - len(total_row_values))

    rows.append(total_row_labels)
    rows.append(total_row_values)

    df = pd.DataFrame(rows, index=list('ABCDEFGH') + ['X', 'Stats',''], columns=[*map(str, range(1, 13)), 'X'])
    df.to_csv(outpath)

# ---------------------------------
# Main
# ---------------------------------

# ---------------------------------
# Process a single FASTQ file
# ---------------------------------

def process_single_file(fastq_file, barcode_csv, outdir, min_length, max_penalty, cpus, flank, adapter_file=None, var_q=10, use_fallback=True, , generate_report=False):
    """Process a single FASTQ file with improved two-tiered barcode matching."""
    os.makedirs(outdir, exist_ok=True)
    row_map, col_map = load_barcodes(barcode_csv)
    
    # Analyze barcode sets (to determine const/var positions and distances)
    print("\\nAnalyzing row barcodes:")
    row_meta = analyze_barcode_set(list(row_map.keys()))
    print("\\nAnalyzing column barcodes:")
    col_meta = analyze_barcode_set(list(col_map.keys()))
    print()
    
    adapters = load_adapters(adapter_file) if adapter_file else []
    pool = Pool(processes=cpus)
    jobs = [(chunk, row_map, col_map, row_meta, col_meta, min_length, max_penalty, flank, adapters, var_q, use_fallback)
            for chunk in chunk_reads(fastq_file)]
    results = pool.map(process_chunk, jobs)
    pool.close()
    pool.join()

    all_stats = merge_dicts([r[0] for r in results])
    all_reads = merge_reads([r[1] for r in results])

    # Save matched reads by well
    for well, reads in all_reads.items():
        path = os.path.join(outdir, f"{well}_reads.fastq")
        SeqIO.write(reads, path, "fastq")

    # Save stats grid
    write_barcode_grid_csv(all_stats, os.path.join(outdir, "barcode_stats.csv"))
    
    # Print method usage statistics
    method_stats = all_stats.get('method', {})
    if method_stats:
        print(f"\n📊 Method Usage Statistics:")
        print(f"  Fast filter (rows): {method_stats.get('row_fast', 0)}")
        print(f"  Alignment fallback (rows): {method_stats.get('row_alignment', 0)}")
        print(f"  Fast filter (cols): {method_stats.get('col_fast', 0)}")
        print(f"  Alignment fallback (cols): {method_stats.get('col_alignment', 0)}")
    
    print(f"✅ Processed {fastq_file}")
    
    # Generate quality report if requested
    if generate_report:
        try:
            print(f"📊 Generating quality report...")
            import subprocess
            # Use sys.executable for portability and construct path relative to current script
            script_dir = os.path.dirname(os.path.abspath(__file__))
            report_script = os.path.join(script_dir, 'generate_quality_report.py')
            subprocess.run([
                sys.executable, report_script,
                outdir, barcode_csv
            ], check=True)
        except Exception as e:
            print(f"⚠️  Warning: Could not generate quality report: {e}")


# ---------------------------------
# Main function - handles both files and directories
# ---------------------------------

def main(fastq_input, barcode_csv, outdir, min_length, max_penalty, cpus, flank, adapter_file=None, var_q=10, use_fallback=True, generate_report=False):
    """
    Main entry point that handles both single files and directories with improved two-tiered barcode matching.
    
    Parameters:
        fastq_input: Path to either a FASTQ file or a directory containing FASTQ files
        barcode_csv: CSV file with barcode-to-well mappings
        outdir: Base output directory (will be adjusted for directory inputs)
        min_length: Minimum read length
        max_penalty: Maximum quality-weighted mismatch penalty
        cpus: Number of CPU cores
        flank: Number of bases at each end to search
        adapter_file: Optional adapter file
        generate_report: Whether to generate quality report after demultiplexing
        var_q: Phred cutoff to call a variable-base match 'confident' (default: 10)
        use_fallback: Use alignment fallback for ambiguous cases (default: True)
    """
    # Check if input is a file or directory
    if os.path.isfile(fastq_input):
        # Single file mode - use outdir as-is or adjust if not specified
        if outdir == "demuxed":
            # Default behavior: extract filename without extension
            basename = get_basename_without_extensions(fastq_input)
            outdir = os.path.join("demplex_data", basename)
        print(f"Processing single file: {fastq_input}")
        process_single_file(fastq_input, barcode_csv, outdir, min_length, max_penalty, cpus, flank, adapter_file, var_q, use_fallback, generate_report)
    elif os.path.isdir(fastq_input):
        # Directory mode - process all FASTQ files
        print(f"Processing directory: {fastq_input}")
        # Find all FASTQ files in the directory
        fastq_files = []
        for ext in ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz']:
            fastq_files.extend(glob.glob(os.path.join(fastq_input, ext)))
        
        if not fastq_files:
            print(f"⚠️  No FASTQ files found in {fastq_input}")
            return
        
        # Get the directory name for output structure
        dir_basename = os.path.basename(os.path.normpath(fastq_input))
        base_outdir = os.path.join("demplex_data", dir_basename)
        
        print(f"Found {len(fastq_files)} FASTQ file(s)")
        for fastq_file in fastq_files:
            # Create output directory for each file
            file_basename = get_basename_without_extensions(fastq_file)
            file_outdir = os.path.join(base_outdir, file_basename)
            print(f"\n📂 Processing: {os.path.basename(fastq_file)}")
            process_single_file(fastq_file, barcode_csv, file_outdir, min_length, max_penalty, cpus, flank, adapter_file, var_q, use_fallback, generate_report)
        
        print(f"\n✅ All files processed. Output in: {base_outdir}")
    else:
        raise ValueError(f"Input path does not exist: {fastq_input}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Demultiplex and analyze barcoded linear PCR reads from single or multiple FASTQ files.",
        epilog="""
Examples:
  # Single file - basic usage (no adapter detection)
  python demux_barcodes.py raw_data/reads.fastq barcodes/251202_primer_well_map_DA.csv
  
  # Single file - with adapter detection
  python demux_barcodes.py raw_data/reads.fastq barcodes/251202_primer_well_map_DA.csv \\
      --adapters barcodes/251205_adapters_DA.py
  
  # Single file - with custom output directory
  python demux_barcodes.py raw_data/reads.fastq barcodes/251202_primer_well_map_DA.csv \\
      --outdir output/ --cpus 4 --flank 100
  
  # Multiple files - process all FASTQ files in a directory
  python demux_barcodes.py raw_data/experiment1/ barcodes/251202_primer_well_map_DA.csv \\
      --adapters barcodes/251205_adapters_DA.py --cpus 4
  
  # Generate quality report after demultiplexing
  python demux_barcodes.py raw_data/reads.fastq barcodes/251202_primer_well_map_DA.csv \\
      --report
      
  Note: When processing a directory, output will be in demplex_data/<dirname>/<filename>/
        When processing a single file, output will be in demplex_data/<filename>/ (unless --outdir is specified)
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("fastq", help="Input FASTQ file or directory containing FASTQ files")
    parser.add_argument("barcodes", help="CSV file with barcode-to-well mappings (e.g., barcodes/251202_primer_well_map_DA.csv)")
    parser.add_argument("--adapters", "-a", dest="adapter_file", default=None,
                        help="Python file defining ADAPTERS list (e.g., barcodes/251205_adapters_DA.py). Optional; when omitted, adapter detection is skipped.")
    parser.add_argument("--outdir", default="demuxed", help="Output directory for matched reads [default: demuxed or auto-generated based on input]")
    parser.add_argument("--min_length", type=int, default=50, help="Minimum read length [default: 50]")
    parser.add_argument("--max_penalty", type=float, default=60, help="Max allowed mismatch penalty (sum of Phred scores) [default: 60]")
    parser.add_argument("--cpus", type=int, default=1, help="Number of CPU cores to use [default: 1]")
    parser.add_argument("--flank", type=int, default=100, help="Number of bases at each end to search for barcodes [default: 100]")
    parser.add_argument("--var_q", type=int, default=10,
                        help="Phred cutoff to call a variable-base match 'confident' [default: 10]")
    parser.add_argument("--no-fallback", dest="use_fallback", action="store_false", default=True,
                        help="Disable alignment fallback for ambiguous cases (use only fast filter) [default: enabled]")
    parser.add_argument("--report", "-r", action="store_true", 
                        help="Generate graphical quality report with MSA, read lengths, and barcode analysis after demultiplexing")
    args = parser.parse_args()

    main(args.fastq, args.barcodes, args.outdir, args.min_length, args.max_penalty, args.cpus, args.flank, args.adapter_file, args.var_q, args.use_fallback, args.report)
