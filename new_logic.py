#!/usr/bin/env python3

"""
Demultiplex nanopore FASTQ by row/column barcodes with:
 - barcode analysis (constant vs variable positions)
 - quality-weighted matching on constant positions
 - variable-region disambiguation tolerant of low-quality bases (but requires uniqueness)
Minimal changes to existing pipeline; multiprocessing-safe (top-level functions).
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from collections import defaultdict
import os
from multiprocessing import Pool
import glob
import itertools

# ---------------------------------
# Utilities / helpers
# ---------------------------------

def revcomp(seq):
    return str(Seq(seq).reverse_complement())

def int_defaultdict():
    return defaultdict(int)

def get_basename_without_extensions(filename):
    basename = os.path.basename(filename)
    for ext in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']:
        if basename.endswith(ext):
            return basename[:-len(ext)]
    return os.path.splitext(basename)[0]

def hamming(a, b):
    """Hamming distance of equal-length strings."""
    if len(a) != len(b):
        raise ValueError("Hamming requires equal lengths")
    return sum(x != y for x, y in zip(a, b))

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
# Matching function (new logic)
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
                # too many/index-constant mismatches; reject this window fast
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
            # - either variable support strictly greater than best other
            # - or confident support strictly greater than best other confident support
            unique_by_var = (target_sup > best_other)
            unique_by_conf = (target_csup > best_other_conf)

            if (unique_by_var or unique_by_conf) and (penalty_const <= max_penalty):
                # accepted: constant region ok AND variable region uniquely identifies barcode
                return True

            # else continue searching windows

    return False

# ---------------------------------
# Load adapters from Python file
# ---------------------------------

def load_adapters(adapter_file):
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
# Load barcodes from CSV (unchanged format)
# ---------------------------------

def load_barcodes(csv_file):
    print(f"Loading barcodes from {csv_file}...")
    df = pd.read_csv(csv_file, engine="python")
    row_map = {}   # full_seq -> row_letter (e.g. 'AATG...': 'A')
    col_map = {}   # full_seq -> col_number (e.g. 'CAAG...': 1)
    for _, r in df.iterrows():
        well = r["Well Position"]
        letter, number = well[0], int(well[1:])
        full = r["Sequence"].strip().upper()
        if "AGAT" not in full:
            # we allow barcodes that don't contain AGAT but warn
            print(f"Warning: 'AGAT' not found in adapter for well {well}: {full[:20]}...")
        # store full sequence as the barcode key
        if r["Sequence Name"].split()[0] == "R":
            row_map[full] = letter
        else:
            col_map[full] = number
    print(f"Loaded {len(row_map)} row primers and {len(col_map)} column primers.")
    return row_map, col_map

# ---------------------------------
# Process a chunk of reads
# ---------------------------------

def process_chunk(args):
    """
    args is a tuple:
      (chunk, row_map, col_map, row_meta, col_meta, min_length, max_penalty, flank, adapters, var_q)
    """
    chunk, row_map, col_map, row_meta, col_meta, min_length, max_penalty, flank, adapters, var_q = args

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
            # try to coerce from ascii qualities if present (rare)
            # but BioPython usually sets phred_quality
            qual = [ord(c) - 33 for c in rec.format("fastq").splitlines()[3].strip()]

        rc_qual = qual[::-1]

        stats['GLOBAL']['total'] += 1
        if len(seq) < min_length:
            stats['GLOBAL']['too_short'] += 1
            stats['unassigned']['too_short'] += 1
            continue

        stats['GLOBAL']['length_ok'] += 1

        # Find row primers (check each barcode in row_map)
        found_rows = [row_map[idx] for idx in row_map
                      if match_barcode_ends_weighted(seq, qual, idx,
                                                     row_group,
                                                     row_meta['const_idx'], row_meta['var_idx'],
                                                     max_penalty, flank, var_q)
                      or match_barcode_ends_weighted(rc_seq, rc_qual, idx,
                                                     row_group,
                                                     row_meta['const_idx'], row_meta['var_idx'],
                                                     max_penalty, flank, var_q)]

        # Find column primers
        found_cols = [col_map[idx] for idx in col_map
                      if match_barcode_ends_weighted(seq, qual, idx,
                                                     col_group,
                                                     col_meta['const_idx'], col_meta['var_idx'],
                                                     max_penalty, flank, var_q)
                      or match_barcode_ends_weighted(rc_seq, rc_qual, idx,
                                                     col_group,
                                                     col_meta['const_idx'], col_meta['var_idx'],
                                                     max_penalty, flank, var_q)]

        # mapping logic with ambiguous tracking
        if len(found_rows) == 1 and len(found_cols) == 1:
            well = f"{found_rows[0]}{found_cols[0]}"
            stats[well]['both'] += 1
            stats[well]['total'] += 1
            stats['GLOBAL']['mapped'] += 1
            grouped_reads[well].append(rec)

        elif len(found_rows) == 1 and len(found_cols) > 1:
            stats['ambiguous']['multiple_cols'] += 1
            stats['unassigned']['ambiguous'] += 1

        elif len(found_cols) == 1 and len(found_rows) > 1:
            stats['ambiguous']['multiple_rows'] += 1
            stats['unassigned']['ambiguous'] += 1

        elif len(found_rows) > 1 and len(found_cols) > 1:
            stats['ambiguous']['multiple_rows_and_cols'] += 1
            stats['unassigned']['ambiguous'] += 1

        elif len(found_rows) == 1:
            stats['row_only'][found_rows[0]] += 1
            stats['GLOBAL']['single'] += 1
            # still "usable" for consensus building — keep in unassigned bucket if needed
            stats['unassigned']['row_only'] += 1

        elif len(found_cols) == 1:
            stats['col_only'][found_cols[0]] += 1
            stats['GLOBAL']['single'] += 1
            stats['unassigned']['col_only'] += 1

        else:
            # no row/col match -> check adapters (if provided)
            if adapters:
                adapter_found = any(
                    match_barcode_ends_weighted(seq, qual, adapter,
                                                adapters,  # adapter group (treat each adapter independently for const/var)
                                                [], [],  # adapter usually all variable but we check all positions (no const/var)
                                                max_penalty,
                                                flank, var_q)
                    or match_barcode_ends_weighted(rc_seq, rc_qual, adapter,
                                                   adapters, [], [], max_penalty, flank, var_q)
                    for adapter in adapters
                )
                if adapter_found:
                    stats['GLOBAL']['adapters_only'] += 1
                else:
                    stats['GLOBAL']['no_match'] += 1
                    stats['unassigned']['neither'] += 1
            else:
                stats['GLOBAL']['no_match'] += 1
                stats['unassigned']['neither'] += 1

    return stats, grouped_reads

# ---------------------------------
# chunking and merging helpers
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
# stats grid writing (updated)
# ---------------------------------

def write_barcode_grid_csv(stats, outpath):
    rows = []
    for letter in "ABCDEFGH":
        row = [stats[f"{letter}{col}"]["both"] for col in range(1, 13)]
        row.append(stats['row_only'][letter])
        rows.append(row)
    row_i = [stats['col_only'][col] for col in range(1, 13)] + [0]
    rows.append(row_i)

    ambiguous = stats.get('ambiguous', {})
    g = stats['GLOBAL']
    total_row_labels = ['total', 'length_ok', 'mapped', 'single', 'adapters_only', 'no_match', 'too_short',
                        'ambiguous_multiple_cols', 'ambiguous_multiple_rows', 'ambiguous_multiple_rows_and_cols']
    total_row_values = [
        g.get('total', 0), g.get('length_ok', 0), g.get('mapped', 0), g.get('single', 0),
        g.get('adapters_only', 0), g.get('no_match', 0), g.get('too_short', 0),
        ambiguous.get('multiple_cols', 0),
        ambiguous.get('multiple_rows', 0),
        ambiguous.get('multiple_rows_and_cols', 0)
    ]
    total_row_values += [''] * (13 - len(total_row_values))

    rows.append(total_row_labels)
    rows.append(total_row_values)

    df = pd.DataFrame(rows, index=list('ABCDEFGH') + ['I', 'Stats', ''], columns=[*map(str, range(1, 13)), 'X'])
    df.to_csv(outpath)

# ---------------------------------
# Main processing functions
# ---------------------------------

def process_single_file(fastq_file, barcode_csv, outdir, min_length, max_penalty, cpus, flank, adapter_file=None, var_q=10):
    os.makedirs(outdir, exist_ok=True)
    row_map, col_map = load_barcodes(barcode_csv)

    # analyze barcode sets (to determine const/var positions and distances)
    row_meta = analyze_barcode_set(list(row_map.keys()))
    col_meta = analyze_barcode_set(list(col_map.keys()))

    adapters = load_adapters(adapter_file) if adapter_file else []

    pool = Pool(processes=cpus)
    jobs = []
    for chunk in chunk_reads(fastq_file):
        jobs.append((chunk, row_map, col_map, row_meta, col_meta, min_length, max_penalty, flank, adapters, var_q))

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
    print(f"✅ Processed {fastq_file}")

def main(fastq_input, barcode_csv, outdir, min_length, max_penalty, cpus, flank, adapter_file=None, var_q=10):
    if os.path.isfile(fastq_input):
        if outdir == "demuxed":
            basename = get_basename_without_extensions(fastq_input)
            outdir = os.path.join("demplex_data", basename)
        print(f"Processing single file: {fastq_input}")
        process_single_file(fastq_input, barcode_csv, outdir, min_length, max_penalty, cpus, flank, adapter_file, var_q)
    elif os.path.isdir(fastq_input):
        print(f"Processing directory: {fastq_input}")
        fastq_files = []
        for ext in ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz']:
            fastq_files.extend(glob.glob(os.path.join(fastq_input, ext)))
        if not fastq_files:
            print(f"⚠️  No FASTQ files found in {fastq_input}")
            return
        dir_basename = os.path.basename(os.path.normpath(fastq_input))
        base_outdir = os.path.join("demplex_data", dir_basename)
        print(f"Found {len(fastq_files)} FASTQ file(s)")
        for fastq_file in fastq_files:
            file_basename = get_basename_without_extensions(fastq_file)
            file_outdir = os.path.join(base_outdir, file_basename)
            print(f"\n📂 Processing: {os.path.basename(fastq_file)}")
            process_single_file(fastq_file, barcode_csv, file_outdir, min_length, max_penalty, cpus, flank, adapter_file, var_q)
        print(f"\n✅ All files processed. Output in: {base_outdir}")
    else:
        raise ValueError(f"Input path does not exist: {fastq_input}")

# ---------------------------------
# CLI
# ---------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Demultiplex and analyze barcoded linear PCR reads from single or multiple FASTQ files.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("fastq", help="Input FASTQ file or directory containing FASTQ files")
    parser.add_argument("barcodes", help="CSV file with barcode-to-well mappings")
    parser.add_argument("--adapters", "-a", dest="adapter_file", default=None,
                        help="Python file defining ADAPTERS list (optional)")
    parser.add_argument("--outdir", default="demuxed", help="Output directory [default: demuxed or auto-generated]")
    parser.add_argument("--min_length", type=int, default=50, help="Minimum read length [default: 50]")
    parser.add_argument("--max_penalty", type=float, default=60, help="Max allowed Phred-sum penalty on constant positions [default: 60]")
    parser.add_argument("--cpus", type=int, default=1, help="Number of CPU cores to use [default: 1]")
    parser.add_argument("--flank", type=int, default=100, help="Number of bases at each end to search for barcodes [default: 100]")
    parser.add_argument("--var_q", type=int, default=10, help="Phred cutoff to call a variable-base match 'confident' [default: 10]")
    args = parser.parse_args()

    main(args.fastq, args.barcodes, args.outdir, args.min_length, args.max_penalty, args.cpus, args.flank, args.adapter_file, args.var_q)
