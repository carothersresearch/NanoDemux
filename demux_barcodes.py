#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from collections import defaultdict
import os
from multiprocessing import Pool



# ---------------------------------
# Helpers
# ---------------------------------

def revcomp(seq):
    return str(Seq(seq).reverse_complement())

def int_defaultdict():
    return defaultdict(int)

def match_barcode_ends_weighted(seq, qual, barcode, max_penalty=60, flank=20):
    """
    Match barcode at either end of the read using quality-weighted mismatches.

    Parameters:
        seq (str): DNA sequence (A/T/C/G).
        qual (List[int]): Phred scores (integers).
        barcode (str): Barcode to match.
        max_penalty (float): Maximum allowed sum of Phred scores at mismatched bases.
        flank (int): How far to search from each end of the read.

    Returns:
        bool: True if barcode found within mismatch quality penalty threshold, else False.
    """
    k = len(barcode)
    ends = [
        (seq[:flank], qual[:flank]),
        (seq[-flank:], qual[-flank:])
    ]

    for region_seq, region_qual in ends:
        for i in range(len(region_seq) - k + 1):
            window = region_seq[i:i + k]
            window_qual = region_qual[i:i + k]
            # Sum Phred scores only where mismatches occur
            penalty = sum(
                q if bc != wc else 0
                for bc, wc, q in zip(barcode, window, window_qual)
            )
            if penalty <= max_penalty:
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
        if "AGAT" not in full:
            raise ValueError(f"No AGAT in adapter: {full}")
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
    chunk, row_map, col_map, min_length, max_penalty, flank, adapters = args
    stats = defaultdict(int_defaultdict)
    grouped_reads = defaultdict(list)

    for rec in chunk:
        seq = str(rec.seq).upper()
        rc_seq = revcomp(seq)
        qual = rec.letter_annotations["phred_quality"]
        rc_qual = qual[::-1]  # reverse complement quality

        stats['GLOBAL']['total'] += 1
        if len(seq) < min_length:
            stats['GLOBAL']['too_short'] += 1
            continue

        stats['GLOBAL']['length_ok'] += 1

        # Find row primers
        found_rows = [row_map[idx] for idx in row_map
                      if match_barcode_ends_weighted(seq, qual, idx, max_penalty, flank)
                      or match_barcode_ends_weighted(rc_seq, rc_qual, idx, max_penalty, flank)]

        # Find column primers
        found_cols = [col_map[idx] for idx in col_map
                      if match_barcode_ends_weighted(seq, qual, idx, max_penalty, flank)
                      or match_barcode_ends_weighted(rc_seq, rc_qual, idx, max_penalty, flank)]

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
                    match_barcode_ends_weighted(seq, qual, adapter, max_penalty, flank) or
                    match_barcode_ends_weighted(rc_seq, rc_qual, adapter, max_penalty, flank)
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
    g = stats['GLOBAL']
    total_row = ['total','length_ok','mapped','single','no_match','too_short'] + [''] * 7
    rows.append(total_row)
    total_row = [g['total'], g['length_ok'], g['mapped'], g['single'], g['no_match'], g['too_short']] + [''] * 7
    rows.append(total_row)
    df = pd.DataFrame(rows, index=list('ABCDEFGH') + ['X', 'Stats',''], columns=[*map(str, range(1, 13)), 'X'])
    df.to_csv(outpath)

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

def main(fastq_file, barcode_csv, outdir, min_length, max_penalty, cpus, flank, adapter_file=None):
    os.makedirs(outdir, exist_ok=True)
    row_map, col_map = load_barcodes(barcode_csv)
    adapters = load_adapters(adapter_file) if adapter_file else []
    pool = Pool(processes=cpus)
    jobs = [(chunk, row_map, col_map, min_length, max_penalty, flank, adapters)
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
    print("✅ Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Demultiplex and analyze barcoded linear PCR reads.",
        epilog="""
Examples:
  # Basic usage (no adapter detection)
  python demux_barcodes.py reads.fastq barcodes/251202_primer_well_map_DA.csv
  
  # With adapter detection from a Python file defining ADAPTERS
  python demux_barcodes.py reads.fastq barcodes/251202_primer_well_map_DA.csv \\
      --adapters barcodes/251205_adapters_DA.py
  
  # With additional options
  python demux_barcodes.py reads.fastq barcodes/251202_primer_well_map_DA.csv \\
      --adapters barcodes/251205_adapters_DA.py \\
      --outdir output/ --cpus 4 --flank 100
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("fastq", help="Input FASTQ file")
    parser.add_argument("barcodes", help="CSV file with barcode-to-well mappings (e.g., barcodes/251202_primer_well_map_DA.csv)")
    parser.add_argument("--adapters", "-a", dest="adapter_file", default=None,
                        help="Python file defining ADAPTERS list (e.g., barcodes/251205_adapters_DA.py). Optional; when omitted, adapter detection is skipped.")
    parser.add_argument("--outdir", default="demuxed", help="Output directory for matched reads [default: demuxed]")
    parser.add_argument("--min_length", type=int, default=50, help="Minimum read length [default: 50]")
    parser.add_argument("--max_penalty", type=float, default=60, help="Max allowed mismatch penalty (sum of Phred scores) [default: 60]")
    parser.add_argument("--cpus", type=int, default=1, help="Number of CPU cores to use [default: 1]")
    parser.add_argument("--flank", type=int, default=100, help="Number of bases at each end to search for barcodes [default: 100]")
    args = parser.parse_args()

    main(args.fastq, args.barcodes, args.outdir, args.min_length, args.max_penalty, args.cpus, args.flank, args.adapter_file)
