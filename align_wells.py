#!/usr/bin/env python3
"""
Perform sequence alignment on demultiplexed well data.

This script:
1. Reads demultiplexed FASTQ files for each well
2. Performs multiple sequence alignment (MSA) using quality-weighted consensus
3. Generates consensus sequences taking quality scores into account
4. Outputs aligned FASTQ files for each well
5. Creates a CSV summary with all wells and their consensus sequences
"""

import argparse
import os
import sys
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import parasail
import numpy as np


def calculate_quality_weighted_consensus(sequences, qualities, min_quality=20, min_coverage=1):
    """
    Calculate consensus sequence from multiple sequences with quality scores.
    
    Parameters:
        sequences: list of sequence strings (aligned or unaligned)
        qualities: list of quality score lists (Phred scores)
        min_quality: minimum quality score to consider a base (default: 20)
        min_coverage: minimum number of reads required at a position
    
    Returns:
        consensus_seq: consensus sequence string
        consensus_qual: consensus quality scores (list)
        alignment_depth: coverage at each position
    """
    if not sequences:
        return "", [], []
    
    # Find the maximum length
    max_len = max(len(seq) for seq in sequences)
    
    # Pad sequences and qualities to max length
    padded_seqs = []
    padded_quals = []
    for seq, qual in zip(sequences, qualities):
        padded_seq = seq + '-' * (max_len - len(seq))
        padded_qual = qual + [0] * (max_len - len(qual))
        padded_seqs.append(padded_seq)
        padded_quals.append(padded_qual)
    
    consensus_seq = []
    consensus_qual = []
    alignment_depth = []
    
    for pos in range(max_len):
        # Collect bases and qualities at this position
        base_qual_pairs = []
        for seq, qual in zip(padded_seqs, padded_quals):
            base = seq[pos]
            q = qual[pos]
            if base != '-' and q >= min_quality:
                base_qual_pairs.append((base, q))
        
        alignment_depth.append(len(base_qual_pairs))
        
        if len(base_qual_pairs) < min_coverage:
            # Not enough coverage - use N with low quality
            consensus_seq.append('N')
            consensus_qual.append(2)
            continue
        
        # Calculate quality-weighted base frequency
        base_weights = defaultdict(float)
        for base, qual in base_qual_pairs:
            # Convert Phred to probability: P_correct = 1 - 10^(-Q/10)
            # Weight by quality score directly (higher quality = more weight)
            weight = qual
            base_weights[base] += weight
        
        # Select base with highest total weight
        best_base = max(base_weights.items(), key=lambda x: x[1])[0]
        consensus_seq.append(best_base)
        
        # Calculate consensus quality as the average quality of matching bases
        matching_quals = [q for b, q in base_qual_pairs if b == best_base]
        avg_qual = int(np.mean(matching_quals)) if matching_quals else 2
        consensus_qual.append(min(avg_qual, 60))  # Cap at Q60
    
    return ''.join(consensus_seq), consensus_qual, alignment_depth


def align_sequences_pairwise(sequences, qualities):
    """
    Align sequences using parasail Smith-Waterman alignment.
    Returns aligned sequences and qualities.
    
    Parameters:
        sequences: list of sequence strings
        qualities: list of quality score lists
    
    Returns:
        aligned_seqs: list of aligned sequence strings
        aligned_quals: list of aligned quality score lists
    """
    if not sequences:
        return [], []
    
    if len(sequences) == 1:
        return sequences, qualities
    
    # Use first sequence as reference
    ref_seq = sequences[0]
    ref_qual = qualities[0]
    
    aligned_seqs = [ref_seq]
    aligned_quals = [ref_qual]
    
    # Create scoring matrix
    matrix = parasail.matrix_create("ACGT", 2, -1)
    
    for seq, qual in zip(sequences[1:], qualities[1:]):
        # Perform Smith-Waterman alignment
        result = parasail.sw_trace_scan_16(ref_seq, seq, 2, 1, matrix)
        
        # Get alignment
        ref_aligned = result.traceback.ref
        query_aligned = result.traceback.query
        
        # Adjust qualities based on alignment
        ref_qual_aligned = []
        query_qual_aligned = []
        
        ref_idx = 0
        query_idx = 0
        
        for i in range(len(ref_aligned)):
            if ref_aligned[i] == '-':
                ref_qual_aligned.append(0)
            else:
                if ref_idx < len(ref_qual):
                    ref_qual_aligned.append(ref_qual[ref_idx])
                else:
                    ref_qual_aligned.append(0)
                ref_idx += 1
            
            if query_aligned[i] == '-':
                query_qual_aligned.append(0)
            else:
                if query_idx < len(qual):
                    query_qual_aligned.append(qual[query_idx])
                else:
                    query_qual_aligned.append(0)
                query_idx += 1
        
        aligned_seqs.append(query_aligned)
        aligned_quals.append(query_qual_aligned)
    
    return aligned_seqs, aligned_quals


def process_well_fastq(fastq_path, max_reads=None, align_sequences=True):
    """
    Process a single well's FASTQ file.
    
    Parameters:
        fastq_path: path to well FASTQ file
        max_reads: maximum number of reads to process (None = all)
        align_sequences: whether to perform MSA (default: True)
    
    Returns:
        dict with:
            - sequences: list of sequences
            - qualities: list of quality score lists
            - records: list of SeqRecord objects
            - consensus_seq: consensus sequence string
            - consensus_qual: consensus quality scores
            - alignment_depth: coverage at each position
    """
    records = []
    sequences = []
    qualities = []
    
    count = 0
    for record in SeqIO.parse(fastq_path, "fastq"):
        if max_reads and count >= max_reads:
            break
        
        records.append(record)
        sequences.append(str(record.seq).upper())
        qualities.append(record.letter_annotations.get("phred_quality", [30] * len(record.seq)))
        count += 1
    
    if not sequences:
        return {
            'sequences': [],
            'qualities': [],
            'records': [],
            'consensus_seq': '',
            'consensus_qual': [],
            'alignment_depth': [],
            'num_reads': 0
        }
    
    # Perform alignment if requested
    if align_sequences and len(sequences) > 1:
        aligned_seqs, aligned_quals = align_sequences_pairwise(sequences, qualities)
    else:
        aligned_seqs = sequences
        aligned_quals = qualities
    
    # Calculate consensus
    consensus_seq, consensus_qual, alignment_depth = calculate_quality_weighted_consensus(
        aligned_seqs, aligned_quals, min_quality=20, min_coverage=1
    )
    
    return {
        'sequences': sequences,
        'qualities': qualities,
        'records': records,
        'aligned_sequences': aligned_seqs,
        'aligned_qualities': aligned_quals,
        'consensus_seq': consensus_seq,
        'consensus_qual': consensus_qual,
        'alignment_depth': alignment_depth,
        'num_reads': len(sequences)
    }


def write_aligned_fastq(well_data, output_path, well_name):
    """
    Write aligned sequences to FASTQ file.
    
    Parameters:
        well_data: dict from process_well_fastq
        output_path: path to output FASTQ file
        well_name: name of the well (e.g., 'A1')
    """
    records_to_write = []
    
    # Add consensus sequence as first record
    if well_data['consensus_seq']:
        consensus_record = SeqRecord(
            Seq(well_data['consensus_seq']),
            id=f"{well_name}_consensus",
            description=f"Quality-weighted consensus sequence for well {well_name} ({well_data['num_reads']} reads)",
            letter_annotations={'phred_quality': well_data['consensus_qual']}
        )
        records_to_write.append(consensus_record)
    
    # Add original reads
    records_to_write.extend(well_data['records'])
    
    SeqIO.write(records_to_write, output_path, "fastq")


def process_demux_directory(demux_dir, output_dir, max_reads=None, align=True):
    """
    Process all well FASTQ files in a demultiplexed directory.
    
    Parameters:
        demux_dir: directory containing demultiplexed FASTQ files
        output_dir: directory for output files
        max_reads: maximum reads per well (None = all)
        align: whether to perform MSA
    
    Returns:
        dict mapping well names to well data
    """
    os.makedirs(output_dir, exist_ok=True)
    
    well_data = {}
    
    # Find all well FASTQ files
    fastq_files = [f for f in os.listdir(demux_dir) if f.endswith('_reads.fastq')]
    
    if not fastq_files:
        print(f"⚠️  No demultiplexed FASTQ files found in {demux_dir}")
        return well_data
    
    print(f"Found {len(fastq_files)} well FASTQ files")
    
    for fastq_file in sorted(fastq_files):
        well_name = fastq_file.replace('_reads.fastq', '')
        fastq_path = os.path.join(demux_dir, fastq_file)
        
        print(f"Processing well {well_name}...", end=' ')
        
        data = process_well_fastq(fastq_path, max_reads=max_reads, align_sequences=align)
        
        if data['num_reads'] == 0:
            print("no reads")
            continue
        
        well_data[well_name] = data
        
        # Write aligned FASTQ
        output_fastq = os.path.join(output_dir, f"{well_name}_aligned.fastq")
        write_aligned_fastq(data, output_fastq, well_name)
        
        print(f"{data['num_reads']} reads, consensus length: {len(data['consensus_seq'])}")
    
    return well_data


def write_consensus_csv(well_data, output_path):
    """
    Write CSV file with consensus sequences for all wells.
    
    Parameters:
        well_data: dict mapping well names to well data
        output_path: path to output CSV file
    """
    rows = []
    
    for well_name in sorted(well_data.keys()):
        data = well_data[well_name]
        
        # Calculate average depth
        avg_depth = np.mean(data['alignment_depth']) if data['alignment_depth'] else 0
        
        # Calculate consensus quality metrics
        avg_qual = np.mean(data['consensus_qual']) if data['consensus_qual'] else 0
        
        rows.append({
            'Well': well_name,
            'Num_Reads': data['num_reads'],
            'Consensus_Length': len(data['consensus_seq']),
            'Consensus_Sequence': data['consensus_seq'],
            'Avg_Coverage': f"{avg_depth:.1f}",
            'Avg_Consensus_Quality': f"{avg_qual:.1f}",
            'Min_Coverage': min(data['alignment_depth']) if data['alignment_depth'] else 0,
            'Max_Coverage': max(data['alignment_depth']) if data['alignment_depth'] else 0
        })
    
    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    
    print(f"\n✅ Consensus CSV saved to {output_path}")
    print(f"   Total wells: {len(rows)}")


def main():
    parser = argparse.ArgumentParser(
        description="Perform sequence alignment on demultiplexed well data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage - align all reads in each well
  python align_wells.py demplex_data/55XPXK_1_P4_323_EG/ aligned_output/
  
  # Limit reads per well and disable alignment
  python align_wells.py demplex_data/experiment/ aligned_output/ --max-reads 1000 --no-align
  
  # Custom output locations
  python align_wells.py demux_dir/ output_dir/ --csv summary.csv

The script will:
  1. Read all *_reads.fastq files from the demux directory
  2. Perform MSA on reads within each well (if --no-align not specified)
  3. Calculate quality-weighted consensus sequences
  4. Write aligned FASTQ files (with consensus as first sequence)
  5. Generate a CSV summary with consensus sequences for all wells
        """
    )
    
    parser.add_argument(
        'demux_dir',
        help='Directory containing demultiplexed FASTQ files (output from demux_barcodes.py)'
    )
    parser.add_argument(
        'output_dir',
        help='Output directory for aligned FASTQ files and CSV summary'
    )
    parser.add_argument(
        '--max-reads',
        type=int,
        default=None,
        help='Maximum number of reads to process per well (default: all reads)'
    )
    parser.add_argument(
        '--min-quality',
        type=int,
        default=20,
        help='Minimum Phred quality score to consider a base for consensus (default: 20)'
    )
    parser.add_argument(
        '--no-align',
        action='store_true',
        help='Skip MSA and calculate consensus from unaligned sequences'
    )
    parser.add_argument(
        '--csv',
        default=None,
        help='Custom path for consensus CSV output (default: <output_dir>/consensus_sequences.csv)'
    )
    
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.isdir(args.demux_dir):
        print(f"❌ Error: Input directory not found: {args.demux_dir}")
        sys.exit(1)
    
    print(f"🧬 NanoDemux Well Alignment Tool")
    print(f"   Input: {args.demux_dir}")
    print(f"   Output: {args.output_dir}")
    print(f"   Alignment: {'disabled' if args.no_align else 'enabled'}")
    if args.max_reads:
        print(f"   Max reads per well: {args.max_reads}")
    print()
    
    # Process all wells
    well_data = process_demux_directory(
        args.demux_dir,
        args.output_dir,
        max_reads=args.max_reads,
        align=not args.no_align
    )
    
    if not well_data:
        print("❌ No wells processed. Check input directory.")
        sys.exit(1)
    
    # Write consensus CSV
    csv_path = args.csv or os.path.join(args.output_dir, 'consensus_sequences.csv')
    write_consensus_csv(well_data, csv_path)
    
    print(f"\n✅ Processing complete!")
    print(f"   Aligned FASTQ files: {args.output_dir}/*_aligned.fastq")
    print(f"   Consensus CSV: {csv_path}")


if __name__ == '__main__':
    main()
