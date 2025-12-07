#!/usr/bin/env python3
"""
Generate graphical data quality report for demultiplexed reads.

This script creates comprehensive visualizations including:
- Multiple Sequence Alignment (MSA) visualization showing relative positions
- Read length distributions
- Barcode presence and position analysis
"""

import argparse
import os
import sys
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import numpy as np
import seaborn as sns
from datetime import datetime
import parasail

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300

# Alignment constants for Smith-Waterman
ALIGNMENT_GAP_OPEN = 2
ALIGNMENT_GAP_EXTEND = 1


def load_barcode_map(barcode_csv):
    """Load barcode sequences from CSV file."""
    df = pd.read_csv(barcode_csv)
    barcodes = {}
    for _, row in df.iterrows():
        seq_name = row['Sequence Name']
        sequence = row['Sequence'].strip().upper()
        well = row['Well Position']
        barcodes[sequence] = {
            'name': seq_name,
            'well': well,
            'type': 'row' if seq_name.startswith('R') else 'col'
        }
    return barcodes


def find_barcode_positions(seq, barcodes, max_mismatches=5):
    """
    Find positions of barcodes in a sequence.
    Returns list of (start, end, barcode_seq, barcode_info) tuples.
    """
    positions = []
    seq_upper = seq.upper()
    
    for barcode_seq, info in barcodes.items():
        # Search for exact matches first
        start = 0
        while True:
            pos = seq_upper.find(barcode_seq, start)
            if pos == -1:
                break
            positions.append((pos, pos + len(barcode_seq), barcode_seq, info))
            start = pos + 1
    
    return positions


def analyze_well_reads(well_fastq, barcodes):
    """Analyze reads from a single well."""
    reads_data = []
    
    for record in SeqIO.parse(well_fastq, "fastq"):
        seq = str(record.seq)
        barcode_positions = find_barcode_positions(seq, barcodes)
        
        reads_data.append({
            'id': record.id,
            'length': len(seq),
            'sequence': seq,
            'barcode_positions': barcode_positions,
            'num_barcodes': len(barcode_positions)
        })
    
    return reads_data


def plot_read_length_distribution(demux_dir, output_file):
    """Plot read length distribution across all wells."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    all_lengths = []
    well_lengths = defaultdict(list)
    
    # Collect all read lengths
    for fastq_file in os.listdir(demux_dir):
        if fastq_file.endswith('_reads.fastq'):
            well = fastq_file.replace('_reads.fastq', '')
            fastq_path = os.path.join(demux_dir, fastq_file)
            
            for record in SeqIO.parse(fastq_path, "fastq"):
                length = len(record.seq)
                all_lengths.append(length)
                well_lengths[well].append(length)
    
    if not all_lengths:
        print("⚠️  No reads found for length analysis")
        plt.close()
        return
    
    # Plot 1: Overall histogram
    ax1.hist(all_lengths, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    ax1.set_xlabel('Read Length (bp)', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title(f'Read Length Distribution\n(n={len(all_lengths)} reads)', fontsize=14, fontweight='bold')
    ax1.axvline(np.median(all_lengths), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(all_lengths):.0f} bp')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Box plot by well (top 10 wells by read count)
    top_wells = sorted(well_lengths.items(), key=lambda x: len(x[1]), reverse=True)[:10]
    if top_wells:
        well_names = [w[0] for w in top_wells]
        well_data = [w[1] for w in top_wells]
        
        bp = ax2.boxplot(well_data, tick_labels=well_names, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightblue')
            patch.set_alpha(0.7)
        
        ax2.set_xlabel('Well', fontsize=12)
        ax2.set_ylabel('Read Length (bp)', fontsize=12)
        ax2.set_title('Read Length by Well (Top 10)', fontsize=14, fontweight='bold')
        ax2.tick_params(axis='x', rotation=45)
        ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Read length distribution saved to {output_file}")
    
    # Return statistics
    return {
        'total_reads': len(all_lengths),
        'mean_length': np.mean(all_lengths),
        'median_length': np.median(all_lengths),
        'min_length': np.min(all_lengths),
        'max_length': np.max(all_lengths),
        'std_length': np.std(all_lengths)
    }


def plot_msa_layout(demux_dir, barcodes, output_file, max_reads=100):
    """
    Create Multiple Sequence Alignment-style visualization showing relative positions.
    Similar to SeqAn ReadLayout visualization.
    """
    fig, ax = plt.subplots(figsize=(16, 10))
    
    # Collect reads from all wells (limit to max_reads for visualization)
    all_reads = []
    well_colors = {}
    color_palette = sns.color_palette("husl", 96)  # 96 well plate
    color_idx = 0
    
    # Count fastq files first for efficient division
    fastq_files = [f for f in os.listdir(demux_dir) if f.endswith('_reads.fastq')]
    reads_per_well = max_reads // len(fastq_files) if fastq_files else max_reads
    
    for fastq_file in sorted(fastq_files):
        well = fastq_file.replace('_reads.fastq', '')
        well_colors[well] = color_palette[color_idx % len(color_palette)]
        color_idx += 1
        
        fastq_path = os.path.join(demux_dir, fastq_file)
        reads_data = analyze_well_reads(fastq_path, barcodes)
        
        for read_data in reads_data[:reads_per_well]:
            read_data['well'] = well
            all_reads.append(read_data)
    
    if not all_reads:
        print("⚠️  No reads found for MSA layout")
        plt.close()
        return
    
    # Sort reads by length for better visualization
    all_reads = sorted(all_reads, key=lambda x: x['length'], reverse=True)[:max_reads]
    
    # Plot each read as a horizontal bar
    y_position = 0
    y_spacing = 1.0
    max_length = max(read['length'] for read in all_reads)
    
    for read in all_reads:
        well = read['well']
        length = read['length']
        
        # Draw main read bar
        rect = Rectangle((0, y_position), length, 0.8, 
                         facecolor=well_colors[well], 
                         edgecolor='black', 
                         linewidth=0.5,
                         alpha=0.6)
        ax.add_patch(rect)
        
        # Draw barcode positions
        for start, end, barcode_seq, info in read['barcode_positions']:
            barcode_color = 'red' if info['type'] == 'row' else 'blue'
            barcode_rect = Rectangle((start, y_position), end - start, 0.8,
                                     facecolor=barcode_color,
                                     edgecolor='black',
                                     linewidth=1,
                                     alpha=0.9)
            ax.add_patch(barcode_rect)
        
        y_position += y_spacing
    
    # Set axis properties
    ax.set_xlim(0, max_length * 1.05)
    ax.set_ylim(-1, y_position + 1)
    ax.set_xlabel('Position (bp)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Reads', fontsize=12, fontweight='bold')
    ax.set_title(f'Multiple Sequence Alignment Layout\n(showing {len(all_reads)} reads with barcode positions)', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor='red', edgecolor='black', label='Row Barcode', alpha=0.9),
        mpatches.Patch(facecolor='blue', edgecolor='black', label='Column Barcode', alpha=0.9),
    ]
    
    # Add top wells to legend
    well_counts = defaultdict(int)
    for read in all_reads:
        well_counts[read['well']] += 1
    
    top_wells = sorted(well_counts.items(), key=lambda x: x[1], reverse=True)[:5]
    for well, count in top_wells:
        legend_elements.append(
            mpatches.Patch(facecolor=well_colors[well], edgecolor='black', 
                          label=f'Well {well} ({count} reads)', alpha=0.6)
        )
    
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ MSA layout visualization saved to {output_file}")


def plot_barcode_position_heatmap(demux_dir, barcodes, output_file, max_reads_per_well=50):
    """Plot heatmap showing barcode positions across reads."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Collect barcode position data
    all_barcode_positions = {'row': [], 'col': []}
    
    for fastq_file in os.listdir(demux_dir):
        if fastq_file.endswith('_reads.fastq'):
            fastq_path = os.path.join(demux_dir, fastq_file)
            reads_data = analyze_well_reads(fastq_path, barcodes)[:max_reads_per_well]
            
            for read_data in reads_data:
                for start, end, barcode_seq, info in read_data['barcode_positions']:
                    barcode_type = info['type']
                    # Normalize position as percentage of read length
                    rel_position = start / read_data['length'] if read_data['length'] > 0 else 0
                    all_barcode_positions[barcode_type].append(rel_position * 100)
    
    # Plot row barcode positions
    if all_barcode_positions['row']:
        ax1.hist(all_barcode_positions['row'], bins=50, color='red', alpha=0.6, edgecolor='black')
        ax1.set_xlabel('Relative Position in Read (%)', fontsize=12)
        ax1.set_ylabel('Count', fontsize=12)
        ax1.set_title(f'Row Barcode Position Distribution\n(n={len(all_barcode_positions["row"])} detections)', 
                     fontsize=13, fontweight='bold')
        ax1.axvline(0, color='green', linestyle='--', linewidth=2, alpha=0.7, label='Read Start')
        ax1.axvline(100, color='green', linestyle='--', linewidth=2, alpha=0.7, label='Read End')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
    
    # Plot column barcode positions
    if all_barcode_positions['col']:
        ax2.hist(all_barcode_positions['col'], bins=50, color='blue', alpha=0.6, edgecolor='black')
        ax2.set_xlabel('Relative Position in Read (%)', fontsize=12)
        ax2.set_ylabel('Count', fontsize=12)
        ax2.set_title(f'Column Barcode Position Distribution\n(n={len(all_barcode_positions["col"])} detections)', 
                     fontsize=13, fontweight='bold')
        ax2.axvline(0, color='green', linestyle='--', linewidth=2, alpha=0.7, label='Read Start')
        ax2.axvline(100, color='green', linestyle='--', linewidth=2, alpha=0.7, label='Read End')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Barcode position heatmap saved to {output_file}")


def plot_barcode_presence_summary(stats_csv, output_file):
    """Plot barcode presence summary from stats CSV."""
    # Read the stats CSV
    df = pd.read_csv(stats_csv, index_col=0)
    
    # Extract the 8x12 grid (rows A-H, columns 1-12)
    grid_data = df.iloc[:8, :12].astype(float)
    
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    # Plot 1: Heatmap of read counts per well
    ax1 = fig.add_subplot(gs[0, :])
    sns.heatmap(grid_data, annot=True, fmt='.0f', cmap='YlOrRd', 
                cbar_kws={'label': 'Read Count'}, ax=ax1,
                linewidths=0.5, linecolor='gray')
    ax1.set_title('Read Distribution Across 96-Well Plate', fontsize=14, fontweight='bold', pad=10)
    ax1.set_xlabel('Column', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Row', fontsize=12, fontweight='bold')
    
    # Plot 2: Bar chart of reads per row
    ax2 = fig.add_subplot(gs[1, 0])
    row_sums = grid_data.sum(axis=1)
    row_sums.plot(kind='bar', ax=ax2, color='steelblue', edgecolor='black')
    ax2.set_title('Reads per Row', fontsize=12, fontweight='bold')
    ax2.set_xlabel('Row', fontsize=11)
    ax2.set_ylabel('Total Reads', fontsize=11)
    ax2.tick_params(axis='x', rotation=0)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Plot 3: Bar chart of reads per column
    ax3 = fig.add_subplot(gs[1, 1])
    col_sums = grid_data.sum(axis=0)
    col_sums.plot(kind='bar', ax=ax3, color='coral', edgecolor='black')
    ax3.set_title('Reads per Column', fontsize=12, fontweight='bold')
    ax3.set_xlabel('Column', fontsize=11)
    ax3.set_ylabel('Total Reads', fontsize=11)
    ax3.tick_params(axis='x', rotation=0)
    ax3.grid(True, alpha=0.3, axis='y')
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Barcode presence summary saved to {output_file}")


def align_read_to_anchor(read_seq, anchor_seq):
    """
    Align a read sequence to an anchor sequence using Smith-Waterman local alignment.
    
    Parameters:
        read_seq (str): Read sequence
        anchor_seq (str): Anchor/reference sequence
    
    Returns:
        dict: Alignment information including position, score, and aligned sequences
    """
    # Create alignment matrix (match=2, mismatch=-1)
    matrix = parasail.matrix_create("ACGT", 2, -1)
    
    # Perform Smith-Waterman alignment with traceback
    result = parasail.sw_trace_scan_16(read_seq, anchor_seq, ALIGNMENT_GAP_OPEN, ALIGNMENT_GAP_EXTEND, matrix)
    
    # Get traceback
    traceback = result.get_traceback()
    
    # Calculate begin positions from end positions and lengths
    # In parasail, end positions are 0-based inclusive
    ref_end = result.end_ref + 1  # Convert to exclusive end
    query_end = result.end_query + 1  # Convert to exclusive end
    ref_begin = max(0, ref_end - result.len_ref)
    query_begin = max(0, query_end - result.len_query)
    
    return {
        'score': result.score,
        'ref_begin': ref_begin,
        'ref_end': ref_end,
        'read_begin': query_begin,
        'read_end': query_end,
        'length': result.len_ref,
        'traceback': traceback
    }


def plot_anchor_msa(demux_dir, anchor_seq, output_file, max_reads=100):
    """
    Create MSA visualization showing alignment of reads to an anchor sequence.
    
    Parameters:
        demux_dir (str): Directory containing demultiplexed FASTQ files
        anchor_seq (str): Anchor/reference sequence to align reads to
        output_file (str): Output file path for the plot
        max_reads (int): Maximum number of reads to visualize
    """
    fig, ax = plt.subplots(figsize=(16, 12))
    
    # Collect reads from all wells
    all_reads = []
    well_colors = {}
    color_palette = sns.color_palette("husl", 96)
    color_idx = 0
    
    # Count fastq files for efficient division
    fastq_files = [f for f in os.listdir(demux_dir) if f.endswith('_reads.fastq')]
    reads_per_well = max(1, max_reads // len(fastq_files)) if fastq_files else max_reads
    
    for fastq_file in sorted(fastq_files):
        well = fastq_file.replace('_reads.fastq', '')
        well_colors[well] = color_palette[color_idx % len(color_palette)]
        color_idx += 1
        
        fastq_path = os.path.join(demux_dir, fastq_file)
        
        read_count = 0
        for record in SeqIO.parse(fastq_path, "fastq"):
            if read_count >= reads_per_well:
                break
            
            seq = str(record.seq).upper()
            
            # Align read to anchor
            alignment_info = align_read_to_anchor(seq, anchor_seq.upper())
            
            all_reads.append({
                'id': record.id,
                'well': well,
                'length': len(seq),
                'alignment': alignment_info
            })
            read_count += 1
    
    if not all_reads:
        print("⚠️  No reads found for anchor MSA")
        plt.close()
        return {
            'num_aligned': 0,
            'mean_score': 0,
            'median_score': 0
        }
    
    # Limit to max_reads
    all_reads = all_reads[:max_reads]
    
    # Sort reads by alignment position on anchor for better visualization
    all_reads.sort(key=lambda x: x['alignment']['ref_begin'])
    
    # Plot anchor sequence as a reference bar at the top
    anchor_length = len(anchor_seq)
    ax.add_patch(Rectangle((0, -2), anchor_length, 1.0, 
                           facecolor='gold', 
                           edgecolor='black', 
                           linewidth=2,
                           alpha=0.8,
                           label='Anchor Sequence'))
    
    # Plot each read aligned to the anchor
    y_position = 0
    y_spacing = 1.0
    max_ref_pos = anchor_length
    
    for read in all_reads:
        well = read['well']
        alignment = read['alignment']
        
        # Only plot if there's a reasonable alignment
        if alignment['score'] > 0:
            # Calculate read position relative to anchor
            ref_start = alignment['ref_begin']
            ref_end = alignment['ref_end']
            
            # Draw aligned portion
            aligned_length = ref_end - ref_start
            rect = Rectangle((ref_start, y_position), aligned_length, 0.8,
                           facecolor=well_colors[well],
                           edgecolor='black',
                           linewidth=0.5,
                           alpha=0.7)
            ax.add_patch(rect)
            
            # Update max position if needed
            max_ref_pos = max(max_ref_pos, ref_end)
        
        y_position += y_spacing
    
    # Set axis properties
    ax.set_xlim(-10, max_ref_pos + 10)
    ax.set_ylim(-3, y_position + 1)
    ax.set_xlabel('Position on Anchor Sequence (bp)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Reads', fontsize=12, fontweight='bold')
    ax.set_title(f'Multiple Sequence Alignment to Anchor\n(showing {len(all_reads)} reads aligned to {len(anchor_seq)} bp anchor)', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor='gold', edgecolor='black', label='Anchor Sequence', alpha=0.8)
    ]
    
    # Add top wells to legend
    well_counts = defaultdict(int)
    for read in all_reads:
        well_counts[read['well']] += 1
    
    top_wells = sorted(well_counts.items(), key=lambda x: x[1], reverse=True)[:5]
    for well, count in top_wells:
        legend_elements.append(
            mpatches.Patch(facecolor=well_colors[well], edgecolor='black',
                          label=f'Well {well} ({count} reads)', alpha=0.7)
        )
    
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3, axis='x')
    
    # Add a horizontal line to separate anchor from reads
    ax.axhline(y=-1, color='black', linestyle='--', linewidth=1, alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Anchor MSA visualization saved to {output_file}")
    
    # Return alignment statistics
    alignment_scores = [r['alignment']['score'] for r in all_reads if r['alignment']['score'] > 0]
    return {
        'num_aligned': len(alignment_scores),
        'mean_score': np.mean(alignment_scores) if alignment_scores else 0,
        'median_score': np.median(alignment_scores) if alignment_scores else 0
    }


def generate_html_report(demux_dir, report_dir, stats, anchor_stats=None):
    """Generate HTML report combining all visualizations."""
    
    # Build anchor section if anchor stats are provided
    anchor_section = ""
    if anchor_stats:
        anchor_section = f"""
        <h2>🎯 Anchor Sequence Alignment</h2>
        <p>Multiple Sequence Alignment of reads to the provided anchor sequence. The gold bar at the top represents the anchor/reference sequence, and each colored bar below shows where a read aligns to the anchor. Reads are sorted by alignment position for clarity.</p>
        
        <div class="stats-grid">
            <div class="stat-card">
                <h3>ALIGNED READS</h3>
                <div class="value">{anchor_stats.get('num_aligned', 0):,}</div>
                <div class="unit">reads</div>
            </div>
            <div class="stat-card">
                <h3>MEAN ALIGNMENT SCORE</h3>
                <div class="value">{anchor_stats.get('mean_score', 0):.1f}</div>
                <div class="unit">score</div>
            </div>
            <div class="stat-card">
                <h3>MEDIAN ALIGNMENT SCORE</h3>
                <div class="value">{anchor_stats.get('median_score', 0):.1f}</div>
                <div class="unit">score</div>
            </div>
        </div>
        
        <img src="anchor_msa.png" alt="Anchor MSA">
        """
    
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>NanoDemux Quality Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            border-left: 4px solid #3498db;
            padding-left: 10px;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        .stat-card h3 {{
            margin: 0 0 10px 0;
            font-size: 14px;
            opacity: 0.9;
        }}
        .stat-card .value {{
            font-size: 32px;
            font-weight: bold;
            margin: 0;
        }}
        .stat-card .unit {{
            font-size: 14px;
            opacity: 0.8;
        }}
        img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            margin: 20px 0;
        }}
        .metadata {{
            background-color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            margin: 20px 0;
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 2px solid #ecf0f1;
            text-align: center;
            color: #7f8c8d;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 NanoDemux Quality Report</h1>
        
        <div class="metadata">
            <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p><strong>Demux Directory:</strong> {os.path.basename(demux_dir)}</p>
        </div>
        
        <h2>📊 Summary Statistics</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <h3>TOTAL READS</h3>
                <div class="value">{stats.get('total_reads', 0):,}</div>
                <div class="unit">reads</div>
            </div>
            <div class="stat-card">
                <h3>MEAN LENGTH</h3>
                <div class="value">{stats.get('mean_length', 0):.0f}</div>
                <div class="unit">bp</div>
            </div>
            <div class="stat-card">
                <h3>MEDIAN LENGTH</h3>
                <div class="value">{stats.get('median_length', 0):.0f}</div>
                <div class="unit">bp</div>
            </div>
            <div class="stat-card">
                <h3>LENGTH RANGE</h3>
                <div class="value">{stats.get('min_length', 0):.0f} - {stats.get('max_length', 0):.0f}</div>
                <div class="unit">bp</div>
            </div>
        </div>
        
        {anchor_section}
        
        <h2>📏 Read Length Distribution</h2>
        <p>Distribution of read lengths across all wells, showing both overall histogram and per-well box plots for the top 10 wells by read count.</p>
        <img src="read_length_distribution.png" alt="Read Length Distribution">
        
        <h2>🧬 Multiple Sequence Alignment Layout</h2>
        <p>Visualization of read positions and lengths, similar to SeqAn ReadLayout. Each horizontal bar represents a read, colored by well. Red and blue regions indicate row and column barcode positions within the reads.</p>
        <img src="msa_layout.png" alt="MSA Layout">
        
        <h2>📍 Barcode Position Analysis</h2>
        <p>Distribution of barcode positions along reads, showing where row and column barcodes are typically found. Positions are shown as percentage of read length (0% = start, 100% = end).</p>
        <img src="barcode_positions.png" alt="Barcode Positions">
        
        <h2>🎯 Barcode Presence Summary</h2>
        <p>Heatmap showing the distribution of successfully demultiplexed reads across the 96-well plate, with row and column summaries.</p>
        <img src="barcode_summary.png" alt="Barcode Summary">
        
        <div class="footer">
            <p>Generated by NanoDemux Quality Report Generator</p>
            <p>Carothers Research Lab</p>
        </div>
    </div>
</body>
</html>
"""
    
    html_file = os.path.join(report_dir, 'quality_report.html')
    with open(html_file, 'w') as f:
        f.write(html_content)
    
    print(f"✅ HTML report saved to {html_file}")
    return html_file


def main():
    parser = argparse.ArgumentParser(
        description="Generate graphical data quality report for demultiplexed reads.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate report for a demuxed directory
  python generate_quality_report.py demplex_data/55XPXK_1_P4_323_EG/ \\
      barcodes/251202_primer_well_map_DA.csv

  # Generate report with custom output directory
  python generate_quality_report.py demplex_data/experiment1/sample1/ \\
      barcodes/primer_well_map.csv --output reports/sample1/
  
  # Generate report with anchor sequence MSA
  python generate_quality_report.py demplex_data/55XPXK_1_P4_323_EG/ \\
      barcodes/251202_primer_well_map_DA.csv \\
      --anchor-seq AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC
        """
    )
    
    parser.add_argument('demux_dir', help='Directory containing demultiplexed FASTQ files and barcode_stats.csv')
    parser.add_argument('barcodes', help='CSV file with barcode definitions (same as used for demultiplexing)')
    parser.add_argument('--output', '-o', default=None, 
                       help='Output directory for report (default: demux_dir/quality_report/)')
    parser.add_argument('--max-reads', type=int, default=100,
                       help='Maximum number of reads to show in MSA layout (default: 100)')
    parser.add_argument('--anchor-seq', default=None,
                       help='Anchor/reference sequence for MSA alignment (optional)')
    
    args = parser.parse_args()
    
    # Validate input
    if not os.path.isdir(args.demux_dir):
        print(f"❌ Error: Directory not found: {args.demux_dir}")
        sys.exit(1)
    
    if not os.path.isfile(args.barcodes):
        print(f"❌ Error: Barcode file not found: {args.barcodes}")
        sys.exit(1)
    
    stats_csv = os.path.join(args.demux_dir, 'barcode_stats.csv')
    if not os.path.isfile(stats_csv):
        print(f"❌ Error: barcode_stats.csv not found in {args.demux_dir}")
        sys.exit(1)
    
    # Set output directory
    if args.output:
        report_dir = args.output
    else:
        report_dir = os.path.join(args.demux_dir, 'quality_report')
    
    os.makedirs(report_dir, exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"NanoDemux Quality Report Generator")
    print(f"{'='*60}")
    print(f"Input directory: {args.demux_dir}")
    print(f"Barcode file: {args.barcodes}")
    print(f"Output directory: {report_dir}")
    if args.anchor_seq:
        print(f"Anchor sequence: {args.anchor_seq[:50]}{'...' if len(args.anchor_seq) > 50 else ''} ({len(args.anchor_seq)} bp)")
    print(f"{'='*60}\n")
    
    # Load barcodes
    barcodes = load_barcode_map(args.barcodes)
    print(f"Loaded {len(barcodes)} barcode sequences\n")
    
    # Generate visualizations
    print("Generating visualizations...")
    
    # 1. Read length distribution
    stats = plot_read_length_distribution(
        args.demux_dir,
        os.path.join(report_dir, 'read_length_distribution.png')
    )
    
    # 2. Anchor MSA (if anchor sequence provided)
    anchor_stats = None
    if args.anchor_seq:
        print("\n🎯 Generating anchor sequence MSA...")
        anchor_stats = plot_anchor_msa(
            args.demux_dir,
            args.anchor_seq,
            os.path.join(report_dir, 'anchor_msa.png'),
            max_reads=args.max_reads
        )
    
    # 3. MSA layout
    plot_msa_layout(
        args.demux_dir,
        barcodes,
        os.path.join(report_dir, 'msa_layout.png'),
        max_reads=args.max_reads
    )
    
    # 4. Barcode position heatmap
    plot_barcode_position_heatmap(
        args.demux_dir,
        barcodes,
        os.path.join(report_dir, 'barcode_positions.png')
    )
    
    # 5. Barcode presence summary
    plot_barcode_presence_summary(
        stats_csv,
        os.path.join(report_dir, 'barcode_summary.png')
    )
    
    # 6. Generate HTML report
    html_file = generate_html_report(args.demux_dir, report_dir, stats, anchor_stats)
    
    print(f"\n{'='*60}")
    print(f"✅ Report generation complete!")
    print(f"{'='*60}")
    print(f"📄 Open the report: {html_file}")
    print(f"{'='*60}\n")


if __name__ == '__main__':
    main()
