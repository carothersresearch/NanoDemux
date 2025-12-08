#!/usr/bin/env python3
"""
Generate graphical quality report for raw (unmultiplexed) FASTQ data.

This script creates comprehensive visualizations for raw nanopore sequencing data including:
- Read length distributions
- Quality score distributions
- Base composition analysis
- N50 and general statistics
"""

import argparse
import os
import sys
from collections import defaultdict, Counter
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from datetime import datetime

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300


def calculate_n50(lengths):
    """Calculate N50 from a list of read lengths."""
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    cumsum = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= total_length / 2:
            return length
    return 0


def analyze_fastq(fastq_file, max_reads=None):
    """
    Analyze raw FASTQ file and extract statistics.
    
    Parameters:
        fastq_file: Path to FASTQ file
        max_reads: Maximum number of reads to analyze (None = all reads)
    
    Returns:
        Dictionary with statistics and data for plotting
    """
    lengths = []
    qualities = []  # Average quality per read
    base_counts = Counter()
    quality_by_position = defaultdict(list)
    
    read_count = 0
    for record in SeqIO.parse(fastq_file, "fastq"):
        read_count += 1
        if max_reads and read_count > max_reads:
            break
        
        seq = str(record.seq)
        lengths.append(len(seq))
        
        # Get quality scores
        if hasattr(record, 'letter_annotations') and 'phred_quality' in record.letter_annotations:
            quals = record.letter_annotations['phred_quality']
            qualities.append(np.mean(quals))
            
            # Quality by position (for the first 1000 bases to keep it manageable)
            for i, q in enumerate(quals[:1000]):
                quality_by_position[i].append(q)
        
        # Base composition
        for base in seq.upper():
            base_counts[base] += 1
    
    # Calculate statistics
    stats = {
        'total_reads': len(lengths),
        'total_bases': sum(lengths),
        'mean_length': np.mean(lengths) if lengths else 0,
        'median_length': np.median(lengths) if lengths else 0,
        'min_length': np.min(lengths) if lengths else 0,
        'max_length': np.max(lengths) if lengths else 0,
        'n50': calculate_n50(lengths),
        'mean_quality': np.mean(qualities) if qualities else 0,
        'median_quality': np.median(qualities) if qualities else 0,
    }
    
    return {
        'stats': stats,
        'lengths': lengths,
        'qualities': qualities,
        'base_counts': base_counts,
        'quality_by_position': quality_by_position
    }


def plot_read_length_distribution(data, output_file):
    """Plot read length distribution."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    lengths = data['lengths']
    stats = data['stats']
    
    if not lengths:
        print("⚠️  No reads found for length analysis")
        plt.close()
        return
    
    # Plot 1: Histogram
    ax1.hist(lengths, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    ax1.set_xlabel('Read Length (bp)', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title(f'Read Length Distribution\n(n={stats["total_reads"]:,} reads)', 
                  fontsize=14, fontweight='bold')
    ax1.axvline(stats['median_length'], color='red', linestyle='--', linewidth=2, 
                label=f'Median: {stats["median_length"]:.0f} bp')
    ax1.axvline(stats['n50'], color='orange', linestyle='--', linewidth=2, 
                label=f'N50: {stats["n50"]:.0f} bp')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Box plot
    bp = ax2.boxplot([lengths], patch_artist=True, orientation='vertical')
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
        patch.set_alpha(0.7)
    
    ax2.set_ylabel('Read Length (bp)', fontsize=12)
    ax2.set_title('Read Length Box Plot', fontsize=14, fontweight='bold')
    ax2.set_xticklabels(['All Reads'])
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Add text annotations
    textstr = f'Mean: {stats["mean_length"]:.0f} bp\n'
    textstr += f'Median: {stats["median_length"]:.0f} bp\n'
    textstr += f'N50: {stats["n50"]:.0f} bp\n'
    textstr += f'Min: {stats["min_length"]:.0f} bp\n'
    textstr += f'Max: {stats["max_length"]:.0f} bp'
    ax2.text(0.02, 0.98, textstr, transform=ax2.transAxes, fontsize=10,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Read length distribution saved to {output_file}")


def plot_quality_distribution(data, output_file):
    """Plot quality score distributions."""
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    qualities = data['qualities']
    quality_by_position = data['quality_by_position']
    stats = data['stats']
    
    if not qualities:
        print("⚠️  No quality data found")
        plt.close()
        return
    
    # Plot 1: Histogram of mean quality per read
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(qualities, bins=50, edgecolor='black', alpha=0.7, color='green')
    ax1.set_xlabel('Mean Quality Score (Phred)', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title(f'Mean Quality Score per Read\n(n={len(qualities):,} reads)', 
                  fontsize=13, fontweight='bold')
    ax1.axvline(stats['mean_quality'], color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {stats["mean_quality"]:.1f}')
    ax1.axvline(stats['median_quality'], color='orange', linestyle='--', linewidth=2, 
                label=f'Median: {stats["median_quality"]:.1f}')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Quality by position
    ax2 = fig.add_subplot(gs[0, 1])
    if quality_by_position:
        positions = sorted(quality_by_position.keys())
        mean_quals = [np.mean(quality_by_position[pos]) for pos in positions]
        q25_quals = [np.percentile(quality_by_position[pos], 25) for pos in positions]
        q75_quals = [np.percentile(quality_by_position[pos], 75) for pos in positions]
        
        ax2.plot(positions, mean_quals, color='blue', linewidth=2, label='Mean')
        ax2.fill_between(positions, q25_quals, q75_quals, alpha=0.3, color='blue', label='IQR')
        ax2.set_xlabel('Position in Read (bp)', fontsize=12)
        ax2.set_ylabel('Quality Score (Phred)', fontsize=12)
        ax2.set_title('Quality Score by Position\n(first 1000 bp)', fontsize=13, fontweight='bold')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.axhline(20, color='red', linestyle='--', linewidth=1, alpha=0.5, label='Q20')
        ax2.axhline(30, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='Q30')
    
    # Plot 3: Box plot of quality
    ax3 = fig.add_subplot(gs[1, 0])
    bp = ax3.boxplot([qualities], patch_artist=True, orientation='vertical')
    for patch in bp['boxes']:
        patch.set_facecolor('lightgreen')
        patch.set_alpha(0.7)
    
    ax3.set_ylabel('Mean Quality Score (Phred)', fontsize=12)
    ax3.set_title('Quality Score Box Plot', fontsize=13, fontweight='bold')
    ax3.set_xticklabels(['All Reads'])
    ax3.grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Quality distribution heatmap
    ax4 = fig.add_subplot(gs[1, 1])
    if quality_by_position:
        # Create heatmap data
        positions_subset = sorted(quality_by_position.keys())[::10]  # Sample every 10th position
        quality_matrix = []
        for pos in positions_subset:
            hist, _ = np.histogram(quality_by_position[pos], bins=range(0, 42, 2))
            quality_matrix.append(hist)
        
        quality_matrix = np.array(quality_matrix).T
        
        im = ax4.imshow(quality_matrix, aspect='auto', cmap='YlGnBu', origin='lower')
        ax4.set_xlabel('Position in Read (bp)', fontsize=12)
        ax4.set_ylabel('Quality Score (Phred)', fontsize=12)
        ax4.set_title('Quality Score Distribution by Position', fontsize=13, fontweight='bold')
        
        # Set tick labels
        ax4.set_yticks(np.arange(0, len(range(0, 42, 2)), 5))
        ax4.set_yticklabels([str(i) for i in range(0, 42, 10)])
        
        fig.colorbar(im, ax=ax4, label='Count')
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Quality distribution saved to {output_file}")


def plot_base_composition(data, output_file):
    """Plot base composition analysis."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    base_counts = data['base_counts']
    stats = data['stats']
    
    if not base_counts:
        print("⚠️  No base composition data found")
        plt.close()
        return
    
    # Get counts for standard bases
    standard_bases = ['A', 'T', 'G', 'C']
    counts = [base_counts.get(base, 0) for base in standard_bases]
    total_standard = sum(counts)
    
    # Calculate GC content
    gc_content = (base_counts.get('G', 0) + base_counts.get('C', 0)) / total_standard * 100 if total_standard > 0 else 0
    
    # Plot 1: Bar chart
    colors = ['#FF6B6B', '#4ECDC4', '#FFD93D', '#6BCB77']
    bars = ax1.bar(standard_bases, counts, color=colors, edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Base', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax1.set_title(f'Base Composition\n(Total bases: {stats["total_bases"]:,})', 
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Add percentage labels on bars
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        percentage = count / total_standard * 100 if total_standard > 0 else 0
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{percentage:.1f}%',
                ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Plot 2: Pie chart
    colors_pie = colors[:len(counts)]
    wedges, texts, autotexts = ax2.pie(counts, labels=standard_bases, colors=colors_pie, 
                                         autopct='%1.1f%%', startangle=90,
                                         explode=[0.05]*len(counts))
    ax2.set_title('Base Composition Proportions', fontsize=14, fontweight='bold')
    
    # Make percentage text bold
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(11)
    
    # Add GC content annotation
    textstr = f'GC Content: {gc_content:.1f}%'
    ax2.text(0.5, -1.3, textstr, transform=ax2.transAxes, fontsize=12,
             horizontalalignment='center', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Base composition saved to {output_file}")


def generate_html_report(fastq_file, report_dir, data):
    """Generate HTML report combining all visualizations."""
    stats = data['stats']
    base_counts = data['base_counts']
    
    # Calculate GC content
    standard_bases = ['A', 'T', 'G', 'C']
    counts = [base_counts.get(base, 0) for base in standard_bases]
    total_standard = sum(counts)
    gc_content = (base_counts.get('G', 0) + base_counts.get('C', 0)) / total_standard * 100 if total_standard > 0 else 0
    
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>NanoDemux Raw Data Quality Report</title>
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
        <h1>🧬 NanoDemux Raw Data Quality Report</h1>
        
        <div class="metadata">
            <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p><strong>Input File:</strong> {os.path.basename(fastq_file)}</p>
        </div>
        
        <h2>📊 Summary Statistics</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <h3>TOTAL READS</h3>
                <div class="value">{stats['total_reads']:,}</div>
                <div class="unit">reads</div>
            </div>
            <div class="stat-card">
                <h3>TOTAL BASES</h3>
                <div class="value">{stats['total_bases'] / 1e6:.1f}</div>
                <div class="unit">Mb</div>
            </div>
            <div class="stat-card">
                <h3>MEAN LENGTH</h3>
                <div class="value">{stats['mean_length']:.0f}</div>
                <div class="unit">bp</div>
            </div>
            <div class="stat-card">
                <h3>MEDIAN LENGTH</h3>
                <div class="value">{stats['median_length']:.0f}</div>
                <div class="unit">bp</div>
            </div>
            <div class="stat-card">
                <h3>N50</h3>
                <div class="value">{stats['n50']:.0f}</div>
                <div class="unit">bp</div>
            </div>
            <div class="stat-card">
                <h3>MEAN QUALITY</h3>
                <div class="value">{stats['mean_quality']:.1f}</div>
                <div class="unit">Phred</div>
            </div>
            <div class="stat-card">
                <h3>GC CONTENT</h3>
                <div class="value">{gc_content:.1f}</div>
                <div class="unit">%</div>
            </div>
            <div class="stat-card">
                <h3>LENGTH RANGE</h3>
                <div class="value">{stats['min_length']:.0f} - {stats['max_length']:.0f}</div>
                <div class="unit">bp</div>
            </div>
        </div>
        
        <h2>📏 Read Length Distribution</h2>
        <p>Distribution of read lengths across all reads in the raw FASTQ file. The histogram shows the frequency distribution, while the box plot provides quartile information. N50 is a weighted median that represents the length at which 50% of the total bases are in reads of that length or longer.</p>
        <img src="read_length_distribution.png" alt="Read Length Distribution">
        
        <h2>📊 Quality Score Analysis</h2>
        <p>Quality score distributions showing mean quality per read and quality by position in reads. Higher Phred scores indicate higher base-calling confidence. Q20 (99% accuracy) and Q30 (99.9% accuracy) thresholds are commonly used benchmarks.</p>
        <img src="quality_distribution.png" alt="Quality Distribution">
        
        <h2>🧬 Base Composition</h2>
        <p>Analysis of nucleotide base composition (A, T, G, C) across all reads. Balanced base composition and appropriate GC content are indicators of good sequencing quality.</p>
        <img src="base_composition.png" alt="Base Composition">
        
        <div class="footer">
            <p>Generated by NanoDemux Raw Data Quality Report Generator</p>
            <p>Carothers Research Lab</p>
        </div>
    </div>
</body>
</html>
"""
    
    html_file = os.path.join(report_dir, 'raw_quality_report.html')
    with open(html_file, 'w') as f:
        f.write(html_content)
    
    print(f"✅ HTML report saved to {html_file}")
    return html_file


def main():
    parser = argparse.ArgumentParser(
        description="Generate graphical quality report for raw (unmultiplexed) FASTQ data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate report for a raw FASTQ file
  python generate_raw_quality_report.py raw_data/sample.fastq

  # Generate report with custom output directory
  python generate_raw_quality_report.py raw_data/sample.fastq --output reports/sample/

  # Analyze only first 10000 reads for faster processing
  python generate_raw_quality_report.py raw_data/large_sample.fastq --max-reads 10000
        """
    )
    
    parser.add_argument('fastq', help='Input FASTQ file (raw, unmultiplexed data)')
    parser.add_argument('--output', '-o', default=None, 
                       help='Output directory for report (default: raw_data_quality_report/)')
    parser.add_argument('--max-reads', type=int, default=None,
                       help='Maximum number of reads to analyze (default: all reads)')
    
    args = parser.parse_args()
    
    # Validate input
    if not os.path.isfile(args.fastq):
        print(f"❌ Error: File not found: {args.fastq}")
        sys.exit(1)
    
    # Set output directory
    if args.output:
        report_dir = args.output
    else:
        report_dir = 'raw_data_quality_report'
    
    os.makedirs(report_dir, exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"NanoDemux Raw Data Quality Report Generator")
    print(f"{'='*60}")
    print(f"Input file: {args.fastq}")
    print(f"Output directory: {report_dir}")
    if args.max_reads:
        print(f"Max reads to analyze: {args.max_reads:,}")
    print(f"{'='*60}\n")
    
    # Analyze FASTQ
    print("Analyzing FASTQ file...")
    data = analyze_fastq(args.fastq, max_reads=args.max_reads)
    
    print(f"Analyzed {data['stats']['total_reads']:,} reads")
    print(f"Total bases: {data['stats']['total_bases']:,}")
    print(f"Mean length: {data['stats']['mean_length']:.0f} bp")
    print(f"N50: {data['stats']['n50']:.0f} bp")
    print(f"Mean quality: {data['stats']['mean_quality']:.1f}\n")
    
    # Generate visualizations
    print("Generating visualizations...")
    
    # 1. Read length distribution
    plot_read_length_distribution(
        data,
        os.path.join(report_dir, 'read_length_distribution.png')
    )
    
    # 2. Quality distribution
    plot_quality_distribution(
        data,
        os.path.join(report_dir, 'quality_distribution.png')
    )
    
    # 3. Base composition
    plot_base_composition(
        data,
        os.path.join(report_dir, 'base_composition.png')
    )
    
    # 4. Generate HTML report
    html_file = generate_html_report(args.fastq, report_dir, data)
    
    print(f"\n{'='*60}")
    print(f"✅ Report generation complete!")
    print(f"{'='*60}")
    print(f"📄 Open the report: {html_file}")
    print(f"{'='*60}\n")


if __name__ == '__main__':
    main()
