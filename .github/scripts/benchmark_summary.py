#!/usr/bin/env python3
"""
Extract benchmark summary for GitHub commit comments.
"""

import json
import sys
from pathlib import Path


def format_benchmark_summary():
    """Generate a concise benchmark summary for commit comments."""
    results_file = Path("benchmarking_data/benchmarks/benchmark_results.json")
    
    if not results_file.exists():
        return "## 📊 Benchmark Results\n\nNo benchmark results found."
    
    with open(results_file, 'r') as f:
        history = json.load(f)
    
    if not history:
        return "## 📊 Benchmark Results\n\nNo benchmark runs completed."
    
    # Get latest runs
    latest = history[-2:] if len(history) >= 2 else history
    
    summary = ["## 📊 Benchmark Results"]
    summary.append("")
    summary.append("### Latest Run")
    summary.append("")
    
    # Format table header
    summary.append("| Dataset | Subset | Mapped | Mapping % | Wells | Time (s) |")
    summary.append("|---------|--------|--------|-----------|-------|----------|")
    
    # Add latest results
    for result in latest:
        stats = result['statistics']
        subset_info = result.get('subset_reads', 'full')
        summary.append(
            f"| {result['name']:<7} | {subset_info:>6} | "
            f"{stats['mapped_reads']:>6,} | {stats['mapping_rate']:>8.2f}% | "
            f"{stats['wells_with_reads']:>5} | {result['execution_time']:>8.2f} |"
        )
    
    # Add comparison if we have history
    if len(history) >= 2:
        first = history[-2]['statistics']
        last = history[-1]['statistics']
        
        mapped_change = last['mapped_reads'] - first['mapped_reads']
        rate_change = last['mapping_rate'] - first['mapping_rate']
        
        summary.append("")
        summary.append("### Change Since Previous Run")
        summary.append(f"- Mapped reads: {mapped_change:+,}")
        summary.append(f"- Mapping rate: {rate_change:+.2f}%")
    
    summary.append("")
    summary.append("---")
    summary.append(f"*Total benchmark runs: {len(history)}*")
    summary.append("")
    summary.append("📦 [Download full report from Actions artifacts](../../actions)")
    
    return "\n".join(summary)


if __name__ == "__main__":
    print(format_benchmark_summary())
