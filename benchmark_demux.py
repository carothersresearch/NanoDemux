#!/usr/bin/env python3
"""
Benchmarking suite for NanoDemux demultiplexing performance.

This script runs the demultiplexing pipeline on sample data and tracks
key performance metrics to measure improvements over time.

Configuration files (benchmark_config.json or benchmark_config_*.json) in each
subfolder define the parameters for running benchmarks on that data.
"""

import os
import sys
import json
import time
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
import pandas as pd
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO


class DemuxBenchmark:
    """Benchmark the demultiplexing pipeline."""
    
    def __init__(self, benchmark_dir: str = "benchmarking_data/benchmarks"):
        self.benchmark_dir = Path(benchmark_dir)
        self.benchmark_dir.mkdir(exist_ok=True, parents=True)
        self.results_file = self.benchmark_dir / "benchmark_results.json"
        self.history = self._load_history()
        self.subset_files = []  # Track subset files for cleanup
    
    def _load_history(self) -> List[Dict]:
        """Load historical benchmark results."""
        if self.results_file.exists():
            with open(self.results_file, 'r') as f:
                return json.load(f)
        return []
    
    def _save_history(self):
        """Save benchmark results to file."""
        with open(self.results_file, 'w') as f:
            json.dump(self.history, f, indent=2)
    
    def create_subset(self, fastq_file: str, num_reads: int) -> str:
        """
        Create a subset of a FASTQ file for faster benchmarking.
        
        Args:
            fastq_file: Path to input FASTQ file
            num_reads: Number of reads to extract
        
        Returns:
            Path to subset FASTQ file
        """
        subset_file = self.benchmark_dir / f"subset_{Path(fastq_file).stem}_{num_reads}.fastq"
        
        # If subset already exists and is recent, reuse it
        if subset_file.exists():
            self.subset_files.append(subset_file)
            return str(subset_file)
        
        print(f"Creating subset: {num_reads} reads from {Path(fastq_file).name}")
        
        with open(fastq_file, 'r') as inf, open(subset_file, 'w') as outf:
            count = 0
            for record in SeqIO.parse(inf, "fastq"):
                SeqIO.write(record, outf, "fastq")
                count += 1
                if count >= num_reads:
                    break
        
        print(f"✓ Created subset with {count} reads")
        self.subset_files.append(subset_file)
        return str(subset_file)
    
    def cleanup_output_fastqs(self, outdir: str):
        """
        Remove FASTQ output files, keeping only the stats CSV.
        
        Args:
            outdir: Output directory containing FASTQ files
        """
        outdir_path = Path(outdir)
        fastq_files = list(outdir_path.glob("*.fastq"))
        
        if fastq_files:
            print(f"Cleaning up {len(fastq_files)} FASTQ files...")
            for fastq_file in fastq_files:
                fastq_file.unlink()
            print(f"✓ Removed FASTQ outputs (kept stats CSV)")
    
    def cleanup_subset_files(self):
        """
        Remove temporary subset FASTQ files created during benchmarking.
        """
        if self.subset_files:
            print(f"\nCleaning up {len(self.subset_files)} subset files...")
            for subset_file in self.subset_files:
                if subset_file.exists():
                    subset_file.unlink()
            print(f"✓ Removed subset FASTQ files")
            self.subset_files = []
    
    def run_demux(self, fastq_file: str, barcode_csv: str, 
                  outdir: str, adapter_file: str = None, **kwargs) -> Tuple[float, str]:
        """
        Run demultiplexing and measure execution time.
        
        Returns:
            Tuple of (execution_time, output_directory)
        """
        # Build command
        cmd = [
            sys.executable, "demux_barcodes.py",
            fastq_file, barcode_csv,
            "--outdir", outdir
        ]
        
        if adapter_file:
            cmd.extend(["--adapters", adapter_file])
        
        # Add optional parameters
        for key, value in kwargs.items():
            cmd.extend([f"--{key}", str(value)])
        
        # Run and time
        print(f"Running: {' '.join(cmd)}")
        start_time = time.time()
        
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True,
            cwd=Path(__file__).parent
        )
        
        elapsed = time.time() - start_time
        
        if result.returncode != 0:
            raise RuntimeError(f"Demux failed: {result.stderr}")
        
        return elapsed, outdir
    
    def extract_stats(self, stats_file: str) -> Dict:
        """
        Extract key statistics from barcode_stats.csv.
        
        Returns dictionary with metrics:
        - total_reads: Total reads processed
        - mapped_reads: Reads assigned to wells (both barcodes)
        - single_barcode: Reads with only one barcode
        - no_match: Reads with no barcodes
        - mapping_rate: Percentage of mapped reads
        - well_distribution: Per-well read counts
        """
        df = pd.read_csv(stats_file, index_col=0)
        
        # Extract global stats (last two rows)
        stats_row = df.iloc[-1]
        
        metrics = {
            'total_reads': int(stats_row['1']),
            'length_ok': int(stats_row['2']),
            'mapped_reads': int(stats_row['3']),
            'single_barcode': int(stats_row['4']),
            'no_match': int(stats_row['6']),
            'too_short': int(stats_row['7']),
        }
        
        # Calculate rates
        if metrics['total_reads'] > 0:
            metrics['mapping_rate'] = (metrics['mapped_reads'] / metrics['total_reads']) * 100
        else:
            metrics['mapping_rate'] = 0.0
        
        # Extract per-well counts (grid A-H, 1-12)
        well_counts = {}
        for row_idx, row_letter in enumerate('ABCDEFGH'):
            for col_num in range(1, 13):
                well = f"{row_letter}{col_num}"
                count = int(df.loc[row_letter, str(col_num)])
                if count > 0:
                    well_counts[well] = count
        
        metrics['well_distribution'] = well_counts
        metrics['wells_with_reads'] = len(well_counts)
        
        # Add ambiguous counts if present
        if len(stats_row) >= 10:
            metrics['adapter_only'] = int(stats_row.get('5', 0))
            metrics['ambiguous_multiple_cols'] = int(stats_row.get('8', 0))
            metrics['ambiguous_multiple_rows'] = int(stats_row.get('9', 0))
            metrics['ambiguous_both'] = int(stats_row.get('10', 0))
        
        return metrics
    
    def benchmark_dataset(self, name: str, fastq_file: str, 
                         barcode_csv: str = None,
                         adapter_file: str = None,
                         subset: int = None,
                         config_name: str = None,
                         **demux_params) -> Dict:
        """
        Run benchmark on a specific dataset.
        
        Args:
            name: Benchmark identifier
            fastq_file: Path to input FASTQ
            barcode_csv: Path to barcode CSV
            adapter_file: Path to adapter Python file
            subset: If specified, only process this many reads (for speed)
            config_name: Name of the configuration used
            **demux_params: Additional parameters for demux_barcodes.py
        
        Returns:
            Dictionary with benchmark results
        """
        print(f"\n{'='*70}")
        print(f"Benchmarking: {name}")
        if config_name:
            print(f"Config: {config_name}")
        print(f"{'='*70}")
        
        # Create subset if requested
        if subset:
            print(f"Using subset of {subset} reads for faster benchmarking")
            fastq_file = self.create_subset(fastq_file, subset)
        
        # Set up output directory
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        outdir = self.benchmark_dir / f"{name}_{timestamp}"
        
        # Run demultiplexing
        exec_time, outdir = self.run_demux(
            fastq_file, barcode_csv, str(outdir), adapter_file, **demux_params
        )
        
        # Extract statistics
        stats_file = Path(outdir) / "barcode_stats.csv"
        if not stats_file.exists():
            raise FileNotFoundError(f"Stats file not found: {stats_file}")
        
        stats = self.extract_stats(str(stats_file))
        
        # Clean up FASTQ outputs (keep only stats CSV)
        self.cleanup_output_fastqs(outdir)
        
        # Compile results
        result = {
            'timestamp': datetime.now().isoformat(),
            'name': name,
            'config_name': config_name,
            'fastq_file': fastq_file,
            'barcode_csv': barcode_csv,
            'adapter_file': adapter_file,
            'subset_reads': subset,
            'parameters': demux_params,
            'execution_time': round(exec_time, 2),
            'statistics': stats,
            'output_dir': str(outdir)
        }
        
        # Add to history and save
        self.history.append(result)
        self._save_history()
        
        # Print summary
        self._print_result(result)
        
        return result
    
    def _print_result(self, result: Dict):
        """Print formatted benchmark result."""
        stats = result['statistics']
        params = result['parameters']
        
        print(f"\n✅ Benchmark Complete: {result['name']}")
        if result.get('config_name'):
            print(f"   Config: {result['config_name']}")
        print(f"   Execution time: {result['execution_time']:.2f}s")
        print(f"\n📊 Statistics:")
        print(f"   Total reads:      {stats['total_reads']:,}")
        print(f"   Mapped reads:     {stats['mapped_reads']:,} ({stats['mapping_rate']:.2f}%)")
        print(f"   Single barcode:   {stats['single_barcode']:,}")
        print(f"   No match:         {stats['no_match']:,}")
        print(f"   Too short:        {stats['too_short']:,}")
        print(f"   Wells populated:  {stats['wells_with_reads']}")
        
        if 'ambiguous_multiple_cols' in stats:
            ambig_total = (stats.get('ambiguous_multiple_cols', 0) + 
                          stats.get('ambiguous_multiple_rows', 0) + 
                          stats.get('ambiguous_both', 0))
            print(f"   Ambiguous:        {ambig_total}")
        
        print(f"\n🔧 Parameters:")
        print(f"   Flank:            {params.get('flank', 'N/A')}")
        print(f"   Max penalty:      {params.get('max_penalty', 'N/A')}")
        print(f"   Min length:       {params.get('min_length', 'N/A')}")
        print(f"   CPUs:             {params.get('cpus', 'N/A')}")
        
        print(f"\n   Output: {result['output_dir']}")

    
    def compare_latest(self, n: int = 5):
        """
        Compare the latest n benchmark runs.
        
        Args:
            n: Number of recent runs to compare
        """
        if len(self.history) < 2:
            print("Need at least 2 benchmark runs to compare")
            return
        
        recent = self.history[-n:]
        
        print(f"\n{'='*70}")
        print(f"Comparison of Last {len(recent)} Runs")
        print(f"{'='*70}\n")
        
        # Create comparison table with parameters
        headers = ['Timestamp', 'Name', 'Config', 'Subset', 'Mapped', 'Map %', 'Wells', 'F/P/L', 'Time']
        rows = []
        
        for r in recent:
            stats = r['statistics']
            params = r.get('parameters', {})
            subset_info = f"{r.get('subset_reads', 'full'):,}" if r.get('subset_reads') else 'full'
            config_name = r.get('config_name', 'N/A')[:12]
            
            # Format parameters as F/P/L (Flank/Penalty/Length)
            param_str = f"{params.get('flank', '?')}/{params.get('max_penalty', '?')}/{params.get('min_length', '?')}"
            
            rows.append([
                r['timestamp'][:19],
                r['name'][:18],
                config_name,
                subset_info,
                f"{stats['mapped_reads']:,}",
                f"{stats['mapping_rate']:.1f}%",
                stats['wells_with_reads'],
                param_str,
                f"{r['execution_time']:.1f}s"
            ])
        
        # Print table
        col_widths = [max(len(str(row[i])) for row in [headers] + rows) 
                      for i in range(len(headers))]
        
        # Header
        header_row = ' | '.join(h.ljust(w) for h, w in zip(headers, col_widths))
        print(header_row)
        print('-' * len(header_row))
        
        # Data rows
        for row in rows:
            print(' | '.join(str(cell).ljust(w) for cell, w in zip(row, col_widths)))
        
        print(f"\nF/P/L = Flank / Max Penalty / Min Length")
        
        # Calculate improvements
        if len(recent) >= 2:
            first = recent[0]['statistics']
            last = recent[-1]['statistics']
            
            mapped_change = last['mapped_reads'] - first['mapped_reads']
            rate_change = last['mapping_rate'] - first['mapping_rate']
            
            print(f"\n📈 Change from first to last:")
            print(f"   Mapped reads: {mapped_change:+,}")
            print(f"   Mapping rate: {rate_change:+.2f}%")
    
    def generate_report(self, output_file: str = None):
        """Generate a detailed markdown report of all benchmarks."""
        if not output_file:
            output_file = self.benchmark_dir / "benchmark_report.md"
        
        with open(output_file, 'w') as f:
            f.write("# NanoDemux Benchmarking Report\n\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"Total Runs: {len(self.history)}\n\n")
            
            f.write("## Summary Statistics\n\n")
            
            if self.history:
                # Calculate averages
                avg_mapping = sum(r['statistics']['mapping_rate'] 
                                 for r in self.history) / len(self.history)
                avg_time = sum(r['execution_time'] for r in self.history) / len(self.history)
                
                f.write(f"- Average mapping rate: {avg_mapping:.2f}%\n")
                f.write(f"- Average execution time: {avg_time:.2f}s\n")
                f.write(f"- Best mapping rate: {max(r['statistics']['mapping_rate'] for r in self.history):.2f}%\n")
                f.write(f"- Worst mapping rate: {min(r['statistics']['mapping_rate'] for r in self.history):.2f}%\n\n")
            
            f.write("## All Benchmark Runs\n\n")
            
            for i, result in enumerate(self.history, 1):
                stats = result['statistics']
                f.write(f"### Run {i}: {result['name']}\n\n")
                f.write(f"**Timestamp:** {result['timestamp']}\n\n")
                f.write(f"**Parameters:**\n")
                for k, v in result['parameters'].items():
                    f.write(f"- {k}: {v}\n")
                f.write("\n")
                
                f.write(f"**Results:**\n")
                f.write(f"- Total reads: {stats['total_reads']:,}\n")
                f.write(f"- Mapped reads: {stats['mapped_reads']:,} ({stats['mapping_rate']:.2f}%)\n")
                f.write(f"- Single barcode: {stats['single_barcode']:,}\n")
                f.write(f"- No match: {stats['no_match']:,}\n")
                f.write(f"- Wells with reads: {stats['wells_with_reads']}\n")
                f.write(f"- Execution time: {result['execution_time']:.2f}s\n\n")
        
        print(f"\n📄 Report saved to: {output_file}")


def load_benchmark_configs(base_dir: str = "benchmarking_data") -> List[Dict]:
    """
    Load benchmark configurations from JSON files in subdirectories.
    
    Each subdirectory should have a benchmark_config.json file containing
    a "configs" array with one or more configuration objects.
    
    Example format:
    {
      "configs": [
        {
          "name": "config1",
          "description": "...",
          "barcodes": "file.csv",
          "adapters": "file.py",
          "flank": 100,
          "max_penalty": 60,
          "cpus": 4,
          "min_length": 50
        }
      ]
    }
    
    Returns:
        List of configuration dictionaries with expanded paths and FASTQ files
    """
    configs = []
    base_path = Path(base_dir)
    
    if not base_path.exists():
        print(f"Warning: {base_dir} not found")
        return configs
    
    # Find all subdirectories
    subdirs = [d for d in base_path.iterdir() if d.is_dir() and d.name != 'benchmarks']
    
    for subdir in sorted(subdirs):
        # Look for benchmark_config.json (single file per directory)
        config_file = subdir / "benchmark_config.json"
        
        if not config_file.exists():
            print(f"⚠️  Skipping {subdir.name}: No benchmark_config.json found")
            continue
        
        # Find all FASTQ files in this directory
        fastq_files = sorted(subdir.glob("*.fastq"))
        
        if not fastq_files:
            print(f"⚠️  Skipping {subdir.name}: No FASTQ files found")
            continue
        
        # Load the config file
        try:
            with open(config_file, 'r') as f:
                config_data = json.load(f)
            
            # Check for configs array
            config_list = config_data.get('configs', [])
            
            # Support legacy single-config format
            if not config_list and 'name' in config_data:
                config_list = [config_data]
            
            if not config_list:
                print(f"⚠️  Warning: No 'configs' array found in {config_file}")
                continue
            
            # Process each config in the array
            for config in config_list:
                config_name = config.get('name', 'unnamed')
                barcodes_file = config.get('barcodes')
                adapters_file = config.get('adapters')
                
                # Make paths absolute relative to subdir
                if barcodes_file:
                    barcodes_path = subdir / barcodes_file
                    if not barcodes_path.exists():
                        print(f"⚠️  Warning: Barcode file not found: {barcodes_path}")
                        barcodes_path = None
                else:
                    barcodes_path = None
                
                if adapters_file:
                    adapters_path = subdir / adapters_file
                    if not adapters_path.exists():
                        print(f"⚠️  Warning: Adapter file not found: {adapters_path}")
                        adapters_path = None
                else:
                    adapters_path = None
                
                # Create a benchmark config entry for each FASTQ file
                for fastq_file in fastq_files:
                    benchmark_config = {
                        'name': f"{subdir.name}_{fastq_file.stem}",
                        'config_name': config_name,
                        'config_file': str(config_file),
                        'fastq': str(fastq_file),
                        'barcodes': str(barcodes_path) if barcodes_path else None,
                        'adapters': str(adapters_path) if adapters_path else None,
                        'subdir': subdir.name,
                        'description': config.get('description', ''),
                        'parameters': {
                            'flank': config.get('flank', 100),
                            'max_penalty': config.get('max_penalty', 60),
                            'cpus': config.get('cpus', 4),
                            'min_length': config.get('min_length', 50)
                        }
                    }
                    configs.append(benchmark_config)
                    
        except json.JSONDecodeError as e:
            print(f"❌ Error parsing {config_file}: {e}")
        except Exception as e:
            print(f"❌ Error loading {config_file}: {e}")
    
    return configs


def discover_datasets(base_dir: str = "benchmarking_data") -> List[Dict]:
    """
    DEPRECATED: Use load_benchmark_configs instead.
    
    Dynamically discover all datasets in subdirectories of benchmarking_data.
    This function is kept for backward compatibility but will use config files
    if available, falling back to auto-detection.
    """
    # Try loading from config files first
    configs = load_benchmark_configs(base_dir)
    if configs:
        return configs
    
    # Fallback to old auto-detection method
    datasets = []
    base_path = Path(base_dir)
    
    if not base_path.exists():
        print(f"Warning: {base_dir} not found")
        return datasets
    
    # Find all subdirectories
    subdirs = [d for d in base_path.iterdir() if d.is_dir() and d.name != 'benchmarks']
    
    for subdir in sorted(subdirs):
        # Find FASTQ files
        fastq_files = list(subdir.glob("*.fastq"))
        
        if not fastq_files:
            continue
        
        # Find barcode CSV (prefer demux_barcodes.csv, fallback to *well_map*.csv)
        barcode_csv = None
        if (subdir / "demux_barcodes.csv").exists():
            barcode_csv = subdir / "demux_barcodes.csv"
        else:
            well_map_files = list(subdir.glob("*well_map*.csv"))
            if well_map_files:
                barcode_csv = well_map_files[0]
        
        # Find adapter Python file
        adapter_file = None
        adapter_files = list(subdir.glob("*adapters*.py"))
        if adapter_files:
            adapter_file = adapter_files[0]
        
        # Create dataset entry for each FASTQ file
        for fastq_file in sorted(fastq_files):
            dataset_name = f"{subdir.name}_{fastq_file.stem}"
            datasets.append({
                'name': dataset_name,
                'fastq': str(fastq_file),
                'barcodes': str(barcode_csv) if barcode_csv else None,
                'adapters': str(adapter_file) if adapter_file else None,
                'subdir': subdir.name,
                'parameters': {
                    'flank': 100,
                    'max_penalty': 60,
                    'cpus': 4,
                    'min_length': 50
                }
            })
    
    return datasets


def main():
    """Run benchmarks on sample datasets."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Benchmark NanoDemux demultiplexing performance",
        epilog="""
Examples:
  # Quick benchmark with 100 reads (fast, ~5s)
  python benchmark_demux.py --subset 100
  
  # Default benchmark with 1000 reads (~20s)
  python benchmark_demux.py
  
  # Full dataset benchmark (slow, several minutes)
  python benchmark_demux.py --full
  
  # Compare last 5 benchmarks
  python benchmark_demux.py --compare 5
  
  # Custom parameters
  python benchmark_demux.py --subset 500 --flank 150 --max-penalty 80
  
  # List available datasets
  python benchmark_demux.py --list
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--compare", "-c", type=int, metavar="N",
        help="Compare last N benchmark runs"
    )
    parser.add_argument(
        "--report", "-r", action="store_true",
        help="Generate markdown report"
    )
    parser.add_argument(
        "--list", "-l", action="store_true",
        help="List all available datasets"
    )
    parser.add_argument(
        "--flank", type=int, default=100,
        help="Flank size parameter (default: 100)"
    )
    parser.add_argument(
        "--max-penalty", type=int, default=60,
        help="Max penalty parameter (default: 60)"
    )
    parser.add_argument(
        "--cpus", type=int, default=4,
        help="Number of CPUs (default: 4)"
    )
    parser.add_argument(
        "--subset", type=int, default=1000,
        help="Number of reads to use for benchmarking (default: 1000, use 0 for all)"
    )
    parser.add_argument(
        "--full", action="store_true",
        help="Run on full dataset (equivalent to --subset 0)"
    )
    parser.add_argument(
        "--data-dir", default="benchmarking_data",
        help="Base directory containing benchmark data subfolders (default: benchmarking_data)"
    )
    parser.add_argument(
        "--override-params", action="store_true",
        help="Override config file parameters with command line arguments"
    )
    
    args = parser.parse_args()
    
    bench = DemuxBenchmark()
    
    # If only comparing or reporting, don't run new benchmarks
    if args.compare:
        bench.compare_latest(args.compare)
        return
    
    if args.report:
        bench.generate_report()
        return
    
    # Load configurations from config files
    configs = load_benchmark_configs(args.data_dir)
    
    if not configs:
        print(f"❌ No benchmark configurations found in {args.data_dir}")
        print(f"   Expected: benchmark_config.json files in subfolders")
        print(f"   Each subfolder should contain:")
        print(f"     - benchmark_config.json (or benchmark_config_*.json)")
        print(f"     - *.fastq files")
        print(f"     - barcode CSV and adapter Python files")
        return
    
    # List datasets if requested
    if args.list:
        print(f"\n📦 Found {len(configs)} benchmark configuration(s) in {args.data_dir}:\n")
        for cfg in configs:
            print(f"  • {cfg['name']}")
            print(f"    Config:   {cfg['config_name']} ({cfg['config_file']})")
            print(f"    FASTQ:    {cfg['fastq']}")
            print(f"    Barcodes: {cfg['barcodes'] or 'NOT FOUND'}")
            print(f"    Adapters: {cfg['adapters'] or 'NOT FOUND'}")
            params = cfg['parameters']
            print(f"    Params:   flank={params['flank']}, max_penalty={params['max_penalty']}, " +
                  f"cpus={params['cpus']}, min_length={params['min_length']}")
            if cfg.get('description'):
                print(f"    Description: {cfg['description']}")
            print()
        return
    
    # Run benchmarks on all configurations
    subset_size = 0 if args.full else args.subset
    
    print(f"\n🔬 Running benchmarks on {len(configs)} configuration(s)")
    print(f"   Subset size: {'FULL' if subset_size == 0 else f'{subset_size:,} reads'}")
    if args.override_params:
        print(f"   ⚠️  Overriding config parameters with CLI args")
        print(f"   Flank: {args.flank}, Max penalty: {args.max_penalty}, CPUs: {args.cpus}")
    else:
        print(f"   Using parameters from config files")
    print()
    
    successful = 0
    failed = 0
    
    for config in configs:
        if not Path(config['fastq']).exists():
            print(f"⚠️  Skipping {config['name']}: FASTQ file not found")
            failed += 1
            continue
        
        if not config['barcodes']:
            print(f"⚠️  Skipping {config['name']}: No barcode CSV found")
            failed += 1
            continue
        
        # Use parameters from config, or override with CLI args
        if args.override_params:
            demux_params = {
                'flank': args.flank,
                'max_penalty': args.max_penalty,
                'cpus': args.cpus,
                'min_length': 50
            }
        else:
            demux_params = config['parameters'].copy()
        
        try:
            bench.benchmark_dataset(
                name=config['name'],
                fastq_file=config['fastq'],
                barcode_csv=config['barcodes'],
                adapter_file=config.get('adapters'),
                subset=subset_size if subset_size > 0 else None,
                config_name=config.get('config_name'),
                **demux_params
            )
            successful += 1
        except Exception as e:
            print(f"❌ Error benchmarking {config['name']}: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print(f"\n{'='*70}")
    print(f"✅ Completed: {successful}/{len(configs)} configurations")
    if failed > 0:
        print(f"❌ Failed: {failed}")
    print(f"{'='*70}")
    
    # Show comparison
    if successful > 0:
        print()
        bench.compare_latest(n=min(10, len(configs) * 2))
        
        # Generate report
        bench.generate_report()
    
    # Clean up subset files
    bench.cleanup_subset_files()


if __name__ == "__main__":
    main()
