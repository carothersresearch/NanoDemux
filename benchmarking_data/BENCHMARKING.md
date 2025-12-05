# Benchmarking Suite Documentation

## Overview

The benchmarking suite tracks demultiplexing performance metrics over time, allowing you to measure improvements as you optimize the `demux_barcodes.py` algorithm. It automatically runs the pipeline on sample datasets, extracts key statistics, and maintains a historical record for comparison.

**Key Features:**
- 🚀 Fast benchmarking with subset sampling (default: 1000 reads)
- 📊 Automatic FASTQ cleanup (keeps only stats CSV)
- 📁 Self-contained in `benchmarking_data/` directory
- 🧭 Flexible barcode/adapter inputs via `--barcodes` and optional `--adapters`
- 🤖 **CI/CD Integration** - Automated benchmarking on every push to main

## CI/CD Integration

Benchmarks run automatically in GitHub Actions on every push to the `main` branch:

### What Happens
1. ✅ Tests run first across all Python versions (3.9-3.12)
2. 📊 Benchmarks run after tests pass (Python 3.11)
3. 🔍 Results are compared with previous runs
4. 💬 Performance summary posted as commit comment
5. 📦 Full results uploaded as GitHub Actions artifacts (90-day retention)

### Viewing Results

**Commit Comments:**
- Check commit comments for benchmark summaries
- Shows key metrics: mapped reads, mapping rate, execution time
- Compares against previous benchmarks

**GitHub Actions Artifacts:**
- Download `benchmark-results` artifact from Actions tab
- Contains `benchmark_results.json` and `benchmark_report.md`
- Retained for 90 days

**Why Quick Benchmarks in CI:**
- Uses `make benchmark-fast` (100 reads, ~5 seconds)
- Keeps CI fast while still tracking performance trends
- Run full benchmarks locally for detailed analysis

## Quick Start

### Run Benchmarks

```bash
# Quick benchmark (100 reads, ~5 seconds)
make benchmark-fast

# Standard benchmark (1000 reads, ~20 seconds)
make benchmark

# Full dataset benchmark (all reads, several minutes)
make benchmark-full

# Or run directly with custom parameters
python benchmark_demux.py --subset 500 --flank 100 --max-penalty 60
```

### Compare Results

```bash
# Compare last 5 benchmark runs
make benchmark-compare

# Or compare custom number
python benchmark_demux.py --compare 10
```

### Generate Report

```bash
# Generate markdown report
python benchmark_demux.py --report
```

## What Gets Measured

The benchmark suite tracks:

### Performance Metrics
- ⏱️ **Execution time** - Total time to process dataset
- 📊 **Mapping rate** - Percentage of reads assigned to wells
- 🎯 **Mapped reads** - Reads with both row + column barcodes
- 📍 **Wells populated** - Number of wells receiving reads
- 📏 **Subset size** - Number of reads processed (for comparison)

### Read Classification
- **Total reads** - All reads in input FASTQ
- **Mapped reads** - Both barcodes found (assigned to well)
- **Single barcode** - Only row OR column barcode found
- **No match** - No barcodes detected
- **Too short** - Reads below minimum length
- **Ambiguous** - Multiple barcode matches detected

### Well Distribution
- Per-well read counts (A1-H12)
- Number of populated wells
- Read distribution across plate

## Directory Structure

All benchmark files are self-contained in `benchmarking_data/`:

```
benchmarking_data/
├── BENCHMARKING.md              # This documentation
├── firstpass/                   # Sample datasets
│   ├── 55XPXK_1_P4_323_EG.fastq
│   ├── VL69M6_1_P4_323_full.fastq
│   └── primer_well_map.csv
└── benchmarks/                  # All benchmark outputs
    ├── benchmark_results.json   # Historical results (JSON)
    ├── benchmark_report.md      # Detailed markdown report
    ├── 55XPXK_20251205_223801/  # Run-specific outputs
    │   └── barcode_stats.csv    # Only stats CSV (no FASTQs)
    └── VL69M6_20251205_223923/
        └── barcode_stats.csv
```

**Note:** 
- FASTQ output files are automatically cleaned up after statistics extraction to save disk space
- Temporary subset FASTQ files are also cleaned up at the end of benchmarking
- Only the `barcode_stats.csv` file is retained for each run

## Output Files

### benchmark_results.json

Complete history of all benchmark runs in JSON format. Each entry contains:
- Timestamp
- Dataset name
- Subset size (number of reads processed)
- Parameters used (flank, max_penalty, cpus, etc.)
- Execution time
- Full statistics
- Output directory path

### benchmark_report.md

Human-readable markdown report with:
- Summary statistics (averages, best/worst)
- Detailed results for each run
- Parameters used
- Trends over time

## Sample Datasets

The benchmark suite uses data from `benchmarking_data/firstpass/`:

| Dataset | File | Reads | Description |
|---------|------|-------|-------------|
| 55XPXK | `55XPXK_1_P4_323_EG.fastq` | ~15,000 | Sample dataset 1 |
| VL69M6 | `VL69M6_1_P4_323_full.fastq` | ~20,000 | Sample dataset 2 |

By default, the script processes a subset of 1000 reads from each file for faster benchmarking.

## Subset Benchmarking

For rapid iteration, the suite supports subset sampling:

- **Default:** 1000 reads (~20 seconds total)
- **Quick:** 100 reads (~5 seconds total) via `--subset 100`
- **Full:** All reads via `--full` flag (several minutes)

Subset files are created temporarily during benchmarking and automatically cleaned up at the end to avoid clutter.

## Command-Line Options

```bash
python benchmark_demux.py [OPTIONS]
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--subset N` | 1000 | Number of reads to process per dataset |
| `--full` | - | Process all reads (overrides --subset) |
| `--compare N` | - | Compare last N benchmark runs (no new benchmarks) |
| `--report` | - | Generate markdown report only (no new benchmarks) |
| `--barcodes PATH` | `benchmarking_data/firstpass/251202_primer_well_map_DA.csv` | Barcode map to use for benchmarking |
| `--adapters PATH` | `benchmarking_data/firstpass/251205_adapters_DA.py` | Adapter file defining `ADAPTERS` list (can be overridden) |
| `--flank` | 100 | Flank size for barcode search |
| `--max-penalty` | 60 | Maximum quality-weighted mismatch penalty |
| `--cpus` | 4 | Number of CPU cores to use |

### Examples

```bash
# Quick test with 100 reads
python benchmark_demux.py --subset 100

# Standard benchmark (1000 reads)
python benchmark_demux.py

# Full dataset with custom parameters
python benchmark_demux.py --full --flank 150 --max-penalty 80

# Compare last 3 runs
python benchmark_demux.py --compare 3

# Generate report without running new benchmarks
python benchmark_demux.py --report
```

## Workflow for Testing Improvements

### 1. Baseline Benchmark

```bash
# Run initial benchmark with current code
make benchmark-fast
```

### 2. Make Code Changes

Edit `demux_barcodes.py` to implement improvements.

### 3. Run New Benchmark

```bash
# Run benchmark with same parameters
make benchmark-fast
```

### 4. Compare Results

```bash
# See if metrics improved
make benchmark-compare
```

### 5. Iterate

Repeat steps 2-4, comparing each iteration to see if changes improve:
- Mapping rate (higher is better)
- Execution time (lower is better)
- Wells populated (higher is better)

## Example Output

```
======================================================================
Benchmarking: 55XPXK
======================================================================
Running: python demux_barcodes.py raw_data/55XPXK_1_P4_323_EG.fastq ...

✅ Benchmark Complete: 55XPXK
   Execution time: 81.31s

📊 Statistics:
   Total reads:      5,000
   Mapped reads:     329 (6.58%)
   Single barcode:   3,264
   No match:         1,052
   Too short:        5
   Wells populated:  50
   Ambiguous:        13

   Output: benchmarks/55XPXK_20251205_223801

======================================================================
Comparison of Last 5 Runs
======================================================================

Timestamp           | Name   | Mapped | Mapping % | Wells | Time (s)
--------------------------------------------------------------------
2025-12-05T22:30:00 | 55XPXK | 89     | 1.78%     | 22    | 75.20
2025-12-05T22:35:00 | 55XPXK | 329    | 6.58%     | 50    | 81.31
2025-12-05T22:40:00 | 55XPXK | 450    | 9.00%     | 65    | 79.45

📈 Change from first to last:
   Mapped reads: +361
   Mapping rate: +7.22%
```

## Interpreting Results

### Good Improvements
- ✅ **Mapping rate increases** - More reads assigned to wells
- ✅ **Wells populated increases** - Better plate coverage
- ✅ **Single barcode decreases** - Fewer partial matches
- ✅ **No match decreases** - Fewer unassigned reads
- ✅ **Execution time stable/decreases** - No performance regression

### Warning Signs
- ⚠️ **Mapping rate decreases** - Algorithm may be too strict
- ⚠️ **Ambiguous increases significantly** - May be too permissive
- ⚠️ **Execution time increases** - Performance regression
- ⚠️ **Wells populated decreases** - Reads concentrating in fewer wells

## Tips for Optimization

### Increasing Mapping Rate

Try adjusting parameters:
```bash
# More permissive matching
python benchmark_demux.py --max-penalty 80 --flank 150

# More strict matching
python benchmark_demux.py --max-penalty 40 --flank 80
```

### Parameter Sweep

Test multiple parameter combinations:
```bash
# Test different flank sizes
for flank in 50 75 100 125 150; do
    python benchmark_demux.py --flank $flank
done

# Compare results
python benchmark_demux.py --compare 5
```

## Integration with Tests

The benchmark suite complements the test suite:

| Test Suite | Purpose | When to Use |
|------------|---------|-------------|
| `make test` | Verify correctness | Before every commit |
| `make benchmark` | Measure performance | When optimizing algorithm |

**Workflow:**
1. Make code changes
2. Run `make test` to ensure correctness
3. Run `make benchmark` to measure performance impact
4. Compare results to verify improvements

## Continuous Monitoring

### Track Over Time

Keep `benchmarks/benchmark_results.json` in version control to:
- Track historical performance
- Identify regressions
- Document optimization progress

### CI/CD Integration

Add to `.github/workflows/tests.yml`:
```yaml
- name: Run benchmarks
  run: |
    python benchmark_demux.py
    python benchmark_demux.py --report
```

## Troubleshooting

### Benchmark Fails

```bash
# Check if data files exist
ls -la raw_data/

# Verify demux script works
python demux_barcodes.py --help

# Run with verbose output
python benchmark_demux.py 2>&1 | tee benchmark.log
```

### Inconsistent Results

- Ensure consistent parameters across runs
- Use same input data
- Check for system load (other processes running)
- Run multiple times and average results

### Missing Data

If benchmark data files are missing:
```bash
# Check data directory
ls -la raw_data/

# Regenerate if needed (or use your own data)
python demux_barcodes.py your_input.fastq primer_well_map.csv --outdir raw_data/
```

## Advanced Usage

### Custom Datasets

Add your own datasets:

```python
datasets = [
    {
        'name': 'my_dataset',
        'fastq': 'raw_data/my_data.fastq',
        'subset': 10000
    }
]
```

### Export Results

```python
from benchmark_demux import DemuxBenchmark

bench = DemuxBenchmark()
results = bench.history

# Export to CSV
import pandas as pd
df = pd.DataFrame([r['statistics'] for r in results])
df.to_csv('benchmark_summary.csv')
```

## Files Modified

The benchmarking suite adds:
- `benchmark_demux.py` - Main benchmarking script
- `benchmarks/` - Output directory (in .gitignore)
- Updated `Makefile` with benchmark targets
- This documentation file

## Next Steps

1. **Run baseline** - `make benchmark` to establish current performance
2. **Make improvements** - Optimize the demux algorithm
3. **Compare** - `make benchmark-compare` to see results
4. **Iterate** - Repeat until satisfied with performance
5. **Document** - Note successful parameter combinations

---

**Status**: Ready to use
**Last Updated**: Dec 5, 2025
**Maintained By**: Benchmarking Suite
