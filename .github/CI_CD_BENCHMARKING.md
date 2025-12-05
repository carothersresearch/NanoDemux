# CI/CD Benchmarking

## Overview

Automated benchmarking is integrated into the GitHub Actions CI/CD pipeline. Benchmarks run automatically on every push to the `main` branch after all tests pass.

## Workflow

```
Push to main → Run Tests (all Python versions) → Run Benchmarks (if tests pass)
                                                 ↓
                                    Download previous results
                                                 ↓
                                    Run benchmark-fast (100 reads)
                                                 ↓
                                    Compare with previous runs
                                                 ↓
                                    Upload results as artifacts
                                                 ↓
                                    Post summary to commit comments
```

## What Gets Benchmarked

- **Trigger**: Every push to `main` branch (after tests pass)
- **Datasets**: 55XPXK and VL69M6 sample datasets
- **Subset Size**: 100 reads per dataset (fast, ~5 seconds total)
- **Python Version**: 3.11
- **Metrics Tracked**:
  - Mapped reads (both barcodes found)
  - Mapping rate (percentage)
  - Wells populated
  - Execution time

## Accessing Results

### 1. Commit Comments

Each benchmarked commit gets an automatic comment with:
- Latest run results (table format)
- Comparison with previous run
- Link to download full artifacts

**Example:**
```
## 📊 Benchmark Results

### Latest Run

| Dataset | Subset | Mapped | Mapping % | Wells | Time (s) |
|---------|--------|--------|-----------|-------|----------|
| 55XPXK  |    100 |      9 |     9.00% |     9 |     2.30 |
| VL69M6  |    100 |     19 |    19.00% |    17 |     2.01 |

### Change Since Previous Run
- Mapped reads: +10
- Mapping rate: +10.00%
```

### 2. GitHub Actions Artifacts

Full benchmark results are uploaded as artifacts:
- Navigate to **Actions** tab → Select workflow run → **Artifacts** section
- Download `benchmark-results` artifact (ZIP file)
- Contains:
  - `benchmark_results.json` - Complete history
  - `benchmark_report.md` - Detailed markdown report
- **Retention**: 90 days

### 3. Actions Logs

View detailed execution logs:
- Actions tab → Select workflow run → `benchmark` job
- See step-by-step benchmark execution
- View comparison tables and timing information

## Files

- **Workflow**: `.github/workflows/tests.yml`
- **Summary Script**: `.github/scripts/benchmark_summary.py`
- **Results**: `benchmarking_data/benchmarks/benchmark_results.json`
- **Reports**: `benchmarking_data/benchmarks/benchmark_report.md`

## Why Fast Benchmarks in CI?

- ⚡ **Speed**: 100 reads takes ~5 seconds (vs. minutes for full datasets)
- 📊 **Trends**: Still tracks performance trends over time
- 💰 **Cost**: Minimizes GitHub Actions minutes usage
- ✅ **Coverage**: Benchmarks both sample datasets
- 🔍 **Local Testing**: Run full benchmarks locally when needed

## Running Locally

To match CI benchmarks:
```bash
make benchmark-fast
```

For more detailed benchmarking:
```bash
make benchmark        # 1000 reads, ~20 seconds
make benchmark-full   # All reads, several minutes
```

## Troubleshooting

### Benchmarks Not Running

- Check that push is to `main` branch
- Verify all tests passed first (benchmarks only run after tests)
- Check Actions tab for any failures

### No Previous Results

- First run won't have comparison data (expected)
- Subsequent runs will download previous artifacts automatically

### Artifacts Not Found

- Artifacts expire after 90 days
- Download and save important results locally
- Benchmark history in `benchmark_results.json` accumulates over time

## Modifying Benchmark Configuration

### Change Subset Size

Edit `.github/workflows/tests.yml`:
```yaml
- name: Run benchmarks
  run: |
    python benchmark_demux.py --subset 500  # Increase to 500 reads
```

### Add Custom Parameters

```yaml
- name: Run benchmarks
  run: |
    python benchmark_demux.py --subset 100 --flank 150 --max-penalty 80
```

### Change Python Version

```yaml
- name: Set up Python 3.11
  uses: actions/setup-python@v4
  with:
    python-version: "3.12"  # Use different version
```

## Benefits

- 📈 **Track Performance**: Automatically track performance across commits
- 🐛 **Catch Regressions**: Identify performance degradations early
- 📊 **Historical Data**: Build comprehensive performance history
- 🤖 **Zero Effort**: Runs automatically, no manual intervention
- 💬 **Visibility**: Results posted directly on commits
- 🔍 **Detailed Analysis**: Full reports available as artifacts

## See Also

- [BENCHMARKING.md](../benchmarking_data/BENCHMARKING.md) - Full benchmarking documentation
- [tests.yml](../workflows/tests.yml) - Complete workflow definition
- [benchmark_summary.py](../scripts/benchmark_summary.py) - Summary generation script
