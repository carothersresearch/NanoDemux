# NanoDemux

[![Tests](https://github.com/carothersresearch/NanoDemux/actions/workflows/tests.yml/badge.svg)](https://github.com/carothersresearch/NanoDemux/actions/workflows/tests.yml)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

NanoDemux is a Python tool for demultiplexing barcoded reads from pooled Oxford Nanopore sequencing experiments. It identifies and separates reads based on dual-indexed barcodes (row and column primers), enabling efficient analysis of multiplexed samples from 96-well plate formats.

## Overview

In pooled nanopore sequencing experiments, multiple samples are combined and sequenced together to increase throughput and reduce costs. Each sample is tagged with unique barcode sequences that allow the reads to be sorted back to their original wells after sequencing. This tool performs that demultiplexing step by:

1. Scanning each read for the presence of row and column barcode sequences
2. Using quality-weighted mismatch scoring to handle sequencing errors
3. Assigning reads to specific wells (e.g., A1, B5, H12) based on barcode combinations
4. Outputting separate FASTQ files for each well

## Key Features

- **Dual-indexing support**: Matches both row (A-H) and column (1-12) barcodes for 96-well plate demultiplexing
- **Quality-aware matching**: Uses Phred quality scores to tolerate sequencing errors intelligently
- **Bidirectional search**: Searches both forward and reverse-complement orientations
- **Multiprocessing**: Parallel processing for faster analysis of large datasets
- **Comprehensive statistics**: Generates detailed mapping statistics and barcode match counts
- **Flexible parameters**: Configurable mismatch tolerance, search regions, and minimum read length

## Installation

### Project Structure

```
NanoDemux/
├── demux_barcodes.py          # Main demultiplexing script
├── benchmark_demux.py         # Benchmarking suite
├── run_tests.py               # Test runner
├── Makefile                   # Build automation and test targets
├── requirements.txt           # Python dependencies
├── primer_well_map.csv        # Example barcode definitions
├── raw_data/                  # Input FASTQ files (add your data here)
│   ├── 55XPXK_1_P4_323_EG.fastq
│   └── VL69M6_1_P4_323_full.fastq
├── benchmarking_data/         # Benchmarking suite data (self-contained)
│   ├── firstpass/             # Sample datasets for benchmarking
│   ├── BENCHMARKING.md        # Benchmarking documentation
│   ├── benchmark_results.json # Historical benchmark results
│   └── */                     # Run-specific outputs (stats only)
├── demplex_data/              # Demultiplexed output examples
│   ├── 55XPXK/
│   └── VL96M6/
├── tests/                     # Comprehensive testing suite
│   ├── test_*.py              # 27 test cases
│   └── *.md                   # Testing documentation
├── .github/workflows/         # CI/CD pipeline
└── DOCKER.md                  # Docker instructions
```

### Testing

This project includes a comprehensive testing suite with 27 test cases covering all functionality:

```bash
# Run all tests (requires pytest)
make test

# Run tests with coverage report
make test-cov

# Show all available test commands
make help
```

See `tests/TESTING_QUICK_REFERENCE.md` for complete testing documentation.

### Benchmarking

Track demultiplexing performance over time as you optimize the algorithm:

```bash
# Quick benchmark (100 reads, ~5 seconds)
make benchmark-fast

# Standard benchmark (1000 reads, ~20 seconds)
make benchmark

# Full dataset benchmark (all reads, several minutes)
make benchmark-full

# Compare recent benchmark results
make benchmark-compare

# View detailed metrics
cat benchmarking_data/benchmarks/benchmark_report.md
```

**Key Features:**
- 🚀 Fast benchmarking with subset sampling
- 📊 Automatic FASTQ cleanup (keeps only stats CSV)
- 📁 Self-contained in `benchmarking_data/` directory
- 🤖 **Automated CI/CD benchmarking** on every push to main
- 🧭 Flexible inputs: pass `--barcodes` to choose a barcode map and optional `--adapters` to test adapter sets

**CI/CD Integration:** Benchmarks run automatically on every push to the main branch. Results are:
- Posted as commit comments
- Uploaded as GitHub Actions artifacts (90-day retention)
- Compared against previous runs to track performance

See `benchmarking_data/BENCHMARKING.md` for complete documentation.

### Option 1: Docker (Recommended)

The easiest way to run NanoDemux is using Docker:

```bash
# Build the Docker image
docker build -t nanodemux .

# Run with your data
docker run --rm -v "$(pwd)":/data nanodemux \
  input.fastq barcodes.csv --outdir output --cpus 4
```

See [DOCKER.md](DOCKER.md) for detailed Docker usage instructions.

### Option 2: Direct Python Installation

Requirements:
- Python 3.6+
- Biopython
- pandas

Install dependencies:

```bash
pip install -r requirements.txt
# or
pip install biopython pandas
```

## Usage

### Basic Command

```bash
python demux_barcodes.py <input.fastq> <barcodes.csv> --outdir <output_directory>
```

### Example

```bash
# Using example barcodes and optional adapters
python demux_barcodes.py raw_data/55XPXK_1_P4_323_EG.fastq \
    barcodes/251202_primer_well_map_DA.csv \
    --adapters barcodes/251205_adapters_DA.py \
    --outdir output/55XPXK --cpus 4 --flank 100

# Using your own data (no adapter detection)
python demux_barcodes.py raw_data/your_data.fastq barcodes/your_barcode_map.csv \
    --outdir demplex_data/your_output --cpus 4
```

### Command-line Options

| Option | Default | Description |
|--------|---------|-------------|
| `fastq` | (required) | Input FASTQ file with pooled reads |
| `barcodes` | (required) | CSV file defining barcode-to-well mappings |
| `--outdir` | `demuxed` | Output directory for demultiplexed FASTQ files |
| `--min_length` | `50` | Minimum read length to process (shorter reads are filtered) |
| `--max_penalty` | `60` | Maximum quality-weighted mismatch penalty (higher = more tolerant) |
| `--cpus` | `1` | Number of CPU cores for parallel processing |
| `--flank` | `100` | Number of bases at each read end to search for barcodes |

## Input Files

### 1. FASTQ File

Standard FASTQ format with nanopore sequencing reads. The tool expects reads to contain adapter sequences with embedded barcodes at either end.

Example:
```
@read_id
AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC...
+
IHHIKLIKHIGHGOKOLJGFFHLNLQPLHIJKJLLEFGFDEEEFHJKMBC...
```

### 2. Barcode CSV File

A CSV file mapping barcode sequences to well positions. Required columns:
- `Well Position`: Well identifier (e.g., A1, B5, H12)
- `Sequence Name`: Barcode identifier (must start with 'R' for row primers, 'F' for column primers)
- `Sequence`: Full adapter sequence containing the barcode

Example format:
```csv
Well Position,Sequence Name,Sequence
A1,R oDA373.D501,AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC
A1,F oDA361.D701,CAAGCAGAAGACGGCATACGAGATATTACTCGGTCTCGTGGGCTCGG
A2,R oDA373.D501,AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC
A2,F oDA362.D702,CAAGCAGAAGACGGCATACGAGATTCCGGAGAGTCTCGTGGGCTCGG
```

💡 **Multiple maps**: Place different barcode maps under `barcodes/` (e.g., `barcodes/<date>_primer_well_map_*.csv`) and pass the desired file on the command line.

### 3. Adapter File (optional)

If you want to detect adapter-only reads, provide a Python file that defines an `ADAPTERS = [...]` list (one adapter sequence per entry). Example: `barcodes/251205_adapters_DA.py`.

- Omit `--adapters` to skip adapter detection.
- Provide any adapter file on the command line to change behavior—no defaults are assumed.

## Output Files

### Demultiplexed FASTQ Files

One FASTQ file per well containing all reads assigned to that well:
- `A1_reads.fastq`, `A2_reads.fastq`, etc.

### Barcode Statistics CSV

`barcode_stats.csv` - A summary grid showing:
- Number of reads assigned to each well
- Row-only and column-only matches
- Total reads, filtered reads, and unmapped reads
- Ambiguous matches (multiple barcodes detected)

## Algorithm

The demultiplexing algorithm works as follows:

1. **Read Processing**: Each read is processed along with its quality scores
2. **Barcode Matching**: 
   - Searches the first and last `--flank` bases of each read
   - Compares against all row and column barcode sequences
   - Calculates quality-weighted mismatch penalty (sum of Phred scores at mismatched positions)
   - A match is called if penalty ≤ `--max_penalty`
3. **Orientation Handling**: Both forward and reverse-complement orientations are checked
4. **Assignment Logic**:
   - **Mapped**: Exactly one row barcode + one column barcode → assigned to that well
   - **Ambiguous**: Multiple barcodes of same type detected → flagged as ambiguous
   - **Single**: Only row or only column barcode found → counted separately
   - **Adapter-only**: Known adapter sequences but no barcodes
   - **No match**: No recognizable sequences found
5. **Output**: Matched reads are written to well-specific FASTQ files

### Quality-Weighted Matching

The quality-weighted mismatch scoring allows the tool to distinguish between likely sequencing errors (low quality scores) and true sequence differences (high quality scores). For example:
- A mismatch at Q10 contributes 10 to the penalty
- A mismatch at Q40 contributes 40 to the penalty
- If the total penalty exceeds `--max_penalty`, the barcode is not considered a match

This approach is more robust than simple edit distance for noisy nanopore data.

## Example Workflow

```bash
# 1. Prepare your barcode mapping file (primer_well_map.csv)
# 2. Run demultiplexing with 4 CPU cores
python demux_barcodes.py my_reads.fastq primer_well_map.csv \
    --outdir demuxed_output \
    --cpus 4 \
    --max_penalty 60 \
    --min_length 50

# 3. Check statistics
cat demuxed_output/barcode_stats.csv

# 4. Process individual well FASTQ files
# e.g., downstream analysis on demuxed_output/A1_reads.fastq
```

## Understanding the Statistics

The tool tracks several categories of reads:
- **total**: All reads in the input file
- **too_short**: Reads below `--min_length`
- **length_ok**: Reads that passed length filter
- **mapped**: Successfully assigned to a well (one row + one column barcode)
- **single**: Only one barcode type found (row-only or column-only)
- **adapter_only**: Adapter sequences detected but no specific barcodes
- **no_match**: No recognizable sequences found
- **ambiguous_***: Reads with multiple barcode matches

The grid in `barcode_stats.csv` shows the distribution of reads across the 96-well plate, with row and column totals.

## Troubleshooting

### Low mapping rates
- Increase `--max_penalty` to allow more mismatches
- Increase `--flank` to search a larger region at read ends
- Check that your barcode CSV matches your library preparation

### Too many ambiguous reads
- Decrease `--max_penalty` to be more stringent
- Verify barcode sequences are sufficiently different

### Memory issues with large files
- Decrease `--cpus` to reduce memory usage
- The tool processes reads in chunks to manage memory

## Citation

If you use NanoDemux in your research, please cite the Carothers Research Lab and this repository.

## License

This project is maintained by the Carothers Research Lab.

## Contact

For questions or issues, please open an issue on the GitHub repository.

