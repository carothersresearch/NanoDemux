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

- **Comprehensive Workflow**: Single command to run complete analysis pipeline with automated HTML reports
- **Dual-indexing support**: Matches both row (A-H) and column (1-12) barcodes for 96-well plate demultiplexing
- **Quality-aware matching**: Uses Phred quality scores to tolerate sequencing errors intelligently
- **Sequence alignment**: Performs MSA and generates quality-weighted consensus sequences for each well
- **Bidirectional search**: Searches both forward and reverse-complement orientations
- **Multiprocessing**: Parallel processing for faster analysis of large datasets
- **Comprehensive statistics**: Generates detailed mapping statistics and barcode match counts
- **Flexible parameters**: Configurable mismatch tolerance, search regions, and minimum read length
- **Quality reports**: Graphical reports for raw and demultiplexed data with MSA visualization

## Installation

### Project Structure

```
NanoDemux/
├── nanodemux                   # 🌟 Comprehensive workflow script (runs all steps)
├── demux_barcodes.py          # Main demultiplexing script (supports single or multi-file)
├── generate_raw_quality_report.py  # Raw data quality analysis
├── generate_quality_report.py      # Demultiplexed data quality reports
├── align_wells.py             # Sequence alignment and consensus generation for wells
├── benchmark_demux.py         # Benchmarking suite
├── run_tests.py               # Test runner
├── Makefile                   # Build automation and test targets
├── requirements.txt           # Python dependencies
├── primer_well_map.csv        # Example barcode definitions
├── raw_data/                  # Input FASTQ files (organize in subfolders)
│   ├── 55XPXK_1_P4_323_EG.fastq          # Single file example
│   ├── VL69M6_1_P4_323_full.fastq        # Single file example
│   └── experiment1/                       # Directory example (batch processing)
│       ├── sample1.fastq
│       └── sample2.fastq
├── benchmarking_data/         # Benchmarking suite data (self-contained)
│   ├── firstpass/             # Sample datasets for benchmarking
│   ├── BENCHMARKING.md        # Benchmarking documentation
│   ├── benchmark_results.json # Historical benchmark results
│   └── */                     # Run-specific outputs (stats only)
├── demplex_data/              # Demultiplexed output (auto-generated)
│   ├── 55XPXK_1_P4_323_EG/    # Single file output
│   │   ├── A1_reads.fastq
│   │   ├── A2_reads.fastq
│   │   └── barcode_stats.csv
│   └── experiment1/           # Multi-file output (batch)
│       ├── sample1/
│       │   ├── A1_reads.fastq
│       │   └── barcode_stats.csv
│       └── sample2/
│           ├── A1_reads.fastq
│           └── barcode_stats.csv
├── tests/                     # Comprehensive testing suite
│   ├── test_*.py              # 30+ test cases
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

### Comprehensive Workflow Script (Recommended)

The `nanodemux` script provides a complete end-to-end workflow that runs all analysis steps and generates a comprehensive HTML report:

```bash
# Basic usage - runs full workflow with default parameters
python nanodemux <input.fastq> <barcodes.csv>

# Custom output directory with more CPUs
python nanodemux raw_data/reads.fastq barcodes/primer_well_map.csv \
    --output results/ --cpus 4

# With all custom parameters
python nanodemux raw_data/reads.fastq barcodes/primer_well_map.csv \
    --output results/ --cpus 4 --flank 150 --max-penalty 80 \
    --anchor-seq AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC \
    --min-quality 25
```

**The workflow automatically:**
1. Generates raw quality report (pre-demultiplexing analysis)
2. Demultiplexes reads by barcodes into wells
3. Generates quality report on demultiplexed data with MSA visualization
4. Performs sequence alignment and generates consensus sequences for each well
5. Creates a comprehensive HTML report combining all results

**Output Structure:**
```
output/
├── 1_raw_quality_report/          # Raw data quality analysis
├── 2_demultiplexed/                # Demultiplexed FASTQ files by well
├── 3_demux_quality_report/         # Quality report on demuxed data
├── 4_aligned/                      # Aligned sequences and consensus
└── comprehensive_report.html       # Summary of all results
```

**Skip Steps:**
```bash
# Skip raw QC and alignment steps
python nanodemux input.fastq barcodes.csv \
    --skip-raw-qc --skip-alignment
```

See `python nanodemux --help` for all available parameters from each workflow step.

### Individual Script Usage

You can also run each step independently:

**Single File Processing:**
```bash
python demux_barcodes.py <input.fastq> <barcodes.csv>
```

**Multiple File Processing (Directory):**
```bash
python demux_barcodes.py <input_directory> <barcodes.csv>
```

### Examples

```bash
# Single file - using example barcodes and optional adapters
python demux_barcodes.py raw_data/55XPXK_1_P4_323_EG.fastq \
    barcodes/251202_primer_well_map_DA.csv \
    --adapters barcodes/251205_adapters_DA.py \
    --cpus 4 --flank 100

# Single file - with custom output directory
python demux_barcodes.py raw_data/your_data.fastq barcodes/your_barcode_map.csv \
    --outdir custom_output/ --cpus 4

# Process all FASTQ files in a directory
python demux_barcodes.py raw_data/experiment1/ barcodes/251202_primer_well_map_DA.csv \
    --adapters barcodes/251205_adapters_DA.py --cpus 4
```

### Output Directory Structure

The tool automatically creates an organized output structure:

- **Single file (default)**: `demplex_data/<filename>/` - Each file gets its own directory
- **Single file (custom)**: Use `--outdir` to specify exact output location
- **Multiple files**: `demplex_data/<directory_name>/<filename>/` - Directory structure is preserved

Example:
```
# Input: raw_data/experiment1/sample1.fastq
# Output: demplex_data/experiment1/sample1/*.fastq + barcode_stats.csv

# Input: raw_data/experiment1/ (containing sample1.fastq, sample2.fastq)
# Output: demplex_data/experiment1/sample1/*.fastq + barcode_stats.csv
#         demplex_data/experiment1/sample2/*.fastq + barcode_stats.csv
```

### Command-line Options

| Option | Default | Description |
|--------|---------|-------------|
| `fastq` | (required) | Input FASTQ file or directory containing FASTQ files |
| `barcodes` | (required) | CSV file defining barcode-to-well mappings |
| `--outdir` | `demuxed` | Output directory (auto-generated for default, see above) |
| `--min_length` | `50` | Minimum read length to process (shorter reads are filtered) |
| `--max_penalty` | `60` | Maximum quality-weighted mismatch penalty (higher = more tolerant) |
| `--cpus` | `1` | Number of CPU cores for parallel processing |
| `--flank` | `100` | Number of bases at each read end to search for barcodes |
| `--var_q` | `10` | Phred quality cutoff to call a variable-base match 'confident' |
| `--no-fallback` | `False` | Disable alignment fallback for ambiguous cases (use only fast filter) |
| `--adapters` | `None` | Optional Python file defining ADAPTERS list for adapter detection |
| `--report` | `False` | Generate graphical quality report with MSA, read lengths, and barcode analysis |
| `--raw-report` | `False` | Generate graphical quality report for raw (unmultiplexed) data |
| `--anchor-seq` | `None` | Optional anchor/reference sequence for MSA alignment in quality report |

## Graphical Quality Reports

NanoDemux can generate two types of quality reports:

### 1. Raw Data Quality Reports

Analyze unmultiplexed FASTQ files before demultiplexing to assess overall data quality:
- **Read Length Distribution**: Histograms and box plots with N50 statistics
- **Quality Score Distribution**: Mean quality per read and quality by position analysis
- **Base Composition**: Nucleotide distribution and GC content
- **General Statistics**: Total reads, total bases, mean/median lengths

#### Generate Raw Data Reports

```bash
# Automatically generate raw report before demultiplexing
python demux_barcodes.py raw_data/reads.fastq barcodes/primer_well_map.csv --raw-report

# Or generate raw report separately
python generate_raw_quality_report.py raw_data/reads.fastq

# Generate raw report with custom output directory
python generate_raw_quality_report.py raw_data/reads.fastq --output reports/raw_qc/

# Analyze only first 10000 reads for faster processing
python generate_raw_quality_report.py raw_data/large_sample.fastq --max-reads 10000
```

### 2. Demultiplexed Data Quality Reports

Analyze demultiplexed data after barcode matching:
- **Multiple Sequence Alignment (MSA) Layout**: Visualization of read positions and lengths, similar to [SeqAn ReadLayout](https://seqan.readthedocs.io/en/seqan-v1.4.2/Tutorial/FragmentStore.html), with barcode positions highlighted
- **Anchor Sequence MSA** (optional): Align all reads to a reference/anchor sequence using Smith-Waterman alignment. Perfect for analyzing how reads align to a known template or reference sequence.
- **Read Length Distribution**: Histograms and box plots showing read length statistics across all wells
- **Barcode Position Analysis**: Heatmaps showing where row and column barcodes are detected within reads
- **Barcode Presence Summary**: 96-well plate heatmap with read distribution statistics

#### Generate Demultiplexed Reports

```bash
# Automatically generate report during demultiplexing
python demux_barcodes.py raw_data/reads.fastq barcodes/primer_well_map.csv --report

# Generate report with anchor sequence MSA during demultiplexing
python demux_barcodes.py raw_data/reads.fastq barcodes/primer_well_map.csv \
    --report --anchor-seq AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC

# Or generate report separately from existing demux output
python generate_quality_report.py demplex_data/55XPXK_1_P4_323_EG/ barcodes/primer_well_map.csv

# Generate report with anchor sequence MSA
python generate_quality_report.py demplex_data/55XPXK_1_P4_323_EG/ \
    barcodes/primer_well_map.csv \
    --anchor-seq AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC
```

Reports are saved as HTML files with embedded visualizations:
- Raw data reports: `raw_data_quality_report/raw_quality_report.html` (default) or in the output directory's `raw_quality_report/` subdirectory
- Demultiplexed reports: `quality_report/` subdirectory of your demux output folder

## Sequence Alignment and Consensus Generation

After demultiplexing, you can perform multiple sequence alignment (MSA) and generate quality-weighted consensus sequences for each well using the `align_wells.py` script.

### Features

- **Quality-Weighted Consensus**: Uses Phred quality scores to calculate consensus sequences
- **Multiple Sequence Alignment**: Performs pairwise Smith-Waterman alignment using parasail
- **FASTQ Output**: Generates aligned FASTQ files with consensus as the first sequence
- **CSV Summary**: Creates a comprehensive CSV with consensus sequences for all wells

### Usage

```bash
# Basic usage - align all reads in each well
python align_wells.py demplex_data/55XPXK_1_P4_323_EG/ aligned_output/

# Limit reads per well for faster processing
python align_wells.py demplex_data/experiment/ aligned_output/ --max-reads 100

# Skip MSA and calculate consensus from unaligned sequences
python align_wells.py demux_dir/ output_dir/ --no-align

# Custom CSV output location
python align_wells.py demux_dir/ output_dir/ --csv consensus_summary.csv
```

### Command-line Options

| Option | Default | Description |
|--------|---------|-------------|
| `demux_dir` | (required) | Directory containing demultiplexed FASTQ files (`*_reads.fastq`) |
| `output_dir` | (required) | Output directory for aligned FASTQ files and CSV summary |
| `--max-reads` | `None` | Maximum number of reads to process per well (all reads by default) |
| `--min-quality` | `20` | Minimum Phred quality score to consider a base for consensus |
| `--no-align` | `False` | Skip MSA and calculate consensus from unaligned sequences |
| `--csv` | `<output_dir>/consensus_sequences.csv` | Custom path for consensus CSV output |

### Output Files

**Aligned FASTQ files** (`<output_dir>/<WELL>_aligned.fastq`):
- First sequence: Quality-weighted consensus sequence for the well
- Remaining sequences: Original reads from the well

**Consensus CSV** (`<output_dir>/consensus_sequences.csv`):
- `Well`: Well identifier (e.g., A1, B5)
- `Num_Reads`: Number of reads in the well
- `Consensus_Length`: Length of the consensus sequence
- `Consensus_Sequence`: The consensus sequence
- `Avg_Coverage`: Average number of reads at each position
- `Avg_Consensus_Quality`: Average Phred quality score of consensus
- `Min_Coverage`, `Max_Coverage`: Coverage depth statistics

### Example Workflow

```bash
# 1. Demultiplex reads
python demux_barcodes.py raw_data/experiment.fastq barcodes/primer_map.csv --cpus 4

# 2. Generate consensus sequences for each well
python align_wells.py demplex_data/experiment/ aligned_wells/ --max-reads 500

# 3. View results
cat aligned_wells/consensus_sequences.csv
head aligned_wells/A1_aligned.fastq
```

## Input Files

### 1. FASTQ File or Directory

You can provide either:
- **Single FASTQ file**: A standard FASTQ format file with nanopore sequencing reads
- **Directory containing FASTQ files**: The tool will automatically find and process all files with extensions `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`

The tool expects reads to contain adapter sequences with embedded barcodes at either end.

Example FASTQ format:
```
@read_id
AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC...
+
IHHIKLIKHIGHGOKOLJGFFHLNLQPLHIJKJLLEFGFDEEEFHJKMBC...
```

Recommended directory structure:
```
raw_data/
├── experiment1/
│   ├── sample1.fastq
│   ├── sample2.fastq
│   └── sample3.fastq
└── experiment2/
    ├── sampleA.fastq
    └── sampleB.fastq
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

## Example Workflows

### Single File Processing

```bash
# 1. Prepare your barcode mapping file (primer_well_map.csv)
# 2. Run demultiplexing on a single FASTQ file with 4 CPU cores
python demux_barcodes.py raw_data/my_reads.fastq barcodes/primer_well_map.csv \
    --cpus 4 \
    --max_penalty 60 \
    --min_length 50

# 3. Check statistics (output will be in demplex_data/my_reads/)
cat demplex_data/my_reads/barcode_stats.csv

# 4. Process individual well FASTQ files
# e.g., downstream analysis on demplex_data/my_reads/A1_reads.fastq
```

### Multiple File Processing (Batch Mode)

```bash
# 1. Organize your FASTQ files in a subdirectory
mkdir -p raw_data/experiment1
mv *.fastq raw_data/experiment1/

# 2. Run demultiplexing on all files in the directory
python demux_barcodes.py raw_data/experiment1/ barcodes/primer_well_map.csv \
    --cpus 4 \
    --max_penalty 60

# 3. Output will be organized in demplex_data/experiment1/
#    - demplex_data/experiment1/sample1/barcode_stats.csv
#    - demplex_data/experiment1/sample1/A1_reads.fastq
#    - demplex_data/experiment1/sample2/barcode_stats.csv
#    - demplex_data/experiment1/sample2/A1_reads.fastq
#    etc.

# 4. Check statistics for each file
for dir in demplex_data/experiment1/*/; do
    echo "Statistics for $(basename $dir):"
    cat "$dir/barcode_stats.csv"
    echo ""
done
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

