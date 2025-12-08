**NanoDemux — User Manual**

**Purpose:** This manual explains, step-by-step, how to use this repository to analyze new FASTQ data: where to place raw files, which supporting files are needed, how to run demultiplexing, how to generate raw and demultiplexed quality reports, and how to run the benchmark suite.

**Prerequisites**
- Python 3.10+ (or the version used in the dev container).
- Install dependencies from `requirements.txt`:

```bash
pip install -r requirements.txt
```

**Quick file map (important paths)**
- Raw FASTQ files: place in `raw_data/` (or anywhere; pass path to scripts explicitly).
- Barcode CSV files: `barcodes/` (examples: `251202_primer_well_map_DA.csv`, `demux_barcodes.csv`).
- Adapter definitions (optional): Python files alongside barcodes (e.g., `barcodes/251205_adapters_DA.py`).
- Demultiplexed output: created by `demux_barcodes.py` in an output directory you specify (default: `demplex_data/<filename>/` or custom via `--outdir`).
- Benchmarks: sample files and configs in `benchmarking_data/`.

**1) Comprehensive Workflow Script (Recommended)**

The `nanodemux` script provides a single-command solution to run the complete NanoDemux analysis pipeline. It automatically executes all steps and generates a comprehensive HTML report.

**Quick Start:**

```bash
python nanodemux /path/to/raw_reads.fastq barcodes/your_barcodes.csv
```

This runs the full workflow:
1. Raw quality analysis (pre-demultiplexing)
2. Demultiplexing by barcodes
3. Quality report on demultiplexed data with MSA visualization
4. Sequence alignment and consensus generation for each well
5. Comprehensive HTML report combining all results

**Output Structure:**
```
nanodemux_output/
├── 1_raw_quality_report/          # Raw data quality analysis
├── 2_demultiplexed/                # Demultiplexed FASTQ files by well
├── 3_demux_quality_report/         # Quality report on demuxed data
├── 4_aligned/                      # Aligned sequences and consensus
└── comprehensive_report.html       # Summary of all results
```

**Common Options:**

```bash
# Custom output directory with more CPUs
python nanodemux raw_data/reads.fastq barcodes/primer_well_map.csv \
    --output results/ --cpus 4

# With adapter detection and custom parameters
python nanodemux raw_data/reads.fastq barcodes/primer_well_map.csv \
    --output results/ --adapters barcodes/adapters.py \
    --cpus 4 --flank 150 --max-penalty 80 --var-q 15 --min-quality 25

# Skip specific steps
python nanodemux raw_data/reads.fastq barcodes/primer_well_map.csv \
    --skip-raw-qc --skip-alignment

# With anchor sequence for MSA visualization
python nanodemux raw_data/reads.fastq barcodes/primer_well_map.csv \
    --anchor-seq AATGATACGGCGACCACCGAGATCTACACTATAGCCTTCGTCGGCAGCGTC
```

**All Available Parameters:**

Raw Quality Analysis:
- `--raw-max-reads`: Maximum reads to analyze in raw QC (default: all)
- `--skip-raw-qc`: Skip raw quality analysis step

Demultiplexing:
- `--adapters` / `-a`: Python file defining ADAPTERS list
- `--min-length`: Minimum read length (default: 50)
- `--max-penalty`: Max mismatch penalty, sum of Phred scores (default: 60)
- `--cpus`: Number of CPU cores for parallel processing (default: 1)
- `--flank`: Bases at each end to search for barcodes (default: 100)
- `--var-q`: Phred cutoff for variable-base match confidence (default: 10)
- `--no-fallback`: Disable alignment fallback for ambiguous cases

Quality Reports:
- `--anchor-seq`: Anchor/reference sequence for MSA alignment
- `--demux-max-reads`: Max reads to show in MSA layout (default: 100)
- `--skip-demux-qc`: Skip demultiplexed quality analysis step

Alignment:
- `--align-max-reads`: Max reads to process per well in alignment (default: all)
- `--min-quality`: Min Phred quality for consensus calling (default: 20)
- `--no-align`: Skip MSA, calculate consensus from unaligned sequences
- `--skip-alignment`: Skip alignment step

General:
- `--output` / `-o`: Output directory for all results (default: nanodemux_output)

For complete help: `python nanodemux --help`

**2) Understanding required barcode CSV format**
The barcode loader expects a CSV with at least these columns (headers must match):
- `Sequence Name` — human readable name for the primer/barcode
- `Sequence` — the DNA sequence of the barcode (no whitespace)
- `Well Position` — the well designation (e.g., `A1`)

Example rows (CSV):

```
Sequence Name,Sequence,Well Position
R01,ACGT... ,A1
C01,TTGG... ,1
```

Note: `generate_quality_report.py` also expects a `barcode_stats.csv` inside the demux directory (this file is produced by `demux_barcodes.py` when you run demultiplexing).

**3) Demultiplexing raw FASTQ files (manual approach)**
Use `demux_barcodes.py` to demultiplex and optionally generate a quality report. (For automatic workflow, use the `nanodemux` script described in section 1.)

Basic usage (single FASTQ file):

```bash
python demux_barcodes.py /path/to/raw_reads.fastq barcodes/your_barcodes.csv \
    --adapters barcodes/your_adapters.py --outdir /path/to/outdir/ --cpus 4
```

Options of interest:
- `--adapters` (`-a`): optional adapter definitions file
- `--outdir`: where to write demultiplexed per-well FASTQ files and statistics
- `--flank`: bases at each end to search for barcodes (default 100)
- `--max_penalty`: alignment/mismatch penalty threshold (default 60)
- `--min_length`: minimum read length to keep (default 50)
- `--var_q`: Phred cutoff for variable-base match confidence (default 10)
- `--cpus`: parallelism
- `--report` (`-r`): automatically generate the demultiplexed quality report after processing

Example (with report):

```bash
python demux_barcodes.py benchmarking_data/251205_AR/KD2X25_3_group3.fastq \
    benchmarking_data/251205_AR/demux_barcodes.csv \
    --adapters benchmarking_data/251205_AR/251205_adapters_AR.py \
    --outdir benchmarking_data/251205_AR/demux_output/ \
    --flank 100 --max_penalty 60 --min_length 50 --cpus 4 --report
```

After completion, the demultiplexed per-well FASTQ files and a `barcode_stats.csv` will be placed in the `--outdir` directory (or an auto-generated folder if you omit `--outdir`). If you provided `--report`, a `quality_report/` folder will be created under the output directory containing PNGs and an HTML report.

**4) Generate a demultiplexed quality report (manual)**
If you already have a demultiplexed directory (contains files named like `WELL_reads.fastq` and `barcode_stats.csv`), run `generate_quality_report.py` to make the visual report.

```bash
python generate_quality_report.py /path/to/demux_dir barcodes/your_barcodes.csv --output /path/to/report_dir --max-reads 200
```

- `demux_dir`: directory with demultiplexed FASTQ files (named `*_reads.fastq`) and a `barcode_stats.csv` file
- `barcodes`: the same barcode CSV used for demultiplexing
- `--output` / `-o`: optional custom report directory (default: `demux_dir/quality_report/`)
- `--max-reads`: how many reads to show in the MSA-style layout visualization

Generated files (typical):
- `read_length_distribution.png`
- `msa_layout.png`
- `barcode_positions.png`
- `barcode_summary.png`
- `quality_report.html`

Open the HTML with your browser or the workspace preview to inspect results.

**5) Generate raw FASTQ quality report**
For raw (un-demultiplexed) FASTQ files, use `generate_raw_quality_report.py`:

```bash
python generate_raw_quality_report.py /path/to/raw.fastq --output /path/to/raw_report_dir --max-reads 5000
```

Defaults: if `--output` is omitted, the script creates `raw_data_quality_report/` in the current directory. Outputs include read length histogram, quality distribution plots, base composition, and an HTML report.

**6) Running the benchmark suite**
A benchmarking script automates creating read subsets and running demultiplexing repeatedly.

Quick run with defaults (uses sample configs in `benchmarking_data/`):

```bash
python benchmark_demux.py
```

Options of interest:
- `--subset N`: number of reads per dataset to use (default 1000). Use 0 for full dataset.
- `--full`: shortcut for `--subset 0` (run on full FASTQ)
- `--flank`, `--max-penalty`, `--cpus` — override default parameters
- `--report`: generate a markdown report of benchmark runs
- `--list`: list available benchmark configs found under `benchmarking_data/`

Benchmark outputs are saved under `benchmarking_data/benchmarks/` and a comparison/markdown report is written (e.g. `benchmark_report.md`).

**7) Common workflows (examples)**
1. Quick inspection of raw data

```bash
python generate_raw_quality_report.py benchmarking_data/251205_AR/KD2X25_3_group3.fastq
# open raw_data_quality_report/raw_quality_report.html
```

2. Complete workflow with single command (recommended)

```bash
python nanodemux raw_data/reads.fastq barcodes/251202_primer_well_map_DA.csv \
    --adapters barcodes/251205_adapters_DA.py --cpus 4
# open nanodemux_output/comprehensive_report.html
```

3. Full demultiplex + demux quality report (manual approach)

```bash
python demux_barcodes.py raw_data/reads.fastq barcodes/251202_primer_well_map_DA.csv \
    --adapters barcodes/251205_adapters_DA.py --outdir demuxed_output/ --report
# open demuxed_output/quality_report/quality_report.html
```

4. Run benchmarks for comparison

```bash
python benchmark_demux.py --subset 100 --cpus 4 --report
# read benchmarking_data/benchmarks/benchmark_report.md
```

**8) Output file conventions**
- Demultiplexed reads: files named like `<WELL>_reads.fastq` (script uses `_reads.fastq` suffix to identify wells)
- `barcode_stats.csv`: tabular matrix of reads per well — used by `generate_quality_report.py` to create plate heatmaps
- Quality reports: saved as PNGs and a single HTML report combining them

**9) Troubleshooting & tips**
- Error: `barcode_stats.csv not found in demux_dir` — run `demux_barcodes.py` (it creates the stats file) or place a valid `barcode_stats.csv` in the folder.
- Barcode CSV load issues — confirm headers include `Sequence Name`, `Sequence`, and `Well Position`.
- If `generate_quality_report.py` finds no `*_reads.fastq` files it will warn and skip some plots — ensure the demux output uses the `_reads.fastq` suffix.
- If scripts error with missing Python packages, run `pip install -r requirements.txt`.
- For faster testing during development use the `--subset` flag in `benchmark_demux.py` to operate on fewer reads.

**10) Contributing / Extending**
- Add new adapter sets in `barcodes/` as `.py` files with an `ADAPTERS` list as required by `demux_barcodes.py`.
- To change plotting styles, edit `generate_quality_report.py` or `generate_raw_quality_report.py` (they use `seaborn`/`matplotlib`).

**11) Contact / Support**
- For issues with this repo, open an issue with a clear description, command used, and the error output.

---
_Last updated: 2025-12-08 - Added comprehensive nanodemux workflow script documentation_
