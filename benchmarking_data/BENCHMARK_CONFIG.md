# Benchmark Configuration Guide

## Overview

The benchmarking system uses JSON configuration files to define parameters for running demultiplexing benchmarks. Each subfolder in `benchmarking_data/` should have a single `benchmark_config.json` file containing one or more configurations.

## Configuration File Format

Each subfolder should have exactly one `benchmark_config.json` file with the following structure:

```json
{
  "configs": [
    {
      "name": "config_name",
      "description": "Description of this configuration",
      "barcodes": "demux_barcodes.csv",
      "adapters": "adapters.py",
      "flank": 100,
      "max_penalty": 60,
      "cpus": 4,
      "min_length": 50
    },
    {
      "name": "another_config",
      "description": "Another configuration with different parameters",
      "barcodes": "demux_barcodes.csv",
      "adapters": "adapters.py",
      "flank": 150,
      "max_penalty": 40,
      "cpus": 4,
      "min_length": 100
    }
  ]
}
```

### Configuration Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `name` | string | Configuration identifier | required |
| `description` | string | Human-readable description | "" |
| `barcodes` | string | Path to barcode CSV file (relative to config) | required |
| `adapters` | string | Path to adapter Python file (relative to config) | optional |
| `flank` | integer | Flank size for barcode matching | 100 |
| `max_penalty` | integer | Maximum penalty for fuzzy matching | 60 |
| `cpus` | integer | Number of CPU cores to use | 4 |
| `min_length` | integer | Minimum read length to process | 50 |

## Directory Structure

```
benchmarking_data/
├── 251205_AR/
│   ├── benchmark_config.json          # Single config file
│   ├── demux_barcodes.csv            # Barcode definitions
│   ├── 251205_adapters_AR.py         # Adapter sequences
│   ├── KD2X25_1_group1.fastq         # Data file 1
│   ├── KD2X25_3_group3.fastq         # Data file 2
│   └── KD2X25_6_group6.fastq         # Data file 3
├── firstpass/
│   ├── benchmark_config.json          # Config with multiple configs inside
│   ├── 251202_primer_well_map_DA.csv
│   ├── 251205_adapters_DA.py
│   ├── 55XPXK_1_P4_323_EG.fastq
│   └── VL69M6_1_P4_323_full.fastq
└── benchmarks/                        # Output directory
```

## Multiple Configurations per Folder

To test different parameter sets on the same data, add multiple configuration objects to the `configs` array:

```json
{
  "configs": [
    {
      "name": "default",
      "description": "Standard parameters",
      "barcodes": "barcodes.csv",
      "adapters": "adapters.py",
      "flank": 100,
      "max_penalty": 60,
      "cpus": 4,
      "min_length": 50
    },
    {
      "name": "high_stringency",
      "description": "Stricter matching parameters",
      "barcodes": "barcodes.csv",
      "adapters": "adapters.py",
      "flank": 150,
      "max_penalty": 40,
      "cpus": 4,
      "min_length": 100
    },
    {
      "name": "relaxed",
      "description": "More permissive matching",
      "barcodes": "barcodes.csv",
      "adapters": "adapters.py",
      "flank": 80,
      "max_penalty": 80,
      "cpus": 4,
      "min_length": 30
    }
  ]
}
```

The benchmark script will run ALL configurations on ALL FASTQ files in the directory.

## Usage Examples

### List Available Configurations
```bash
python benchmark_demux.py --list
```

### Run Benchmarks with Config Parameters
```bash
# Use parameters from config files (default behavior)
python benchmark_demux.py --subset 100

# Full dataset
python benchmark_demux.py --full
```

### Override Config Parameters
```bash
# Override ALL config parameters with command line arguments
python benchmark_demux.py --subset 100 --override-params --flank 150 --max-penalty 40
```

### Compare Results
```bash
# Compare last 10 benchmark runs
python benchmark_demux.py --compare 10

# Generate detailed report
python benchmark_demux.py --report
```

## Output

Each benchmark run creates:
- A timestamped output directory in `benchmarking_data/benchmarks/`
- A `barcode_stats.csv` file with demultiplexing statistics
- Entry in `benchmark_results.json` with full run details

## Best Practices

1. **One file per folder**: Keep all configs in a single `benchmark_config.json` per subfolder
2. **Descriptive names**: Use clear config names like `default`, `high_stringency`, `low_penalty`
3. **Document parameters**: Use the `description` field to explain what each config tests
4. **Systematic testing**: Start with default params, then create variants to test sensitivity
5. **Version control**: Commit config files to track parameter changes over time

## Interpreting Results

The comparison table shows:
- **Config**: Which configuration was used
- **Map %**: Percentage of reads successfully mapped to wells
- **F/P/L**: Flank/Max Penalty/Min Length parameters
- **Wells**: Number of wells with reads

Higher mapping rates generally indicate better parameter tuning for your data.

## Adding New Data

To benchmark new data:

1. Create subfolder in `benchmarking_data/`
2. Add FASTQ files
3. Add barcode CSV and adapter Python files
4. Create `benchmark_config.json` with one or more configs
5. Run `python benchmark_demux.py --list` to verify
6. Run `python benchmark_demux.py --subset 100` to test
