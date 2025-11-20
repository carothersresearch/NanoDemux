# Docker Usage Guide for NanoDemux

This guide explains how to build and run NanoDemux using Docker.

## Prerequisites

- Docker installed on your system
- FASTQ input files and barcode CSV file

## Building the Docker Image

From the repository root directory, build the Docker image:

```bash
docker build -t nanodemux .
```

This will create a Docker image named `nanodemux` with all required dependencies (Python, Biopython, Pandas).

## Running NanoDemux with Docker

### Basic Usage

Mount your data directory and run the demultiplexing:

```bash
docker run --rm -v /path/to/your/data:/data nanodemux \
  input.fastq barcodes.csv --outdir output
```

### With Multiple CPUs

Use the `--cpus` flag to enable parallel processing:

```bash
docker run --rm -v /path/to/your/data:/data nanodemux \
  input.fastq barcodes.csv --outdir output --cpus 4
```

### Full Example with All Parameters

```bash
docker run --rm -v /path/to/your/data:/data nanodemux \
  reads.fastq primer_well_map.csv \
  --outdir demuxed_output \
  --min_length 50 \
  --max_penalty 60 \
  --cpus 4 \
  --flank 50
```

### Example with Repository Sample Data

From the repository directory:

```bash
docker run --rm -v "$(pwd)":/data nanodemux \
  55XPXK_1_P4_323_EG.fastq primer_well_map.csv \
  --outdir output --cpus 2
```

## Volume Mounting

The Docker container expects data files to be mounted at `/data`. The `-v` flag mounts your host directory:

- `-v /host/path:/data` - Mounts `/host/path` from your host system to `/data` in the container
- Use `$(pwd)` to mount the current directory
- On Windows (PowerShell), use `${PWD}` instead of `$(pwd)`

## Output Files

Output files will be written to the specified output directory within your mounted volume:

- `<outdir>/<well>_reads.fastq` - One FASTQ file per well (e.g., A1_reads.fastq)
- `<outdir>/barcode_stats.csv` - Statistics summary grid

## Permissions

Note: Files created by the Docker container will be owned by the root user. You may need to change ownership after running:

```bash
sudo chown -R $USER:$USER output/
```

Or run Docker with user permissions:

```bash
docker run --rm -u $(id -u):$(id -g) -v "$(pwd)":/data nanodemux \
  input.fastq barcodes.csv --outdir output
```

## Getting Help

View the help message:

```bash
docker run --rm nanodemux --help
```

## Troubleshooting

### Permission Denied Errors

If you encounter permission errors accessing output files, use the `-u` flag as shown above to run as your user.

### File Not Found Errors

Make sure your input files are in the directory you're mounting, and use the filename as it appears in the mounted directory (not the full host path).

### Out of Memory

If processing large files, you may need to allocate more memory to Docker:

```bash
docker run --rm -m 8g -v "$(pwd)":/data nanodemux \
  large_file.fastq barcodes.csv --outdir output
```

## Advanced: Building with a Different Python Version

To use a different Python version, edit the Dockerfile's first line:

```dockerfile
FROM python:3.12-slim
```

Then rebuild the image.
