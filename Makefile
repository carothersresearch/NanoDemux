.PHONY: help test test-fast test-cov test-verbose test-unittest test-pytest install lint clean benchmark benchmark-fast benchmark-full benchmark-compare report

help:
	@echo "NanoDemux Makefile targets:"
	@echo ""
	@echo "  make install          Install all dependencies including dev/test dependencies"
	@echo "  make test             Run full test suite (unittest + pytest)"
	@echo "  make test-fast        Run tests with minimal output"
	@echo "  make test-verbose     Run tests with verbose output"
	@echo "  make test-cov         Run tests with coverage report (HTML)"
	@echo "  make test-unittest    Run only unittest tests"
	@echo "  make test-pytest      Run only pytest tests"
	@echo "  make benchmark        Run benchmarking suite on sample data (1000 reads)"
	@echo "  make benchmark-fast   Quick benchmark (100 reads, ~5 seconds)"
	@echo "  make benchmark-full   Full dataset benchmark (slow, several minutes)"
	@echo "  make benchmark-compare Compare recent benchmark results"
	@echo "  make report           Generate quality report for example demux output"
	@echo "  make lint             Run code quality checks"
	@echo "  make clean            Remove temporary files and cache"
	@echo ""

install:
	pip install -r requirements.txt

test:
	python nanodemux_scripts/run_tests.py

test-fast:
	python -m pytest tests/ -q

test-verbose:
	python nanodemux_scripts/run_tests.py --verbose

test-cov:
	python nanodemux_scripts/run_tests.py --cov
	@echo ""
	@echo "Coverage report generated in htmlcov/index.html"

test-unittest:
	python nanodemux_scripts/run_tests.py --unittest

test-pytest:
	python nanodemux_scripts/run_tests.py --pytest

lint:
	python -m flake8 nanodemux_scripts/demux_barcodes.py tests/ --max-line-length=100 || true
	python -m pylint nanodemux_scripts/demux_barcodes.py --disable=all --enable=E,F || true

benchmark:
	python nanodemux_scripts/benchmark_demux.py

benchmark-fast:
	python nanodemux_scripts/benchmark_demux.py --subset 100

benchmark-full:
	python nanodemux_scripts/benchmark_demux.py --full

benchmark-compare:
	python nanodemux_scripts/benchmark_demux.py --compare 5

report:
	@echo "Generating quality report for example demux output..."
	@if [ -d "demplex_data/55XPXK_1_P4_323_EG" ]; then \
		python nanodemux_scripts/generate_quality_report.py demplex_data/55XPXK_1_P4_323_EG/ barcodes/251202_primer_well_map_DA.csv; \
	else \
		echo "Error: Example demux output not found. Run demultiplexing first:"; \
		echo "  python nanodemux_scripts/demux_barcodes.py raw_data/55XPXK_1_P4_323_EG.fastq barcodes/251202_primer_well_map_DA.csv --cpus 2"; \
	fi

clean:
	find . -type d -name __pycache__ -exec rm -rf {} + || true
	find . -type f -name "*.pyc" -delete
	rm -rf .pytest_cache
	rm -rf htmlcov
	rm -rf .coverage
	rm -rf *.egg-info
	rm -rf dist build

.DEFAULT_GOAL := help
