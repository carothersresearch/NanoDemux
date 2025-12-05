.PHONY: help test test-fast test-cov test-verbose test-unittest test-pytest install lint clean benchmark benchmark-fast benchmark-full benchmark-compare

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
	@echo "  make lint             Run code quality checks"
	@echo "  make clean            Remove temporary files and cache"
	@echo ""

install:
	pip install -r requirements.txt

test:
	python run_tests.py

test-fast:
	python -m pytest tests/ -q

test-verbose:
	python run_tests.py --verbose

test-cov:
	python run_tests.py --cov
	@echo ""
	@echo "Coverage report generated in htmlcov/index.html"

test-unittest:
	python run_tests.py --unittest

test-pytest:
	python run_tests.py --pytest

lint:
	python -m flake8 demux_barcodes.py tests/ --max-line-length=100 || true
	python -m pylint demux_barcodes.py --disable=all --enable=E,F || true

benchmark:
	python benchmark_demux.py

benchmark-fast:
	python benchmark_demux.py --subset 100

benchmark-full:
	python benchmark_demux.py --full

benchmark-compare:
	python benchmark_demux.py --compare 5

clean:
	find . -type d -name __pycache__ -exec rm -rf {} + || true
	find . -type f -name "*.pyc" -delete
	rm -rf .pytest_cache
	rm -rf htmlcov
	rm -rf .coverage
	rm -rf *.egg-info
	rm -rf dist build

.DEFAULT_GOAL := help
