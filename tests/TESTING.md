# Testing Suite Documentation

## Overview

This comprehensive testing suite provides multiple layers of testing for the NanoDemux barcode demultiplexing pipeline:

- **27 test cases** covering unit tests, integration tests, and end-to-end pipeline tests
- **Dual framework support**: unittest and pytest for flexibility
- **100% test passage** with all edge cases covered
- **CI/CD ready** with GitHub Actions workflow included

## Quick Start

### 1. Run All Tests (Recommended)
```bash
make test
# or
python run_tests.py
```

### 2. Run Tests with Coverage
```bash
make test-cov
# Generates HTML report in htmlcov/index.html
```

### 3. Run Quick Tests
```bash
make test-fast
# Minimal output, fastest execution
```

## Test Organization

### Directory Structure
```
tests/
├── __init__.py                    # Package marker
├── conftest.py                    # Pytest fixtures and configuration
├── test_demux.py                  # Unit tests (22 test cases)
├── test_integration.py            # Integration tests (unittest-based)
├── test_pytest_integration.py      # Integration tests (pytest-based)
└── README.md                       # Detailed test documentation
```

### Test Categories

#### Unit Tests (test_demux.py) - 17 tests
- **Helper Functions** (3 tests)
  - Reverse complement calculation
  - Default dictionary initialization
  
- **Barcode Matching** (5 tests)
  - Exact sequence matching
  - Mismatch handling
  - Reverse complement detection
  - Flank size constraints
  - Quality-weighted scoring
  
- **Barcode Loading** (2 tests)
  - CSV parsing
  - Well position mapping
  
- **Merging Functions** (2 tests)
  - Statistics dictionary merging
  - Read set merging
  
- **FASTQ Handling** (3 tests)
  - File I/O validation
  - Sequence length verification
  - Quality score parsing

#### Integration Tests - 10 tests

**unittest-based (7 tests)**
- Full pipeline execution
- Output directory creation
- Stats file generation and validation
- FASTQ output format validation
- Barcode loading consistency
- Chunk processing

**pytest-based (3 tests)**
- Full pipeline with fixtures
- Stats file structure validation
- Multiple flank size configurations
- FASTQ output format verification

## Available Commands

### Using Make (Recommended)
```bash
make test              # Run all tests
make test-fast         # Quick test run
make test-verbose      # Detailed output
make test-cov          # Generate coverage report
make test-unittest     # Only unittest tests
make test-pytest       # Only pytest tests
make clean             # Remove temp files
make help              # Show all targets
```

### Using Python Directly
```bash
# Run all tests
python run_tests.py

# With options
python run_tests.py --verbose              # Verbose output
python run_tests.py --cov                  # With coverage
python run_tests.py --unittest             # Only unittest
python run_tests.py --pytest               # Only pytest

# With unittest directly
python -m unittest discover -s tests -p "test_*.py" -v

# With pytest directly
pytest tests/ -v
pytest tests/ --cov=. --cov-report=html
```

## Test Coverage

The test suite validates:

1. ✅ **Core Algorithm**
   - Quality-weighted barcode matching
   - Reverse complement detection
   - Mismatch penalty calculation

2. ✅ **Data Loading**
   - CSV parsing and validation
   - Barcode-to-position mapping
   - File I/O operations

3. ✅ **Pipeline Execution**
   - Multi-chunk processing
   - Statistics aggregation
   - Output file generation
   - Well-based read organization

4. ✅ **Edge Cases**
   - Short reads (< min_length)
   - Reads with no barcodes
   - Reads with single barcode
   - Multiple barcode ambiguity
   - Quality score variations

5. ✅ **Parameter Validation**
   - Different flank sizes
   - Varying penalty thresholds
   - CPU core configurations

## Expected Test Results

All tests should pass with output similar to:

```
============================== 27 passed in 0.29s ==============================
unittest             ✅ PASSED
pytest               ✅ PASSED

✅ All tests passed!
```

## Continuous Integration

The repository includes GitHub Actions workflow (`.github/workflows/tests.yml`) that:

- Runs tests on Python 3.9, 3.10, 3.11, 3.12
- Executes on every push to main/develop
- Runs on all pull requests to main
- Generates coverage reports
- Uploads to Codecov (optional)

## Adding New Tests

When adding new functionality:

1. **Write a test first** (Test-Driven Development)
2. **Add to appropriate test file** or create new one
3. **Use descriptive names** following existing patterns
4. **Include docstrings** explaining what is tested
5. **Run full suite** before committing

### Test Template

```python
def test_new_feature(self):
    """Test description of what is being validated."""
    # Setup: Initialize test data
    input_data = prepare_test_input()
    
    # Execute: Run the code being tested
    result = function_under_test(input_data)
    
    # Assert: Verify the result is correct
    self.assertEqual(result, expected_value)
```

## Debugging Failed Tests

### Run with verbose output
```bash
python -m pytest tests/ -v -s
```

### Stop on first failure
```bash
python -m pytest tests/ -x
```

### Run specific test
```bash
python -m unittest tests.test_demux.TestBarcodeMatching.test_exact_match
```

### Debug with pdb
```bash
python -m pytest tests/ --pdb
```

### Show local variables on failure
```bash
python -m pytest tests/ -l
```

## Performance Notes

- Full test suite runs in **~0.3 seconds**
- No external network dependencies
- Tests use temporary directories automatically cleaned up
- Suitable for CI/CD pipelines

## Dependencies

The test suite requires:
- `pytest>=7.0.0` - For pytest-based tests
- `pytest-cov>=3.0.0` - For coverage reports
- `biopython>=1.79` - For sequence handling
- `pandas>=1.3.0` - For CSV operations

Install with:
```bash
pip install -r requirements.txt
```

## Future Enhancements

- [ ] Add performance benchmarks
- [ ] Property-based testing (hypothesis)
- [ ] Integration with real large FASTQ files
- [ ] Memory profiling tests
- [ ] Mutation testing for code quality
- [ ] Load testing with large barcode sets

## Support

For test failures or questions:
1. Check the detailed test output with `-v` flag
2. Review test documentation in `tests/README.md`
3. Inspect test implementation in relevant test file
4. Check the main code comments for implementation details
