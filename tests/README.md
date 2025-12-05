# NanoDemux Testing Suite

This directory contains comprehensive unit and integration tests for the NanoDemux barcode demultiplexing pipeline.

## Test Structure

### Test Files

- **`test_demux.py`** - Unit tests for individual functions
  - `TestHelperFunctions` - Tests for utility functions (revcomp, defaultdict)
  - `TestBarcodeMatching` - Tests for barcode matching algorithm
  - `TestBarcodeLoading` - Tests for loading barcodes from CSV
  - `TestMergingFunctions` - Tests for dict/read merging
  - `TestFASTQHandling` - Tests for FASTQ I/O

- **`test_integration.py`** - Integration tests (unittest-based)
  - `TestDemultiplexingPipeline` - Full pipeline tests
  - `TestProcessChunk` - Tests for chunked processing

- **`test_pytest_integration.py`** - Integration tests (pytest-based)
  - `TestDemultiplexingWithPytest` - Fixture-based pipeline tests

- **`conftest.py`** - Pytest fixtures and configuration
  - `temp_directory` - Provides temporary test directories
  - `sample_barcodes_csv` - Creates sample barcode CSV
  - `sample_fastq` - Creates sample FASTQ file
  - `output_directory` - Provides output directory

## Running Tests

### Option 1: Using the test runner script

```bash
# Run all tests
python run_tests.py

# Run with verbose output
python run_tests.py --verbose

# Run with coverage report
python run_tests.py --cov

# Run only unittest tests
python run_tests.py --unittest

# Run only pytest tests
python run_tests.py --pytest
```

### Option 2: Run unittest directly

```bash
# Run all unit tests
python -m unittest discover -s tests -p "test_*.py" -v

# Run specific test file
python -m unittest tests.test_demux -v

# Run specific test class
python -m unittest tests.test_demux.TestBarcodeMatching -v

# Run specific test
python -m unittest tests.test_demux.TestBarcodeMatching.test_exact_match -v
```

### Option 3: Run pytest directly

```bash
# Run all pytest tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=. --cov-report=html

# Run specific test file
pytest tests/test_pytest_integration.py -v

# Run specific test
pytest tests/test_pytest_integration.py::TestDemultiplexingWithPytest::test_full_pipeline -v
```

## Test Coverage

The test suite covers:

1. **Helper Functions**
   - Reverse complement calculation
   - Default dictionary initialization

2. **Barcode Matching Algorithm**
   - Exact sequence matching
   - Quality-weighted mismatch scoring
   - Reverse complement matching
   - Flank size constraints

3. **Barcode CSV Loading**
   - CSV parsing
   - Well position mapping
   - Row/column barcode separation

4. **Dictionary Merging**
   - Merging statistics from multiple chunks
   - Merging read sets

5. **FASTQ Handling**
   - File reading/writing
   - Sequence length validation
   - Quality score parsing

6. **Full Pipeline**
   - End-to-end demultiplexing
   - Output file generation
   - Statistics reporting
   - Different parameter configurations

## Test Fixtures

The test suite uses fixtures (temporary files/directories) to ensure tests are:
- **Isolated** - Each test has its own working directory
- **Reproducible** - Same results every time
- **Clean** - Temporary files are cleaned up automatically

## Continuous Integration

To add this to CI/CD (GitHub Actions, GitLab CI, etc.):

```yaml
# Example GitHub Actions workflow
- name: Run tests
  run: |
    pip install -r requirements.txt
    python run_tests.py --cov
    
- name: Upload coverage
  uses: codecov/codecov-action@v3
```

## Adding New Tests

When adding new functionality:

1. **Write tests first** (Test-Driven Development)
2. **Add unit tests** to `test_demux.py` for individual functions
3. **Add integration tests** to test the full pipeline
4. **Use fixtures** in pytest for shared resources
5. **Run full test suite** before committing

Example test template:

```python
def test_new_feature(self):
    """Test description."""
    # Setup
    input_data = "test_input"
    
    # Execute
    result = function_under_test(input_data)
    
    # Assert
    self.assertEqual(result, expected_value)
```

## Debugging Tests

```bash
# Run with print statements visible
python -m pytest tests/ -s -v

# Run with PDB debugger
python -m pytest tests/ --pdb

# Stop on first failure
python -m pytest tests/ -x

# Show local variables on failure
python -m pytest tests/ -l
```

## Coverage Reports

After running `python run_tests.py --cov`, view the HTML report:

```bash
open htmlcov/index.html  # macOS
xdg-open htmlcov/index.html  # Linux
start htmlcov/index.html  # Windows
```

Target: Maintain >80% code coverage for critical functions.

## Known Issues / Future Improvements

- [ ] Add performance benchmarks
- [ ] Add memory usage tests
- [ ] Test with large real FASTQ files
- [ ] Add mutation testing
- [ ] Add property-based tests (hypothesis)
