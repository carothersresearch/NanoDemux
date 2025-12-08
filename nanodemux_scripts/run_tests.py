#!/usr/bin/env python3
"""
Test runner script - run all tests with various options.
"""

import subprocess
import sys
import argparse
from pathlib import Path


def run_unittest(verbose=False):
    """Run unittest-based tests."""
    print("\n" + "=" * 70)
    print("Running unittest-based tests...")
    print("=" * 70)
    
    cmd = [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"]
    if verbose:
        cmd.append("-v")
    
    result = subprocess.run(cmd)
    return result.returncode == 0


def run_pytest(verbose=False, cov=False):
    """Run pytest-based tests."""
    print("\n" + "=" * 70)
    print("Running pytest-based tests...")
    print("=" * 70)
    
    cmd = [sys.executable, "-m", "pytest", "tests/", "-v" if verbose else "-q"]
    if cov:
        cmd.extend(["--cov=.", "--cov-report=html"])
    
    result = subprocess.run(cmd)
    return result.returncode == 0


def run_all_tests(verbose=False, cov=False):
    """Run all test suites."""
    results = {
        "unittest": run_unittest(verbose),
        "pytest": run_pytest(verbose, cov),
    }
    
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)
    
    for suite, passed in results.items():
        status = "✅ PASSED" if passed else "❌ FAILED"
        print(f"{suite:20} {status}")
    
    all_passed = all(results.values())
    
    if all_passed:
        print("\n✅ All tests passed!")
        return 0
    else:
        print("\n❌ Some tests failed!")
        return 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run test suite for NanoDemux")
    parser.add_argument("--unittest", action="store_true", help="Run only unittest tests")
    parser.add_argument("--pytest", action="store_true", help="Run only pytest tests")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--cov", action="store_true", help="Generate coverage report")
    
    args = parser.parse_args()
    
    # Determine which tests to run
    if args.unittest:
        exit_code = 0 if run_unittest(args.verbose) else 1
    elif args.pytest:
        exit_code = 0 if run_pytest(args.verbose, args.cov) else 1
    else:
        exit_code = run_all_tests(args.verbose, args.cov)
    
    sys.exit(exit_code)
