"""
Unit and regression test for the benchmarkpl package.
"""

# Import package, test suite, and other packages as needed
import benchmarkpl
import pytest
import sys

def test_benchmarkpl_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "benchmarkpl" in sys.modules
