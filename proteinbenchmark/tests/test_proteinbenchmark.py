"""
Unit and regression test for the proteinbenchmark package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import proteinbenchmark


def test_proteinbenchmark_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "proteinbenchmark" in sys.modules
