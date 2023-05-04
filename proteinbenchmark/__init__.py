"""Benchmarks for OpenFF protein force fields"""

# Add imports here
# Handle versioneer
from ._version import get_versions
from .proteinbenchmark import *

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
