"""
Dark Matter Index (DMI) Calculator Package

A toolkit for computing the fraction of k-mer hashes in metagenomic samples
that do not appear in reference databases.

Main classes:
    DMICalculator: Core DMI computation engine
    DMIResult: Container for DMI results
    DMIBatchProcessor: Batch processing for many samples

Functions:
    compute_dmi_for_sample: Convenience function for single samples
    load_hashes_from_duckdb: Load sample hashes from DuckDB
    batch_compute_dmi_parallel: Parallel batch processing

Example:
    from dark_matter_index import DMICalculator
    
    calc = DMICalculator('/path/to/reference_hashes_k31.bin')
    result = calc.compute_dmi_from_signature('sample.sig')
    print(f"DMI: {result.dmi:.4f}")
"""

from .dmi_calculator import (
    DMICalculator,
    DMIResult,
    compute_dmi_for_sample,
    load_hashes_from_duckdb
)

__version__ = "1.0.0"
__author__ = "David Koslicki Lab"
