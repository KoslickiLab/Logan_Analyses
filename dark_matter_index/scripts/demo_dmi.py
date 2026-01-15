#!/usr/bin/env python3
"""
Demonstration script for Dark Matter Index calculation.

This script shows how to use the DMI calculator with:
1. A sourmash signature file (.sig)
2. Raw hash arrays
3. Hashes fetched from DuckDB

Run this after building the reference hash set.
"""

import numpy as np
import json
from pathlib import Path
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from dmi_calculator import DMICalculator, DMIResult, load_hashes_from_duckdb


def demo_with_signature_file():
    """Demonstrate DMI calculation from a signature file."""
    print("="*60)
    print("DEMO: Computing DMI from a Sourmash Signature File")
    print("="*60)
    
    # Example signature file (adjust path as needed)
    sig_file = "/mnt/user-data/uploads/c8c92aa83f79218d31da2e9151860fd4.sig"
    
    # Load and display signature info
    with open(sig_file, 'r') as f:
        sig_data = json.load(f)
    
    if isinstance(sig_data, list):
        sig_data = sig_data[0]
    
    sample_name = sig_data.get('name', 'Unknown')
    hashes = sig_data['signatures'][0]['mins']
    ksize = sig_data['signatures'][0]['ksize']
    
    print(f"Sample name: {sample_name}")
    print(f"K-mer size: {ksize}")
    print(f"Number of hashes: {len(hashes):,}")
    print(f"First 5 hashes: {hashes[:5]}")
    print(f"Hash value range: {min(hashes):,} to {max(hashes):,}")
    print()
    
    # If reference exists, compute DMI
    reference_path = Path("/scratch/dmk333_new/dmi_reference/reference_hashes_k31.bin")
    
    if reference_path.exists():
        print("Computing DMI...")
        calculator = DMICalculator(str(reference_path), use_mmap=True)
        result = calculator.compute_dmi_from_signature(sig_file, ksize=31)
        
        print(f"\nResults for {result.sample_name}:")
        print(f"  DMI: {result.dmi:.4f} ({result.dmi*100:.2f}% dark matter)")
        print(f"  Total hashes: {result.total_hashes:,}")
        print(f"  Mapped to reference: {result.mapped_hashes:,}")
        print(f"  Unmapped (dark matter): {result.unmapped_hashes:,}")
    else:
        print(f"Reference file not found at {reference_path}")
        print("Run build_reference_hashes.py first to create the reference set.")
        print()
        print("To simulate DMI calculation, here's what would happen:")
        print("1. Load reference hashes (~160GB) as memory-mapped array")
        print("2. For each sample hash, binary search in reference")
        print("3. Count matches vs non-matches")
        print("4. DMI = unmapped / total")


def demo_with_raw_hashes():
    """Demonstrate DMI calculation from raw hash arrays."""
    print("\n" + "="*60)
    print("DEMO: Computing DMI from Raw Hash Arrays")
    print("="*60)
    
    # Create some example hashes
    np.random.seed(42)
    sample_hashes = np.random.randint(0, 2**63, size=1000, dtype=np.uint64)
    
    print(f"Generated {len(sample_hashes):,} random hashes")
    
    reference_path = Path("/scratch/dmk333_new/dmi_reference/reference_hashes_k31.bin")
    
    if reference_path.exists():
        calculator = DMICalculator(str(reference_path), use_mmap=True)
        result = calculator.compute_dmi(sample_hashes, sample_name="random_test")
        
        print(f"\nResults:")
        print(f"  DMI: {result.dmi:.4f}")
        print(f"  (Random hashes should have DMI ≈ 1.0 since they're unlikely to match)")
    else:
        print("Reference file not found. Skipping actual computation.")


def demo_architecture():
    """Explain the DMI computation architecture."""
    print("\n" + "="*60)
    print("ARCHITECTURE OVERVIEW")
    print("="*60)
    
    print("""
    The Dark Matter Index (DMI) measures the fraction of k-mer hashes in a 
    metagenomic sample that do NOT appear in any reference database.
    
    DMI = |sample_hashes - reference_hashes| / |sample_hashes|
    
    DATA STRUCTURES:
    ================
    
    Reference Hash Set:
    - ~20 billion unique 64-bit hashes (160 GB)
    - Stored as sorted binary file (numpy uint64 array)
    - Memory-mapped for efficient access
    - Binary search: O(log n) per query ≈ 34 comparisons
    
    DATABASES INCLUDED:
    ==================
    1. GTDB (Genome Taxonomy Database)
    2. GenBank WGS (Whole Genome Shotgun)
    3. GenBank Genomes (Complete genomes)
    4. AllTheBacteria
    5. SILVA (ribosomal RNA)
    6. GenBank TLS (Targeted Locus Study)
    7. GenBank TSA (Transcriptome Shotgun Assembly)
    8. Logan Obelisks (RNA elements)
    9. Logan Plasmids
    10. Serratus Viruses
    
    COMPUTATION FLOW:
    =================
    
    For a single sample:
    1. Load sample hashes (from .sig file or DuckDB)
    2. Memory-map reference array (no RAM needed for 160GB file)
    3. np.searchsorted() for vectorized binary search
    4. Count matches vs non-matches
    5. Return DMI = unmapped / total
    
    For batch processing:
    1. Load input parquet with sample metadata
    2. Parallel workers each memory-map the reference
    3. Each worker processes chunk of samples
    4. Results merged and saved to output parquet
    
    PERFORMANCE:
    ============
    
    - Reference loading: O(1) with memory mapping
    - Per-hash lookup: O(log 20B) ≈ 34 comparisons
    - Per-sample (100K hashes): ~10-50 ms
    - Batch (1M samples): ~3-10 hours with 32 workers
    
    INTERPRETATION:
    ===============
    
    DMI = 0.0: All hashes match known references
    DMI = 0.1: 10% of sample is "dark matter"
    DMI = 0.5: Half the sample is uncharacterized
    DMI = 1.0: Nothing matches (100% dark matter)
    
    High DMI suggests:
    - Novel organisms not in databases
    - Unusual environment with uncharacterized microbes
    - Possible sequencing artifacts (check with other QC)
    """)


def show_usage_examples():
    """Show command-line usage examples."""
    print("\n" + "="*60)
    print("COMMAND-LINE USAGE EXAMPLES")
    print("="*60)
    
    print("""
    # Step 1: Build the reference hash set (one-time, ~1-2 hours)
    # ============================================================
    
    python build_reference_hashes.py \\
        DB_info.json \\
        /scratch/dmk333_new/dmi_reference/ \\
        --chunk-size 1.0 \\
        --ksize 31
    
    # Step 2a: Compute DMI for a single signature file
    # =================================================
    
    python dmi_calculator.py \\
        /scratch/dmk333_new/dmi_reference/reference_hashes_k31.bin \\
        sample.sig
    
    # Step 2b: Batch compute DMI for all SRA samples
    # ===============================================
    
    python batch_compute_dmi.py \\
        --reference /scratch/dmk333_new/dmi_reference/reference_hashes_k31.bin \\
        --input /path/to/filtered_data.parquet \\
        --output /path/to/filtered_data_with_dmi.parquet \\
        --database /scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db \\
        --workers 32 \\
        --chunk-size 1000
    
    # Step 3: Analyze results
    # =======================
    
    python -c "
    import pandas as pd
    df = pd.read_parquet('filtered_data_with_dmi.parquet')
    print(df.groupby('organism')['dmi'].describe())
    "
    """)


if __name__ == "__main__":
    demo_with_signature_file()
    demo_with_raw_hashes()
    demo_architecture()
    show_usage_examples()
