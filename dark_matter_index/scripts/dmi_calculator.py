#!/usr/bin/env python3
"""
Dark Matter Index (DMI) Calculator

This module provides functionality to compute the Dark Matter Index for metagenomic
samples. The DMI represents the fraction of k-mer hashes in a sample that do NOT
appear in any reference database, indicating potentially novel or uncharacterized
genetic material.

DMI = |sample_hashes - reference_hashes| / |sample_hashes|

A DMI of 0.0 means all hashes are explained by references.
A DMI of 1.0 means no hashes match any reference (100% "dark matter").

Usage:
    # Initialize calculator with reference hash file
    calculator = DMICalculator('/path/to/reference_hashes.bin')
    
    # Compute DMI for a sample
    dmi, stats = calculator.compute_dmi(sample_hashes)
    
    # Or from a sourmash signature file
    dmi, stats = calculator.compute_dmi_from_signature('/path/to/sample.sig')

Author: David Koslicki lab
Date: 2026
"""

import numpy as np
import json
import gzip
from pathlib import Path
from typing import Union, Tuple, Dict, List, Optional
from dataclasses import dataclass
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class DMIResult:
    """Container for DMI computation results."""
    dmi: float                      # Dark Matter Index (0.0 to 1.0)
    total_hashes: int               # Total hashes in sample
    unmapped_hashes: int            # Hashes not in reference
    mapped_hashes: int              # Hashes found in reference
    sample_name: Optional[str] = None
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for easy serialization."""
        return {
            'dmi': self.dmi,
            'total_hashes': self.total_hashes,
            'unmapped_hashes': self.unmapped_hashes,
            'mapped_hashes': self.mapped_hashes,
            'mapped_fraction': 1.0 - self.dmi,
            'sample_name': self.sample_name
        }


class DMICalculator:
    """
    Calculator for Dark Matter Index using a sorted reference hash set.
    
    The reference hashes are stored as a sorted numpy array of uint64 values,
    enabling efficient O(log n) membership testing via binary search.
    
    For 20 billion hashes:
    - File size: ~160 GB
    - Memory (if loaded): ~160 GB
    - Memory-mapped mode: minimal RAM footprint
    - Lookup time per hash: O(log 20B) ≈ 34 comparisons
    
    Attributes:
        reference_hashes: numpy array of sorted uint64 reference hashes
        n_reference: number of reference hashes
        use_mmap: whether using memory-mapped mode
    """
    
    def __init__(
        self, 
        reference_path: Union[str, Path],
        use_mmap: bool = True,
        preload: bool = False
    ):
        """
        Initialize the DMI calculator with a reference hash file.
        
        Args:
            reference_path: Path to the binary file containing sorted uint64 hashes
            use_mmap: If True, memory-map the file (low RAM, slightly slower)
                     If False, load entirely into RAM (faster, needs ~160GB+ RAM)
            preload: If True with mmap, preload pages into memory for faster access
        """
        self.reference_path = Path(reference_path)
        self.use_mmap = use_mmap
        
        if not self.reference_path.exists():
            raise FileNotFoundError(f"Reference file not found: {reference_path}")
        
        # Determine file size and number of hashes
        file_size = self.reference_path.stat().st_size
        if file_size % 8 != 0:
            raise ValueError(f"Invalid reference file: size {file_size} not divisible by 8")
        
        self.n_reference = file_size // 8
        logger.info(f"Loading reference with {self.n_reference:,} hashes ({file_size/1e9:.2f} GB)")
        
        # Load or memory-map the reference hashes
        if use_mmap:
            # Memory-mapped mode: minimal RAM usage
            self.reference_hashes = np.memmap(
                self.reference_path,
                dtype=np.uint64,
                mode='r',
                shape=(self.n_reference,)
            )
            if preload:
                logger.info("Preloading memory-mapped pages...")
                # Touch all pages to bring into memory
                _ = self.reference_hashes.sum()
            logger.info("Reference loaded in memory-mapped mode")
        else:
            # Full load mode: faster lookups but needs RAM
            logger.info("Loading reference fully into RAM...")
            self.reference_hashes = np.fromfile(
                self.reference_path,
                dtype=np.uint64
            )
            logger.info("Reference loaded into RAM")
        
        # Verify sorted order (check first and last few elements)
        self._verify_sorted()
    
    def _verify_sorted(self, n_check: int = 1000):
        """Verify the reference array is sorted (spot check)."""
        if self.n_reference < n_check:
            check_indices = np.arange(self.n_reference - 1)
        else:
            check_indices = np.linspace(0, self.n_reference - 2, n_check, dtype=int)
        
        for i in check_indices:
            if self.reference_hashes[i] > self.reference_hashes[i + 1]:
                raise ValueError(f"Reference hashes not sorted at index {i}")
        
        logger.debug("Reference sort order verified")
    
    def hash_in_reference(self, hash_value: int) -> bool:
        """
        Check if a single hash exists in the reference set.
        
        Uses binary search: O(log n) time complexity.
        
        Args:
            hash_value: 64-bit hash to check
            
        Returns:
            True if hash exists in reference, False otherwise
        """
        idx = np.searchsorted(self.reference_hashes, hash_value)
        return idx < self.n_reference and self.reference_hashes[idx] == hash_value
    
    def count_mapped_hashes(self, sample_hashes: np.ndarray) -> int:
        """
        Count how many sample hashes exist in the reference set.
        
        Uses vectorized binary search for efficiency.
        
        Args:
            sample_hashes: numpy array of uint64 sample hashes
            
        Returns:
            Number of sample hashes found in reference
        """
        # Ensure proper dtype
        if sample_hashes.dtype != np.uint64:
            sample_hashes = sample_hashes.astype(np.uint64)
        
        # Use searchsorted to find insertion points
        # This is O(m log n) where m = sample size, n = reference size
        indices = np.searchsorted(self.reference_hashes, sample_hashes)
        
        # Check which insertions land on actual matches
        # Need to handle edge case where index == n_reference
        valid_indices = indices < self.n_reference
        
        # Create boolean mask for matches
        matches = np.zeros(len(sample_hashes), dtype=bool)
        matches[valid_indices] = (
            self.reference_hashes[indices[valid_indices]] == sample_hashes[valid_indices]
        )
        
        return int(np.sum(matches))
    
    def get_unmapped_hashes(self, sample_hashes: np.ndarray) -> np.ndarray:
        """
        Return the subset of sample hashes NOT in the reference set.
        
        Args:
            sample_hashes: numpy array of uint64 sample hashes
            
        Returns:
            numpy array of hashes not found in reference
        """
        if sample_hashes.dtype != np.uint64:
            sample_hashes = sample_hashes.astype(np.uint64)
        
        indices = np.searchsorted(self.reference_hashes, sample_hashes)
        valid_indices = indices < self.n_reference
        
        matches = np.zeros(len(sample_hashes), dtype=bool)
        matches[valid_indices] = (
            self.reference_hashes[indices[valid_indices]] == sample_hashes[valid_indices]
        )
        
        return sample_hashes[~matches]
    
    def compute_dmi(
        self, 
        sample_hashes: Union[np.ndarray, List[int]],
        sample_name: Optional[str] = None
    ) -> DMIResult:
        """
        Compute the Dark Matter Index for a sample.
        
        Args:
            sample_hashes: Array or list of 64-bit hash values from the sample
            sample_name: Optional name/identifier for the sample
            
        Returns:
            DMIResult with DMI value and statistics
        """
        # Convert to numpy array if needed
        if isinstance(sample_hashes, list):
            sample_hashes = np.array(sample_hashes, dtype=np.uint64)
        elif sample_hashes.dtype != np.uint64:
            sample_hashes = sample_hashes.astype(np.uint64)
        
        total = len(sample_hashes)
        
        if total == 0:
            logger.warning(f"Sample {sample_name} has no hashes")
            return DMIResult(
                dmi=np.nan,
                total_hashes=0,
                unmapped_hashes=0,
                mapped_hashes=0,
                sample_name=sample_name
            )
        
        # Count mapped hashes
        mapped = self.count_mapped_hashes(sample_hashes)
        unmapped = total - mapped
        dmi = unmapped / total
        
        return DMIResult(
            dmi=dmi,
            total_hashes=total,
            unmapped_hashes=unmapped,
            mapped_hashes=mapped,
            sample_name=sample_name
        )
    
    def compute_dmi_from_signature(
        self, 
        sig_path: Union[str, Path],
        ksize: int = 31
    ) -> DMIResult:
        """
        Compute DMI from a sourmash signature file (.sig or .sig.gz).
        
        Args:
            sig_path: Path to the signature file
            ksize: k-mer size to use (default 31)
            
        Returns:
            DMIResult with DMI value and statistics
        """
        sig_path = Path(sig_path)
        
        # Load signature file (handle gzip if needed)
        if sig_path.suffix == '.gz' or str(sig_path).endswith('.sig.gz'):
            with gzip.open(sig_path, 'rt') as f:
                sig_data = json.load(f)
        else:
            with open(sig_path, 'r') as f:
                sig_data = json.load(f)
        
        # Handle both list format and single signature
        if isinstance(sig_data, list):
            sig_data = sig_data[0]
        
        # Extract sample name
        sample_name = sig_data.get('name', sig_path.stem)
        
        # Find the signature with matching ksize
        signatures = sig_data.get('signatures', [])
        target_sig = None
        
        for sig in signatures:
            if sig.get('ksize') == ksize:
                target_sig = sig
                break
        
        if target_sig is None:
            raise ValueError(f"No signature with ksize={ksize} found in {sig_path}")
        
        # Extract hash values
        hashes = np.array(target_sig['mins'], dtype=np.uint64)
        
        return self.compute_dmi(hashes, sample_name=sample_name)


def load_hashes_from_duckdb(
    db_path: str,
    sample_id: str,
    table: str = "sigs_dna.signature_mins",
    column: str = "min_hash"
) -> np.ndarray:
    """
    Load hash values for a sample from a DuckDB database.
    
    Args:
        db_path: Path to DuckDB database
        sample_id: Sample accession ID
        table: Table name containing the hashes
        column: Column name for the hash values
        
    Returns:
        numpy array of uint64 hash values
    """
    import duckdb
    
    conn = duckdb.connect(db_path, read_only=True)
    
    query = f"""
        SELECT {column}
        FROM {table}
        WHERE sample_id = ?
    """
    
    result = conn.execute(query, [sample_id]).fetchall()
    conn.close()
    
    if not result:
        logger.warning(f"No hashes found for sample {sample_id}")
        return np.array([], dtype=np.uint64)
    
    hashes = np.array([row[0] for row in result], dtype=np.uint64)
    return hashes


# Convenience function for command-line usage
def compute_dmi_for_sample(
    reference_path: str,
    sample_hashes: Union[str, np.ndarray, List[int]],
    ksize: int = 31
) -> DMIResult:
    """
    Convenience function to compute DMI for a single sample.
    
    Args:
        reference_path: Path to reference hash binary file
        sample_hashes: Either path to .sig file or array of hash values
        ksize: k-mer size (only used if sample_hashes is a path)
        
    Returns:
        DMIResult
    """
    calculator = DMICalculator(reference_path, use_mmap=True)
    
    if isinstance(sample_hashes, (str, Path)):
        return calculator.compute_dmi_from_signature(sample_hashes, ksize=ksize)
    else:
        return calculator.compute_dmi(sample_hashes)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Compute Dark Matter Index for a sample")
    parser.add_argument("reference", help="Path to reference hash binary file")
    parser.add_argument("sample", help="Path to sample signature file (.sig or .sig.gz)")
    parser.add_argument("--ksize", type=int, default=31, help="k-mer size (default: 31)")
    
    args = parser.parse_args()
    
    result = compute_dmi_for_sample(args.reference, args.sample, args.ksize)
    
    print(f"Sample: {result.sample_name}")
    print(f"DMI: {result.dmi:.4f} ({result.dmi*100:.2f}% dark matter)")
    print(f"Total hashes: {result.total_hashes:,}")
    print(f"Mapped to reference: {result.mapped_hashes:,} ({result.mapped_hashes/result.total_hashes*100:.2f}%)")
    print(f"Unmapped (dark matter): {result.unmapped_hashes:,}")
