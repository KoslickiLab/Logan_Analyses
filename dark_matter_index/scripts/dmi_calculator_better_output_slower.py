#!/usr/bin/env python3
"""
compute_dmi.py — Compute Dark Matter Index for a metagenome sample

This is the user-facing script for computing the Dark Matter Index (DMI) from 
a sourmash signature file (.sig or .sig.gz). The DMI measures the fraction of 
k-mer hashes in your sample that do NOT appear in any reference database.

DMI = |sample_hashes - reference_hashes| / |sample_hashes|

Interpretation:
    DMI = 0.0  → All hashes are explained by known references
    DMI = 0.5  → Half of your sample is uncharacterized "dark matter"
    DMI = 1.0  → Nothing matches any reference database

Usage:
    # Basic usage
    python compute_dmi.py --reference reference_hashes_k31.bin --input sample.sig --ksize 31
    
    # With gzipped signature
    python compute_dmi.py -r reference_hashes_k31.bin -i sample.sig.gz -k 31
    
    # Output to JSON file
    python compute_dmi.py -r reference_hashes_k31.bin -i sample.sig -k 31 -o result.json
    
    # Process multiple files
    python compute_dmi.py -r reference_hashes_k31.bin -i *.sig -k 31 -o results.json
    
    # Quiet mode (just print DMI value)
    python compute_dmi.py -r reference_hashes_k31.bin -i sample.sig -k 31 --quiet

Author: David Koslicki Lab, Penn State University
"""

import numpy as np
import json
import gzip
import sys
import argparse
from pathlib import Path
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, asdict
import time


@dataclass
class DMIResult:
    """Container for DMI computation results."""
    sample_name: str
    dmi: float
    dark_matter_percent: float
    total_hashes: int
    unmapped_hashes: int
    mapped_hashes: int
    input_file: str
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class DMICalculator:
    """
    Calculator for Dark Matter Index using a sorted reference hash set.
    
    The reference hashes are stored as a sorted numpy array of uint64 values,
    enabling efficient O(log n) membership testing via binary search.
    """
    
    def __init__(self, reference_path: str, use_mmap: bool = True, verbose: bool = True):
        """
        Initialize the DMI calculator.
        
        Args:
            reference_path: Path to the sorted binary reference hash file
            use_mmap: If True, memory-map the file (recommended for large files)
            verbose: If True, print loading progress
        """
        self.reference_path = Path(reference_path)
        self.verbose = verbose
        
        if not self.reference_path.exists():
            raise FileNotFoundError(f"Reference file not found: {reference_path}")
        
        # Calculate number of hashes from file size
        file_size = self.reference_path.stat().st_size
        if file_size % 8 != 0:
            raise ValueError(f"Invalid reference file: size {file_size} not divisible by 8")
        
        self.n_reference = file_size // 8
        
        if verbose:
            print(f"Loading reference: {self.n_reference:,} hashes ({file_size/1e9:.2f} GB)", 
                  file=sys.stderr)
        
        # Memory-map the reference file
        if use_mmap:
            self.reference = np.memmap(
                self.reference_path,
                dtype=np.uint64,
                mode='r',
                shape=(self.n_reference,)
            )
        else:
            self.reference = np.fromfile(self.reference_path, dtype=np.uint64)
        
        if verbose:
            print(f"Reference loaded successfully", file=sys.stderr)
    
    def compute_dmi(self, sample_hashes: np.ndarray) -> tuple:
        """
        Compute DMI for a set of sample hashes.
        
        Args:
            sample_hashes: numpy array of uint64 hash values
            
        Returns:
            Tuple of (dmi, mapped_count, unmapped_count)
        """
        if len(sample_hashes) == 0:
            return np.nan, 0, 0
        
        # Ensure proper dtype
        if sample_hashes.dtype != np.uint64:
            sample_hashes = sample_hashes.astype(np.uint64)
        
        # Binary search for each hash in the sorted reference
        indices = np.searchsorted(self.reference, sample_hashes)
        
        # Check which insertions land on actual matches
        valid = indices < self.n_reference
        matches = np.zeros(len(sample_hashes), dtype=bool)
        matches[valid] = self.reference[indices[valid]] == sample_hashes[valid]
        
        mapped = int(np.sum(matches))
        unmapped = len(sample_hashes) - mapped
        dmi = unmapped / len(sample_hashes)
        
        return dmi, mapped, unmapped
    
    def compute_dmi_from_signature(self, sig_path: str, ksize: int = 31) -> DMIResult:
        """
        Compute DMI from a sourmash signature file.
        
        Args:
            sig_path: Path to .sig or .sig.gz file
            ksize: k-mer size to use (default: 31)
            
        Returns:
            DMIResult object with all statistics
        """
        sig_path = Path(sig_path)
        
        # Load signature (handle gzip)
        if sig_path.suffix == '.gz' or str(sig_path).endswith('.sig.gz'):
            with gzip.open(sig_path, 'rt') as f:
                sig_data = json.load(f)
        else:
            with open(sig_path, 'r') as f:
                sig_data = json.load(f)
        
        # Handle list format
        if isinstance(sig_data, list):
            sig_data = sig_data[0]
        
        # Extract sample name
        sample_name = sig_data.get('name', sig_path.stem)
        
        # Find signature with matching ksize
        signatures = sig_data.get('signatures', [])
        target_sig = None
        
        for sig in signatures:
            if sig.get('ksize') == ksize:
                target_sig = sig
                break
        
        if target_sig is None:
            available_ksizes = [s.get('ksize') for s in signatures]
            raise ValueError(
                f"No signature with ksize={ksize} found in {sig_path}. "
                f"Available ksizes: {available_ksizes}"
            )
        
        # Extract hashes
        hashes = np.array(target_sig['mins'], dtype=np.uint64)
        
        # Compute DMI
        dmi, mapped, unmapped = self.compute_dmi(hashes)
        
        return DMIResult(
            sample_name=sample_name,
            dmi=dmi,
            dark_matter_percent=dmi * 100,
            total_hashes=len(hashes),
            unmapped_hashes=unmapped,
            mapped_hashes=mapped,
            input_file=str(sig_path)
        )


def print_result(result: DMIResult, quiet: bool = False):
    """Print DMI result to stdout."""
    if quiet:
        print(f"{result.dmi:.6f}")
    else:
        print(f"\n{'='*60}")
        print(f"Sample: {result.sample_name}")
        print(f"{'='*60}")
        print(f"Dark Matter Index (DMI): {result.dmi:.4f}")
        print(f"Dark Matter:             {result.dark_matter_percent:.2f}%")
        print(f"")
        print(f"Total hashes:            {result.total_hashes:,}")
        print(f"Mapped to reference:     {result.mapped_hashes:,} ({100*result.mapped_hashes/result.total_hashes:.2f}%)")
        print(f"Unmapped (dark matter):  {result.unmapped_hashes:,} ({result.dark_matter_percent:.2f}%)")
        print(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Compute Dark Matter Index for metagenome samples",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python compute_dmi.py -r reference_hashes_k31.bin -i sample.sig -k 31
    
    # Process multiple files
    python compute_dmi.py -r reference_hashes_k31.bin -i *.sig -k 31 -o results.json
    
    # Quiet mode (just DMI values)
    python compute_dmi.py -r reference_hashes_k31.bin -i sample.sig -k 31 --quiet

Interpretation:
    DMI = 0.0  → All k-mers are in reference databases (well-characterized sample)
    DMI = 0.3  → 30% of k-mers are "dark matter" (moderately novel)
    DMI = 0.7  → 70% of k-mers are uncharacterized (highly novel sample!)
        """
    )
    
    parser.add_argument(
        "--reference", "-r",
        required=True,
        help="Path to reference hash binary file (reference_hashes_k31.bin)"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        nargs='+',
        help="Input signature file(s) (.sig or .sig.gz)"
    )
    parser.add_argument(
        "--output", "-o",
        help="Output JSON file for results (optional)"
    )
    parser.add_argument(
        "--ksize", "-k",
        type=int,
        required=True,
        help="k-mer size (REQUIRED, e.g., 31)"
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Quiet mode: only print DMI values"
    )
    
    args = parser.parse_args()
    
    # Initialize calculator
    calculator = DMICalculator(
        args.reference, 
        use_mmap=True, 
        verbose=not args.quiet
    )
    
    # Process each input file
    results = []
    
    for sig_path in args.input:
        try:
            if not args.quiet:
                print(f"\nProcessing: {sig_path}", file=sys.stderr)
            
            start = time.time()
            result = calculator.compute_dmi_from_signature(sig_path, ksize=args.ksize)
            elapsed = time.time() - start
            
            results.append(result)
            print_result(result, quiet=args.quiet)
            
            if not args.quiet:
                print(f"Time: {elapsed*1000:.1f} ms", file=sys.stderr)
                
        except Exception as e:
            print(f"Error processing {sig_path}: {e}", file=sys.stderr)
            if not args.quiet:
                raise
    
    # Save results to JSON if requested
    if args.output and results:
        output_data = {
            'results': [r.to_dict() for r in results],
            'summary': {
                'n_samples': len(results),
                'mean_dmi': np.mean([r.dmi for r in results]),
                'median_dmi': np.median([r.dmi for r in results]),
                'min_dmi': min(r.dmi for r in results),
                'max_dmi': max(r.dmi for r in results),
            }
        }
        
        with open(args.output, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        if not args.quiet:
            print(f"\nResults saved to: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
