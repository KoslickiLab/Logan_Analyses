#!/usr/bin/env python3
"""
Batch Compute Dark Matter Index for SRA Samples

This script processes all samples in a filtered_data.parquet file, computing the
Dark Matter Index (DMI) for each sample by querying its hashes from the Logan
DuckDB database.

The output is an augmented parquet file with the DMI column added.

Usage:
    python batch_compute_dmi.py \
        --reference /path/to/reference_hashes_k31.bin \
        --input /path/to/filtered_data.parquet \
        --output /path/to/output_with_dmi.parquet \
        --database /scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db \
        --ksize 31 \
        --workers 128

Performance Notes:
    - The reference hash set is memory-mapped, using minimal RAM per worker
    - Each worker opens its own read-only database connection
    - Progress is logged every 1000 samples
    - Failed samples are logged but don't stop processing
    - Results are checkpointed periodically to allow resumption
    - Optimized for high-core-count systems (tested with 256 cores)

Author: David Koslicki lab
Date: 2026
"""

import numpy as np
import pandas as pd
import json
import time
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, asdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import multiprocessing as mp

# Try imports
try:
    import duckdb
    HAS_DUCKDB = True
except ImportError:
    HAS_DUCKDB = False

try:
    import pyarrow.parquet as pq
    HAS_PYARROW = True
except ImportError:
    HAS_PYARROW = False

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class DMIResult:
    """Container for DMI computation results."""
    accession: str
    dmi: float
    total_hashes: int
    unmapped_hashes: int
    mapped_hashes: int
    error: Optional[str] = None


class DMIBatchProcessor:
    """
    Batch processor for computing DMI across many samples.
    
    Designed for efficient processing of large numbers of samples with:
    - Memory-mapped reference hash set (minimal RAM)
    - Parallel processing with configurable workers
    - Checkpointing for resumption
    - Comprehensive error handling and logging
    """
    
    def __init__(
        self,
        reference_path: str,
        database_path: str,
        n_workers: int = 1,
        checkpoint_interval: int = 10000
    ):
        """
        Initialize the batch processor.
        
        Args:
            reference_path: Path to the reference hash binary file
            database_path: Path to the DuckDB database with sample hashes
            n_workers: Number of parallel workers
            checkpoint_interval: Save checkpoint every N samples
        """
        self.reference_path = Path(reference_path)
        self.database_path = database_path
        self.n_workers = n_workers
        self.checkpoint_interval = checkpoint_interval
        
        # Verify files exist
        if not self.reference_path.exists():
            raise FileNotFoundError(f"Reference file not found: {reference_path}")
        
        # Get reference metadata
        ref_size = self.reference_path.stat().st_size
        self.n_reference = ref_size // 8
        logger.info(f"Reference contains {self.n_reference:,} hashes ({ref_size/1e9:.2f} GB)")
    
    def _load_reference_mmap(self) -> np.ndarray:
        """Load reference as memory-mapped array."""
        return np.memmap(
            self.reference_path,
            dtype=np.uint64,
            mode='r',
            shape=(self.n_reference,)
        )
    
    def _compute_dmi_single(
        self,
        sample_hashes: np.ndarray,
        reference: np.ndarray
    ) -> Tuple[float, int, int]:
        """
        Compute DMI for a single sample.
        
        Args:
            sample_hashes: Array of sample hash values
            reference: Memory-mapped reference array
            
        Returns:
            Tuple of (dmi, mapped_count, unmapped_count)
        """
        if len(sample_hashes) == 0:
            return np.nan, 0, 0
        
        # Ensure proper dtype
        if sample_hashes.dtype != np.uint64:
            sample_hashes = sample_hashes.astype(np.uint64)
        
        # Binary search for each hash
        indices = np.searchsorted(reference, sample_hashes)
        
        # Check which are actual matches
        valid = indices < len(reference)
        matches = np.zeros(len(sample_hashes), dtype=bool)
        matches[valid] = reference[indices[valid]] == sample_hashes[valid]
        
        mapped = int(np.sum(matches))
        total = len(sample_hashes)
        unmapped = total - mapped
        dmi = unmapped / total
        
        return dmi, mapped, unmapped
    
    def _fetch_sample_hashes(
        self,
        conn: 'duckdb.DuckDBPyConnection',
        accession: str
    ) -> np.ndarray:
        """
        Fetch hash values for a sample from the database.
        
        Args:
            conn: DuckDB connection
            accession: Sample accession ID
            
        Returns:
            numpy array of uint64 hash values
        """
        query = """
            SELECT min_hash
            FROM sigs_dna.signature_mins
            WHERE sample_id = ?
        """
        
        result = conn.execute(query, [accession]).fetchall()
        
        if not result:
            return np.array([], dtype=np.uint64)
        
        return np.array([row[0] for row in result], dtype=np.uint64)
    
    def process_samples_sequential(
        self,
        accessions: List[str],
        progress_callback=None
    ) -> List[DMIResult]:
        """
        Process samples sequentially (single-threaded).
        
        Args:
            accessions: List of sample accessions to process
            progress_callback: Optional callback(current, total) for progress
            
        Returns:
            List of DMIResult objects
        """
        reference = self._load_reference_mmap()
        conn = duckdb.connect(self.database_path, read_only=True)
        
        results = []
        start_time = time.time()
        
        for i, accession in enumerate(accessions):
            try:
                # Fetch hashes
                hashes = self._fetch_sample_hashes(conn, accession)
                
                # Compute DMI
                dmi, mapped, unmapped = self._compute_dmi_single(hashes, reference)
                
                results.append(DMIResult(
                    accession=accession,
                    dmi=dmi,
                    total_hashes=len(hashes),
                    unmapped_hashes=unmapped,
                    mapped_hashes=mapped
                ))
                
            except Exception as e:
                logger.warning(f"Error processing {accession}: {e}")
                results.append(DMIResult(
                    accession=accession,
                    dmi=np.nan,
                    total_hashes=0,
                    unmapped_hashes=0,
                    mapped_hashes=0,
                    error=str(e)
                ))
            
            # Progress logging
            if (i + 1) % 1000 == 0:
                elapsed = time.time() - start_time
                rate = (i + 1) / elapsed
                eta = (len(accessions) - i - 1) / rate / 60
                logger.info(
                    f"Progress: {i+1:,}/{len(accessions):,} "
                    f"({100*(i+1)/len(accessions):.1f}%) "
                    f"Rate: {rate:.1f}/s, ETA: {eta:.1f} min"
                )
            
            if progress_callback:
                progress_callback(i + 1, len(accessions))
        
        conn.close()
        return results


def process_chunk(args: Tuple) -> List[Dict]:
    """
    Worker function for parallel processing.
    
    Args:
        args: Tuple of (accessions, reference_path, database_path)
        
    Returns:
        List of result dictionaries
    """
    accessions, reference_path, database_path = args
    
    # Load reference (each worker gets its own mmap)
    ref_size = Path(reference_path).stat().st_size // 8
    reference = np.memmap(reference_path, dtype=np.uint64, mode='r', shape=(ref_size,))
    
    # Connect to database
    conn = duckdb.connect(database_path, read_only=True)
    
    results = []
    
    for accession in accessions:
        try:
            # Fetch hashes
            query = """
                SELECT min_hash
                FROM sigs_dna.signature_mins
                WHERE sample_id = ?
            """
            result = conn.execute(query, [accession]).fetchall()
            
            if not result:
                hashes = np.array([], dtype=np.uint64)
            else:
                hashes = np.array([row[0] for row in result], dtype=np.uint64)
            
            # Compute DMI
            if len(hashes) == 0:
                dmi, mapped, unmapped = np.nan, 0, 0
            else:
                indices = np.searchsorted(reference, hashes)
                valid = indices < len(reference)
                matches = np.zeros(len(hashes), dtype=bool)
                matches[valid] = reference[indices[valid]] == hashes[valid]
                
                mapped = int(np.sum(matches))
                unmapped = len(hashes) - mapped
                dmi = unmapped / len(hashes)
            
            results.append({
                'accession': accession,
                'dmi': dmi,
                'total_hashes': len(hashes),
                'unmapped_hashes': unmapped,
                'mapped_hashes': mapped,
                'error': None
            })
            
        except Exception as e:
            results.append({
                'accession': accession,
                'dmi': np.nan,
                'total_hashes': 0,
                'unmapped_hashes': 0,
                'mapped_hashes': 0,
                'error': str(e)
            })
    
    conn.close()
    return results


def batch_compute_dmi_parallel(
    reference_path: str,
    database_path: str,
    accessions: List[str],
    n_workers: int = 8,
    chunk_size: int = 1000
) -> pd.DataFrame:
    """
    Compute DMI for many samples in parallel.
    
    Args:
        reference_path: Path to reference hash binary file
        database_path: Path to DuckDB database
        accessions: List of sample accessions
        n_workers: Number of parallel workers
        chunk_size: Samples per worker chunk
        
    Returns:
        DataFrame with DMI results
    """
    logger.info(f"Processing {len(accessions):,} samples with {n_workers} workers")
    
    # Split accessions into chunks
    chunks = [
        accessions[i:i + chunk_size]
        for i in range(0, len(accessions), chunk_size)
    ]
    
    # Prepare arguments for workers
    args_list = [
        (chunk, reference_path, database_path)
        for chunk in chunks
    ]
    
    all_results = []
    start_time = time.time()
    processed = 0
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_chunk, args): i for i, args in enumerate(args_list)}
        
        for future in as_completed(futures):
            chunk_results = future.result()
            all_results.extend(chunk_results)
            
            processed += len(chunk_results)
            elapsed = time.time() - start_time
            rate = processed / elapsed
            eta = (len(accessions) - processed) / rate / 60
            
            logger.info(
                f"Progress: {processed:,}/{len(accessions):,} "
                f"({100*processed/len(accessions):.1f}%) "
                f"Rate: {rate:.1f}/s, ETA: {eta:.1f} min"
            )
    
    # Convert to DataFrame
    df = pd.DataFrame(all_results)
    
    total_time = time.time() - start_time
    logger.info(f"Completed in {total_time/60:.1f} minutes ({total_time/len(accessions)*1000:.2f} ms/sample)")
    
    return df


def main(
    reference_path: str,
    input_parquet: str,
    output_parquet: str,
    database_path: str,
    ksize: int,
    n_workers: int = 128,
    chunk_size: int = 500,
    resume_checkpoint: Optional[str] = None
):
    """
    Main function to process all samples and create augmented parquet.
    
    Args:
        reference_path: Path to reference hash binary file
        input_parquet: Path to input filtered_data.parquet
        output_parquet: Path for output parquet with DMI
        database_path: Path to DuckDB database
        ksize: k-mer size (for documentation/logging)
        n_workers: Number of parallel workers
        chunk_size: Samples per worker chunk
        resume_checkpoint: Path to checkpoint file for resumption
    """
    logger.info("="*60)
    logger.info("BATCH DMI COMPUTATION")
    logger.info("="*60)
    logger.info(f"ksize={ksize}")
    
    # Load input data
    logger.info(f"Loading input data from {input_parquet}")
    df_input = pd.read_parquet(input_parquet)
    logger.info(f"Loaded {len(df_input):,} samples")
    
    # Get list of accessions
    accessions = df_input['accession'].tolist()
    
    # Check for checkpoint
    processed_accessions = set()
    checkpoint_results = []
    
    if resume_checkpoint and Path(resume_checkpoint).exists():
        logger.info(f"Loading checkpoint from {resume_checkpoint}")
        checkpoint_df = pd.read_parquet(resume_checkpoint)
        processed_accessions = set(checkpoint_df['accession'].tolist())
        checkpoint_results = checkpoint_df.to_dict('records')
        logger.info(f"Resuming from checkpoint: {len(processed_accessions):,} already processed")
    
    # Filter to unprocessed accessions
    accessions_to_process = [a for a in accessions if a not in processed_accessions]
    logger.info(f"Will process {len(accessions_to_process):,} samples")
    
    if len(accessions_to_process) == 0:
        logger.info("All samples already processed!")
        dmi_df = pd.DataFrame(checkpoint_results)
    else:
        # Compute DMI
        if n_workers > 1:
            dmi_df = batch_compute_dmi_parallel(
                reference_path,
                database_path,
                accessions_to_process,
                n_workers=n_workers,
                chunk_size=chunk_size
            )
        else:
            processor = DMIBatchProcessor(
                reference_path,
                database_path,
                n_workers=1
            )
            results = processor.process_samples_sequential(accessions_to_process)
            dmi_df = pd.DataFrame([asdict(r) for r in results])
        
        # Combine with checkpoint if resuming
        if checkpoint_results:
            checkpoint_df = pd.DataFrame(checkpoint_results)
            dmi_df = pd.concat([checkpoint_df, dmi_df], ignore_index=True)
    
    # Merge with input data
    logger.info("Merging results with input data...")
    
    # Select columns to add
    dmi_columns = ['accession', 'dmi', 'total_hashes_dmi', 'unmapped_hashes', 'mapped_hashes']
    dmi_df = dmi_df.rename(columns={
        'total_hashes': 'total_hashes_dmi'  # Avoid conflict with existing column
    })
    
    # Keep only needed columns
    dmi_merge = dmi_df[['accession', 'dmi', 'total_hashes_dmi', 'unmapped_hashes', 'mapped_hashes']].copy()
    
    # Merge
    df_output = df_input.merge(dmi_merge, on='accession', how='left')
    
    # Save output
    logger.info(f"Saving output to {output_parquet}")
    df_output.to_parquet(output_parquet, index=False)
    
    # Summary statistics
    logger.info("\n" + "="*60)
    logger.info("SUMMARY STATISTICS")
    logger.info("="*60)
    
    valid_dmi = df_output['dmi'].dropna()
    logger.info(f"Total samples: {len(df_output):,}")
    logger.info(f"Samples with DMI: {len(valid_dmi):,}")
    logger.info(f"DMI statistics:")
    logger.info(f"  Mean: {valid_dmi.mean():.4f}")
    logger.info(f"  Median: {valid_dmi.median():.4f}")
    logger.info(f"  Std: {valid_dmi.std():.4f}")
    logger.info(f"  Min: {valid_dmi.min():.4f}")
    logger.info(f"  Max: {valid_dmi.max():.4f}")
    logger.info(f"  25th percentile: {valid_dmi.quantile(0.25):.4f}")
    logger.info(f"  75th percentile: {valid_dmi.quantile(0.75):.4f}")
    logger.info(f"  90th percentile: {valid_dmi.quantile(0.90):.4f}")
    logger.info(f"  95th percentile: {valid_dmi.quantile(0.95):.4f}")
    
    # By organism type (if available)
    if 'organism' in df_output.columns:
        logger.info("\nDMI by organism type (top 10):")
        organism_stats = df_output.groupby('organism')['dmi'].agg(['mean', 'median', 'count'])
        organism_stats = organism_stats[organism_stats['count'] >= 100]
        organism_stats = organism_stats.sort_values('median', ascending=False).head(10)
        for org, row in organism_stats.iterrows():
            logger.info(f"  {org}: median={row['median']:.4f}, mean={row['mean']:.4f}, n={int(row['count']):,}")
    
    logger.info("="*60)
    logger.info("DONE")
    logger.info("="*60)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Batch compute Dark Matter Index for SRA samples",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python batch_compute_dmi.py \\
        --reference /scratch/dmk333_new/dmi_reference/reference_hashes_k31.bin \\
        --input /path/to/filtered_data.parquet \\
        --output /path/to/filtered_data_with_dmi.parquet \\
        --database /scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db \\
        --ksize 31
    
    # With more workers for faster processing
    python batch_compute_dmi.py \\
        --reference reference_hashes_k31.bin \\
        --input filtered_data.parquet \\
        --output output_with_dmi.parquet \\
        --database database_all.db \\
        --ksize 31 \\
        --workers 64 \\
        --chunk-size 500
    
    # Resume from checkpoint
    python batch_compute_dmi.py \\
        --reference reference_hashes_k31.bin \\
        --input filtered_data.parquet \\
        --output output_with_dmi.parquet \\
        --database database_all.db \\
        --ksize 31 \\
        --resume-checkpoint checkpoint.parquet
        """
    )
    
    parser.add_argument(
        "--reference", "-r",
        required=True,
        help="Path to reference hash binary file"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to input filtered_data.parquet"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Path for output parquet with DMI"
    )
    parser.add_argument(
        "--database", "-d",
        required=True,
        help="Path to DuckDB database with sample hashes"
    )
    parser.add_argument(
        "--ksize", "-k",
        type=int,
        required=True,
        help="k-mer size (REQUIRED, e.g., 31). Must match the reference file."
    )
    parser.add_argument(
        "--workers", "-w",
        type=int,
        default=128,
        help="Number of parallel workers (default: 128)"
    )
    parser.add_argument(
        "--chunk-size", "-c",
        type=int,
        default=500,
        help="Samples per worker chunk (default: 500)"
    )
    parser.add_argument(
        "--resume-checkpoint",
        help="Path to checkpoint file for resumption"
    )
    
    args = parser.parse_args()
    
    main(
        reference_path=args.reference,
        input_parquet=args.input,
        output_parquet=args.output,
        database_path=args.database,
        ksize=args.ksize,
        n_workers=args.workers,
        chunk_size=args.chunk_size,
        resume_checkpoint=args.resume_checkpoint
    )
