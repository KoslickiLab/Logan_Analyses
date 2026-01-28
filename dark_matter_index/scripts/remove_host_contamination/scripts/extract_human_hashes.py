#!/usr/bin/env python3
"""
Extract Human Genome Hashes from Sourmash Signatures to DuckDB

This script reads sourmash signature files containing human genome sketches
and extracts all unique hashes into a DuckDB database for efficient lookup
during host contamination removal.

The output database can be used with compute_dmi_host_filtered.py to exclude
human-origin hashes from metagenome DMI calculations.

Usage:
    python extract_human_hashes.py \
        --input /path/to/human_genome_sketches.sig.zip \
        --output /path/to/human_hashes.db \
        --ksize 31

Prerequisites:
    pip install sourmash duckdb pandas

Author: David Koslicki lab
Date: 2026
"""

import argparse
import logging
import time
from pathlib import Path
from typing import Set, Iterator
import sys

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

try:
    import duckdb
    import sourmash
    from sourmash import load_file_as_signatures
except ImportError as e:
    logger.error(f"Missing dependency: {e}")
    logger.error("Install with: pip install duckdb sourmash")
    sys.exit(1)


def load_signatures_from_path(sig_path: str, ksize: int) -> Iterator[sourmash.SourmashSignature]:
    """
    Load signatures from a sourmash signature file.
    
    Handles .sig, .sig.gz, .sig.zip, and directory formats.
    Filters by ksize if specified.
    
    Args:
        sig_path: Path to signature file or directory
        ksize: k-mer size to filter for (None = all)
        
    Yields:
        sourmash.SourmashSignature objects
    """
    logger.info(f"Loading signatures from: {sig_path}")
    
    try:
        # Use sourmash's load_file_as_signatures which handles all formats
        sigs = load_file_as_signatures(sig_path, ksize=ksize)
        for sig in sigs:
            yield sig
    except Exception as e:
        logger.error(f"Error loading signatures: {e}")
        raise


def extract_hashes_from_signatures(
    sig_path: str,
    ksize: int,
    scaled: int = None
) -> tuple[Set[int], int]:
    """
    Extract all unique hashes from sourmash signatures.
    
    Args:
        sig_path: Path to signature file(s)
        ksize: k-mer size to use
        scaled: Expected scaled value (for validation)
        
    Returns:
        Tuple of (set of unique hashes, count of signatures processed)
    """
    all_hashes = set()
    sig_count = 0
    
    for sig in load_signatures_from_path(sig_path, ksize):
        # Get the minhash object
        mh = sig.minhash
        
        # Validate parameters
        if mh.ksize != ksize:
            logger.warning(f"Signature {sig.name} has ksize={mh.ksize}, expected {ksize}, skipping")
            continue
            
        if scaled is not None and mh.scaled != scaled:
            logger.warning(f"Signature {sig.name} has scaled={mh.scaled}, expected {scaled}")
        
        # Get hashes
        sig_hashes = mh.hashes.keys() if hasattr(mh.hashes, 'keys') else mh.hashes
        n_hashes = len(sig_hashes)
        
        # Add to set
        all_hashes.update(sig_hashes)
        sig_count += 1
        
        if sig_count % 50 == 0:
            logger.info(f"  Processed {sig_count} signatures, "
                       f"{len(all_hashes):,} unique hashes so far")
    
    return all_hashes, sig_count


def create_hash_database(
    hashes: Set[int],
    output_path: str,
    source_path: str,
    ksize: int,
    n_signatures: int,
    chunk_size: int = 10_000_000,
    threads: int = 64,
    memory: str = "128GB"
):
    """
    Create a DuckDB database containing the hash set.
    
    The database contains:
    - human_hashes table: All unique hashes
    - metadata table: Information about the source
    
    Args:
        hashes: Set of hash values
        output_path: Path for output database
        source_path: Original signature file path (for metadata)
        ksize: k-mer size
        n_signatures: Number of signatures processed
        chunk_size: Number of hashes to insert per batch
        threads: DuckDB threads
        memory: DuckDB memory limit
    """
    logger.info(f"Creating database at: {output_path}")
    logger.info(f"Total unique hashes: {len(hashes):,}")
    
    # Create database
    conn = duckdb.connect(output_path)
    conn.execute(f"SET threads TO {threads}")
    conn.execute(f"SET memory_limit = '{memory}'")
    
    try:
        # Create tables
        conn.execute("""
            CREATE TABLE IF NOT EXISTS human_hashes (
                hash UBIGINT PRIMARY KEY
            )
        """)
        
        conn.execute("""
            CREATE TABLE IF NOT EXISTS metadata (
                key VARCHAR PRIMARY KEY,
                value VARCHAR
            )
        """)
        
        # Insert metadata
        import json
        from datetime import datetime
        
        metadata = {
            'source_path': source_path,
            'ksize': str(ksize),
            'n_signatures': str(n_signatures),
            'n_unique_hashes': str(len(hashes)),
            'created_at': datetime.now().isoformat(),
            'sourmash_version': sourmash.__version__
        }
        
        for key, value in metadata.items():
            conn.execute(
                "INSERT INTO metadata (key, value) VALUES (?, ?)",
                [key, value]
            )
        
        # Convert set to sorted list for deterministic ordering
        logger.info("Sorting hashes...")
        hash_list = sorted(hashes)
        
        # Insert hashes in chunks
        logger.info(f"Inserting hashes in chunks of {chunk_size:,}...")
        start_time = time.time()
        
        import pandas as pd
        
        for i in range(0, len(hash_list), chunk_size):
            chunk_end = min(i + chunk_size, len(hash_list))
            chunk = hash_list[i:chunk_end]
            
            # Create DataFrame and insert
            chunk_df = pd.DataFrame({'hash': chunk})
            conn.register('chunk_data', chunk_df)
            
            conn.execute("""
                INSERT INTO human_hashes (hash)
                SELECT hash FROM chunk_data
            """)
            
            conn.unregister('chunk_data')
            
            elapsed = time.time() - start_time
            rate = chunk_end / elapsed
            eta = (len(hash_list) - chunk_end) / rate / 60 if rate > 0 else 0
            
            logger.info(f"  Inserted {chunk_end:,}/{len(hash_list):,} hashes "
                       f"({100*chunk_end/len(hash_list):.1f}%), ETA: {eta:.1f} min")
        
        # Create index for fast lookups
        logger.info("Creating index on hash column...")
        idx_start = time.time()
        conn.execute("CREATE INDEX IF NOT EXISTS idx_human_hash ON human_hashes(hash)")
        logger.info(f"  Index created in {time.time() - idx_start:.1f}s")
        
        # Run ANALYZE
        logger.info("Analyzing table...")
        conn.execute("ANALYZE human_hashes")
        
        # Print summary
        count = conn.execute("SELECT COUNT(*) FROM human_hashes").fetchone()[0]
        logger.info(f"\nDatabase created successfully:")
        logger.info(f"  Total hashes: {count:,}")
        logger.info(f"  File: {output_path}")
        
        # Get file size
        db_size = Path(output_path).stat().st_size
        logger.info(f"  Size: {db_size/1e9:.2f} GB")
        
    finally:
        conn.close()


def verify_database(db_path: str):
    """
    Verify the created database.
    """
    logger.info(f"\nVerifying database: {db_path}")
    
    conn = duckdb.connect(db_path, read_only=True)
    
    try:
        # Check table exists
        tables = conn.execute("""
            SELECT table_name 
            FROM information_schema.tables 
            WHERE table_schema = 'main'
        """).fetchall()
        table_names = [t[0] for t in tables]
        
        assert 'human_hashes' in table_names, "Missing human_hashes table"
        assert 'metadata' in table_names, "Missing metadata table"
        
        # Get counts
        count = conn.execute("SELECT COUNT(*) FROM human_hashes").fetchone()[0]
        
        # Get metadata
        metadata = dict(conn.execute("SELECT key, value FROM metadata").fetchall())
        
        logger.info("  Verification passed!")
        logger.info(f"  Hashes in database: {count:,}")
        logger.info(f"  Source: {metadata.get('source_path', 'unknown')}")
        logger.info(f"  k-mer size: {metadata.get('ksize', 'unknown')}")
        logger.info(f"  Signatures processed: {metadata.get('n_signatures', 'unknown')}")
        
        # Sample some hashes
        sample = conn.execute("SELECT hash FROM human_hashes LIMIT 5").fetchall()
        logger.info(f"  Sample hashes: {[h[0] for h in sample]}")
        
    finally:
        conn.close()


def main():
    parser = argparse.ArgumentParser(
        description="Extract human genome hashes from sourmash signatures to DuckDB",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script extracts FracMinHash signatures from human genome assemblies
and stores them in a DuckDB database for efficient lookup during
host contamination removal.

The output database can be used with compute_dmi_host_filtered.py.

Examples:
    # Basic usage
    python extract_human_hashes.py \\
        --input human_genome_sketches_all.sig.zip \\
        --output human_hashes.db \\
        --ksize 31

    # With explicit scaled validation
    python extract_human_hashes.py \\
        --input human_genome_sketches_all.sig.zip \\
        --output human_hashes.db \\
        --ksize 31 \\
        --scaled 1000
        """
    )
    
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to sourmash signature file(s) (.sig, .sig.gz, .sig.zip)"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Path for output DuckDB database"
    )
    parser.add_argument(
        "--ksize", "-k",
        type=int,
        required=True,
        help="k-mer size to extract (e.g., 31)"
    )
    parser.add_argument(
        "--scaled",
        type=int,
        default=None,
        help="Expected scaled value for validation (default: no validation)"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=64,
        help="DuckDB threads (default: 64)"
    )
    parser.add_argument(
        "--memory",
        type=str,
        default="128GB",
        help="DuckDB memory limit (default: 128GB)"
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10_000_000,
        help="Hashes to insert per batch (default: 10M)"
    )
    
    args = parser.parse_args()
    
    total_start = time.time()
    
    logger.info("="*60)
    logger.info("EXTRACT HUMAN GENOME HASHES")
    logger.info("="*60)
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    logger.info(f"k-mer size: {args.ksize}")
    logger.info(f"Scaled: {args.scaled or 'any'}")
    
    # Check input exists
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    # Check output doesn't exist
    output_path = Path(args.output)
    if output_path.exists():
        logger.warning(f"Output file exists and will be overwritten: {args.output}")
        output_path.unlink()
    
    # Step 1: Extract hashes from signatures
    logger.info("\n" + "="*60)
    logger.info("STEP 1: Extracting hashes from signatures")
    logger.info("="*60)
    
    extract_start = time.time()
    hashes, n_sigs = extract_hashes_from_signatures(
        args.input,
        args.ksize,
        args.scaled
    )
    extract_time = time.time() - extract_start
    
    logger.info(f"\nExtraction complete:")
    logger.info(f"  Signatures processed: {n_sigs}")
    logger.info(f"  Unique hashes: {len(hashes):,}")
    logger.info(f"  Time: {extract_time/60:.1f} minutes")
    
    if len(hashes) == 0:
        logger.error("No hashes extracted! Check input file and ksize.")
        sys.exit(1)
    
    # Step 2: Create database
    logger.info("\n" + "="*60)
    logger.info("STEP 2: Creating DuckDB database")
    logger.info("="*60)
    
    db_start = time.time()
    create_hash_database(
        hashes,
        args.output,
        args.input,
        args.ksize,
        n_sigs,
        chunk_size=args.chunk_size,
        threads=args.threads,
        memory=args.memory
    )
    db_time = time.time() - db_start
    
    logger.info(f"\nDatabase creation complete:")
    logger.info(f"  Time: {db_time/60:.1f} minutes")
    
    # Step 3: Verify
    verify_database(args.output)
    
    # Summary
    total_time = time.time() - total_start
    
    logger.info("\n" + "="*60)
    logger.info("COMPLETE")
    logger.info("="*60)
    logger.info(f"Total time: {total_time/60:.1f} minutes")
    logger.info(f"Output: {args.output}")
    logger.info(f"\nNext step: Use with compute_dmi_host_filtered.py")


if __name__ == "__main__":
    main()
