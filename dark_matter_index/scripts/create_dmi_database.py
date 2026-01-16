#!/usr/bin/env python3
"""
Create DMI Database - Prepare a specialized DuckDB database for DMI computation

This script creates a new, compact DuckDB database containing:
1. Sample hashes: Only the accessions you care about, at the specified ksize
2. Reference hashes: The entire reference set, loaded and indexed

This approach leverages DuckDB's native optimizations for joins and aggregations,
making DMI computation orders of magnitude faster than Python-based approaches.

Usage:
    python create_dmi_database.py \\
        --source /path/to/database_all.db \\
        --samples /path/to/filtered_data.parquet \\
        --reference /path/to/reference_hashes_k31.bin \\
        --output /path/to/dmi_database.db \\
        --ksize 31

Time estimates:
    - Extracting sample hashes: 10-30 min (depends on number of samples)
    - Loading reference hashes: 20-60 min (for 10B hashes)
    - Creating indexes: 10-30 min
    - Total: ~1-2 hours

Author: David Koslicki lab
Date: 2026
"""

import argparse
import logging
import time
import numpy as np
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

try:
    import duckdb
    import pandas as pd
except ImportError as e:
    logger.error(f"Missing dependency: {e}")
    logger.error("Install with: pip install duckdb pandas pyarrow")
    exit(1)


def extract_sample_hashes(
    source_db: str,
    samples_parquet: str,
    output_db: duckdb.DuckDBPyConnection,
    ksize: int,
    batch_size: int = 100000
):
    """
    Extract sample hashes from source database for specified accessions.
    
    Uses a chunked approach to avoid memory issues with large sample sets.
    """
    logger.info("="*60)
    logger.info("STEP 1: Extracting sample hashes")
    logger.info("="*60)
    
    # Load accessions from parquet
    logger.info(f"Loading accessions from {samples_parquet}")
    df_samples = pd.read_parquet(samples_parquet)
    accessions = df_samples['accession'].tolist()
    logger.info(f"Found {len(accessions):,} accessions to extract")
    
    # Create sample_hashes table
    logger.info("Creating sample_hashes table...")
    output_db.execute("""
        CREATE TABLE IF NOT EXISTS sample_hashes (
            sample_id VARCHAR,
            min_hash UBIGINT
        )
    """)
    
    # Attach source database
    logger.info(f"Attaching source database: {source_db}")
    output_db.execute(f"ATTACH '{source_db}' AS source (READ_ONLY)")
    
    # Process in batches
    total_hashes = 0
    start_time = time.time()
    n_batches = (len(accessions) + batch_size - 1) // batch_size
    
    for batch_idx in range(n_batches):
        batch_start = batch_idx * batch_size
        batch_end = min((batch_idx + 1) * batch_size, len(accessions))
        batch_accessions = accessions[batch_start:batch_end]
        
        logger.info(f"Batch {batch_idx + 1}/{n_batches}: "
                   f"accessions {batch_start:,} to {batch_end:,}")
        
        # Register batch as temporary table
        batch_df = pd.DataFrame({'sample_id': batch_accessions})
        output_db.register('batch_accessions', batch_df)
        
        # Insert hashes for this batch
        batch_start_time = time.time()
        result = output_db.execute(f"""
            INSERT INTO sample_hashes
            SELECT s.sample_id, s.min_hash
            FROM source.sigs_dna.signature_mins s
            INNER JOIN batch_accessions b ON s.sample_id = b.sample_id
            WHERE s.ksize = {ksize}
        """)
        
        # Get count of inserted rows
        count_result = output_db.execute("""
            SELECT COUNT(*) FROM sample_hashes
        """).fetchone()[0]
        
        batch_hashes = count_result - total_hashes
        total_hashes = count_result
        
        batch_time = time.time() - batch_start_time
        elapsed = time.time() - start_time
        rate = batch_end / elapsed
        eta = (len(accessions) - batch_end) / rate / 60 if rate > 0 else 0
        
        logger.info(f"  Inserted {batch_hashes:,} hashes in {batch_time:.1f}s "
                   f"(total: {total_hashes:,}, ETA: {eta:.1f} min)")
        
        output_db.unregister('batch_accessions')
    
    # Detach source
    output_db.execute("DETACH source")
    
    total_time = time.time() - start_time
    logger.info(f"\nSample hash extraction complete:")
    logger.info(f"  Total hashes: {total_hashes:,}")
    logger.info(f"  Time: {total_time/60:.1f} minutes")
    
    return total_hashes


def load_reference_hashes(
    reference_path: str,
    output_db: duckdb.DuckDBPyConnection,
    chunk_size: int = 100_000_000
):
    """
    Load reference hashes from binary file into DuckDB table.
    
    The binary file contains sorted uint64 values. We load them in chunks
    to manage memory usage.
    """
    logger.info("="*60)
    logger.info("STEP 2: Loading reference hashes")
    logger.info("="*60)
    
    ref_path = Path(reference_path)
    file_size = ref_path.stat().st_size
    n_hashes = file_size // 8
    
    logger.info(f"Reference file: {reference_path}")
    logger.info(f"Total hashes: {n_hashes:,} ({file_size/1e9:.2f} GB)")
    
    # Create reference_hashes table
    logger.info("Creating reference_hashes table...")
    output_db.execute("""
        CREATE TABLE IF NOT EXISTS reference_hashes (
            hash UBIGINT
        )
    """)
    
    # Memory-map the file and load in chunks
    logger.info(f"Loading in chunks of {chunk_size:,} hashes...")
    
    reference = np.memmap(ref_path, dtype=np.uint64, mode='r', shape=(n_hashes,))
    
    start_time = time.time()
    total_loaded = 0
    
    for i in range(0, n_hashes, chunk_size):
        chunk_end = min(i + chunk_size, n_hashes)
        chunk = reference[i:chunk_end].copy()  # Copy to get contiguous array
        
        # Register as temporary table and insert
        chunk_df = pd.DataFrame({'hash': chunk})
        output_db.register('chunk_data', chunk_df)
        
        output_db.execute("""
            INSERT INTO reference_hashes
            SELECT hash FROM chunk_data
        """)
        
        output_db.unregister('chunk_data')
        del chunk_df, chunk
        
        total_loaded = chunk_end
        elapsed = time.time() - start_time
        rate = total_loaded / elapsed
        eta = (n_hashes - total_loaded) / rate / 60 if rate > 0 else 0
        
        logger.info(f"  Loaded {total_loaded:,}/{n_hashes:,} hashes "
                   f"({100*total_loaded/n_hashes:.1f}%), ETA: {eta:.1f} min")
    
    del reference
    
    total_time = time.time() - start_time
    logger.info(f"\nReference hash loading complete:")
    logger.info(f"  Total hashes: {total_loaded:,}")
    logger.info(f"  Time: {total_time/60:.1f} minutes")
    
    return total_loaded


def create_indexes(output_db: duckdb.DuckDBPyConnection):
    """
    Create indexes for fast joins and lookups.
    """
    logger.info("="*60)
    logger.info("STEP 3: Creating indexes")
    logger.info("="*60)
    
    # Index on reference_hashes.hash for fast lookups during JOIN
    logger.info("Creating index on reference_hashes(hash)...")
    start = time.time()
    output_db.execute("""
        CREATE INDEX IF NOT EXISTS idx_reference_hash 
        ON reference_hashes(hash)
    """)
    logger.info(f"  Done in {time.time()-start:.1f}s")
    
    # Index on sample_hashes.sample_id for fast GROUP BY
    logger.info("Creating index on sample_hashes(sample_id)...")
    start = time.time()
    output_db.execute("""
        CREATE INDEX IF NOT EXISTS idx_sample_id 
        ON sample_hashes(sample_id)
    """)
    logger.info(f"  Done in {time.time()-start:.1f}s")
    
    # Index on sample_hashes.min_hash for fast JOIN
    logger.info("Creating index on sample_hashes(min_hash)...")
    start = time.time()
    output_db.execute("""
        CREATE INDEX IF NOT EXISTS idx_sample_hash 
        ON sample_hashes(min_hash)
    """)
    logger.info(f"  Done in {time.time()-start:.1f}s")
    
    logger.info("Index creation complete")


def analyze_database(output_db: duckdb.DuckDBPyConnection):
    """
    Run ANALYZE to update statistics for query optimization.
    """
    logger.info("="*60)
    logger.info("STEP 4: Analyzing tables")
    logger.info("="*60)
    
    logger.info("Running ANALYZE on sample_hashes...")
    start = time.time()
    output_db.execute("ANALYZE sample_hashes")
    logger.info(f"  Done in {time.time()-start:.1f}s")
    
    logger.info("Running ANALYZE on reference_hashes...")
    start = time.time()
    output_db.execute("ANALYZE reference_hashes")
    logger.info(f"  Done in {time.time()-start:.1f}s")


def print_summary(output_db: duckdb.DuckDBPyConnection, output_path: str):
    """
    Print summary of the created database.
    """
    logger.info("="*60)
    logger.info("DATABASE SUMMARY")
    logger.info("="*60)
    
    # Table sizes
    sample_count = output_db.execute("SELECT COUNT(*) FROM sample_hashes").fetchone()[0]
    ref_count = output_db.execute("SELECT COUNT(*) FROM reference_hashes").fetchone()[0]
    unique_samples = output_db.execute("SELECT COUNT(DISTINCT sample_id) FROM sample_hashes").fetchone()[0]
    
    logger.info(f"sample_hashes table:")
    logger.info(f"  Total rows: {sample_count:,}")
    logger.info(f"  Unique samples: {unique_samples:,}")
    logger.info(f"  Avg hashes/sample: {sample_count/unique_samples:,.0f}")
    
    logger.info(f"\nreference_hashes table:")
    logger.info(f"  Total rows: {ref_count:,}")
    
    # File size
    db_size = Path(output_path).stat().st_size
    logger.info(f"\nDatabase file size: {db_size/1e9:.2f} GB")
    
    logger.info(f"\nOutput: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Create specialized DuckDB database for DMI computation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script creates a compact DuckDB database optimized for computing the
Dark Matter Index. The resulting database contains:
  - sample_hashes: Hashes for your samples of interest
  - reference_hashes: The complete reference hash set

With proper indexes, DMI computation becomes a simple, fast SQL query.

Example:
    python create_dmi_database.py \\
        --source /scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db \\
        --samples /path/to/filtered_data.parquet \\
        --reference /path/to/reference_hashes_k31.bin \\
        --output /path/to/dmi_database.db \\
        --ksize 31
        """
    )
    
    parser.add_argument(
        "--source", "-s",
        required=True,
        help="Path to source DuckDB database (database_all.db)"
    )
    parser.add_argument(
        "--samples",
        required=True,
        help="Path to parquet file with 'accession' column"
    )
    parser.add_argument(
        "--reference", "-r",
        required=True,
        help="Path to reference hashes binary file"
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
        help="k-mer size (REQUIRED, e.g., 31)"
    )
    parser.add_argument(
        "--batch-size", "-b",
        type=int,
        default=50000,
        help="Batch size for sample extraction (default: 50000)"
    )
    parser.add_argument(
        "--skip-reference",
        action="store_true",
        help="Skip loading reference hashes (if already loaded)"
    )
    parser.add_argument(
        "--skip-samples",
        action="store_true",
        help="Skip extracting sample hashes (if already loaded)"
    )
    
    args = parser.parse_args()
    
    total_start = time.time()
    
    logger.info("="*60)
    logger.info("CREATE DMI DATABASE")
    logger.info("="*60)
    logger.info(f"Source database: {args.source}")
    logger.info(f"Samples parquet: {args.samples}")
    logger.info(f"Reference hashes: {args.reference}")
    logger.info(f"Output database: {args.output}")
    logger.info(f"ksize: {args.ksize}")
    
    # Create/open output database
    output_db = duckdb.connect(args.output)
    
    # Set DuckDB to use all available threads
    output_db.execute("SET threads TO 200")
    output_db.execute("SET memory_limit = '3000GB'")
    
    try:
        # Step 1: Extract sample hashes
        if not args.skip_samples:
            extract_sample_hashes(
                args.source,
                args.samples,
                output_db,
                args.ksize,
                batch_size=args.batch_size
            )
        else:
            logger.info("Skipping sample hash extraction (--skip-samples)")
        
        # Step 2: Load reference hashes
        if not args.skip_reference:
            load_reference_hashes(
                args.reference,
                output_db
            )
        else:
            logger.info("Skipping reference hash loading (--skip-reference)")
        
        # Step 3: Create indexes
        create_indexes(output_db)
        
        # Step 4: Analyze tables
        analyze_database(output_db)
        
        # Summary
        print_summary(output_db, args.output)
        
    finally:
        output_db.close()
    
    total_time = time.time() - total_start
    
    logger.info("="*60)
    logger.info(f"COMPLETE - Total time: {total_time/60:.1f} minutes")
    logger.info("="*60)
    logger.info(f"\nNext step: Run compute_dmi_native.py with this database")


if __name__ == "__main__":
    main()
