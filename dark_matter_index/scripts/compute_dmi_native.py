#!/usr/bin/env python3
"""
Compute DMI using Native DuckDB Operations

This script computes the Dark Matter Index using pure DuckDB SQL operations,
leveraging DuckDB's optimized JOIN and GROUP BY implementations.

The approach:
1. LEFT JOIN sample_hashes with reference_hashes
2. GROUP BY sample_id and count matched/unmatched hashes
3. Calculate DMI = unmapped / total

This is dramatically faster than Python-based approaches because:
- DuckDB parallelizes the JOIN and aggregation automatically
- Indexes are used for efficient hash lookups
- No Python interpreter overhead for the heavy computation

Prerequisites:
    Run create_dmi_database.py first to create the indexed database.

Usage:
    python compute_dmi_native.py \\
        --database /path/to/dmi_database.db \\
        --input /path/to/filtered_data.parquet \\
        --output /path/to/filtered_data_with_dmi.parquet

Expected time: 10-60 minutes depending on data size and hardware.

Author: David Koslicki lab
Date: 2026
"""

import argparse
import logging
import time
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


def verify_database(conn: duckdb.DuckDBPyConnection):
    """
    Verify the database has the required tables.
    """
    logger.info("Verifying database structure...")
    
    # Check tables exist
    tables = conn.execute("""
        SELECT table_name 
        FROM information_schema.tables 
        WHERE table_schema = 'main'
    """).fetchall()
    table_names = [t[0] for t in tables]
    
    if 'sample_hashes' not in table_names:
        raise ValueError("Database missing 'sample_hashes' table. Run create_dmi_database.py first.")
    if 'reference_hashes' not in table_names:
        raise ValueError("Database missing 'reference_hashes' table. Run create_dmi_database.py first.")
    
    # Get counts
    sample_count = conn.execute("SELECT COUNT(*) FROM sample_hashes").fetchone()[0]
    ref_count = conn.execute("SELECT COUNT(*) FROM reference_hashes").fetchone()[0]
    unique_samples = conn.execute("SELECT COUNT(DISTINCT sample_id) FROM sample_hashes").fetchone()[0]
    
    logger.info(f"  sample_hashes: {sample_count:,} rows, {unique_samples:,} unique samples")
    logger.info(f"  reference_hashes: {ref_count:,} rows")
    
    return unique_samples


def compute_dmi_sql(conn: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """
    Compute DMI for all samples using a single SQL query.
    
    The query uses LEFT JOIN to find which hashes match the reference,
    then aggregates by sample to compute DMI.
    """
    logger.info("="*60)
    logger.info("Computing DMI via SQL")
    logger.info("="*60)
    
    # The main DMI computation query
    # Uses LEFT JOIN: sample hashes that don't match reference will have NULL
    query = """
        SELECT 
            s.sample_id AS accession,
            COUNT(*) AS total_hashes,
            SUM(CASE WHEN r.hash IS NULL THEN 1 ELSE 0 END) AS unmapped_hashes,
            SUM(CASE WHEN r.hash IS NOT NULL THEN 1 ELSE 0 END) AS mapped_hashes,
            SUM(CASE WHEN r.hash IS NULL THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS dmi
        FROM sample_hashes s
        LEFT JOIN reference_hashes r ON s.min_hash = r.hash
        GROUP BY s.sample_id
    """
    
    logger.info("Executing DMI query...")
    logger.info("This may take 10-60 minutes depending on data size...")
    
    start_time = time.time()
    
    # Execute and fetch results
    result_df = conn.execute(query).fetchdf()
    
    elapsed = time.time() - start_time
    
    logger.info(f"\nQuery complete in {elapsed/60:.1f} minutes")
    logger.info(f"Computed DMI for {len(result_df):,} samples")
    
    return result_df


def compute_dmi_chunked(conn: duckdb.DuckDBPyConnection, chunk_size: int = 10000) -> pd.DataFrame:
    """
    Compute DMI in chunks for very large datasets.
    
    This approach processes samples in batches to provide progress feedback
    and reduce memory usage.
    """
    logger.info("="*60)
    logger.info("Computing DMI via SQL (chunked)")
    logger.info("="*60)
    
    # Get list of unique samples
    logger.info("Getting list of samples...")
    samples = conn.execute("""
        SELECT DISTINCT sample_id FROM sample_hashes ORDER BY sample_id
    """).fetchall()
    sample_ids = [s[0] for s in samples]
    
    logger.info(f"Found {len(sample_ids):,} unique samples")
    
    all_results = []
    start_time = time.time()
    n_chunks = (len(sample_ids) + chunk_size - 1) // chunk_size
    
    for chunk_idx in range(n_chunks):
        chunk_start = chunk_idx * chunk_size
        chunk_end = min((chunk_idx + 1) * chunk_size, len(sample_ids))
        chunk_samples = sample_ids[chunk_start:chunk_end]
        
        # Register chunk as temp table
        chunk_df = pd.DataFrame({'sample_id': chunk_samples})
        conn.register('chunk_samples', chunk_df)
        
        # Compute DMI for this chunk
        query = """
            SELECT 
                s.sample_id AS accession,
                COUNT(*) AS total_hashes,
                SUM(CASE WHEN r.hash IS NULL THEN 1 ELSE 0 END) AS unmapped_hashes,
                SUM(CASE WHEN r.hash IS NOT NULL THEN 1 ELSE 0 END) AS mapped_hashes,
                SUM(CASE WHEN r.hash IS NULL THEN 1 ELSE 0 END) * 1.0 / COUNT(*) AS dmi
            FROM sample_hashes s
            INNER JOIN chunk_samples c ON s.sample_id = c.sample_id
            LEFT JOIN reference_hashes r ON s.min_hash = r.hash
            GROUP BY s.sample_id
        """
        
        chunk_result = conn.execute(query).fetchdf()
        all_results.append(chunk_result)
        
        conn.unregister('chunk_samples')
        
        elapsed = time.time() - start_time
        rate = chunk_end / elapsed
        eta = (len(sample_ids) - chunk_end) / rate / 60 if rate > 0 else 0
        
        logger.info(f"Chunk {chunk_idx + 1}/{n_chunks}: "
                   f"processed {chunk_end:,}/{len(sample_ids):,} samples "
                   f"({100*chunk_end/len(sample_ids):.1f}%), ETA: {eta:.1f} min")
    
    # Combine all chunks
    result_df = pd.concat(all_results, ignore_index=True)
    
    total_time = time.time() - start_time
    logger.info(f"\nComplete in {total_time/60:.1f} minutes")
    
    return result_df


def main():
    parser = argparse.ArgumentParser(
        description="Compute DMI using native DuckDB operations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script computes the Dark Matter Index using pure SQL, which is much
faster than Python-based approaches.

Prerequisites:
    First run create_dmi_database.py to create the indexed database.

Examples:
    # Basic usage
    python compute_dmi_native.py \\
        --database dmi_database.db \\
        --input filtered_data.parquet \\
        --output filtered_data_with_dmi.parquet
    
    # Chunked processing (shows progress, lower memory)
    python compute_dmi_native.py \\
        --database dmi_database.db \\
        --input filtered_data.parquet \\
        --output filtered_data_with_dmi.parquet \\
        --chunked --chunk-size 5000
        """
    )
    
    parser.add_argument(
        "--database", "-d",
        required=True,
        help="Path to DMI database (from create_dmi_database.py)"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to input parquet file with 'accession' column"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Path for output parquet file with DMI"
    )
    parser.add_argument(
        "--chunked",
        action="store_true",
        help="Process in chunks (shows progress, uses less memory)"
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10000,
        help="Number of samples per chunk (default: 10000)"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=128,
        help="Number of DuckDB threads (default: 128)"
    )
    parser.add_argument(
        "--memory",
        type=str,
        default="256GB",
        help="DuckDB memory limit (default: 256GB)"
    )
    
    args = parser.parse_args()
    
    total_start = time.time()
    
    logger.info("="*60)
    logger.info("COMPUTE DMI (Native DuckDB)")
    logger.info("="*60)
    logger.info(f"Database: {args.database}")
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    logger.info(f"Threads: {args.threads}")
    logger.info(f"Memory limit: {args.memory}")
    
    # Connect to database (read-only since we're just querying)
    conn = duckdb.connect(args.database, read_only=True)
    conn.execute(f"SET threads TO {args.threads}")
    conn.execute(f"SET memory_limit = '{args.memory}'")
    
    try:
        # Verify database
        n_samples = verify_database(conn)
        
        # Compute DMI
        if args.chunked:
            dmi_results = compute_dmi_chunked(conn, chunk_size=args.chunk_size)
        else:
            dmi_results = compute_dmi_sql(conn)
        
        # Load input parquet
        logger.info(f"\nLoading input data from {args.input}")
        df_input = pd.read_parquet(args.input)
        logger.info(f"Input has {len(df_input):,} rows")
        
        # Merge DMI results
        logger.info("Merging DMI results with input data...")
        
        # Rename total_hashes to avoid conflicts
        dmi_results = dmi_results.rename(columns={'total_hashes': 'total_hashes_dmi'})
        
        df_output = df_input.merge(
            dmi_results[['accession', 'dmi', 'total_hashes_dmi', 'unmapped_hashes', 'mapped_hashes']],
            on='accession',
            how='left'
        )
        
        # Check for samples without DMI
        missing_dmi = df_output['dmi'].isna().sum()
        if missing_dmi > 0:
            logger.warning(f"{missing_dmi:,} samples have no DMI (not in database)")
        
        # Save output
        logger.info(f"Saving output to {args.output}")
        df_output.to_parquet(args.output, index=False)
        
        # Summary statistics
        logger.info("\n" + "="*60)
        logger.info("DMI SUMMARY")
        logger.info("="*60)
        
        valid_dmi = df_output['dmi'].dropna()
        if len(valid_dmi) > 0:
            logger.info(f"Samples with DMI: {len(valid_dmi):,}")
            logger.info(f"Mean DMI:   {valid_dmi.mean():.4f}")
            logger.info(f"Median DMI: {valid_dmi.median():.4f}")
            logger.info(f"Std DMI:    {valid_dmi.std():.4f}")
            logger.info(f"Min DMI:    {valid_dmi.min():.4f}")
            logger.info(f"Max DMI:    {valid_dmi.max():.4f}")
            
            # Distribution
            logger.info("\nDMI Distribution:")
            logger.info(f"  DMI < 0.1:  {(valid_dmi < 0.1).sum():,} ({100*(valid_dmi < 0.1).mean():.1f}%)")
            logger.info(f"  0.1-0.3:    {((valid_dmi >= 0.1) & (valid_dmi < 0.3)).sum():,} ({100*((valid_dmi >= 0.1) & (valid_dmi < 0.3)).mean():.1f}%)")
            logger.info(f"  0.3-0.5:    {((valid_dmi >= 0.3) & (valid_dmi < 0.5)).sum():,} ({100*((valid_dmi >= 0.3) & (valid_dmi < 0.5)).mean():.1f}%)")
            logger.info(f"  0.5-0.7:    {((valid_dmi >= 0.5) & (valid_dmi < 0.7)).sum():,} ({100*((valid_dmi >= 0.5) & (valid_dmi < 0.7)).mean():.1f}%)")
            logger.info(f"  DMI >= 0.7: {(valid_dmi >= 0.7).sum():,} ({100*(valid_dmi >= 0.7).mean():.1f}%)")
        
        # By organism if available
        if 'organism' in df_output.columns and len(valid_dmi) > 0:
            logger.info("\nTop 10 organisms by median DMI (n >= 100):")
            organism_stats = df_output.groupby('organism')['dmi'].agg(['median', 'mean', 'count'])
            organism_stats = organism_stats[organism_stats['count'] >= 100]
            organism_stats = organism_stats.sort_values('median', ascending=False).head(10)
            for org, row in organism_stats.iterrows():
                logger.info(f"  {org}: median={row['median']:.4f}, mean={row['mean']:.4f}, n={int(row['count']):,}")
        
    finally:
        conn.close()
    
    total_time = time.time() - total_start
    
    logger.info("\n" + "="*60)
    logger.info(f"COMPLETE - Total time: {total_time/60:.1f} minutes")
    logger.info(f"Output: {args.output}")
    logger.info("="*60)


if __name__ == "__main__":
    main()
