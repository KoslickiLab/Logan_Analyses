#!/usr/bin/env python3
"""
Compute DMI with Host Contamination Removal

This script computes the Dark Matter Index while excluding host-derived hashes
from the calculation. It's designed for human-associated metagenome samples
where host contamination would artificially inflate the mapped hash count.

Key features:
1. Filter input samples by organism type (e.g., "human gut metagenome")
2. Exclude host hashes from DMI calculation
3. Output compatible with existing analysis scripts

The modified DMI calculation:
    filtered_sample_hashes = sample_hashes - host_hashes
    DMI = |filtered_sample_hashes - reference_hashes| / |filtered_sample_hashes|

Prerequisites:
    1. Run extract_human_hashes.py to create the host hash database
    2. Run create_dmi_database.py to create the sample/reference database

Usage:
    python compute_dmi_host_filtered.py \\
        --database /path/to/dmi_database.db \\
        --host-hashes /path/to/human_hashes.db \\
        --input /path/to/filtered_data.parquet \\
        --output /path/to/filtered_data_with_dmi_host_removed.parquet \\
        --organisms "human gut metagenome" "human feces metagenome"

Author: David Koslicki lab
Date: 2026
"""

import argparse
import logging
import time
from pathlib import Path
from typing import List, Optional
import sys

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

try:
    import duckdb
    import pandas as pd
    import pyarrow.parquet as pq
except ImportError as e:
    logger.error(f"Missing dependency: {e}")
    logger.error("Install with: pip install duckdb pandas pyarrow")
    sys.exit(1)


def verify_database(conn: duckdb.DuckDBPyConnection):
    """
    Verify the main DMI database has the required tables.
    """
    logger.info("Verifying main database structure...")
    
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
    
    sample_count = conn.execute("SELECT COUNT(*) FROM sample_hashes").fetchone()[0]
    ref_count = conn.execute("SELECT COUNT(*) FROM reference_hashes").fetchone()[0]
    unique_samples = conn.execute("SELECT COUNT(DISTINCT sample_id) FROM sample_hashes").fetchone()[0]
    
    logger.info(f"  sample_hashes: {sample_count:,} rows, {unique_samples:,} unique samples")
    logger.info(f"  reference_hashes: {ref_count:,} rows")
    
    return unique_samples


def verify_host_database(host_db_path: str) -> int:
    """
    Verify the host hash database and return hash count.
    """
    logger.info(f"Verifying host database: {host_db_path}")
    
    conn = duckdb.connect(host_db_path, read_only=True)
    
    try:
        tables = conn.execute("""
            SELECT table_name 
            FROM information_schema.tables 
            WHERE table_schema = 'main'
        """).fetchall()
        table_names = [t[0] for t in tables]
        
        if 'human_hashes' not in table_names:
            raise ValueError("Host database missing 'human_hashes' table. "
                           "Run extract_human_hashes.py first.")
        
        count = conn.execute("SELECT COUNT(*) FROM human_hashes").fetchone()[0]
        logger.info(f"  Host hashes: {count:,}")
        
        # Get metadata if available
        if 'metadata' in table_names:
            metadata = dict(conn.execute("SELECT key, value FROM metadata").fetchall())
            logger.info(f"  Source: {metadata.get('source_path', 'unknown')}")
            logger.info(f"  k-mer size: {metadata.get('ksize', 'unknown')}")
        
        return count
        
    finally:
        conn.close()


def load_organisms_from_file(filepath: str) -> List[str]:
    """
    Load organism names from a file (one per line).
    
    Args:
        filepath: Path to file with organism names
        
    Returns:
        List of organism names (stripped, non-empty lines)
    """
    organisms = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):  # Skip empty lines and comments
                organisms.append(line)
    return organisms


def parse_organisms_arg(organisms_arg: Optional[List[str]], organisms_file: Optional[str]) -> Optional[List[str]]:
    """
    Parse organisms from either command line or file.
    
    Handles comma-separated values if passed as a single string.
    
    Args:
        organisms_arg: List from --organisms argument
        organisms_file: Path from --organisms-file argument
        
    Returns:
        List of organism names, or None if no filter
    """
    organisms = []
    
    # Load from file if specified
    if organisms_file:
        if not Path(organisms_file).exists():
            logger.error(f"Organisms file not found: {organisms_file}")
            sys.exit(1)
        file_organisms = load_organisms_from_file(organisms_file)
        logger.info(f"Loaded {len(file_organisms)} organisms from {organisms_file}")
        organisms.extend(file_organisms)
    
    # Process command line arguments
    if organisms_arg:
        for org in organisms_arg:
            # Handle comma-separated values
            if ',' in org:
                organisms.extend([o.strip() for o in org.split(',') if o.strip()])
            else:
                organisms.append(org)
    
    if not organisms:
        return None
    
    # Deduplicate while preserving order
    seen = set()
    unique_organisms = []
    for org in organisms:
        org_lower = org.lower()
        if org_lower not in seen:
            seen.add(org_lower)
            unique_organisms.append(org)
    
    return unique_organisms


def filter_samples_by_organism(
    input_parquet: str,
    organisms: Optional[List[str]] = None,
    organism_column: str = 'organism'
) -> pd.DataFrame:
    """
    Load and filter input parquet by organism type.
    
    Args:
        input_parquet: Path to input parquet file
        organisms: List of organism names to include (None = all)
        organism_column: Name of the organism column
        
    Returns:
        Filtered DataFrame
    """
    logger.info(f"Loading input data from {input_parquet}")
    df = pd.read_parquet(input_parquet)
    logger.info(f"Total samples: {len(df):,}")
    
    if organisms is None or len(organisms) == 0:
        logger.info("No organism filter specified, using all samples")
        return df
    
    if organism_column not in df.columns:
        logger.warning(f"Column '{organism_column}' not found in input. "
                      f"Available columns: {list(df.columns)}")
        logger.warning("Proceeding without organism filtering")
        return df
    
    # Filter by organism
    logger.info(f"Filtering by {len(organisms)} organism type(s):")
    for org in organisms:
        logger.info(f"  - {org}")
    
    # Create case-insensitive filter
    organisms_lower = [o.lower() for o in organisms]
    mask = df[organism_column].str.lower().isin(organisms_lower)
    
    df_filtered = df[mask].copy()
    
    logger.info(f"Samples after organism filter: {len(df_filtered):,}")
    
    # Show breakdown by organism
    if len(df_filtered) > 0:
        logger.info("Samples per organism type:")
        for org in organisms:
            count = (df_filtered[organism_column].str.lower() == org.lower()).sum()
            if count > 0:
                logger.info(f"  - {org}: {count:,}")
    
    if len(df_filtered) == 0:
        logger.warning("No samples match the organism filter!")
        logger.info(f"Available organisms (top 20): "
                   f"{df[organism_column].value_counts().head(20).index.tolist()}")
    
    return df_filtered


def compute_dmi_with_host_removal(
    conn: duckdb.DuckDBPyConnection,
    host_db_path: str,
    sample_ids: List[str],
    chunk_size: int = 5000
) -> pd.DataFrame:
    """
    Compute DMI for samples, excluding host hashes from the calculation.
    
    The algorithm:
    1. For each sample, get its hashes
    2. Exclude any hashes that appear in the host database
    3. Of the remaining hashes, count how many are in the reference
    4. DMI = unmapped_non_host / total_non_host
    
    Also tracks:
    - Original hash counts (before host removal)
    - Host hash counts (removed)
    - Non-host hash counts (used for DMI)
    
    Args:
        conn: Connection to main DMI database
        host_db_path: Path to host hash database
        sample_ids: List of sample IDs to process
        chunk_size: Number of samples per batch
        
    Returns:
        DataFrame with DMI results
    """
    logger.info("="*60)
    logger.info("Computing DMI with host contamination removal")
    logger.info("="*60)
    
    # Attach host database
    logger.info(f"Attaching host database: {host_db_path}")
    conn.execute(f"ATTACH '{host_db_path}' AS host (READ_ONLY)")
    
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
        
        # Complex query that:
        # 1. Gets sample hashes
        # 2. Marks which are host hashes
        # 3. For non-host hashes, marks which are in reference
        # 4. Aggregates by sample
        query = """
            WITH sample_hash_status AS (
                SELECT 
                    s.sample_id,
                    s.min_hash,
                    CASE WHEN h.hash IS NOT NULL THEN 1 ELSE 0 END AS is_host,
                    CASE WHEN r.hash IS NOT NULL THEN 1 ELSE 0 END AS is_reference
                FROM sample_hashes s
                INNER JOIN chunk_samples c ON s.sample_id = c.sample_id
                LEFT JOIN host.human_hashes h ON s.min_hash = h.hash
                LEFT JOIN reference_hashes r ON s.min_hash = r.hash
            )
            SELECT 
                sample_id AS accession,
                
                -- Original counts (before host removal)
                COUNT(*) AS total_hashes_original,
                
                -- Host contamination
                SUM(is_host) AS host_hashes,
                
                -- Non-host counts (used for DMI)
                SUM(CASE WHEN is_host = 0 THEN 1 ELSE 0 END) AS total_hashes,
                
                -- Reference matching (non-host only)
                SUM(CASE WHEN is_host = 0 AND is_reference = 1 THEN 1 ELSE 0 END) AS mapped_hashes,
                SUM(CASE WHEN is_host = 0 AND is_reference = 0 THEN 1 ELSE 0 END) AS unmapped_hashes,
                
                -- DMI calculation (based on non-host hashes only)
                CASE 
                    WHEN SUM(CASE WHEN is_host = 0 THEN 1 ELSE 0 END) > 0 
                    THEN SUM(CASE WHEN is_host = 0 AND is_reference = 0 THEN 1 ELSE 0 END) * 1.0 
                         / SUM(CASE WHEN is_host = 0 THEN 1 ELSE 0 END)
                    ELSE NULL 
                END AS dmi,
                
                -- Host contamination fraction
                SUM(is_host) * 1.0 / COUNT(*) AS host_fraction
                
            FROM sample_hash_status
            GROUP BY sample_id
        """
        
        chunk_result = conn.execute(query).fetchdf()
        all_results.append(chunk_result)
        
        conn.unregister('chunk_samples')
        
        elapsed = time.time() - start_time
        rate = chunk_end / elapsed if elapsed > 0 else 0
        eta = (len(sample_ids) - chunk_end) / rate / 60 if rate > 0 else 0
        
        logger.info(f"Chunk {chunk_idx + 1}/{n_chunks}: "
                   f"processed {chunk_end:,}/{len(sample_ids):,} samples "
                   f"({100*chunk_end/len(sample_ids):.1f}%), ETA: {eta:.1f} min")
    
    # Detach host database
    conn.execute("DETACH host")
    
    # Combine all chunks
    result_df = pd.concat(all_results, ignore_index=True)
    
    total_time = time.time() - start_time
    logger.info(f"\nDMI computation complete in {total_time/60:.1f} minutes")
    
    return result_df


def compute_dmi_sql_with_host_removal(
    conn: duckdb.DuckDBPyConnection,
    host_db_path: str
) -> pd.DataFrame:
    """
    Compute DMI for ALL samples at once (non-chunked version).
    
    Use this for smaller datasets where progress tracking isn't critical.
    """
    logger.info("="*60)
    logger.info("Computing DMI with host removal (single query)")
    logger.info("="*60)
    
    # Attach host database
    logger.info(f"Attaching host database: {host_db_path}")
    conn.execute(f"ATTACH '{host_db_path}' AS host (READ_ONLY)")
    
    logger.info("Executing DMI query (this may take a while)...")
    start_time = time.time()
    
    query = """
        WITH sample_hash_status AS (
            SELECT 
                s.sample_id,
                s.min_hash,
                CASE WHEN h.hash IS NOT NULL THEN 1 ELSE 0 END AS is_host,
                CASE WHEN r.hash IS NOT NULL THEN 1 ELSE 0 END AS is_reference
            FROM sample_hashes s
            LEFT JOIN host.human_hashes h ON s.min_hash = h.hash
            LEFT JOIN reference_hashes r ON s.min_hash = r.hash
        )
        SELECT 
            sample_id AS accession,
            COUNT(*) AS total_hashes_original,
            SUM(is_host) AS host_hashes,
            SUM(CASE WHEN is_host = 0 THEN 1 ELSE 0 END) AS total_hashes,
            SUM(CASE WHEN is_host = 0 AND is_reference = 1 THEN 1 ELSE 0 END) AS mapped_hashes,
            SUM(CASE WHEN is_host = 0 AND is_reference = 0 THEN 1 ELSE 0 END) AS unmapped_hashes,
            CASE 
                WHEN SUM(CASE WHEN is_host = 0 THEN 1 ELSE 0 END) > 0 
                THEN SUM(CASE WHEN is_host = 0 AND is_reference = 0 THEN 1 ELSE 0 END) * 1.0 
                     / SUM(CASE WHEN is_host = 0 THEN 1 ELSE 0 END)
                ELSE NULL 
            END AS dmi,
            SUM(is_host) * 1.0 / COUNT(*) AS host_fraction
        FROM sample_hash_status
        GROUP BY sample_id
    """
    
    result_df = conn.execute(query).fetchdf()
    
    elapsed = time.time() - start_time
    logger.info(f"Query complete in {elapsed/60:.1f} minutes")
    
    conn.execute("DETACH host")
    
    return result_df


def print_summary_statistics(df: pd.DataFrame):
    """
    Print comprehensive summary statistics.
    
    Args:
        df: DataFrame with DMI results already merged in (df_output)
    """
    logger.info("\n" + "="*60)
    logger.info("DMI SUMMARY (Host Contamination Removed)")
    logger.info("="*60)
    
    if 'dmi' not in df.columns:
        logger.warning("No 'dmi' column found in output - skipping summary")
        return
    
    valid_dmi = df['dmi'].dropna()
    
    if len(valid_dmi) == 0:
        logger.warning("No valid DMI values computed!")
        return
    
    logger.info(f"Samples with valid DMI: {len(valid_dmi):,}")
    
    # DMI statistics
    logger.info(f"\nDMI Statistics:")
    logger.info(f"  Mean:   {valid_dmi.mean():.4f}")
    logger.info(f"  Median: {valid_dmi.median():.4f}")
    logger.info(f"  Std:    {valid_dmi.std():.4f}")
    logger.info(f"  Min:    {valid_dmi.min():.4f}")
    logger.info(f"  Max:    {valid_dmi.max():.4f}")
    
    # DMI distribution
    logger.info("\nDMI Distribution:")
    logger.info(f"  DMI < 0.1:  {(valid_dmi < 0.1).sum():,} ({100*(valid_dmi < 0.1).mean():.1f}%)")
    logger.info(f"  0.1-0.3:    {((valid_dmi >= 0.1) & (valid_dmi < 0.3)).sum():,} ({100*((valid_dmi >= 0.1) & (valid_dmi < 0.3)).mean():.1f}%)")
    logger.info(f"  0.3-0.5:    {((valid_dmi >= 0.3) & (valid_dmi < 0.5)).sum():,} ({100*((valid_dmi >= 0.3) & (valid_dmi < 0.5)).mean():.1f}%)")
    logger.info(f"  0.5-0.7:    {((valid_dmi >= 0.5) & (valid_dmi < 0.7)).sum():,} ({100*((valid_dmi >= 0.5) & (valid_dmi < 0.7)).mean():.1f}%)")
    logger.info(f"  DMI >= 0.7: {(valid_dmi >= 0.7).sum():,} ({100*(valid_dmi >= 0.7).mean():.1f}%)")
    
    # Host contamination statistics
    if 'host_fraction' in df.columns:
        valid_host = df['host_fraction'].dropna()
        if len(valid_host) > 0:
            logger.info(f"\nHost Contamination Statistics:")
            logger.info(f"  Mean host fraction:   {valid_host.mean():.4f} ({100*valid_host.mean():.2f}%)")
            logger.info(f"  Median host fraction: {valid_host.median():.4f} ({100*valid_host.median():.2f}%)")
            logger.info(f"  Max host fraction:    {valid_host.max():.4f} ({100*valid_host.max():.2f}%)")
            
            # Host contamination distribution
            logger.info("\nHost Contamination Distribution:")
            logger.info(f"  < 1%:  {(valid_host < 0.01).sum():,} ({100*(valid_host < 0.01).mean():.1f}%)")
            logger.info(f"  1-5%:  {((valid_host >= 0.01) & (valid_host < 0.05)).sum():,} ({100*((valid_host >= 0.01) & (valid_host < 0.05)).mean():.1f}%)")
            logger.info(f"  5-10%: {((valid_host >= 0.05) & (valid_host < 0.10)).sum():,} ({100*((valid_host >= 0.05) & (valid_host < 0.10)).mean():.1f}%)")
            logger.info(f"  > 10%: {(valid_host >= 0.10).sum():,} ({100*(valid_host >= 0.10).mean():.1f}%)")
    
    # Hash statistics
    if 'total_hashes_original' in df.columns and 'total_hashes_dmi' in df.columns:
        orig_hashes = df['total_hashes_original'].sum()
        final_hashes = df['total_hashes_dmi'].sum()
        removed_hashes = orig_hashes - final_hashes
        
        if orig_hashes > 0:
            logger.info(f"\nHash Statistics:")
            logger.info(f"  Original total hashes: {orig_hashes:,}")
            logger.info(f"  Host hashes removed:   {removed_hashes:,} ({100*removed_hashes/orig_hashes:.2f}%)")
            logger.info(f"  Final hashes for DMI:  {final_hashes:,}")
    
    # By organism if available
    if 'organism' in df.columns:
        # Filter to rows with valid DMI
        df_valid = df[df['dmi'].notna()].copy()
        
        if len(df_valid) > 0:
            logger.info("\nTop 10 organisms by sample count:")
            
            # Build aggregation based on available columns
            agg_dict = {'dmi': ['median', 'mean', 'count']}
            if 'host_fraction' in df_valid.columns:
                agg_dict['host_fraction'] = 'median'
            
            try:
                organism_stats = df_valid.groupby('organism').agg(agg_dict)
                
                # Flatten multi-level column names
                organism_stats.columns = ['_'.join(col).strip('_') if isinstance(col, tuple) else col 
                                          for col in organism_stats.columns]
                
                organism_stats = organism_stats.sort_values('dmi_count', ascending=False).head(10)
                
                for org, row in organism_stats.iterrows():
                    count = int(row['dmi_count'])
                    dmi_med = row['dmi_median']
                    
                    if 'host_fraction_median' in row.index and pd.notna(row['host_fraction_median']):
                        logger.info(f"  {org}: n={count:,}, "
                                   f"DMI median={dmi_med:.4f}, "
                                   f"host median={row['host_fraction_median']:.4f}")
                    else:
                        logger.info(f"  {org}: n={count:,}, DMI median={dmi_med:.4f}")
            except Exception as e:
                logger.warning(f"Could not compute organism statistics: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute DMI with host contamination removal",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script computes the Dark Matter Index while excluding host-derived hashes.
Designed for human-associated metagenomes where host contamination would
artificially inflate the mapped hash count.

Prerequisites:
    1. Run extract_human_hashes.py to create the host hash database
    2. Run create_dmi_database.py to create the sample/reference database

Output columns added to parquet:
    - dmi: Dark Matter Index (computed on non-host hashes)
    - total_hashes: Non-host hashes (used for DMI calculation)
    - total_hashes_original: Original hash count (before host removal)
    - host_hashes: Number of host hashes removed
    - host_fraction: Fraction of hashes that were host-derived
    - mapped_hashes: Non-host hashes found in reference
    - unmapped_hashes: Non-host hashes not in reference

Examples:
    # RECOMMENDED: Use a file with organism names (one per line)
    python compute_dmi_host_filtered.py \\
        --database dmi_database.db \\
        --host-hashes human_hashes.db \\
        --input filtered_data.parquet \\
        --output filtered_data_with_dmi_host_removed.parquet \\
        --organisms-file human_organisms.txt

    # Alternative: Comma-separated format (avoids shell quoting issues)
    python compute_dmi_host_filtered.py \\
        --database dmi_database.db \\
        --host-hashes human_hashes.db \\
        --input filtered_data.parquet \\
        --output filtered_data_with_dmi_host_removed.parquet \\
        --organisms "human gut metagenome,human feces metagenome,human skin metagenome"

    # All samples (no organism filter)
    python compute_dmi_host_filtered.py \\
        --database dmi_database.db \\
        --host-hashes human_hashes.db \\
        --input filtered_data.parquet \\
        --output filtered_data_with_dmi_host_removed.parquet

    # List available organisms in data
    python compute_dmi_host_filtered.py \\
        --database dmi_database.db \\
        --host-hashes human_hashes.db \\
        --input filtered_data.parquet \\
        --output out.parquet \\
        --list-organisms
        """
    )
    
    parser.add_argument(
        "--database", "-d",
        required=True,
        help="Path to main DMI database (from create_dmi_database.py)"
    )
    parser.add_argument(
        "--host-hashes", "-H",
        required=True,
        help="Path to host hash database (from extract_human_hashes.py)"
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
        "--organisms",
        nargs='+',
        default=None,
        help="Organism types to include. Use quotes around multi-word names or "
             "use comma-separated format: 'human gut metagenome,human feces metagenome'. "
             "Alternatively, use --organisms-file for a file-based list."
    )
    parser.add_argument(
        "--organisms-file",
        default=None,
        help="Path to file with organism names (one per line). "
             "Lines starting with # are treated as comments. "
             "Can be combined with --organisms."
    )
    parser.add_argument(
        "--organism-column",
        default='organism',
        help="Name of organism column in input parquet (default: 'organism')"
    )
    parser.add_argument(
        "--chunked",
        action="store_true",
        help="Process in chunks (shows progress, uses less memory)"
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=5000,
        help="Number of samples per chunk (default: 5000)"
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
    parser.add_argument(
        "--list-organisms",
        action="store_true",
        help="List available organisms in input parquet and exit"
    )
    
    args = parser.parse_args()
    
    total_start = time.time()
    
    logger.info("="*60)
    logger.info("COMPUTE DMI WITH HOST CONTAMINATION REMOVAL")
    logger.info("="*60)
    logger.info(f"Main database: {args.database}")
    logger.info(f"Host hashes: {args.host_hashes}")
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    logger.info(f"Threads: {args.threads}")
    logger.info(f"Memory: {args.memory}")
    
    # Parse organisms from args and/or file
    organisms = parse_organisms_arg(args.organisms, args.organisms_file)
    if organisms:
        logger.info(f"Organism filter: {len(organisms)} type(s)")
    else:
        logger.info("Organism filter: none (all samples)")
    
    # Check files exist
    for path, name in [(args.database, "Main database"), 
                       (args.host_hashes, "Host hashes"),
                       (args.input, "Input parquet")]:
        if not Path(path).exists():
            logger.error(f"{name} not found: {path}")
            sys.exit(1)
    
    # Handle --list-organisms
    if args.list_organisms:
        df = pd.read_parquet(args.input)
        if args.organism_column in df.columns:
            logger.info(f"\nAvailable organisms in {args.input}:")
            counts = df[args.organism_column].value_counts()
            for org, count in counts.head(50).items():
                logger.info(f"  {org}: {count:,}")
            if len(counts) > 50:
                logger.info(f"  ... and {len(counts) - 50} more")
        else:
            logger.error(f"Column '{args.organism_column}' not found")
        sys.exit(0)
    
    # Connect to main database
    conn = duckdb.connect(args.database, read_only=True)
    conn.execute(f"SET threads TO {args.threads}")
    conn.execute(f"SET memory_limit = '{args.memory}'")
    
    try:
        # Verify databases
        verify_database(conn)
        verify_host_database(args.host_hashes)
        
        # Load and filter input data
        df_input = filter_samples_by_organism(
            args.input,
            organisms,  # Use parsed organisms list
            args.organism_column
        )
        
        if len(df_input) == 0:
            logger.error("No samples to process after filtering!")
            sys.exit(1)
        
        # Get sample IDs that are actually in the database
        logger.info("Finding samples in database...")
        sample_ids = df_input['accession'].tolist()
        
        # Check which samples exist in database
        sample_df = pd.DataFrame({'sample_id': sample_ids})
        conn.register('requested_samples', sample_df)
        
        existing = conn.execute("""
            SELECT DISTINCT s.sample_id 
            FROM sample_hashes s
            INNER JOIN requested_samples r ON s.sample_id = r.sample_id
        """).fetchall()
        
        conn.unregister('requested_samples')
        
        existing_ids = [e[0] for e in existing]
        
        if len(existing_ids) < len(sample_ids):
            logger.warning(f"{len(sample_ids) - len(existing_ids):,} samples not found in database")
        
        logger.info(f"Processing {len(existing_ids):,} samples")
        
        # Compute DMI
        if args.chunked or len(existing_ids) > 10000:
            dmi_results = compute_dmi_with_host_removal(
                conn,
                args.host_hashes,
                existing_ids,
                chunk_size=args.chunk_size
            )
        else:
            dmi_results = compute_dmi_sql_with_host_removal(conn, args.host_hashes)
        
        # Rename total_hashes_dmi for compatibility with existing scripts
        dmi_results = dmi_results.rename(columns={'total_hashes': 'total_hashes_dmi'})
        
        # Merge results with input
        logger.info("Merging results with input data...")
        
        # Select columns to merge
        merge_cols = ['accession', 'dmi', 'total_hashes_dmi', 'total_hashes_original',
                      'host_hashes', 'host_fraction', 'unmapped_hashes', 'mapped_hashes']
        merge_cols = [c for c in merge_cols if c in dmi_results.columns]
        
        df_output = df_input.merge(
            dmi_results[merge_cols],
            on='accession',
            how='left'
        )
        
        # Check for missing DMI
        missing_dmi = df_output['dmi'].isna().sum()
        if missing_dmi > 0:
            logger.warning(f"{missing_dmi:,} samples have no DMI (not in database)")
        
        # Print summary
        print_summary_statistics(df_output)
        
        # Save output
        logger.info(f"\nSaving output to {args.output}")
        df_output.to_parquet(args.output, index=False)
        
        # Output file info
        output_size = Path(args.output).stat().st_size
        logger.info(f"Output file size: {output_size/1e6:.1f} MB")
        
    finally:
        conn.close()
    
    total_time = time.time() - total_start
    
    logger.info("\n" + "="*60)
    logger.info(f"COMPLETE - Total time: {total_time/60:.1f} minutes")
    logger.info(f"Output: {args.output}")
    logger.info("="*60)


if __name__ == "__main__":
    main()
