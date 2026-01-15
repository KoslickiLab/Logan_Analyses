#!/usr/bin/env python3
"""
Build Reference Hash Set for Dark Matter Index Calculation

This script creates a unified, sorted binary file containing all unique hashes
from the reference databases (GTDB, GenBank, AllTheBacteria, SILVA, etc.).

The output is a binary file of sorted uint64 values that can be efficiently
queried using binary search.

Strategy for 20 billion hashes:
1. Query each database individually using DuckDB
2. Stream results directly to disk in sorted chunks (external merge sort)
3. Merge all chunks into final sorted, deduplicated file

Memory usage: Configurable chunk size (default ~10GB per chunk)
Disk usage: ~160GB for final file + temporary space during build

Usage:
    python build_reference_hashes.py /path/to/DB_info.json /path/to/output_dir
    
    # With custom chunk size (in billions of hashes)
    python build_reference_hashes.py DB_info.json output_dir --chunk-size 1.0

Author: David Koslicki lab
Date: 2026
"""

import numpy as np
import json
import heapq
import tempfile
import shutil
from pathlib import Path
from typing import List, Dict, Iterator, Tuple
from dataclasses import dataclass
import logging
import time
import argparse

# Try to import duckdb - will be used on HPC
try:
    import duckdb
    HAS_DUCKDB = True
except ImportError:
    HAS_DUCKDB = False
    print("Warning: duckdb not available. This script requires duckdb for database access.")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class DatabaseConfig:
    """Configuration for a single hash database."""
    name: str
    alias: str
    path: str
    table: str
    column: str
    set_id: str
    ksize: str  # column name for ksize filter, or null


def load_database_configs(config_path: str, set_id_filter: str = "right") -> List[DatabaseConfig]:
    """
    Load database configurations from JSON file.
    
    Args:
        config_path: Path to DB_info.json
        set_id_filter: Only include databases with this set_id
        
    Returns:
        List of DatabaseConfig objects
    """
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    databases = []
    for db in config['databases']:
        if db['set_id'] == set_id_filter:
            databases.append(DatabaseConfig(
                name=db['name'],
                alias=db['alias'],
                path=db['path'],
                table=db['table'],
                column=db['column'],
                set_id=db['set_id'],
                ksize=db.get('ksize')
            ))
    
    return databases


def stream_hashes_from_database(
    db_config: DatabaseConfig,
    ksize_value: int,
    batch_size: int = 10_000_000
) -> Iterator[np.ndarray]:
    """
    Stream hash values from a DuckDB database in batches.
    
    Uses cursor-based pagination to handle very large tables without
    loading everything into memory.
    
    NOTE: Individual databases are already deduplicated, so we use SELECT
    (not SELECT DISTINCT) here. Deduplication happens only on the final union.
    
    Args:
        db_config: Database configuration
        ksize_value: k-mer size to filter on (if applicable)
        batch_size: Number of rows to fetch per batch
        
    Yields:
        numpy arrays of uint64 hash values
    """
    logger.info(f"Streaming from {db_config.name} ({db_config.path})")
    
    # CRITICAL: Always connect in read-only mode to prevent any modifications
    conn = duckdb.connect(db_config.path, read_only=True)
    
    # Build query with optional ksize filter
    # NOTE: No DISTINCT here since individual DBs are already deduplicated
    if db_config.ksize and db_config.ksize != "null":
        query = f"""
            SELECT {db_config.column}
            FROM {db_config.table}
            WHERE {db_config.ksize} = {ksize_value}
        """
    else:
        query = f"""
            SELECT {db_config.column}
            FROM {db_config.table}
        """
    
    logger.info(f"Executing query: {query[:100]}...")
    
    # Execute query and stream results
    result = conn.execute(query)
    
    total_fetched = 0
    while True:
        batch = result.fetchmany(batch_size)
        if not batch:
            break
        
        # Convert to numpy array
        hashes = np.array([row[0] for row in batch], dtype=np.uint64)
        total_fetched += len(hashes)
        
        if total_fetched % (batch_size * 10) == 0:
            logger.info(f"  Fetched {total_fetched:,} hashes from {db_config.name}")
        
        yield hashes
    
    conn.close()
    logger.info(f"  Completed {db_config.name}: {total_fetched:,} hashes")


def count_hashes_in_database(db_config: DatabaseConfig, ksize_value: int) -> int:
    """
    Count the number of hashes in a database.
    
    NOTE: Individual databases are already deduplicated, so we count all rows.
    
    Args:
        db_config: Database configuration
        ksize_value: k-mer size to filter on
        
    Returns:
        Count of hashes
    """
    # CRITICAL: Always connect in read-only mode
    conn = duckdb.connect(db_config.path, read_only=True)
    
    if db_config.ksize and db_config.ksize != "null":
        query = f"""
            SELECT COUNT({db_config.column})
            FROM {db_config.table}
            WHERE {db_config.ksize} = {ksize_value}
        """
    else:
        query = f"""
            SELECT COUNT({db_config.column})
            FROM {db_config.table}
        """
    
    count = conn.execute(query).fetchone()[0]
    conn.close()
    
    return count


def write_sorted_chunk(hashes: np.ndarray, chunk_path: Path) -> int:
    """
    Sort hashes and write to a binary file.
    
    Args:
        hashes: numpy array of uint64 hashes
        chunk_path: Output file path
        
    Returns:
        Number of unique hashes written
    """
    # Sort in place
    hashes.sort()
    
    # Remove duplicates
    unique_hashes = np.unique(hashes)
    
    # Write to file
    unique_hashes.tofile(chunk_path)
    
    return len(unique_hashes)


def merge_sorted_files(
    chunk_paths: List[Path],
    output_path: Path,
    buffer_size: int = 100_000_000  # ~800MB buffer per file
) -> int:
    """
    Merge multiple sorted hash files into a single sorted, deduplicated file.
    
    Uses a min-heap for efficient k-way merge.
    
    Args:
        chunk_paths: List of paths to sorted chunk files
        output_path: Path for merged output file
        buffer_size: Number of hashes to buffer per input file
        
    Returns:
        Total number of unique hashes in merged file
    """
    logger.info(f"Merging {len(chunk_paths)} sorted chunks...")
    
    # Open all chunk files as memory-mapped arrays
    chunk_arrays = []
    chunk_indices = []
    chunk_sizes = []
    
    for path in chunk_paths:
        size = path.stat().st_size // 8
        chunk_sizes.append(size)
        chunk_arrays.append(np.memmap(path, dtype=np.uint64, mode='r', shape=(size,)))
        chunk_indices.append(0)  # Current position in each chunk
    
    # Initialize min-heap with first element from each chunk
    # Heap items: (hash_value, chunk_index)
    heap = []
    for i, arr in enumerate(chunk_arrays):
        if len(arr) > 0:
            heapq.heappush(heap, (arr[0], i))
            chunk_indices[i] = 1
    
    # Merge using output buffer
    output_buffer = np.zeros(buffer_size, dtype=np.uint64)
    buffer_pos = 0
    total_written = 0
    last_value = None
    
    # Open output file for appending
    with open(output_path, 'wb') as out_file:
        while heap:
            # Get smallest value
            value, chunk_idx = heapq.heappop(heap)
            
            # Skip duplicates
            if last_value is None or value != last_value:
                output_buffer[buffer_pos] = value
                buffer_pos += 1
                last_value = value
                
                # Flush buffer if full
                if buffer_pos >= buffer_size:
                    output_buffer[:buffer_pos].tofile(out_file)
                    total_written += buffer_pos
                    buffer_pos = 0
                    
                    if total_written % 1_000_000_000 == 0:
                        logger.info(f"  Written {total_written/1e9:.2f}B hashes")
            
            # Add next element from this chunk to heap
            if chunk_indices[chunk_idx] < chunk_sizes[chunk_idx]:
                next_val = chunk_arrays[chunk_idx][chunk_indices[chunk_idx]]
                heapq.heappush(heap, (next_val, chunk_idx))
                chunk_indices[chunk_idx] += 1
        
        # Write remaining buffer
        if buffer_pos > 0:
            output_buffer[:buffer_pos].tofile(out_file)
            total_written += buffer_pos
    
    logger.info(f"Merge complete: {total_written:,} unique hashes")
    return total_written


def build_reference_hashes(
    config_path: str,
    output_dir: str,
    chunk_size_billions: float,
    ksize: int
) -> Path:
    """
    Build the unified reference hash file from all databases using external merge sort.
    
    This method uses chunked processing with external merge sort to handle
    datasets larger than available RAM.
    
    Individual databases are already deduplicated. This function:
    1. Streams hashes from each database (no DISTINCT needed per-DB)
    2. Sorts and deduplicates within each chunk
    3. Merges all chunks, deduplicating across the union
    
    Args:
        config_path: Path to DB_info.json
        output_dir: Directory for output files
        chunk_size_billions: Size of sorting chunks in billions of hashes
        ksize: k-mer size to use (REQUIRED)
        
    Returns:
        Path to the final reference hash file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create temp directory for chunks
    temp_dir = output_dir / "temp_chunks"
    temp_dir.mkdir(exist_ok=True)
    
    # Load database configurations
    databases = load_database_configs(config_path, set_id_filter="right")
    logger.info(f"Found {len(databases)} reference databases")
    logger.info(f"ksize={ksize}")
    
    for db in databases:
        logger.info(f"  - {db.name}: {db.path}")
    
    # Calculate chunk size in number of hashes
    chunk_size = int(chunk_size_billions * 1e9)
    logger.info(f"Using chunk size of {chunk_size:,} hashes (~{chunk_size * 8 / 1e9:.1f} GB)")
    
    # Process each database and create sorted chunks
    chunk_paths = []
    current_chunk = np.zeros(chunk_size, dtype=np.uint64)
    chunk_pos = 0
    chunk_num = 0
    
    start_time = time.time()
    
    for db in databases:
        db_start = time.time()
        
        # Stream hashes (no DISTINCT - individual DBs already deduplicated)
        for batch in stream_hashes_from_database(db, ksize_value=ksize):
            # Add batch to current chunk
            batch_len = len(batch)
            
            while batch_len > 0:
                space_in_chunk = chunk_size - chunk_pos
                to_copy = min(batch_len, space_in_chunk)
                
                current_chunk[chunk_pos:chunk_pos + to_copy] = batch[:to_copy]
                chunk_pos += to_copy
                batch = batch[to_copy:]
                batch_len = len(batch)
                
                # Write chunk if full (sort and deduplicate within chunk)
                if chunk_pos >= chunk_size:
                    chunk_path = temp_dir / f"chunk_{chunk_num:04d}.bin"
                    n_unique = write_sorted_chunk(current_chunk[:chunk_pos], chunk_path)
                    chunk_paths.append(chunk_path)
                    logger.info(f"Wrote chunk {chunk_num}: {n_unique:,} unique hashes")
                    
                    chunk_num += 1
                    chunk_pos = 0
        
        db_time = time.time() - db_start
        logger.info(f"Processed {db.name} in {db_time:.1f}s")
    
    # Write final partial chunk
    if chunk_pos > 0:
        chunk_path = temp_dir / f"chunk_{chunk_num:04d}.bin"
        n_unique = write_sorted_chunk(current_chunk[:chunk_pos], chunk_path)
        chunk_paths.append(chunk_path)
        logger.info(f"Wrote final chunk {chunk_num}: {n_unique:,} unique hashes")
    
    # Merge all chunks (this also deduplicates across the union)
    output_path = output_dir / f"reference_hashes_k{ksize}.bin"
    total_unique = merge_sorted_files(chunk_paths, output_path)
    
    # Clean up temp files
    logger.info("Cleaning up temporary files...")
    shutil.rmtree(temp_dir)
    
    total_time = time.time() - start_time
    file_size = output_path.stat().st_size
    
    logger.info(f"\n{'='*60}")
    logger.info(f"BUILD COMPLETE (chunked method)")
    logger.info(f"{'='*60}")
    logger.info(f"Output file: {output_path}")
    logger.info(f"Total unique hashes: {total_unique:,}")
    logger.info(f"File size: {file_size / 1e9:.2f} GB")
    logger.info(f"Total time: {total_time / 60:.1f} minutes")
    logger.info(f"ksize: {ksize}")
    logger.info(f"{'='*60}")
    
    # Write metadata file
    metadata = {
        'n_hashes': total_unique,
        'file_size_bytes': file_size,
        'ksize': ksize,
        'databases': [db.name for db in databases],
        'build_time_seconds': total_time,
        'method': 'chunked'
    }
    
    with open(output_dir / "reference_metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)
    
    return output_path


def build_reference_hashes_memory_efficient(
    config_path: str,
    output_dir: str,
    ksize: int
) -> Path:
    """
    Alternative build method using DuckDB's UNION for cross-database deduplication.
    
    This approach uses DuckDB's query engine to handle the merge and deduplication
    across databases. Individual databases are already deduplicated, so we only
    need to deduplicate the union.
    
    Args:
        config_path: Path to DB_info.json
        output_dir: Directory for output files
        ksize: k-mer size to use (REQUIRED)
        
    Returns:
        Path to the final reference hash file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    databases = load_database_configs(config_path, set_id_filter="right")
    logger.info(f"Found {len(databases)} reference databases")
    
    output_path = output_dir / f"reference_hashes_k{ksize}.bin"
    
    # Build UNION query across all databases
    # Individual DBs are already deduplicated, so no DISTINCT within each SELECT
    # We use UNION (not UNION ALL) to deduplicate across databases
    
    # CRITICAL: Use read-only in-memory connection
    conn = duckdb.connect(read_only=False)  # In-memory DB, not modifying source files
    
    # Attach each database in READ-ONLY mode
    union_parts = []
    for i, db in enumerate(databases):
        alias = f"db{i}"
        # CRITICAL: Attach in read-only mode to prevent any modifications
        conn.execute(f"ATTACH '{db.path}' AS {alias} (READ_ONLY)")
        
        if db.ksize and db.ksize != "null":
            # No DISTINCT here - individual DBs are already deduplicated
            union_parts.append(f"""
                SELECT {db.column} AS hash 
                FROM {alias}.{db.table} 
                WHERE {db.ksize} = {ksize}
            """)
        else:
            union_parts.append(f"""
                SELECT {db.column} AS hash 
                FROM {alias}.{db.table}
            """)
    
    # UNION (not UNION ALL) to deduplicate across databases
    union_query = " UNION ".join(union_parts)
    final_query = f"""
        SELECT hash 
        FROM ({union_query}) 
        ORDER BY hash
    """
    
    logger.info("Executing unified query across all databases...")
    logger.info(f"Using ksize={ksize}")
    logger.info("This may take a while for 20B+ hashes...")
    
    # Stream results to binary file
    start_time = time.time()
    result = conn.execute(final_query)
    
    buffer = np.zeros(100_000_000, dtype=np.uint64)
    buffer_pos = 0
    total_written = 0
    
    with open(output_path, 'wb') as f:
        while True:
            batch = result.fetchmany(10_000_000)
            if not batch:
                break
            
            for row in batch:
                buffer[buffer_pos] = row[0]
                buffer_pos += 1
                
                if buffer_pos >= len(buffer):
                    buffer.tofile(f)
                    total_written += buffer_pos
                    buffer_pos = 0
                    
                    if total_written % 1_000_000_000 == 0:
                        logger.info(f"Written {total_written/1e9:.1f}B hashes")
        
        # Write remaining
        if buffer_pos > 0:
            buffer[:buffer_pos].tofile(f)
            total_written += buffer_pos
    
    conn.close()
    
    total_time = time.time() - start_time
    file_size = output_path.stat().st_size
    
    logger.info(f"\n{'='*60}")
    logger.info(f"BUILD COMPLETE (duckdb method)")
    logger.info(f"{'='*60}")
    logger.info(f"Output file: {output_path}")
    logger.info(f"Total unique hashes: {total_written:,}")
    logger.info(f"File size: {file_size / 1e9:.2f} GB")
    logger.info(f"Total time: {total_time / 60:.1f} minutes")
    logger.info(f"ksize: {ksize}")
    
    # Write metadata
    metadata = {
        'n_hashes': total_written,
        'file_size_bytes': file_size,
        'ksize': ksize,
        'databases': [db.name for db in databases],
        'build_time_seconds': total_time,
        'method': 'duckdb'
    }
    
    with open(output_dir / "reference_metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)
    
    return output_path


def build_reference_hashes_inmemory(
    config_path: str,
    output_dir: str,
    ksize: int
) -> Path:
    """
    Build reference by loading all hashes into RAM, sorting, and deduplicating.
    
    This is the FASTEST method but requires enough RAM to hold all hashes.
    For 20B hashes, you need ~160GB RAM minimum, plus overhead (~200-250GB total).
    
    Individual databases are already deduplicated, so we only deduplicate the
    union at the end using np.unique().
    
    Recommended for systems with 500GB+ RAM.
    
    Args:
        config_path: Path to DB_info.json
        output_dir: Directory for output files
        ksize: k-mer size to use (REQUIRED)
        
    Returns:
        Path to the final reference hash file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    databases = load_database_configs(config_path, set_id_filter="right")
    logger.info(f"Found {len(databases)} reference databases")
    logger.info(f"Using IN-MEMORY build (fastest, requires ~250GB+ RAM)")
    logger.info(f"ksize={ksize}")
    
    # Collect all hashes into a list of arrays
    all_hashes = []
    total_count = 0
    start_time = time.time()
    
    for db in databases:
        db_start = time.time()
        db_hashes = []
        
        # Individual DBs are already deduplicated, no need for DISTINCT
        for batch in stream_hashes_from_database(db, ksize_value=ksize, batch_size=50_000_000):
            db_hashes.append(batch)
            total_count += len(batch)
        
        if db_hashes:
            combined = np.concatenate(db_hashes)
            all_hashes.append(combined)
            logger.info(f"  {db.name}: {len(combined):,} hashes in {time.time()-db_start:.1f}s")
        
        del db_hashes  # Free memory
    
    logger.info(f"Total hashes collected (with cross-DB duplicates): {total_count:,}")
    logger.info("Concatenating all arrays...")
    
    # Concatenate all
    all_combined = np.concatenate(all_hashes)
    del all_hashes  # Free memory
    
    logger.info(f"Sorting {len(all_combined):,} hashes...")
    all_combined.sort()
    
    # np.unique() deduplicates the UNION of all databases
    logger.info("Deduplicating union of all databases...")
    unique_hashes = np.unique(all_combined)
    del all_combined  # Free memory
    
    logger.info(f"Unique hashes after union deduplication: {len(unique_hashes):,}")
    
    # Write to file
    output_path = output_dir / f"reference_hashes_k{ksize}.bin"
    logger.info(f"Writing to {output_path}...")
    unique_hashes.tofile(output_path)
    
    total_time = time.time() - start_time
    file_size = output_path.stat().st_size
    
    logger.info(f"\n{'='*60}")
    logger.info(f"BUILD COMPLETE (inmemory method)")
    logger.info(f"{'='*60}")
    logger.info(f"Output file: {output_path}")
    logger.info(f"Total unique hashes: {len(unique_hashes):,}")
    logger.info(f"File size: {file_size / 1e9:.2f} GB")
    logger.info(f"Total time: {total_time / 60:.1f} minutes")
    logger.info(f"ksize: {ksize}")
    
    # Write metadata
    metadata = {
        'n_hashes': int(len(unique_hashes)),
        'file_size_bytes': file_size,
        'ksize': ksize,
        'databases': [db.name for db in databases],
        'build_time_seconds': total_time,
        'method': 'inmemory'
    }
    
    with open(output_dir / "reference_metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)
    
    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build reference hash set for Dark Matter Index calculation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Standard build with default settings
    python build_reference_hashes.py DB_info.json /scratch/dmk333_new/dmi_reference/ --ksize 31
    
    # With smaller chunks (lower memory usage)
    python build_reference_hashes.py DB_info.json output/ --ksize 31 --chunk-size 0.5
    
    # Memory-efficient mode using DuckDB's query engine
    python build_reference_hashes.py DB_info.json output/ --ksize 31 --method duckdb
        """
    )
    
    parser.add_argument("config", help="Path to DB_info.json")
    parser.add_argument("output_dir", help="Directory for output files")
    parser.add_argument(
        "--chunk-size", 
        type=float, 
        default=2.0,
        help="Chunk size in billions of hashes (default: 2.0, ~16GB per chunk)"
    )
    parser.add_argument(
        "--ksize", "-k",
        type=int,
        required=True,
        help="k-mer size (REQUIRED, e.g., 31)"
    )
    parser.add_argument(
        "--method",
        choices=["chunked", "duckdb", "inmemory"],
        default="inmemory",
        help="Build method: 'inmemory' (fastest, needs 250GB+ RAM), "
             "'chunked' (external merge sort), or 'duckdb' (in-database). "
             "Default: inmemory"
    )
    
    args = parser.parse_args()
    
    if args.method == "chunked":
        build_reference_hashes(
            args.config,
            args.output_dir,
            chunk_size_billions=args.chunk_size,
            ksize=args.ksize
        )
    elif args.method == "duckdb":
        build_reference_hashes_memory_efficient(
            args.config,
            args.output_dir,
            ksize=args.ksize
        )
    else:  # inmemory
        build_reference_hashes_inmemory(
            args.config,
            args.output_dir,
            ksize=args.ksize
        )
