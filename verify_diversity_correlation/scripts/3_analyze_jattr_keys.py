#!/usr/bin/env python3
"""
Analyze jattr JSON column in metadata_geo_joined table.
For each key appearing in >1000 samples, extract the top 20 most frequent values.
"""

import duckdb
import os
from pathlib import Path

# Configuration
DB_PATH = "/scratch/shared_data_new/Logan_yacht_data/metadata/aws_sra_metadata/metadata_geo_joined_5M.duckdb"
OUTPUT_DIR = "/scratch/shared_data_new/Logan_yacht_data/metadata/aws_sra_metadata/jattr_analysis"
MIN_KEY_COUNT = 5000
TOP_N_VALUES = 20

# Create output directory
Path(OUTPUT_DIR).mkdir(exist_ok=True)

# Connect to database
print(f"Connecting to database: {DB_PATH}")
con = duckdb.connect(DB_PATH, read_only=True)
con.execute("SET threads TO 200")

# Define the filter condition for samples
SAMPLE_FILTER = """
    jattr IS NOT NULL 
    AND assay_type='WGS' 
    AND ((organism ILIKE '%metageno%') OR (librarysource IN ('METAGENOMIC', 'METATRANSCRIPTOMIC')))
"""

# Step 1: Get all keys with more than MIN_KEY_COUNT occurrences
print(f"\nFinding keys that appear in more than {MIN_KEY_COUNT:,} samples...")
keys_query = f"""
    SELECT
        k AS key,
        COUNT(*) AS rows_with_key
    FROM metadata_geo_joined
    CROSS JOIN UNNEST(json_keys(jattr)) AS t(k)
    WHERE {SAMPLE_FILTER}
    GROUP BY k
    HAVING COUNT(*) > {MIN_KEY_COUNT}
    ORDER BY rows_with_key DESC, key
"""

keys_result = con.execute(keys_query).fetchall()
print(f"Found {len(keys_result)} keys with more than {MIN_KEY_COUNT:,} samples")

# Step 2: For each key, get top 20 values
for idx, (key, total_count) in enumerate(keys_result, 1):
    print(f"\n[{idx}/{len(keys_result)}] Processing key: '{key}' ({total_count:,} samples)")
    
    # Query to get top 20 values for this key
    values_query = f"""
        SELECT
            json_extract_string(jattr, '$.' || '{key}') AS value,
            COUNT(*) AS count
        FROM metadata_geo_joined
        WHERE {SAMPLE_FILTER}
            AND json_extract_string(jattr, '$.' || '{key}') IS NOT NULL
        GROUP BY value
        ORDER BY count DESC
        LIMIT {TOP_N_VALUES}
    """
    
    try:
        values_result = con.execute(values_query).fetchall()
        
        # Create safe filename (replace problematic characters)
        safe_key = key.replace('/', '_').replace('\\', '_').replace(' ', '_')
        output_file = os.path.join(OUTPUT_DIR, f"jattr_key_{safe_key}.tsv")
        
        # Write results to file
        with open(output_file, 'w') as f:
            # Write header with metadata
            f.write(f"# Key: {key}\n")
            f.write(f"# Total samples with this key: {total_count:,}\n")
            f.write(f"# Top {TOP_N_VALUES} values:\n")
            f.write("#\n")
            f.write("value\tcount\tpercent\n")
            
            # Write data
            for value, count in values_result:
                percent = (count / total_count) * 100
                # Handle potential None or weird values
                value_str = str(value) if value is not None else "<NULL>"
                f.write(f"{value_str}\t{count}\t{percent:.2f}\n")
        
        print(f"  → Wrote {len(values_result)} values to {output_file}")
        
    except Exception as e:
        print(f"  ✗ Error processing key '{key}': {e}")
        continue

# Summary
print(f"\n{'='*60}")
print(f"Analysis complete!")
print(f"Results written to: {OUTPUT_DIR}")
print(f"Total keys analyzed: {len(keys_result)}")
print(f"{'='*60}")

con.close()
