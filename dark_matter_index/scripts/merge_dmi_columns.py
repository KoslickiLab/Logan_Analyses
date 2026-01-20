#!/usr/bin/env python3
"""
Merge DMI Columns into Metadata Parquet

This script merges DMI columns (dmi, total_hashes_dmi, unmapped_hashes, mapped_hashes)
from a parquet file that has DMI computed into another parquet file that has
additional metadata.

Usage:
    python merge_dmi_columns.py \
        --parquet-with-dmi /path/to/hash_diversity_data_with_dmi.parquet \
        --parquet-without-dmi /path/to/filtered_data.parquet \
        --output /path/to/filtered_data_with_dmi.parquet

Author: David Koslicki lab
Date: 2026
"""

import argparse
import logging
import pandas as pd
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Merge DMI columns from one parquet file into another",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python merge_dmi_columns.py \\
        --parquet-with-dmi hash_diversity_data_with_dmi.parquet \\
        --parquet-without-dmi filtered_data.parquet \\
        --output filtered_data_with_dmi.parquet
        """
    )
    
    parser.add_argument(
        "--parquet-with-dmi",
        required=True,
        help="Path to parquet file containing DMI columns"
    )
    parser.add_argument(
        "--parquet-without-dmi",
        required=True,
        help="Path to parquet file with metadata (missing DMI columns)"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Path for output parquet file"
    )
    
    args = parser.parse_args()
    
    logger.info("="*60)
    logger.info("MERGE DMI COLUMNS")
    logger.info("="*60)
    logger.info(f"Parquet with DMI: {args.parquet_with_dmi}")
    logger.info(f"Parquet without DMI: {args.parquet_without_dmi}")
    logger.info(f"Output: {args.output}")
    
    # Load both parquet files
    logger.info("\nLoading parquet files...")
    df_with_dmi = pd.read_parquet(args.parquet_with_dmi)
    df_without_dmi = pd.read_parquet(args.parquet_without_dmi)
    
    logger.info(f"  Parquet with DMI: {len(df_with_dmi):,} rows, {len(df_with_dmi.columns)} columns")
    logger.info(f"  Parquet without DMI: {len(df_without_dmi):,} rows, {len(df_without_dmi.columns)} columns")
    
    # Check that required columns exist
    dmi_columns = ['dmi', 'total_hashes_dmi', 'unmapped_hashes', 'mapped_hashes']
    missing_cols = [c for c in dmi_columns if c not in df_with_dmi.columns]
    if missing_cols:
        raise ValueError(f"Parquet with DMI is missing columns: {missing_cols}")
    
    if 'accession' not in df_with_dmi.columns:
        raise ValueError("Parquet with DMI is missing 'accession' column")
    if 'accession' not in df_without_dmi.columns:
        raise ValueError("Parquet without DMI is missing 'accession' column")
    
    # Check accession overlap
    logger.info("\nChecking accession overlap...")
    accessions_with_dmi = set(df_with_dmi['accession'])
    accessions_without_dmi = set(df_without_dmi['accession'])
    
    # Accessions in metadata but not in DMI file
    missing_from_dmi = accessions_without_dmi - accessions_with_dmi
    # Accessions in DMI file but not in metadata
    extra_in_dmi = accessions_with_dmi - accessions_without_dmi
    # Accessions in both
    common_accessions = accessions_with_dmi & accessions_without_dmi
    
    logger.info(f"  Accessions in parquet with DMI: {len(accessions_with_dmi):,}")
    logger.info(f"  Accessions in parquet without DMI: {len(accessions_without_dmi):,}")
    logger.info(f"  Common accessions: {len(common_accessions):,}")
    
    if missing_from_dmi:
        logger.warning(f"  ⚠️  {len(missing_from_dmi):,} accessions in metadata file are MISSING from DMI file!")
        logger.warning(f"     First 10: {list(missing_from_dmi)[:10]}")
    else:
        logger.info(f"  ✓ All accessions in metadata file have DMI values")
    
    if extra_in_dmi:
        logger.info(f"  Note: {len(extra_in_dmi):,} accessions in DMI file are not in metadata file (expected)")
    
    # Extract only the DMI columns for merging
    df_dmi_only = df_with_dmi[['accession'] + dmi_columns].copy()
    
    # Merge
    logger.info("\nMerging DMI columns...")
    df_output = df_without_dmi.merge(
        df_dmi_only,
        on='accession',
        how='left'
    )
    
    # Report merge success
    n_with_dmi = df_output['dmi'].notna().sum()
    n_total = len(df_output)
    pct_with_dmi = 100 * n_with_dmi / n_total
    
    logger.info(f"\nMerge results:")
    logger.info(f"  Total rows in output: {n_total:,}")
    logger.info(f"  Rows with DMI: {n_with_dmi:,} ({pct_with_dmi:.2f}%)")
    logger.info(f"  Rows without DMI: {n_total - n_with_dmi:,}")
    
    # DMI Statistics
    logger.info("\n" + "="*60)
    logger.info("DMI STATISTICS")
    logger.info("="*60)
    
    valid_dmi = df_output['dmi'].dropna()
    if len(valid_dmi) > 0:
        logger.info(f"\nOverall DMI distribution (n={len(valid_dmi):,}):")
        logger.info(f"  Mean:   {valid_dmi.mean():.4f}")
        logger.info(f"  Median: {valid_dmi.median():.4f}")
        logger.info(f"  Std:    {valid_dmi.std():.4f}")
        logger.info(f"  Min:    {valid_dmi.min():.4f}")
        logger.info(f"  Max:    {valid_dmi.max():.4f}")
        logger.info(f"  25th:   {valid_dmi.quantile(0.25):.4f}")
        logger.info(f"  75th:   {valid_dmi.quantile(0.75):.4f}")
        logger.info(f"  90th:   {valid_dmi.quantile(0.90):.4f}")
        
        # Distribution buckets
        logger.info("\nDMI Distribution:")
        logger.info(f"  DMI < 0.1:  {(valid_dmi < 0.1).sum():,} ({100*(valid_dmi < 0.1).mean():.1f}%)")
        logger.info(f"  0.1-0.3:    {((valid_dmi >= 0.1) & (valid_dmi < 0.3)).sum():,} ({100*((valid_dmi >= 0.1) & (valid_dmi < 0.3)).mean():.1f}%)")
        logger.info(f"  0.3-0.5:    {((valid_dmi >= 0.3) & (valid_dmi < 0.5)).sum():,} ({100*((valid_dmi >= 0.3) & (valid_dmi < 0.5)).mean():.1f}%)")
        logger.info(f"  0.5-0.7:    {((valid_dmi >= 0.5) & (valid_dmi < 0.7)).sum():,} ({100*((valid_dmi >= 0.5) & (valid_dmi < 0.7)).mean():.1f}%)")
        logger.info(f"  DMI >= 0.7: {(valid_dmi >= 0.7).sum():,} ({100*(valid_dmi >= 0.7).mean():.1f}%)")
    
    # Top organisms by DMI
    if 'organism' in df_output.columns and len(valid_dmi) > 0:
        logger.info("\n" + "-"*60)
        logger.info("TOP 10 ORGANISMS BY MEDIAN DMI (n >= 100)")
        logger.info("-"*60)
        
        organism_stats = df_output.groupby('organism').agg(
            median_dmi=('dmi', 'median'),
            mean_dmi=('dmi', 'mean'),
            std_dmi=('dmi', 'std'),
            count=('dmi', 'count')
        ).reset_index()
        
        # Filter to organisms with >= 100 samples
        organism_stats_filtered = organism_stats[organism_stats['count'] >= 100]
        organism_stats_filtered = organism_stats_filtered.sort_values('median_dmi', ascending=False)
        
        logger.info(f"\nOrganisms with >= 100 samples: {len(organism_stats_filtered)}")
        logger.info("\nTop 10 by median DMI:")
        for i, row in organism_stats_filtered.head(10).iterrows():
            logger.info(f"  {row['organism']}")
            logger.info(f"    median={row['median_dmi']:.4f}, mean={row['mean_dmi']:.4f}, "
                       f"std={row['std_dmi']:.4f}, n={int(row['count']):,}")
        
        logger.info("\n" + "-"*60)
        logger.info("BOTTOM 10 ORGANISMS BY MEDIAN DMI (n >= 100)")
        logger.info("-"*60)
        logger.info("\nBottom 10 by median DMI:")
        for i, row in organism_stats_filtered.tail(10).iterrows():
            logger.info(f"  {row['organism']}")
            logger.info(f"    median={row['median_dmi']:.4f}, mean={row['mean_dmi']:.4f}, "
                       f"std={row['std_dmi']:.4f}, n={int(row['count']):,}")
    
    # By biome if available
    if 'biome' in df_output.columns and len(valid_dmi) > 0:
        logger.info("\n" + "-"*60)
        logger.info("DMI BY BIOME (n >= 100)")
        logger.info("-"*60)
        
        biome_stats = df_output.groupby('biome').agg(
            median_dmi=('dmi', 'median'),
            mean_dmi=('dmi', 'mean'),
            count=('dmi', 'count')
        ).reset_index()
        
        biome_stats_filtered = biome_stats[biome_stats['count'] >= 100]
        biome_stats_filtered = biome_stats_filtered.sort_values('median_dmi', ascending=False)
        
        logger.info(f"\nBiomes with >= 100 samples: {len(biome_stats_filtered)}")
        for i, row in biome_stats_filtered.iterrows():
            logger.info(f"  {row['biome']}: median={row['median_dmi']:.4f}, "
                       f"mean={row['mean_dmi']:.4f}, n={int(row['count']):,}")
    
    # By platform if available
    if 'platform' in df_output.columns and len(valid_dmi) > 0:
        logger.info("\n" + "-"*60)
        logger.info("DMI BY PLATFORM")
        logger.info("-"*60)
        
        platform_stats = df_output.groupby('platform').agg(
            median_dmi=('dmi', 'median'),
            mean_dmi=('dmi', 'mean'),
            count=('dmi', 'count')
        ).reset_index()
        
        platform_stats = platform_stats.sort_values('count', ascending=False)
        
        for i, row in platform_stats.iterrows():
            logger.info(f"  {row['platform']}: median={row['median_dmi']:.4f}, "
                       f"mean={row['mean_dmi']:.4f}, n={int(row['count']):,}")
    
    # Save output
    logger.info("\n" + "="*60)
    logger.info("SAVING OUTPUT")
    logger.info("="*60)
    
    df_output.to_parquet(args.output, index=False)
    logger.info(f"Saved to: {args.output}")
    logger.info(f"Output has {len(df_output):,} rows and {len(df_output.columns)} columns")
    
    # List columns
    logger.info(f"\nColumns in output:")
    for col in df_output.columns:
        logger.info(f"  - {col}")
    
    logger.info("\n" + "="*60)
    logger.info("DONE")
    logger.info("="*60)


if __name__ == "__main__":
    main()
