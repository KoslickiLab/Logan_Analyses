#!/usr/bin/env python3
"""
Metadata Enrichment Analysis for Duplicated Sample Components

This script analyzes metadata patterns in clustered groups of duplicated samples
to identify:
1. Metadata fields that explain component structure
2. Highly homogeneous components (samples sharing key properties)
3. Outliers/surprising patterns within otherwise homogeneous groups
"""

import duckdb
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
from scipy.stats import chi2_contingency, fisher_exact
import json
from pathlib import Path
import argparse
import struct


def parse_ewkb_point(ewkb_hex):
    """
    Parse EWKB (Extended Well-Known Binary) format to extract lat/lon.
    
    Args:
        ewkb_hex: Hex string of EWKB point (e.g., '0101000020E6100000...')
    
    Returns:
        (longitude, latitude) tuple or (None, None) if parsing fails
    """
    if not ewkb_hex or pd.isna(ewkb_hex):
        return None, None
    
    try:
        # Convert hex to bytes
        ewkb_bytes = bytes.fromhex(ewkb_hex)
        
        # EWKB Point structure:
        # byte order (1 byte) + geometry type (4 bytes) + SRID (4 bytes) + X (8 bytes) + Y (8 bytes)
        # For WGS84: SRID = 4326 (0xE6100000 in little endian)
        
        # Skip byte order (1) + type (4) + SRID (4) = 9 bytes
        # Then read two doubles (8 bytes each) for X (longitude) and Y (latitude)
        if len(ewkb_bytes) >= 25:  # Minimum size for point with SRID
            # Read coordinates (little endian doubles)
            lon, lat = struct.unpack('<dd', ewkb_bytes[9:25])
            return lon, lat
    except Exception as e:
        # If parsing fails, return None
        return None, None
    
    return None, None

# Configure pandas display
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 50)


def load_and_merge_data(component_file, accession_file, db_path):
    """Load all data sources and merge into single dataframe"""
    print("Loading component membership...")
    
    # Load JSON format
    with open(component_file, 'r') as f:
        components_json = json.load(f)
    
    # Detect format: check if values are integers or strings
    first_component = next(iter(components_json.values()))
    if len(first_component) > 0:
        first_value = first_component[0]
        is_index_format = isinstance(first_value, int)
    else:
        # Empty component, assume index format for safety
        is_index_format = True
    
    if is_index_format:
        print("  Detected index-based format (integers)")
        print(f"  Loading accessions from {accession_file}...")
        with open(accession_file, 'r') as f:
            accessions = [line.strip() for line in f]
        print(f"  Loaded {len(accessions):,} accessions")
    else:
        print("  Detected accession-based format (strings)")
        accessions = None  # Not needed
    
    # Convert JSON to dataframe format
    component_data = []
    for component_name, sample_ids in components_json.items():
        # Extract component ID from name (e.g., "component_0" -> 0)
        component_id = int(component_name.split('_')[1])
        component_size = len(sample_ids)
        
        for sample_id in sample_ids:
            if is_index_format:
                # Map integer index to accession
                accession = accessions[sample_id]
                component_data.append({
                    'sample_id': sample_id,
                    'component_id': component_id,
                    'component_size': component_size,
                    'accession': accession
                })
            else:
                # Use string directly as accession
                component_data.append({
                    'sample_id': sample_id,  # String in this case
                    'component_id': component_id,
                    'component_size': component_size,
                    'accession': sample_id  # sample_id IS the accession
                })
    
    components_df = pd.DataFrame(component_data)
    print(f"  Loaded {len(components_df):,} sample-component mappings")
    print(f"  Unique components: {components_df['component_id'].nunique():,}")
    
    print("\nLoading metadata from DuckDB...")
    conn = duckdb.connect(db_path, read_only=True)
    
    # Get list of accessions we actually need
    needed_accessions = components_df['accession'].unique().tolist()
    print(f"  Need metadata for {len(needed_accessions):,} unique accessions")
    
    # Load only relevant metadata records (much more efficient!)
    print("  Loading metadata table (filtered)...")
    metadata_df = conn.execute("""
        SELECT * FROM metadata 
        WHERE acc IN (SELECT unnest(?::VARCHAR[]))
    """, [needed_accessions]).df()
    print(f"  Loaded {len(metadata_df):,} metadata records")
    
    # For location metadata, we need biosamples from the metadata we just loaded
    print("  Loading location metadata (filtered)...")
    biosamples = metadata_df['biosample'].dropna().unique().tolist()
    
    if len(biosamples) > 0:
        # Filter to lat_lon attribute and relevant biosamples
        location_df = conn.execute("""
            SELECT 
                accession,
                lat_lon,
                country,
                biome,
                elevation,
                confidence
            FROM location_metadata 
            WHERE attribute_name = 'lat_lon'
            AND accession IN (SELECT unnest(?::VARCHAR[]))
        """, [biosamples]).df()
        print(f"  Loaded {len(location_df):,} location records")
    else:
        print("  No biosamples found - skipping location metadata")
        location_df = pd.DataFrame()
    
    conn.close()
    
    # Merge everything
    print("\nMerging data...")
    # First merge components with metadata (on accession)
    merged = components_df.merge(metadata_df, left_on='accession', right_on='acc', how='left')
    
    # Then merge with location data (on biosample, not accession!)
    if len(location_df) > 0:
        # Parse EWKB lat/lon coordinates
        print("  Parsing lat/lon coordinates from EWKB format...")
        location_df[['longitude', 'latitude']] = location_df['lat_lon'].apply(
            lambda x: pd.Series(parse_ewkb_point(x))
        )
        
        parsed_count = location_df[['longitude', 'latitude']].notna().all(axis=1).sum()
        print(f"  Successfully parsed {parsed_count:,} coordinates")
        
        # Merge on biosample
        merged = merged.merge(
            location_df, 
            left_on='biosample',  # Join on biosample from metadata
            right_on='accession',  # To accession in location_metadata
            how='left', 
            suffixes=('', '_loc')
        )
    else:
        # Add empty columns if no location data
        merged['longitude'] = None
        merged['latitude'] = None
    
    print(f"  Final merged dataset: {len(merged):,} rows")
    
    return merged


def analyze_component_sizes(df):
    """Analyze distribution of component sizes"""
    print("\n" + "="*80)
    print("COMPONENT SIZE DISTRIBUTION")
    print("="*80)
    
    component_sizes = df.groupby('component_id')['component_size'].first()
    
    print(f"\nTotal components: {len(component_sizes):,}")
    print(f"Total duplicated samples: {len(df):,}")
    print(f"\nSize statistics:")
    print(component_sizes.describe())
    
    # Size distribution
    print(f"\nSize distribution:")
    size_counts = component_sizes.value_counts().sort_index()
    for size in [2, 3, 4, 5, 10, 20, 50, 100, 500, 1000, 5000]:
        if size in size_counts.index:
            print(f"  Size {size:>5}: {size_counts[size]:>6} components")
        count_ge = (component_sizes >= size).sum()
        if count_ge > 0:
            print(f"  Size ≥{size:>4}: {count_ge:>6} components")
    
    # Largest components
    print(f"\nLargest 20 components:")
    largest = component_sizes.nlargest(20)
    for comp_id, size in largest.items():
        print(f"  Component {comp_id:>6}: {size:>8,} samples")


def calculate_field_entropy(df, field, component_id):
    """Calculate entropy of a field within a component (measure of homogeneity)"""
    values = df[df['component_id'] == component_id][field].dropna()
    if len(values) == 0:
        return None, None, None
    
    # Handle list-type fields (from JSON)
    if values.dtype == object:
        # Try to handle lists
        expanded = []
        for val in values:
            if isinstance(val, list):
                expanded.extend(val)
            else:
                expanded.append(val)
        values = pd.Series(expanded)
    
    counts = values.value_counts()
    total = len(values)
    
    # Calculate entropy
    if len(counts) == 1:
        entropy = 0.0  # Perfectly homogeneous
    else:
        probs = counts / total
        entropy = -np.sum(probs * np.log2(probs))
    
    # Get most common value and its frequency
    most_common = counts.index[0]
    most_common_freq = counts.iloc[0] / total
    
    return entropy, most_common, most_common_freq


def analyze_field_enrichment(df, field, min_component_size=10, top_n=20):
    """
    Analyze how homogeneous each component is for a given metadata field.
    Returns components ranked by homogeneity (low entropy).
    """
    results = []
    
    components = df[df['component_size'] >= min_component_size]['component_id'].unique()
    
    for comp_id in components:
        comp_df = df[df['component_id'] == comp_id]
        comp_size = len(comp_df)
        
        entropy, most_common, freq = calculate_field_entropy(df, field, comp_id)
        
        if entropy is not None:
            results.append({
                'component_id': comp_id,
                'component_size': comp_size,
                'field': field,
                'entropy': entropy,
                'most_common_value': most_common,
                'frequency': freq,
                'n_unique_values': comp_df[field].nunique()
            })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # Sort by frequency (homogeneity), then component size
        results_df = results_df.sort_values(['frequency', 'component_size'], 
                                           ascending=[False, False])
    
    return results_df.head(top_n)


def find_outliers_in_component(df, component_id, key_fields):
    """
    Find outliers within a component based on key metadata fields.
    Returns samples that differ from the dominant pattern.
    """
    comp_df = df[df['component_id'] == component_id].copy()
    
    # Determine dominant value for each field
    dominant_values = {}
    for field in key_fields:
        if field in comp_df.columns:
            mode_val = comp_df[field].mode()
            if len(mode_val) > 0:
                dominant_values[field] = mode_val.iloc[0]
    
    # Find samples that differ from dominant pattern
    outliers = []
    for idx, row in comp_df.iterrows():
        differences = []
        for field, dominant_val in dominant_values.items():
            if pd.notna(row[field]) and row[field] != dominant_val:
                differences.append(f"{field}: {row[field]} (expected {dominant_val})")
        
        if differences:
            outliers.append({
                'accession': row['accession'],
                'sample_id': row['sample_id'],
                'differences': '; '.join(differences),
                'n_differences': len(differences)
            })
    
    return pd.DataFrame(outliers), dominant_values


def analyze_combination_patterns(df, fields, min_component_size=10, top_n=20):
    """
    Analyze combinations of metadata fields to find multi-field signatures
    that characterize components.
    """
    results = []
    
    components = df[df['component_size'] >= min_component_size]['component_id'].unique()
    
    for comp_id in components:
        comp_df = df[df['component_id'] == comp_id]
        comp_size = len(comp_df)
        
        # Create combination signature
        signatures = []
        for _, row in comp_df.iterrows():
            sig_parts = []
            for field in fields:
                if field in comp_df.columns and pd.notna(row[field]):
                    val = str(row[field])[:50]  # Truncate long values
                    sig_parts.append(f"{field}={val}")
            signatures.append('|'.join(sig_parts))
        
        # Count unique signatures
        sig_counts = Counter(signatures)
        most_common_sig, most_common_count = sig_counts.most_common(1)[0]
        
        results.append({
            'component_id': comp_id,
            'component_size': comp_size,
            'n_unique_signatures': len(sig_counts),
            'most_common_signature': most_common_sig,
            'most_common_count': most_common_count,
            'homogeneity': most_common_count / comp_size
        })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        results_df = results_df.sort_values(['homogeneity', 'component_size'], 
                                           ascending=[False, False])
    
    return results_df.head(top_n)


def main(component_file="components.json", 
         accession_file="accessions_mbases_geq_10.txt",
         db_path="metagenome_metadata_with_geo.duckdb",
         output_dir="./"):
    """
    Main analysis function
    
    Args:
        component_file: Path to component membership JSON file
        accession_file: Path to accessions list
        db_path: Path to metadata database
        output_dir: Directory for output files
    """
    # Create output directory if needed
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Load and merge data
    df = load_and_merge_data(component_file, accession_file, db_path)
    
    # Analyze component sizes
    analyze_component_sizes(df)
    
    # Key metadata fields to analyze
    single_value_fields = [
        'center_name',
        'bioproject', 
        'sra_study',
        'organism',
        'instrument',
        'platform',
        'librarylayout',
        'libraryselection',
        'geo_loc_name_country_calc',
        'geo_loc_name_country_continent_calc',
        'country',
        'biome',
    ]
    
    print("\n" + "="*80)
    print("SINGLE FIELD ENRICHMENT ANALYSIS")
    print("="*80)
    
    enrichment_results = {}
    for field in single_value_fields:
        if field in df.columns:
            print(f"\n--- {field} ---")
            result = analyze_field_enrichment(df, field, min_component_size=10, top_n=10)
            if len(result) > 0:
                enrichment_results[field] = result
                print(result.to_string(index=False))
            else:
                print("  No results")
    
    # Multi-field combination analysis
    print("\n" + "="*80)
    print("MULTI-FIELD COMBINATION PATTERNS")
    print("="*80)
    
    # Common combinations to check
    field_combinations = [
        ['center_name', 'bioproject'],
        ['center_name', 'bioproject', 'sra_study'],
        ['bioproject', 'organism'],
        ['center_name', 'instrument'],
        ['geo_loc_name_country_calc', 'bioproject'],
    ]
    
    for fields in field_combinations:
        print(f"\n--- Combination: {' + '.join(fields)} ---")
        result = analyze_combination_patterns(df, fields, min_component_size=10, top_n=10)
        if len(result) > 0:
            print(result.to_string(index=False))
    
    # Outlier detection for highly homogeneous components
    print("\n" + "="*80)
    print("OUTLIER DETECTION IN HOMOGENEOUS COMPONENTS")
    print("="*80)
    
    # Focus on center_name as it's often very informative
    if 'center_name' in enrichment_results:
        top_center_components = enrichment_results['center_name'].head(10)
        
        for _, row in top_center_components.iterrows():
            comp_id = row['component_id']
            comp_size = row['component_size']
            
            print(f"\n--- Component {comp_id} (size={comp_size}) ---")
            print(f"Dominant center_name: {row['most_common_value']} ({row['frequency']:.1%})")
            
            # Check for outliers across multiple fields
            key_fields = ['center_name', 'bioproject', 'geo_loc_name_country_calc', 
                         'organism', 'instrument']
            outliers_df, dominant_vals = find_outliers_in_component(df, comp_id, key_fields)
            
            print(f"Dominant pattern:")
            for field, val in dominant_vals.items():
                print(f"  {field}: {val}")
            
            if len(outliers_df) > 0:
                print(f"\nFound {len(outliers_df)} outliers:")
                print(outliers_df.to_string(index=False))
            else:
                print("\nNo outliers found - perfectly homogeneous!")
    
    # Save detailed results
    print("\n" + "="*80)
    print("SAVING DETAILED RESULTS")
    print("="*80)
    
    # Save component metadata summary
    print("\nGenerating component metadata summaries...")
    component_summaries = []
    for comp_id in df['component_id'].unique():
        comp_df = df[df['component_id'] == comp_id]
        summary = {
            'component_id': comp_id,
            'size': len(comp_df),
        }
        
        # Add summary stats for key fields
        for field in ['center_name', 'bioproject', 'geo_loc_name_country_calc', 'organism']:
            if field in comp_df.columns:
                value_counts = comp_df[field].value_counts()
                if len(value_counts) > 0:
                    summary[f'{field}_top'] = value_counts.index[0]
                    summary[f'{field}_top_freq'] = value_counts.iloc[0] / len(comp_df)
                    summary[f'{field}_n_unique'] = len(value_counts)
        
        component_summaries.append(summary)
    
    summary_df = pd.DataFrame(component_summaries)
    summary_df = summary_df.sort_values('size', ascending=False)
    output_file = output_path / 'component_metadata_summary.csv'
    summary_df.to_csv(output_file, index=False)
    print(f"  Saved: {output_file}")
    
    # Save full merged data for interactive exploration
    print("\nSaving full merged dataset...")
    output_file = output_path / 'merged_components_metadata.csv'
    df.to_csv(output_file, index=False)
    print(f"  Saved: {output_file}")
    
    print("\nAnalysis complete!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Comprehensive metadata analysis for duplicated sample components',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--components', default='components.json',
                       help='Component membership file (JSON format)')
    parser.add_argument('--accessions', default='accessions_mbases_geq_10.txt',
                       help='Accessions list file (only needed for index-based JSON format)')
    parser.add_argument('--database', default='metagenome_metadata_with_geo.duckdb',
                       help='Metadata database file')
    parser.add_argument('--output-dir', default='./',
                       help='Output directory for results')
    
    args = parser.parse_args()
    
    main(
        component_file=args.components,
        accession_file=args.accessions,
        db_path=args.database,
        output_dir=args.output_dir
    )
