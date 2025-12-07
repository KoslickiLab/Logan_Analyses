#!/usr/bin/env python3
"""
Quick Insights Summary

Provides a rapid, high-level overview of component-metadata relationships.
Run this first to get oriented before diving into detailed analyses.
"""

import pandas as pd
import numpy as np
from collections import Counter
import duckdb
import json
import argparse
from pathlib import Path


def quick_summary(component_file='components.json',
                 accession_file='accessions_mbases_geq_10.txt',
                 db_path='metagenome_metadata_with_geo.duckdb'):
    """Generate quick insights about component metadata patterns"""
    
    print("="*80)
    print("QUICK INSIGHTS SUMMARY")
    print("="*80)
    
    # Load basic data
    print("\n1. Loading data...")
    
    # Load JSON format
    with open(component_file, 'r') as f:
        components_json = json.load(f)
    
    # Detect format: check if values are integers or strings
    first_component = next(iter(components_json.values()))
    if len(first_component) > 0:
        first_value = first_component[0]
        is_index_format = isinstance(first_value, int)
    else:
        is_index_format = True
    
    if is_index_format:
        # Load accessions file for index-based format
        with open(accession_file, 'r') as f:
            accessions = [line.strip() for line in f]
    else:
        # Don't need accessions file for string-based format
        accessions = None
    
    # Convert JSON to dataframe format
    component_data = []
    for component_name, sample_ids in components_json.items():
        component_id = int(component_name.split('_')[1])
        component_size = len(sample_ids)
        
        for sample_id in sample_ids:
            if is_index_format:
                accession = accessions[sample_id]
                component_data.append({
                    'sample_id': sample_id,
                    'component_id': component_id,
                    'component_size': component_size,
                    'accession': accession
                })
            else:
                component_data.append({
                    'sample_id': sample_id,
                    'component_id': component_id,
                    'component_size': component_size,
                    'accession': sample_id  # sample_id IS the accession
                })
    
    components_df = pd.DataFrame(component_data)
    
    conn = duckdb.connect(db_path, read_only=True)
    
    # Get list of accessions we need
    needed_accessions = components_df['accession'].unique().tolist()
    
    # Load only relevant metadata (much more efficient!)
    metadata_df = conn.execute("""
        SELECT * FROM metadata 
        WHERE acc IN (SELECT unnest(?::VARCHAR[]))
    """, [needed_accessions]).df()
    
    # Load location metadata for relevant biosamples
    biosamples = metadata_df['biosample'].dropna().unique().tolist()
    
    if len(biosamples) > 0:
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
    else:
        location_df = pd.DataFrame()
    
    conn.close()
    
    # Merge metadata
    merged = components_df.merge(metadata_df, left_on='accession', right_on='acc', how='left')
    
    # Merge location metadata (on biosample)
    if len(location_df) > 0:
        merged = merged.merge(location_df, left_on='biosample', right_on='accession', 
                             how='left', suffixes=('', '_loc'))
    
    print(f"   Total duplicated samples: {len(merged):,}")
    print(f"   Total components: {components_df['component_id'].nunique():,}")
    
    # Component size overview
    print("\n2. Component sizes:")
    sizes = components_df.groupby('component_id')['component_size'].first()
    print(f"   Median size: {sizes.median():.0f}")
    print(f"   Mean size: {sizes.mean():.1f}")
    print(f"   Largest component: {sizes.max():,} samples")
    print(f"   Components with ≥100 samples: {(sizes >= 100).sum():,}")
    print(f"   Components with ≥1000 samples: {(sizes >= 1000).sum():,}")
    
    # Quick metadata scan - which fields are most homogeneous?
    print("\n3. Most predictive metadata fields (by homogeneity):")
    
    test_fields = ['center_name', 'bioproject', 'organism', 'instrument', 
                   'geo_loc_name_country_calc']
    
    field_scores = []
    
    for field in test_fields:
        if field not in merged.columns:
            continue
        
        # For components ≥20 samples, calculate average homogeneity
        homogeneities = []
        
        for comp_id in merged[merged['component_size'] >= 20]['component_id'].unique():
            comp_values = merged[merged['component_id'] == comp_id][field].dropna()
            if len(comp_values) >= 10:  # Need enough data
                # Homogeneity = fraction belonging to most common value
                mode_freq = comp_values.value_counts().iloc[0] / len(comp_values)
                homogeneities.append(mode_freq)
        
        if len(homogeneities) > 0:
            avg_homogeneity = np.mean(homogeneities)
            field_scores.append((field, avg_homogeneity))
    
    field_scores.sort(key=lambda x: x[1], reverse=True)
    
    for field, score in field_scores:
        stars = '★' * int(score * 10)
        print(f"   {field:30s} {score:.1%} {stars}")
    
    # Geographic diversity
    print("\n4. Geographic diversity:")
    country_dist = merged['geo_loc_name_country_calc'].value_counts().head(10)
    print(f"   Samples have data from {merged['geo_loc_name_country_calc'].nunique()} countries")
    print(f"   Top countries:")
    for country, count in country_dist.items():
        pct = 100 * count / len(merged)
        print(f"      {country:30s} {count:>8,} ({pct:>5.1f}%)")
    
    # Bioproject diversity
    print("\n5. Study (bioproject) diversity:")
    n_bioprojects = merged['bioproject'].nunique()
    print(f"   Samples span {n_bioprojects:,} bioprojects")
    
    # How many components span multiple bioprojects?
    multi_bp_components = 0
    for comp_id in merged['component_id'].unique():
        comp_bp = merged[merged['component_id'] == comp_id]['bioproject'].dropna()
        if len(comp_bp) > 0 and comp_bp.nunique() > 1:
            multi_bp_components += 1
    
    print(f"   Components with multiple bioprojects: {multi_bp_components:,} "
          f"({100*multi_bp_components/merged['component_id'].nunique():.1f}%)")
    
    # Sequencing center concentration
    print("\n6. Sequencing center concentration:")
    center_dist = merged['center_name'].value_counts().head(10)
    print(f"   Samples from {merged['center_name'].nunique()} centers")
    print(f"   Top centers:")
    for center, count in center_dist.items():
        pct = 100 * count / len(merged)
        print(f"      {center:30s} {count:>8,} ({pct:>5.1f}%)")
    
    # Temporal distribution
    print("\n7. Temporal distribution:")
    merged['release_year'] = pd.to_datetime(merged['releasedate'], errors='coerce').dt.year
    year_dist = merged['release_year'].value_counts().sort_index()
    
    if len(year_dist) > 0:
        print(f"   Data spans {year_dist.index.min():.0f} to {year_dist.index.max():.0f}")
        print(f"   Peak year: {year_dist.idxmax():.0f} ({year_dist.max():,} samples)")
    
    # Red flags / interesting patterns
    print("\n8. Preliminary red flags:")
    
    # Very large homogeneous components
    for comp_id in merged[merged['component_size'] >= 1000]['component_id'].unique()[:5]:
        comp_df = merged[merged['component_id'] == comp_id]
        center_mode = comp_df['center_name'].mode()
        bp_mode = comp_df['bioproject'].mode()
        
        if len(center_mode) > 0 and len(bp_mode) > 0:
            center_freq = (comp_df['center_name'] == center_mode.iloc[0]).sum() / len(comp_df)
            
            if center_freq > 0.95:
                print(f"   ⚠ Component {comp_id}: {len(comp_df):,} samples, "
                      f"{center_freq:.0%} from {center_mode.iloc[0]}")
    
    # Geographic mixing
    geo_mixed = 0
    for comp_id in merged[merged['component_size'] >= 20]['component_id'].unique()[:100]:
        comp_countries = merged[merged['component_id'] == comp_id]['geo_loc_name_country_calc'].dropna()
        if len(comp_countries) >= 10 and comp_countries.nunique() >= 3:
            geo_mixed += 1
    
    if geo_mixed > 0:
        print(f"   ⚠ Found {geo_mixed} large components with samples from ≥3 countries")
    
    # Recommendations
    print("\n9. Recommended next steps:")
    
    top_field = field_scores[0][0] if field_scores else 'center_name'
    
    print(f"   → Run full analysis: python analyze_component_metadata.py")
    print(f"   → Most predictive field appears to be: {top_field}")
    print(f"   → Check for cross-bioproject duplications: {multi_bp_components} found")
    
    if geo_mixed > 0:
        print(f"   → Investigate geographic mixing: {geo_mixed} components with diverse locations")
    
    print("\n" + "="*80)
    print("See ANALYSIS_GUIDE.md for detailed interpretation guidance")
    print("="*80)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Quick insights summary for component metadata analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--components', default='components.json',
                       help='Component membership file (JSON format)')
    parser.add_argument('--accessions', default='accessions_mbases_geq_10.txt',
                       help='Accessions list file (only needed for index-based JSON format)')
    parser.add_argument('--database', default='metagenome_metadata_with_geo.duckdb',
                       help='Metadata database file')
    
    args = parser.parse_args()
    
    quick_summary(
        component_file=args.components,
        accession_file=args.accessions,
        db_path=args.database
    )
