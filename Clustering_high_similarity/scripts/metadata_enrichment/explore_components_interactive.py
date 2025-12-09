#!/usr/bin/env python3
"""
Interactive Component Explorer

This script provides functions to explore specific components in detail,
visualize metadata patterns, and investigate surprising relationships.
"""

import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
from math import radians, sin, cos, sqrt, atan2


def haversine_distance(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points on Earth.
    
    Args:
        lat1, lon1: Latitude and longitude of first point (in degrees)
        lat2, lon2: Latitude and longitude of second point (in degrees)
    
    Returns:
        Distance in kilometers
    """
    # Convert to radians
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    
    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    
    # Earth radius in kilometers
    R = 6371.0
    
    return R * c

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)


def load_merged_data(filepath='merged_components_metadata.csv'):
    """Load the merged dataset"""
    print(f"Loading {filepath}...")
    df = pd.read_csv(filepath, low_memory=False)
    print(f"Loaded {len(df):,} rows, {df['component_id'].nunique():,} components")
    return df


def explore_component(df, component_id):
    """
    Deep dive into a specific component.
    Shows all metadata fields and their distributions.
    """
    comp_df = df[df['component_id'] == component_id].copy()
    
    print("\n" + "="*80)
    print(f"COMPONENT {component_id} ANALYSIS")
    print("="*80)
    print(f"Size: {len(comp_df):,} samples")
    
    # Key metadata fields
    key_fields = [
        'center_name', 'bioproject', 'sra_study', 'organism',
        'instrument', 'platform', 'librarylayout', 'libraryselection',
        'geo_loc_name_country_calc', 'geo_loc_name_country_continent_calc',
        'country', 'biome', 'releasedate', 'collection_date_sam'
    ]
    
    print("\nMetadata field distributions:")
    print("-" * 80)
    
    for field in key_fields:
        if field in comp_df.columns:
            values = comp_df[field].dropna()
            if len(values) > 0:
                n_unique = values.nunique()
                print(f"\n{field}:")
                print(f"  Unique values: {n_unique}")
                
                if n_unique <= 20:
                    # Show all values if not too many
                    counts = values.value_counts()
                    for val, count in counts.items():
                        pct = 100 * count / len(values)
                        print(f"    {val}: {count} ({pct:.1f}%)")
                else:
                    # Show top 10
                    counts = values.value_counts().head(10)
                    print(f"  Top 10 values:")
                    for val, count in counts.items():
                        pct = 100 * count / len(values)
                        print(f"    {val}: {count} ({pct:.1f}%)")
    
    # Show sample accessions
    print("\n" + "-" * 80)
    print("Sample accessions (first 20):")
    for acc in comp_df['accession'].head(20):
        print(f"  {acc}")
    
    if len(comp_df) > 20:
        print(f"  ... and {len(comp_df) - 20} more")
    
    return comp_df


def compare_components(df, comp_id1, comp_id2):
    """
    Compare metadata distributions between two components.
    Useful for understanding what distinguishes different groups.
    """
    comp1 = df[df['component_id'] == comp_id1]
    comp2 = df[df['component_id'] == comp_id2]
    
    print("\n" + "="*80)
    print(f"COMPARING COMPONENTS {comp_id1} vs {comp_id2}")
    print("="*80)
    print(f"Component {comp_id1}: {len(comp1):,} samples")
    print(f"Component {comp_id2}: {len(comp2):,} samples")
    
    # Compare key fields
    comparison_fields = [
        'center_name', 'bioproject', 'organism', 'instrument',
        'geo_loc_name_country_calc'
    ]
    
    for field in comparison_fields:
        if field in df.columns:
            print(f"\n{field}:")
            
            # Top values in each component
            vals1 = comp1[field].value_counts().head(5)
            vals2 = comp2[field].value_counts().head(5)
            
            print(f"  Component {comp_id1}:")
            for val, count in vals1.items():
                pct = 100 * count / len(comp1)
                print(f"    {val}: {count} ({pct:.1f}%)")
            
            print(f"  Component {comp_id2}:")
            for val, count in vals2.items():
                pct = 100 * count / len(comp2)
                print(f"    {val}: {count} ({pct:.1f}%)")


def find_geographic_outliers(df, min_component_size=20, output_dir='./'):
    """
    Find components where most samples are from one location but some are not.
    This is one of the "surprising patterns" the user mentioned.
    
    Args:
        df: DataFrame with merged component and metadata
        min_component_size: Minimum component size to analyze
        output_dir: Directory for output files
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "="*80)
    print("GEOGRAPHIC OUTLIERS")
    print("="*80)
    print("Finding components with samples mostly from one country but with exceptions...\n")
    
    results = []
    
    components = df[df['component_size'] >= min_component_size]['component_id'].unique()
    
    for comp_id in components:
        comp_df = df[df['component_id'] == comp_id]
        
        # Look at country distribution
        if 'geo_loc_name_country_calc' in comp_df.columns:
            countries = comp_df['geo_loc_name_country_calc'].dropna()
            if len(countries) >= min_component_size * 0.8:  # At least 80% have country info
                country_counts = countries.value_counts()
                
                if len(country_counts) > 1:  # Multiple countries
                    top_country = country_counts.index[0]
                    top_count = country_counts.iloc[0]
                    top_pct = top_count / len(countries)
                    
                    # Look for cases where dominant country is >70% but not 100%
                    if 0.7 < top_pct < 0.99:
                        # Get the outlier countries
                        outlier_countries = country_counts[1:]
                        
                        # Calculate geographic distances if lat/lon available
                        max_distance_km = None
                        if 'latitude' in comp_df.columns and 'longitude' in comp_df.columns:
                            coords = comp_df[['latitude', 'longitude']].dropna()
                            if len(coords) >= 2:
                                # Calculate pairwise distances
                                distances = []
                                coords_list = coords.values.tolist()
                                for i in range(len(coords_list)):
                                    for j in range(i+1, len(coords_list)):
                                        lat1, lon1 = coords_list[i]
                                        lat2, lon2 = coords_list[j]
                                        dist = haversine_distance(lat1, lon1, lat2, lon2)
                                        distances.append(dist)
                                
                                if distances:
                                    max_distance_km = max(distances)
                        
                        results.append({
                            'component_id': comp_id,
                            'size': len(comp_df),
                            'dominant_country': top_country,
                            'dominant_pct': top_pct,
                            'n_countries': len(country_counts),
                            'outlier_countries': ', '.join([f"{c}({n})" for c, n in outlier_countries.items()]),
                            'max_distance_km': max_distance_km if max_distance_km else 'N/A'
                        })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        results_df = results_df.sort_values('size', ascending=False)
        print(results_df.to_string(index=False))
        
        # Save detailed outlier info
        print("\n\nSaving geographic outlier details...")
        output_file = output_path / 'geographic_outliers.csv'
        results_df.to_csv(output_file, index=False)
        print(f"  Saved: {output_file}")
    else:
        print("No geographic outliers found with current criteria.")
    
    return results_df


def analyze_geographic_distances(df, min_component_size=10, output_dir='./'):
    """
    Analyze geographic distances within components using precise lat/lon coordinates.
    Identifies components with large geographic spread (samples far apart).
    
    Args:
        df: DataFrame with merged component and metadata
        min_component_size: Minimum component size to analyze
        output_dir: Directory for output files
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "="*80)
    print("GEOGRAPHIC DISTANCE ANALYSIS")
    print("="*80)
    print("Finding components with large geographic spread...\n")
    
    results = []
    
    # Filter to components with lat/lon data
    df_with_coords = df[df['latitude'].notna() & df['longitude'].notna()]
    
    components = df_with_coords[df_with_coords['component_size'] >= min_component_size]['component_id'].unique()
    
    for comp_id in components:
        comp_df = df_with_coords[df_with_coords['component_id'] == comp_id]
        
        if len(comp_df) < 2:
            continue
        
        # Calculate all pairwise distances
        coords = comp_df[['latitude', 'longitude']].values
        distances = []
        
        for i in range(len(coords)):
            for j in range(i+1, len(coords)):
                lat1, lon1 = coords[i]
                lat2, lon2 = coords[j]
                dist = haversine_distance(lat1, lon1, lat2, lon2)
                distances.append(dist)
        
        if distances:
            max_dist = max(distances)
            mean_dist = np.mean(distances)
            median_dist = np.median(distances)
            
            # Get country info
            countries = comp_df['geo_loc_name_country_calc'].dropna().unique()
            n_countries = len(countries)
            country_str = ', '.join(countries[:5])  # First 5
            if n_countries > 5:
                country_str += f' (+{n_countries-5} more)'
            
            results.append({
                'component_id': comp_id,
                'size': len(comp_df),
                'n_coords': len(coords),
                'max_distance_km': max_dist,
                'mean_distance_km': mean_dist,
                'median_distance_km': median_dist,
                'n_countries': n_countries,
                'countries': country_str
            })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # Sort by max distance
        results_df = results_df.sort_values('max_distance_km', ascending=False)
        
        print(f"Analyzed {len(results_df)} components with coordinate data")
        print(f"\nTop 20 components by maximum distance:")
        print(results_df.head(20).to_string(index=False))
        
        # Save results
        output_file = output_path / 'geographic_distances.csv'
        results_df.to_csv(output_file, index=False)
        print(f"\n  Saved: {output_file}")
        
        # Summary statistics
        print(f"\nSummary statistics:")
        print(f"  Components with >1000 km spread: {(results_df['max_distance_km'] > 1000).sum()}")
        print(f"  Components with >5000 km spread: {(results_df['max_distance_km'] > 5000).sum()}")
        print(f"  Components with >10000 km spread: {(results_df['max_distance_km'] > 10000).sum()}")
    else:
        print("No components with sufficient coordinate data found.")
    
    return results_df


def find_bioproject_mixtures(df, min_component_size=20, output_dir='./'):
    """
    Find components with samples from multiple bioprojects.
    This could indicate data reuse or contamination.
    
    Args:
        df: DataFrame with merged component and metadata
        min_component_size: Minimum component size to analyze
        output_dir: Directory for output files
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "="*80)
    print("BIOPROJECT MIXTURES")
    print("="*80)
    print("Finding components spanning multiple bioprojects...\n")
    
    results = []
    
    components = df[df['component_size'] >= min_component_size]['component_id'].unique()
    
    for comp_id in components:
        comp_df = df[df['component_id'] == comp_id]
        
        if 'bioproject' in comp_df.columns:
            bioprojects = comp_df['bioproject'].dropna()
            if len(bioprojects) >= min_component_size * 0.8:
                bp_counts = bioprojects.value_counts()
                
                if len(bp_counts) > 1:  # Multiple bioprojects
                    results.append({
                        'component_id': comp_id,
                        'size': len(comp_df),
                        'n_bioprojects': len(bp_counts),
                        'bioprojects': ', '.join([f"{bp}({n})" for bp, n in bp_counts.items()]),
                        'largest_bioproject': bp_counts.index[0],
                        'largest_bp_count': bp_counts.iloc[0]
                    })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        results_df = results_df.sort_values(['n_bioprojects', 'size'], ascending=[False, False])
        print(results_df.to_string(index=False))
        
        print("\n\nSaving bioproject mixture details...")
        output_file = output_path / 'bioproject_mixtures.csv'
        results_df.to_csv(output_file, index=False)
        print(f"  Saved: {output_file}")
    else:
        print("No bioproject mixtures found.")
    
    return results_df


def plot_component_metadata_heatmap(df, component_ids, field, output_file='component_heatmap.png'):
    """
    Create a heatmap showing metadata field distribution across components.
    """
    # Create a matrix: components x field_values
    data_for_heatmap = []
    
    for comp_id in component_ids:
        comp_df = df[df['component_id'] == comp_id]
        values = comp_df[field].value_counts()
        data_for_heatmap.append(values)
    
    # Convert to dataframe
    heatmap_df = pd.DataFrame(data_for_heatmap)
    heatmap_df.index = [f"Comp_{cid}" for cid in component_ids]
    heatmap_df = heatmap_df.fillna(0)
    
    # Create heatmap
    plt.figure(figsize=(max(12, len(heatmap_df.columns) * 0.5), max(8, len(component_ids) * 0.3)))
    sns.heatmap(heatmap_df, cmap='YlOrRd', annot=True, fmt='g', cbar_kws={'label': 'Count'})
    plt.title(f'Distribution of {field} across components')
    plt.xlabel(field)
    plt.ylabel('Component')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved heatmap: {output_file}")
    plt.close()


def analyze_temporal_patterns(df, min_component_size=20, output_dir='./'):
    """
    Analyze temporal patterns in component formation.
    Do components tend to group samples from similar time periods?
    
    Args:
        df: DataFrame with merged component and metadata
        min_component_size: Minimum component size to analyze
        output_dir: Directory for output files
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "="*80)
    print("TEMPORAL CLUSTERING ANALYSIS")
    print("="*80)
    
    # Parse release dates
    df['release_year'] = pd.to_datetime(df['releasedate'], errors='coerce').dt.year
    
    results = []
    components = df[df['component_size'] >= min_component_size]['component_id'].unique()
    
    for comp_id in components:
        comp_df = df[df['component_id'] == comp_id]
        years = comp_df['release_year'].dropna()
        
        if len(years) >= min_component_size * 0.5:
            year_range = years.max() - years.min()
            year_mode = years.mode()
            year_mode_val = year_mode.iloc[0] if len(year_mode) > 0 else None
            year_mode_freq = (years == year_mode_val).sum() / len(years) if year_mode_val else 0
            
            results.append({
                'component_id': comp_id,
                'size': len(comp_df),
                'year_range': year_range,
                'min_year': years.min(),
                'max_year': years.max(),
                'mode_year': year_mode_val,
                'mode_year_freq': year_mode_freq,
                'n_years': years.nunique()
            })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # Components tightly clustered in time
        print("\nComponents tightly clustered in time (year_range <= 1):")
        tight_temporal = results_df[results_df['year_range'] <= 1].sort_values('size', ascending=False)
        print(tight_temporal.head(20).to_string(index=False))
        
        # Components spanning many years
        print("\n\nComponents spanning many years (year_range >= 5):")
        wide_temporal = results_df[results_df['year_range'] >= 5].sort_values('year_range', ascending=False)
        print(wide_temporal.head(20).to_string(index=False))
        
        # Save
        output_file = output_path / 'temporal_patterns.csv'
        results_df.to_csv(output_file, index=False)
        print(f"\n  Saved: {output_file}")
    
    return results_df


def analyze_single_component(df, component_id, output_dir='./'):
    """
    Detailed analysis of a single component
    
    Args:
        df: DataFrame with merged component and metadata
        component_id: Component ID to analyze
        output_dir: Directory for output files
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Filter to this component
    comp_df = df[df['component_id'] == component_id]
    
    if len(comp_df) == 0:
        print(f"\n❌ ERROR: Component {component_id} not found in dataset!")
        print(f"\nAvailable components: {sorted(df['component_id'].unique())[:20]}")
        if df['component_id'].nunique() > 20:
            print(f"... and {df['component_id'].nunique() - 20} more")
        return
    
    print("\n" + "="*80)
    print(f"DETAILED ANALYSIS: COMPONENT {component_id}")
    print("="*80)
    
    # Basic stats
    print(f"\nComponent Size: {len(comp_df):,} samples")
    print(f"Component ID: {component_id}")
    
    # Key metadata fields to analyze
    analysis_fields = [
        'center_name',
        'bioproject',
        'sra_study',
        'organism',
        'instrument',
        'platform',
        'librarylayout',
        'libraryselection',
        'geo_loc_name_country_calc',
        'country',
        'biome'
    ]
    
    # Show distribution for each field
    print("\n" + "-"*80)
    print("METADATA DISTRIBUTIONS")
    print("-"*80)
    
    for field in analysis_fields:
        if field in comp_df.columns:
            values = comp_df[field].dropna()
            if len(values) > 0:
                value_counts = values.value_counts()
                total = len(values)
                
                print(f"\n{field}:")
                print(f"  Total with data: {total:,} / {len(comp_df):,} ({100*total/len(comp_df):.1f}%)")
                print(f"  Unique values: {len(value_counts)}")
                
                if len(value_counts) <= 10:
                    # Show all values if <=10
                    for val, count in value_counts.items():
                        pct = 100 * count / total
                        print(f"    {val}: {count:,} ({pct:.1f}%)")
                else:
                    # Show top 10
                    print("  Top 10 values:")
                    for val, count in value_counts.head(10).items():
                        pct = 100 * count / total
                        print(f"    {val}: {count:,} ({pct:.1f}%)")
                    print(f"  ... and {len(value_counts) - 10} more values")
    
    # Geographic analysis
    print("\n" + "-"*80)
    print("GEOGRAPHIC ANALYSIS")
    print("-"*80)
    
    if 'geo_loc_name_country_calc' in comp_df.columns:
        countries = comp_df['geo_loc_name_country_calc'].dropna()
        if len(countries) > 0:
            country_counts = countries.value_counts()
            print(f"\nCountries represented: {len(country_counts)}")
            for country, count in country_counts.items():
                pct = 100 * count / len(countries)
                print(f"  {country}: {count:,} ({pct:.1f}%)")
    
    # Calculate geographic distances if coordinates available
    if 'latitude' in comp_df.columns and 'longitude' in comp_df.columns:
        coords = comp_df[['latitude', 'longitude']].dropna()
        if len(coords) >= 2:
            print(f"\nSamples with coordinates: {len(coords):,}")
            
            # Calculate pairwise distances
            coords_list = coords.values.tolist()
            distances = []
            for i in range(len(coords_list)):
                for j in range(i+1, len(coords_list)):
                    lat1, lon1 = coords_list[i]
                    lat2, lon2 = coords_list[j]
                    dist = haversine_distance(lat1, lon1, lat2, lon2)
                    distances.append(dist)
            
            if distances:
                print(f"\nGeographic spread:")
                print(f"  Max distance: {max(distances):,.1f} km")
                print(f"  Mean distance: {np.mean(distances):,.1f} km")
                print(f"  Median distance: {np.median(distances):,.1f} km")
    
    # Temporal analysis
    print("\n" + "-"*80)
    print("TEMPORAL ANALYSIS")
    print("-"*80)
    
    # Extract years from any date fields
    date_fields = ['releasedate', 'published', 'updated']
    years = []
    
    for field in date_fields:
        if field in comp_df.columns:
            dates = comp_df[field].dropna()
            for date in dates:
                try:
                    year = int(str(date)[:4])
                    if 2000 <= year <= 2030:
                        years.append(year)
                except:
                    pass
    
    if years:
        years = sorted(years)
        print(f"\nYear range: {min(years)} - {max(years)}")
        print(f"Span: {max(years) - min(years)} years")
        
        # Year distribution
        from collections import Counter
        year_counts = Counter(years)
        print("\nYear distribution:")
        for year in sorted(year_counts.keys()):
            count = year_counts[year]
            pct = 100 * count / len(years)
            print(f"  {year}: {count:,} ({pct:.1f}%)")
    
    # List all accessions
    print("\n" + "-"*80)
    print("SAMPLE ACCESSIONS")
    print("-"*80)
    
    accessions = comp_df['accession'].tolist()
    print(f"\nTotal samples: {len(accessions):,}")
    print(f"\nFirst 20 accessions:")
    for i, acc in enumerate(accessions[:20], 1):
        print(f"  {i:3d}. {acc}")
    
    if len(accessions) > 20:
        print(f"\n... and {len(accessions) - 20:,} more")
    
    # Save detailed report
    print("\n" + "-"*80)
    print("SAVING DETAILED REPORT")
    print("-"*80)
    
    report_file = output_path / f'component_{component_id}_detailed_report.txt'
    
    with open(report_file, 'w') as f:
        f.write(f"DETAILED ANALYSIS REPORT: COMPONENT {component_id}\n")
        f.write("="*80 + "\n\n")
        f.write(f"Component Size: {len(comp_df):,} samples\n")
        f.write(f"Analysis Date: {pd.Timestamp.now()}\n\n")
        
        f.write("METADATA DISTRIBUTIONS\n")
        f.write("-"*80 + "\n\n")
        
        for field in analysis_fields:
            if field in comp_df.columns:
                values = comp_df[field].dropna()
                if len(values) > 0:
                    value_counts = values.value_counts()
                    f.write(f"{field}:\n")
                    f.write(f"  Coverage: {len(values):,} / {len(comp_df):,} ({100*len(values)/len(comp_df):.1f}%)\n")
                    f.write(f"  Unique values: {len(value_counts)}\n")
                    for val, count in value_counts.items():
                        pct = 100 * count / len(values)
                        f.write(f"    {val}: {count:,} ({pct:.1f}%)\n")
                    f.write("\n")
        
        f.write("\nALL ACCESSIONS\n")
        f.write("-"*80 + "\n\n")
        for acc in accessions:
            f.write(f"{acc}\n")
    
    print(f"\n✓ Detailed report saved to: {report_file}")
    
    # Also save the component data as CSV
    csv_file = output_path / f'component_{component_id}_samples.csv'
    comp_df.to_csv(csv_file, index=False)
    print(f"✓ Sample data saved to: {csv_file}")
    
    print("\n" + "="*80)
    print(f"COMPONENT {component_id} ANALYSIS COMPLETE")
    print("="*80)


def main(input_file='merged_components_metadata.csv', output_dir='./', component_id=None):
    """
    Main analysis function
    
    Args:
        input_file: Path to merged components metadata CSV
        output_dir: Directory for output files
        component_id: Optional component ID to analyze in detail (if None, run all analyses)
    """
    # Load data
    df = load_merged_data(input_file)
    
    # If specific component requested, analyze just that one
    if component_id is not None:
        print(f"\n🔍 Focusing on Component {component_id}...")
        analyze_single_component(df, component_id, output_dir)
        return
    
    # Otherwise run all outlier analyses
    print("\nRunning outlier detection analyses...")
    
    geo_outliers = find_geographic_outliers(df, min_component_size=20, output_dir=output_dir)
    geo_distances = analyze_geographic_distances(df, min_component_size=10, output_dir=output_dir)
    bp_mixtures = find_bioproject_mixtures(df, min_component_size=20, output_dir=output_dir)
    temporal_patterns = analyze_temporal_patterns(df, min_component_size=20, output_dir=output_dir)
    
    print("\n" + "="*80)
    print("ALL ANALYSES COMPLETE")
    print("="*80)
    print("\nTo analyze a specific component in detail, run:")
    print(f"  python explore_components_interactive.py --component-id <ID>")
    print("\nExample:")
    print(f"  python explore_components_interactive.py --component-id 6603")
    print("\nOr to explore interactively in Python:")
    print("  from explore_components_interactive import load_merged_data, analyze_single_component")
    print("  df = load_merged_data('merged_components_metadata.csv')")
    print("  analyze_single_component(df, component_id=6603)")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Interactive component exploration and outlier detection',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Examples:
  # Analyze all components (default behavior)
  python explore_components_interactive.py --input merged_components_metadata.csv
  
  # Analyze specific component in detail
  python explore_components_interactive.py --input merged_components_metadata.csv --component-id 6603
  
  # Analyze component 42 with custom output directory
  python explore_components_interactive.py --component-id 42 --output-dir results/component_42/
        """
    )
    parser.add_argument('--input', default='merged_components_metadata.csv',
                       help='Input merged metadata CSV file')
    parser.add_argument('--output-dir', default='./',
                       help='Output directory for results')
    parser.add_argument('--component-id', type=int, default=None,
                       help='Specific component ID to analyze in detail (if not provided, analyzes all components)')
    
    args = parser.parse_args()
    
    main(
        input_file=args.input,
        output_dir=args.output_dir,
        component_id=args.component_id
    )
