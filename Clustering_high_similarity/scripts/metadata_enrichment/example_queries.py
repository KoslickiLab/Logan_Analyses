#!/usr/bin/env python3
"""
Example Queries for Component Metadata Investigation

This file contains ready-to-use pandas queries for exploring components
and their metadata patterns. Use in IPython/Jupyter or modify for batch processing.
"""

import pandas as pd
import argparse


# ============================================================================
# SETUP
# ============================================================================

def load_data(filepath='merged_components_metadata.csv'):
    """Load the merged dataset"""
    print(f"Loading {filepath}...")
    df = pd.read_csv(filepath, low_memory=False)
    print(f"Loaded {len(df):,} rows")
    return df


# ============================================================================
# FINDING INTERESTING COMPONENTS
# ============================================================================

def find_large_single_center_components(df, min_size=100, min_freq=0.95):
    """
    Find large components dominated by a single sequencing center.
    These often represent systematic technical duplications.
    """
    print(f"\nLarge components (≥{min_size}) with ≥{min_freq:.0%} from one center:")
    
    results = []
    for comp_id in df[df['component_size'] >= min_size]['component_id'].unique():
        comp_df = df[df['component_id'] == comp_id]
        center_counts = comp_df['center_name'].value_counts()
        
        if len(center_counts) > 0:
            top_center = center_counts.index[0]
            top_freq = center_counts.iloc[0] / len(comp_df)
            
            if top_freq >= min_freq:
                results.append({
                    'component_id': comp_id,
                    'size': len(comp_df),
                    'center': top_center,
                    'center_freq': top_freq
                })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        results_df = results_df.sort_values('size', ascending=False)
        print(results_df.to_string(index=False))
    else:
        print("None found with current criteria")
    
    return results_df


def find_geographically_diverse_components(df, min_size=20, min_countries=3):
    """
    Find components with samples from multiple countries.
    Unexpected geographic diversity is often interesting.
    """
    print(f"\nComponents (≥{min_size}) with samples from ≥{min_countries} countries:")
    
    results = []
    for comp_id in df[df['component_size'] >= min_size]['component_id'].unique():
        comp_df = df[df['component_id'] == comp_id]
        countries = comp_df['geo_loc_name_country_calc'].dropna()
        
        if len(countries) >= min_size * 0.8:  # Most samples have country data
            n_countries = countries.nunique()
            
            if n_countries >= min_countries:
                country_list = countries.value_counts().head(5)
                country_str = ', '.join([f"{c}({n})" for c, n in country_list.items()])
                
                results.append({
                    'component_id': comp_id,
                    'size': len(comp_df),
                    'n_countries': n_countries,
                    'countries': country_str
                })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        results_df = results_df.sort_values('n_countries', ascending=False)
        print(results_df.to_string(index=False))
    else:
        print("None found with current criteria")
    
    return results_df


def find_cross_study_components(df, min_size=10, min_studies=2):
    """
    Find components spanning multiple bioprojects/studies.
    Could indicate sample reuse or data contamination.
    """
    print(f"\nComponents (≥{min_size}) spanning ≥{min_studies} bioprojects:")
    
    results = []
    for comp_id in df[df['component_size'] >= min_size]['component_id'].unique():
        comp_df = df[df['component_id'] == comp_id]
        studies = comp_df['bioproject'].dropna()
        
        if len(studies) >= min_size * 0.8:
            n_studies = studies.nunique()
            
            if n_studies >= min_studies:
                study_list = studies.value_counts().head(3)
                study_str = ', '.join([f"{s}({n})" for s, n in study_list.items()])
                
                results.append({
                    'component_id': comp_id,
                    'size': len(comp_df),
                    'n_studies': n_studies,
                    'studies': study_str
                })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        results_df = results_df.sort_values('n_studies', ascending=False)
        print(results_df.to_string(index=False))
    else:
        print("None found with current criteria")
    
    return results_df


def find_temporally_spread_components(df, min_size=20, min_year_range=5):
    """
    Find components with samples spanning many years.
    Long temporal spread might indicate reference samples or ongoing duplications.
    """
    print(f"\nComponents (≥{min_size}) spanning ≥{min_year_range} years:")
    
    df['release_year'] = pd.to_datetime(df['releasedate'], errors='coerce').dt.year
    
    results = []
    for comp_id in df[df['component_size'] >= min_size]['component_id'].unique():
        comp_df = df[df['component_id'] == comp_id]
        years = comp_df['release_year'].dropna()
        
        if len(years) >= min_size * 0.5:
            year_range = years.max() - years.min()
            
            if year_range >= min_year_range:
                results.append({
                    'component_id': comp_id,
                    'size': len(comp_df),
                    'min_year': int(years.min()),
                    'max_year': int(years.max()),
                    'year_range': int(year_range)
                })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        results_df = results_df.sort_values('year_range', ascending=False)
        print(results_df.to_string(index=False))
    else:
        print("None found with current criteria")
    
    return results_df


def find_mixed_organism_components(df, min_size=10):
    """
    Find components with samples from multiple organisms.
    Usually suspicious - suggests mislabeling or contamination.
    """
    print(f"\nComponents (≥{min_size}) with multiple organism types:")
    
    results = []
    for comp_id in df[df['component_size'] >= min_size]['component_id'].unique():
        comp_df = df[df['component_id'] == comp_id]
        organisms = comp_df['organism'].dropna()
        
        if len(organisms) >= min_size * 0.8:
            n_organisms = organisms.nunique()
            
            if n_organisms > 1:
                org_list = organisms.value_counts()
                org_str = ', '.join([f"{o}({n})" for o, n in org_list.items()])
                
                results.append({
                    'component_id': comp_id,
                    'size': len(comp_df),
                    'n_organisms': n_organisms,
                    'organisms': org_str
                })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        results_df = results_df.sort_values('n_organisms', ascending=False)
        print(results_df.to_string(index=False))
    else:
        print("None found (this is good!)")
    
    return results_df


# ============================================================================
# INVESTIGATING SPECIFIC COMPONENTS
# ============================================================================

def investigate_component(df, component_id):
    """
    Detailed investigation of a specific component.
    Shows all key metadata and identifies patterns/outliers.
    """
    comp_df = df[df['component_id'] == component_id].copy()
    
    print(f"\n{'='*80}")
    print(f"COMPONENT {component_id} - DETAILED INVESTIGATION")
    print(f"{'='*80}")
    print(f"Size: {len(comp_df):,} samples\n")
    
    # Key metadata fields
    fields_to_check = {
        'center_name': 'Sequencing Center',
        'bioproject': 'BioProject',
        'sra_study': 'SRA Study',
        'organism': 'Organism',
        'instrument': 'Instrument',
        'geo_loc_name_country_calc': 'Country',
        'releasedate': 'Release Date'
    }
    
    for field, label in fields_to_check.items():
        if field in comp_df.columns:
            values = comp_df[field].dropna()
            
            if len(values) > 0:
                print(f"{label}:")
                n_unique = values.nunique()
                
                if n_unique == 1:
                    print(f"  ✓ Homogeneous: {values.iloc[0]}")
                else:
                    counts = values.value_counts()
                    top_freq = counts.iloc[0] / len(values)
                    
                    if top_freq > 0.9:
                        print(f"  ~ Nearly homogeneous ({top_freq:.0%}): {counts.index[0]}")
                        print(f"    Exceptions: {', '.join([str(v) for v in counts.index[1:6]])}")
                    else:
                        print(f"  ✗ Diverse ({n_unique} unique values):")
                        for val, count in counts.head(5).items():
                            print(f"    {val}: {count} ({100*count/len(values):.1f}%)")
                
                print()
    
    # Sample list
    print("Sample accessions (first 20):")
    for i, acc in enumerate(comp_df['accession'].head(20), 1):
        print(f"  {i:2d}. {acc}")
    
    if len(comp_df) > 20:
        print(f"  ... and {len(comp_df)-20} more")
    
    return comp_df


def compare_two_components(df, comp_id1, comp_id2):
    """
    Side-by-side comparison of two components.
    Useful for understanding what distinguishes different groups.
    """
    comp1 = df[df['component_id'] == comp_id1]
    comp2 = df[df['component_id'] == comp_id2]
    
    print(f"\n{'='*80}")
    print(f"COMPARING COMPONENT {comp_id1} vs COMPONENT {comp_id2}")
    print(f"{'='*80}")
    print(f"Component {comp_id1}: {len(comp1):,} samples")
    print(f"Component {comp_id2}: {len(comp2):,} samples\n")
    
    fields = ['center_name', 'bioproject', 'organism', 'geo_loc_name_country_calc', 'instrument']
    
    for field in fields:
        if field in df.columns:
            print(f"{field}:")
            
            vals1 = comp1[field].value_counts().head(3)
            vals2 = comp2[field].value_counts().head(3)
            
            print(f"  Component {comp_id1}:")
            for val, count in vals1.items():
                print(f"    {val}: {count} ({100*count/len(comp1):.1f}%)")
            
            print(f"  Component {comp_id2}:")
            for val, count in vals2.items():
                print(f"    {val}: {count} ({100*count/len(comp2):.1f}%)")
            
            print()


# ============================================================================
# BATCH QUERIES
# ============================================================================

def summary_by_field(df, field, top_n=20):
    """
    Show which values of a field appear most in components.
    Useful for understanding dominant patterns.
    """
    print(f"\nComponents by {field} (top {top_n}):")
    
    # Count samples per field value
    value_counts = df[field].value_counts().head(top_n)
    
    print(f"\nSamples per {field}:")
    for val, count in value_counts.items():
        pct = 100 * count / len(df)
        print(f"  {val:40s} {count:>8,} ({pct:>5.1f}%)")
    
    # Count components per field value
    print(f"\nComponents per {field}:")
    comp_counts = df.groupby(field)['component_id'].nunique().sort_values(ascending=False).head(top_n)
    
    for val, count in comp_counts.items():
        print(f"  {val:40s} {count:>6,} components")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main(input_file='merged_components_metadata.csv'):
    """
    Run a standard set of queries to find interesting patterns
    
    Args:
        input_file: Path to merged components metadata CSV
    """
    
    df = load_data(input_file)
    
    print("\n" + "="*80)
    print("RUNNING STANDARD QUERY SET")
    print("="*80)
    
    # Find different types of interesting components
    print("\n[1/6] Large single-center components...")
    find_large_single_center_components(df, min_size=100, min_freq=0.95)
    
    print("\n[2/6] Geographically diverse components...")
    find_geographically_diverse_components(df, min_size=20, min_countries=3)
    
    print("\n[3/6] Cross-study duplications...")
    find_cross_study_components(df, min_size=10, min_studies=2)
    
    print("\n[4/6] Temporally spread components...")
    find_temporally_spread_components(df, min_size=20, min_year_range=5)
    
    print("\n[5/6] Mixed organism components...")
    find_mixed_organism_components(df, min_size=10)
    
    print("\n[6/6] Summary by sequencing center...")
    summary_by_field(df, 'center_name', top_n=15)
    
    print("\n" + "="*80)
    print("EXAMPLE QUERIES COMPLETE")
    print("="*80)
    print("\nTo investigate specific components, use:")
    print("  investigate_component(df, component_id=XXX)")
    print("  compare_two_components(df, comp_id1=XXX, comp_id2=YYY)")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Example queries for component metadata investigation',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', default='merged_components_metadata.csv',
                       help='Input merged metadata CSV file')
    
    args = parser.parse_args()
    
    main(input_file=args.input)
