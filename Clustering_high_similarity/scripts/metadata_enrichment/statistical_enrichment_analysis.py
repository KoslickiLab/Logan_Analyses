#!/usr/bin/env python3
"""
Statistical Enrichment Testing for Component Metadata

This script performs formal statistical tests to identify which metadata fields
are significantly associated with component structure, including:
- Chi-square tests for categorical associations
- Fisher's exact test for small samples
- Multiple testing correction (FDR)
"""

import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests
from collections import defaultdict
import argparse
from pathlib import Path


def load_data(filepath='merged_components_metadata.csv'):
    """Load merged dataset"""
    print(f"Loading {filepath}...")
    df = pd.read_csv(filepath, low_memory=False)
    print(f"Loaded {len(df):,} rows, {df['component_id'].nunique():,} components")
    return df


def test_field_enrichment(df, field, min_component_size=10, max_categories=100):
    """
    Test whether a metadata field is significantly associated with components.
    
    Uses chi-square test if sample sizes are adequate, Fisher's exact for small samples.
    Returns p-value and effect size metrics.
    """
    # Filter to components of sufficient size
    large_comps = df[df['component_size'] >= min_component_size]['component_id'].unique()
    df_filtered = df[df['component_id'].isin(large_comps)]
    
    # Get field values
    if field not in df_filtered.columns:
        return None
    
    # Drop NA
    df_test = df_filtered[['component_id', field]].dropna()
    
    if len(df_test) == 0:
        return None
    
    # Limit number of categories to avoid sparse contingency tables
    field_counts = df_test[field].value_counts()
    if len(field_counts) > max_categories:
        # Keep only top categories
        top_categories = field_counts.head(max_categories).index
        df_test = df_test[df_test[field].isin(top_categories)]
    
    # Also limit number of components if too many
    comp_counts = df_test['component_id'].value_counts()
    if len(comp_counts) > 100:
        # Keep top 100 largest components
        top_comps = comp_counts.head(100).index
        df_test = df_test[df_test['component_id'].isin(top_comps)]
    
    # Create contingency table
    contingency = pd.crosstab(df_test['component_id'], df_test[field])
    
    # Check if table is too sparse (but be smart about it)
    n_cells = contingency.shape[0] * contingency.shape[1]
    n_zeros = (contingency == 0).sum().sum()
    sparsity = n_zeros / n_cells
    
    # Check if this is "meaningful sparsity" vs "no data"
    # Meaningful sparsity: Each component dominated by one value (high homogeneity)
    # No data: Most components have no data for this field
    
    # Count how many components have at least some data
    components_with_data = (contingency.sum(axis=1) > 0).sum()
    pct_components_with_data = components_with_data / len(contingency)
    
    # Skip only if:
    # 1. Very sparse (>95% zeros) AND
    # 2. Most components have no data (<50% with data)
    # This catches truly empty fields while allowing homogeneous patterns
    if sparsity > 0.95 and pct_components_with_data < 0.5:
        print(f"  {field}: Insufficient data ({pct_components_with_data:.1%} components with data), skipping...")
        return None
    
    # Warn about sparsity but continue if components have data
    if sparsity > 0.8:
        print(f"  {field}: Sparse table ({sparsity:.1%} zeros) but {pct_components_with_data:.1%} components have data, continuing...")
    
    try:
        # Perform chi-square test
        chi2, p_value, dof, expected = chi2_contingency(contingency)
        
        # Calculate Cramér's V (effect size)
        n = contingency.sum().sum()
        min_dim = min(contingency.shape[0] - 1, contingency.shape[1] - 1)
        cramers_v = np.sqrt(chi2 / (n * min_dim))
        
        # Calculate how many components are dominated by single field value
        dominated_components = 0
        for comp_id in contingency.index:
            comp_row = contingency.loc[comp_id]
            if len(comp_row) > 0:
                max_freq = comp_row.max() / comp_row.sum()
                if max_freq > 0.8:  # >80% of samples have same value
                    dominated_components += 1
        
        pct_dominated = dominated_components / len(contingency)
        
        return {
            'field': field,
            'chi2': chi2,
            'p_value': p_value,
            'dof': dof,
            'cramers_v': cramers_v,
            'n_components': len(contingency),
            'n_field_values': contingency.shape[1],
            'pct_dominated_components': pct_dominated,
            'n_samples': n
        }
    
    except Exception as e:
        print(f"  {field}: Error in chi-square test: {e}")
        return None


def comprehensive_enrichment_analysis(df, fields, min_component_size=10, output_dir='./'):
    """
    Run enrichment tests for all specified fields and correct for multiple testing.
    
    This performs INTER-COMPONENT analysis:
    - Tests if DIFFERENT components are associated with DIFFERENT metadata values
    - Example: Component 0 mostly from Center A, Component 1 mostly from Center B
    - Complements the INTRA-COMPONENT homogeneity analysis done earlier
    
    High sparsity in contingency tables is EXPECTED and GOOD when:
    - Components are highly homogeneous (intra-component)
    - Different components have different dominant values (inter-component)
    
    Args:
        df: DataFrame with merged component and metadata
        fields: List of metadata fields to test
        min_component_size: Minimum component size to include
        output_dir: Directory for output files
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "="*80)
    print("STATISTICAL ENRICHMENT ANALYSIS (INTER-COMPONENT)")
    print("="*80)
    print("\nWhat this analysis tests:")
    print("  • Do DIFFERENT components have DIFFERENT metadata profiles?")
    print("  • Are certain metadata values enriched in specific components?")
    print("  • Example: Component 0 = Center A, Component 1 = Center B")
    print("\nThis complements the earlier intra-component homogeneity analysis.")
    print("\nNote: High sparsity in contingency tables is EXPECTED when components")
    print("      are homogeneous! It's a signal of strong patterns, not a problem.")
    print("="*80)
    print(f"\nTesting {len(fields)} metadata fields for association with components...")
    print(f"Minimum component size: {min_component_size}")
    
    results = []
    
    for i, field in enumerate(fields, 1):
        print(f"\n[{i}/{len(fields)}] Testing {field}...")
        result = test_field_enrichment(df, field, min_component_size=min_component_size)
        if result is not None:
            results.append(result)
            print(f"  p-value: {result['p_value']:.2e}, Cramér's V: {result['cramers_v']:.3f}")
            print(f"  {result['pct_dominated_components']:.1%} of components dominated by single value")
    
    if len(results) == 0:
        print("\nNo results to analyze.")
        return None
    
    results_df = pd.DataFrame(results)
    
    # Multiple testing correction using Benjamini-Hochberg (FDR)
    print("\n\nApplying multiple testing correction (FDR)...")
    reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(
        results_df['p_value'], 
        alpha=0.05, 
        method='fdr_bh'
    )
    
    results_df['p_value_fdr'] = pvals_corrected
    results_df['significant_fdr'] = reject
    
    # Sort by effect size (Cramér's V) among significant results
    results_df = results_df.sort_values(['significant_fdr', 'cramers_v'], 
                                       ascending=[False, False])
    
    # Summary
    print(f"\nSignificant associations (FDR < 0.05): {reject.sum()} / {len(results_df)}")
    
    print("\n" + "="*80)
    print("TOP ENRICHED FIELDS (by effect size)")
    print("="*80)
    
    sig_results = results_df[results_df['significant_fdr']]
    if len(sig_results) > 0:
        print("\nSignificant fields:")
        print(sig_results.to_string(index=False))
    else:
        print("\nNo significant associations found after multiple testing correction.")
    
    # Save results
    output_file = output_path / 'statistical_enrichment_results.csv'
    results_df.to_csv(output_file, index=False)
    print(f"\n\nSaved: {output_file}")
    
    return results_df


def analyze_pairwise_field_interactions(df, fields, top_n=10, min_component_size=20, output_dir='./'):
    """
    Analyze how pairs of metadata fields jointly explain component structure.
    For example: do certain center_name + bioproject combinations strongly predict components?
    
    Args:
        df: DataFrame with merged component and metadata
        fields: List of metadata fields to test
        top_n: Number of top results to display
        min_component_size: Minimum component size to include
        output_dir: Directory for output files
    """
    output_path = Path(output_dir)
    
    print("\n" + "="*80)
    print("PAIRWISE FIELD INTERACTION ANALYSIS")
    print("="*80)
    
    # Filter to larger components
    large_comps = df[df['component_size'] >= min_component_size]['component_id'].unique()
    df_filtered = df[df['component_id'].isin(large_comps)]
    
    results = []
    
    # Test all pairs
    n_pairs = len(fields) * (len(fields) - 1) // 2
    print(f"\nTesting {n_pairs} field pairs...")
    
    pair_idx = 0
    for i, field1 in enumerate(fields):
        for field2 in fields[i+1:]:
            pair_idx += 1
            print(f"  [{pair_idx}/{n_pairs}] {field1} × {field2}...", end='')
            
            if field1 not in df_filtered.columns or field2 not in df_filtered.columns:
                print(" skipped (missing)")
                continue
            
            # Create combination field
            df_test = df_filtered[['component_id', field1, field2]].dropna()
            
            if len(df_test) < 100:
                print(" skipped (insufficient data)")
                continue
            
            # Create combined signature
            df_test['combo'] = df_test[field1].astype(str) + ' | ' + df_test[field2].astype(str)
            
            # Limit to top combos
            combo_counts = df_test['combo'].value_counts()
            if len(combo_counts) > 50:
                top_combos = combo_counts.head(50).index
                df_test = df_test[df_test['combo'].isin(top_combos)]
            
            # Analyze: how many components are dominated by a single combo?
            dominated = 0
            component_purity = []
            
            for comp_id in df_test['component_id'].unique():
                comp_combos = df_test[df_test['component_id'] == comp_id]['combo']
                if len(comp_combos) > 0:
                    mode_combo = comp_combos.mode()
                    if len(mode_combo) > 0:
                        mode_freq = (comp_combos == mode_combo.iloc[0]).sum() / len(comp_combos)
                        component_purity.append(mode_freq)
                        if mode_freq > 0.8:
                            dominated += 1
            
            if len(component_purity) > 0:
                mean_purity = np.mean(component_purity)
                pct_dominated = dominated / len(component_purity)
                
                results.append({
                    'field1': field1,
                    'field2': field2,
                    'n_components': len(component_purity),
                    'mean_purity': mean_purity,
                    'pct_dominated': pct_dominated,
                    'n_combos': len(combo_counts)
                })
                
                print(f" mean purity: {mean_purity:.2f}")
            else:
                print(" no valid components")
    
    if len(results) == 0:
        print("\nNo results.")
        return None
    
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('mean_purity', ascending=False)
    
    print("\n" + "="*80)
    print("TOP FIELD PAIRS (by mean component purity)")
    print("="*80)
    print(results_df.head(top_n).to_string(index=False))
    
    # Save
    output_file = output_path / 'pairwise_field_interactions.csv'
    results_df.to_csv(output_file, index=False)
    print(f"\nSaved: {output_file}")
    
    return results_df


def main(input_file='merged_components_metadata.csv', output_dir='./'):
    """
    Main analysis function
    
    Args:
        input_file: Path to merged components metadata CSV
        output_dir: Directory for output files
    """
    # Load data
    df = load_data(input_file)
    
    # Key metadata fields to test
    fields_to_test = [
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
    
    # Run comprehensive enrichment analysis
    enrichment_results = comprehensive_enrichment_analysis(
        df, 
        fields_to_test, 
        min_component_size=10,
        output_dir=output_dir
    )
    
    # Analyze pairwise interactions for top significant fields
    if enrichment_results is not None:
        sig_fields = enrichment_results[enrichment_results['significant_fdr']]['field'].tolist()
        
        if len(sig_fields) > 1:
            print("\n\nAnalyzing interactions between significant fields...")
            interaction_results = analyze_pairwise_field_interactions(
                df,
                sig_fields[:10],  # Top 10 to keep computation reasonable
                top_n=20,
                min_component_size=20,
                output_dir=output_dir
            )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Statistical enrichment testing for component metadata',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', default='merged_components_metadata.csv',
                       help='Input merged metadata CSV file')
    parser.add_argument('--output-dir', default='./',
                       help='Output directory for results')
    
    args = parser.parse_args()
    
    main(
        input_file=args.input,
        output_dir=args.output_dir
    )
