#!/usr/bin/env python3
"""
Sensitivity Analysis for Hash-Diversity Correlation
====================================================
Tests how the correlation between hashes and alpha diversity varies
with different coverage thresholds.

This helps determine if the observed relationship is robust across
different parameter choices.

Usage:
    python hash_diversity_sensitivity.py --output-dir results_sensitivity --n-samples 5000 --n-jobs 64
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from tqdm import tqdm

# Import functions from main script
sys.path.insert(0, str(Path(__file__).parent))
from hash_diversity_correlation import (
    get_wgs_samples, process_sample_batch, Config
)

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.3)


def analyze_at_coverage(samples_df: pd.DataFrame, coverage: float, n_jobs: int) -> Dict:
    """
    Run correlation analysis at a specific coverage threshold.
    
    Returns:
        Dictionary with analysis results
    """
    # Split samples into batches
    sample_ids = samples_df['acc'].tolist()
    batch_size = max(1, len(sample_ids) // (n_jobs * 4))
    batches = [sample_ids[i:i+batch_size] for i in range(0, len(sample_ids), batch_size)]
    
    # Process batches in parallel
    all_results = []
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = {
            executor.submit(process_sample_batch, batch, coverage): i 
            for i, batch in enumerate(batches)
        }
        
        for future in as_completed(futures):
            try:
                batch_results = future.result()
                all_results.extend(batch_results)
            except Exception as e:
                continue
    
    if not all_results:
        return None
    
    # Convert to DataFrame
    results_df = pd.DataFrame(all_results)
    results_df = results_df.merge(
        samples_df[['acc', 'mbases']], 
        left_on='sample_id', 
        right_on='acc', 
        how='left'
    ).drop('acc', axis=1)
    
    # Calculate normalized metrics
    results_df['hashes_per_mb'] = results_df['num_hashes'] / results_df['mbases']
    results_df['diversity_per_mb'] = results_df['alpha_diversity'] / results_df['mbases']
    
    # Remove missing data
    df_clean = results_df.dropna(subset=['hashes_per_mb', 'diversity_per_mb'])
    
    if len(df_clean) < 10:  # Need minimum samples for meaningful correlation
        return None
    
    # Calculate correlations
    pearson_r, pearson_p = pearsonr(df_clean['hashes_per_mb'], df_clean['diversity_per_mb'])
    spearman_r, spearman_p = spearmanr(df_clean['hashes_per_mb'], df_clean['diversity_per_mb'])
    
    # Fit linear regression
    X = df_clean['hashes_per_mb'].values.reshape(-1, 1)
    y = df_clean['diversity_per_mb'].values
    
    model = LinearRegression()
    model.fit(X, y)
    y_pred = model.predict(X)
    r2 = r2_score(y, y_pred)
    
    return {
        'coverage': coverage,
        'n_samples': len(df_clean),
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p,
        'r2': r2,
        'slope': model.coef_[0],
        'intercept': model.intercept_,
        'mean_hashes_per_mb': df_clean['hashes_per_mb'].mean(),
        'mean_diversity_per_mb': df_clean['diversity_per_mb'].mean(),
        'median_hashes_per_mb': df_clean['hashes_per_mb'].median(),
        'median_diversity_per_mb': df_clean['diversity_per_mb'].median(),
    }


def run_sensitivity_analysis(config: Config, coverages: List[float]) -> pd.DataFrame:
    """
    Run correlation analysis across multiple coverage thresholds.
    """
    print("\n" + "="*70)
    print("SENSITIVITY ANALYSIS: Testing multiple coverage thresholds")
    print("="*70)
    print(f"Coverage values to test: {coverages}")
    print(f"Using {config.n_jobs} parallel workers per coverage")
    
    # Get samples once
    samples_df = get_wgs_samples(config)
    
    # Analyze each coverage threshold
    results = []
    for coverage in tqdm(coverages, desc="Testing coverages"):
        result = analyze_at_coverage(samples_df, coverage, config.n_jobs)
        if result:
            results.append(result)
            print(f"\nCoverage {coverage}: r={result['pearson_r']:.4f}, "
                  f"p={result['pearson_p']:.2e}, n={result['n_samples']:,}")
        else:
            print(f"\nCoverage {coverage}: Insufficient data")
    
    results_df = pd.DataFrame(results)
    
    # Save results
    output_file = config.output_dir / "data" / "sensitivity_results.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\nSaved sensitivity results to: {output_file}")
    
    return results_df


def create_sensitivity_plots(results_df: pd.DataFrame, config: Config):
    """
    Create visualizations showing how correlations vary with coverage.
    """
    print("\n" + "="*70)
    print("Creating sensitivity plots")
    print("="*70)
    
    # 1. Correlation coefficient vs coverage
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Pearson r
    axes[0, 0].plot(results_df['coverage'], results_df['pearson_r'], 
                    'o-', linewidth=2, markersize=8, color='steelblue')
    axes[0, 0].axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    axes[0, 0].set_xlabel('Coverage threshold', fontweight='bold')
    axes[0, 0].set_ylabel('Pearson r', fontweight='bold')
    axes[0, 0].set_title('Correlation Strength vs Coverage', fontweight='bold')
    axes[0, 0].grid(True, alpha=0.3)
    
    # R²
    axes[0, 1].plot(results_df['coverage'], results_df['r2'], 
                    'o-', linewidth=2, markersize=8, color='darkgreen')
    axes[0, 1].set_xlabel('Coverage threshold', fontweight='bold')
    axes[0, 1].set_ylabel('R²', fontweight='bold')
    axes[0, 1].set_title('Variance Explained vs Coverage', fontweight='bold')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Sample count
    axes[1, 0].plot(results_df['coverage'], results_df['n_samples'], 
                    'o-', linewidth=2, markersize=8, color='darkred')
    axes[1, 0].set_xlabel('Coverage threshold', fontweight='bold')
    axes[1, 0].set_ylabel('Number of samples', fontweight='bold')
    axes[1, 0].set_title('Sample Count vs Coverage', fontweight='bold')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].ticklabel_format(style='plain', axis='y')
    
    # P-value (log scale)
    axes[1, 1].semilogy(results_df['coverage'], results_df['pearson_p'], 
                        'o-', linewidth=2, markersize=8, color='darkorange')
    axes[1, 1].axhline(y=0.05, color='red', linestyle='--', 
                       label='p=0.05', alpha=0.7)
    axes[1, 1].axhline(y=0.01, color='red', linestyle='--', 
                       label='p=0.01', alpha=0.7)
    axes[1, 1].set_xlabel('Coverage threshold', fontweight='bold')
    axes[1, 1].set_ylabel('p-value (log scale)', fontweight='bold')
    axes[1, 1].set_title('Statistical Significance vs Coverage', fontweight='bold')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(config.output_dir / "plots" / "sensitivity_analysis.png", 
                dpi=config.dpi, bbox_inches='tight')
    plt.close()
    print("Saved sensitivity analysis plot")
    
    # 2. Mean values vs coverage
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    axes[0].plot(results_df['coverage'], results_df['mean_hashes_per_mb'], 
                 'o-', linewidth=2, markersize=8, label='Mean')
    axes[0].plot(results_df['coverage'], results_df['median_hashes_per_mb'], 
                 's-', linewidth=2, markersize=8, label='Median')
    axes[0].set_xlabel('Coverage threshold', fontweight='bold')
    axes[0].set_ylabel('Hashes per Mb', fontweight='bold')
    axes[0].set_title('Hash Density vs Coverage', fontweight='bold')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(results_df['coverage'], results_df['mean_diversity_per_mb'], 
                 'o-', linewidth=2, markersize=8, label='Mean')
    axes[1].plot(results_df['coverage'], results_df['median_diversity_per_mb'], 
                 's-', linewidth=2, markersize=8, label='Median')
    axes[1].set_xlabel('Coverage threshold', fontweight='bold')
    axes[1].set_ylabel('Diversity per Mb', fontweight='bold')
    axes[1].set_title('Alpha Diversity vs Coverage', fontweight='bold')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(config.output_dir / "plots" / "mean_values_vs_coverage.png", 
                dpi=config.dpi, bbox_inches='tight')
    plt.close()
    print("Saved mean values plot")
    
    # 3. Heatmap-style summary
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Prepare data for heatmap
    metrics = ['pearson_r', 'spearman_r', 'r2', 'slope']
    metric_labels = ['Pearson r', 'Spearman ρ', 'R²', 'Slope']
    
    heatmap_data = results_df[metrics].T
    heatmap_data.columns = [f"{c:.4f}" for c in results_df['coverage']]
    heatmap_data.index = metric_labels
    
    sns.heatmap(heatmap_data, annot=True, fmt='.4f', cmap='RdYlGn', 
                center=0, ax=ax, cbar_kws={'label': 'Value'})
    ax.set_xlabel('Coverage threshold', fontweight='bold')
    ax.set_ylabel('Metric', fontweight='bold')
    ax.set_title('Correlation Metrics Across Coverage Thresholds', 
                 fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(config.output_dir / "plots" / "metrics_heatmap.png", 
                dpi=config.dpi, bbox_inches='tight')
    plt.close()
    print("Saved metrics heatmap")


def generate_sensitivity_report(results_df: pd.DataFrame, config: Config):
    """
    Generate a report summarizing the sensitivity analysis.
    """
    report_file = config.output_dir / "reports" / "sensitivity_report.txt"
    
    with open(report_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("SENSITIVITY ANALYSIS REPORT\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Coverage thresholds tested: {len(results_df)}\n")
        f.write(f"Range: {results_df['coverage'].min()} - {results_df['coverage'].max()}\n\n")
        
        f.write("-"*70 + "\n")
        f.write("SUMMARY OF RESULTS\n")
        f.write("-"*70 + "\n\n")
        
        # Find best and worst correlations
        best_idx = results_df['pearson_r'].abs().idxmax()
        worst_idx = results_df['pearson_r'].abs().idxmin()
        
        f.write("Strongest correlation:\n")
        f.write(f"  Coverage: {results_df.loc[best_idx, 'coverage']}\n")
        f.write(f"  Pearson r: {results_df.loc[best_idx, 'pearson_r']:.6f}\n")
        f.write(f"  R²: {results_df.loc[best_idx, 'r2']:.6f}\n")
        f.write(f"  p-value: {results_df.loc[best_idx, 'pearson_p']:.4e}\n")
        f.write(f"  Samples: {results_df.loc[best_idx, 'n_samples']:,}\n\n")
        
        f.write("Weakest correlation:\n")
        f.write(f"  Coverage: {results_df.loc[worst_idx, 'coverage']}\n")
        f.write(f"  Pearson r: {results_df.loc[worst_idx, 'pearson_r']:.6f}\n")
        f.write(f"  R²: {results_df.loc[worst_idx, 'r2']:.6f}\n")
        f.write(f"  p-value: {results_df.loc[worst_idx, 'pearson_p']:.4e}\n")
        f.write(f"  Samples: {results_df.loc[worst_idx, 'n_samples']:,}\n\n")
        
        f.write("-"*70 + "\n")
        f.write("TREND ANALYSIS\n")
        f.write("-"*70 + "\n\n")
        
        # Analyze trends
        f.write("Correlation strength across coverage thresholds:\n")
        r_range = results_df['pearson_r'].max() - results_df['pearson_r'].min()
        f.write(f"  Range: {r_range:.6f}\n")
        f.write(f"  Mean: {results_df['pearson_r'].mean():.6f}\n")
        f.write(f"  Std: {results_df['pearson_r'].std():.6f}\n\n")
        
        # Count significant results
        sig_count = (results_df['pearson_p'] < 0.05).sum()
        f.write(f"Statistical significance:\n")
        f.write(f"  Results with p < 0.05: {sig_count}/{len(results_df)}\n")
        f.write(f"  Results with p < 0.01: {(results_df['pearson_p'] < 0.01).sum()}/{len(results_df)}\n")
        f.write(f"  Results with p < 0.001: {(results_df['pearson_p'] < 0.001).sum()}/{len(results_df)}\n\n")
        
        f.write("Sample count variation:\n")
        f.write(f"  Min samples: {results_df['n_samples'].min():,} (coverage={results_df.loc[results_df['n_samples'].idxmin(), 'coverage']})\n")
        f.write(f"  Max samples: {results_df['n_samples'].max():,} (coverage={results_df.loc[results_df['n_samples'].idxmax(), 'coverage']})\n\n")
        
        f.write("-"*70 + "\n")
        f.write("INTERPRETATION\n")
        f.write("-"*70 + "\n\n")
        
        # Assess robustness
        if r_range < 0.1:
            robustness = "highly robust"
        elif r_range < 0.2:
            robustness = "robust"
        elif r_range < 0.3:
            robustness = "moderately robust"
        else:
            robustness = "variable"
        
        f.write(f"The correlation between hashes and alpha diversity is {robustness}\n")
        f.write(f"across the tested coverage thresholds.\n\n")
        
        if sig_count == len(results_df):
            f.write("All tested coverage thresholds show statistically significant correlations,\n")
            f.write("indicating that the relationship is not dependent on this parameter choice.\n\n")
        elif sig_count > len(results_df) * 0.8:
            f.write("Most tested coverage thresholds show statistically significant correlations,\n")
            f.write("suggesting the relationship is generally robust.\n\n")
        else:
            f.write("The statistical significance varies with coverage threshold,\n")
            f.write("suggesting the relationship may be sensitive to this parameter.\n\n")
        
        f.write("-"*70 + "\n")
        f.write("DETAILED RESULTS TABLE\n")
        f.write("-"*70 + "\n\n")
        
        # Write full results table
        f.write(results_df.to_string(index=False))
        f.write("\n\n")
        
        f.write("="*70 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*70 + "\n")
    
    print(f"Saved sensitivity report to: {report_file}")


def main():
    """Main sensitivity analysis pipeline."""
    parser = argparse.ArgumentParser(
        description="Sensitivity analysis for hash-diversity correlation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--output-dir', '-o',
        type=str,
        required=True,
        help='Output directory for results'
    )
    parser.add_argument(
        '--n-samples', '-n',
        type=int,
        default=None,
        help='Number of samples to analyze (default: all available samples)'
    )
    parser.add_argument(
        '--coverages',
        type=str,
        default='0.015625,0.03125,0.0625,0.125,0.25,0.5,1.0',
        help='Comma-separated list of coverage thresholds to test'
    )
    parser.add_argument(
        '--n-jobs', '-j',
        type=int,
        default=64,
        help='Number of parallel workers'
    )
    parser.add_argument(
        '--min-mbases', '-m',
        type=float,
        default=100.0,
        help='Minimum megabases for sample inclusion'
    )
    parser.add_argument(
        '--random-seed',
        type=int,
        default=42,
        help='Random seed for reproducibility'
    )
    
    args = parser.parse_args()
    
    # Parse coverages
    coverages = [float(c) for c in args.coverages.split(',')]
    
    # Create configuration
    config = Config(
        output_dir=args.output_dir,
        n_samples=args.n_samples,
        coverage=coverages[0],  # Default, not used in sensitivity analysis
        n_jobs=args.n_jobs,
        min_mbases=args.min_mbases,
        random_seed=args.random_seed
    )
    
    print("\n" + "="*70)
    print("HASH-DIVERSITY CORRELATION SENSITIVITY ANALYSIS")
    print("="*70)
    print(f"Configuration:")
    print(f"  Output directory: {config.output_dir}")
    print(f"  Number of samples: {'all' if config.n_samples is None else f'{config.n_samples:,}'}")
    print(f"  Coverage thresholds: {coverages}")
    print(f"  Parallel workers: {config.n_jobs}")
    print("="*70)
    
    start_time = time.time()
    
    # Run sensitivity analysis
    results_df = run_sensitivity_analysis(config, coverages)
    
    # Create visualizations
    create_sensitivity_plots(results_df, config)
    
    # Generate report
    generate_sensitivity_report(results_df, config)
    
    elapsed = time.time() - start_time
    
    print("\n" + "="*70)
    print("SENSITIVITY ANALYSIS COMPLETE")
    print("="*70)
    print(f"Total time: {elapsed/60:.1f} minutes")
    print(f"Coverage thresholds tested: {len(results_df)}")
    print(f"\nOutputs saved to: {config.output_dir}/")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()
