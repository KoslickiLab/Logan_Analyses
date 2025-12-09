#!/usr/bin/env python3
"""
Hash-Diversity Correlation Analysis
====================================
Tests the hypothesis that distinct hashes per basepair correlates with alpha diversity
in WGS metagenomic samples.

This script:
1. Queries WGS metagenomic samples from metadata database
2. Extracts hash counts and taxonomic data from YACHT database
3. Calculates normalized metrics (per million bases)
4. Performs correlation analysis
5. Generates publication-quality visualizations

Usage:
    python hash_diversity_correlation.py --output-dir results --n-samples 10000 --coverage 0.0625 --n-jobs 128
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
import warnings

import duckdb
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# Set publication-quality style
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.5)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

# Database paths
YACHT_DB = "/scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db"
METADATA_DB = "/scratch/shared_data_new/Logan_yacht_data/metadata/aws_sra_metadata/metadata_geo_joined.duckdb"


@dataclass
class Config:
    """Configuration for the analysis."""
    output_dir: Path
    n_samples: Optional[int] = None  # None = all samples
    coverage: float = 0.0625
    n_jobs: int = 64
    random_seed: int = 42
    min_mbases: float = 100.0
    normalization_factor: float = 1_000_000.0  # Per million bases
    figure_size: Tuple[float, float] = (10, 8)
    dpi: int = 300
    
    def __post_init__(self):
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "plots").mkdir(exist_ok=True)
        (self.output_dir / "data").mkdir(exist_ok=True)
        (self.output_dir / "reports").mkdir(exist_ok=True)


def get_wgs_samples(config: Config) -> pd.DataFrame:
    """
    Query metadata database to get WGS metagenomic sample accessions.
    
    Returns:
        DataFrame with columns: acc, mbases, reads
    """
    print("\n" + "="*70)
    print("STEP 1: Querying WGS metagenomic samples")
    print("="*70)
    
    query = f"""
    SELECT 
        acc,
        mbases
    FROM metadata_geo_joined 
    WHERE assay_type = 'WGS' 
        AND libraryselection = 'RANDOM' 
        AND mbases > {config.min_mbases}
    ORDER BY acc
    """
    
    print(f"Connecting to {METADATA_DB}...")
    conn = duckdb.connect(METADATA_DB, read_only=True, config={'threads': 1})
    
    try:
        print("Executing query...")
        df = conn.execute(query).fetchdf()
        print(f"Found {len(df):,} WGS metagenomic samples")
        
        if config.n_samples and len(df) > config.n_samples:
            print(f"Randomly sampling {config.n_samples:,} samples (seed={config.random_seed})")
            df = df.sample(n=config.n_samples, random_state=config.random_seed)
            df = df.sort_values('acc').reset_index(drop=True)
        
        print(f"\nSample statistics:")
        print(f"  Mbases range: {df['mbases'].min():.1f} - {df['mbases'].max():.1f}")
        print(f"  Mbases mean: {df['mbases'].mean():.1f}")
        print(f"  Mbases median: {df['mbases'].median():.1f}")
        
        # Save sample list
        sample_file = config.output_dir / "data" / "selected_samples.csv"
        df.to_csv(sample_file, index=False)
        print(f"\nSaved sample list to: {sample_file}")
        
        return df
        
    finally:
        conn.close()


def process_sample_batch(sample_ids: List[str], coverage: float) -> List[Dict]:
    """
    Process a batch of samples to extract hash counts and alpha diversity.
    
    This function is designed to be run in parallel.
    
    Returns:
        List of dicts with keys: sample_id, num_hashes, alpha_diversity, mbases
    """
    conn = duckdb.connect(YACHT_DB, read_only=True, config={'threads': 1})
    results = []
    
    try:
        for sample_id in sample_ids:
            try:
                # Query taxa profiles for this sample at specified coverage
                query = f"""
                SELECT 
                    sample_id,
                    organism_id,
                    tax_id,
                    num_total_kmers_in_sample_sketch
                FROM taxa_profiles.profiles 
                WHERE sample_id = '{sample_id}' 
                    AND min_coverage = {coverage}
                    AND in_sample_est = true
                """
                
                df = conn.execute(query).fetchdf()
                
                if df.empty:
                    # No data for this sample at this coverage
                    continue
                
                # Extract number of distinct hashes (should be same for all rows)
                num_hashes = df['num_total_kmers_in_sample_sketch'].iloc[0]
                
                # Calculate alpha diversity as number of unique taxa
                # Use organism_id as primary identifier, fall back to tax_id if needed
                unique_organisms = df['organism_id'].nunique()
                unique_taxa = df['tax_id'].nunique()
                alpha_diversity = max(unique_organisms, unique_taxa)
                
                results.append({
                    'sample_id': sample_id,
                    'num_hashes': num_hashes,
                    'alpha_diversity': alpha_diversity,
                    'num_records': len(df)
                })
                
            except Exception as e:
                # Skip samples that cause errors
                print(f"Warning: Error processing {sample_id}: {e}", file=sys.stderr)
                continue
                
    finally:
        conn.close()
    
    return results


def extract_hash_and_diversity_data(samples_df: pd.DataFrame, config: Config) -> pd.DataFrame:
    """
    Extract hash counts and alpha diversity for all samples using parallel processing.
    
    Args:
        samples_df: DataFrame with sample accessions and mbases
        config: Configuration object
        
    Returns:
        DataFrame with columns: sample_id, num_hashes, alpha_diversity, mbases
    """
    print("\n" + "="*70)
    print("STEP 2: Extracting hash counts and alpha diversity")
    print("="*70)
    print(f"Coverage threshold: {config.coverage}")
    print(f"Using {config.n_jobs} parallel workers")
    
    # Split samples into batches for parallel processing
    sample_ids = samples_df['acc'].tolist()
    batch_size = max(1, len(sample_ids) // (config.n_jobs * 4))  # 4 batches per worker
    batches = [sample_ids[i:i+batch_size] for i in range(0, len(sample_ids), batch_size)]
    
    print(f"Processing {len(sample_ids):,} samples in {len(batches)} batches")
    print(f"Batch size: ~{batch_size} samples")
    
    # Process batches in parallel
    all_results = []
    start_time = time.time()
    
    with ProcessPoolExecutor(max_workers=config.n_jobs) as executor:
        # Submit all jobs
        futures = {
            executor.submit(process_sample_batch, batch, config.coverage): i 
            for i, batch in enumerate(batches)
        }
        
        # Collect results with progress bar
        with tqdm(total=len(batches), desc="Processing batches") as pbar:
            for future in as_completed(futures):
                try:
                    batch_results = future.result()
                    all_results.extend(batch_results)
                except Exception as e:
                    print(f"\nError in batch: {e}", file=sys.stderr)
                pbar.update(1)
    
    elapsed = time.time() - start_time
    print(f"\nProcessing complete in {elapsed:.1f}s ({len(sample_ids)/elapsed:.1f} samples/s)")
    
    # Convert to DataFrame
    results_df = pd.DataFrame(all_results)
    
    if results_df.empty:
        print("ERROR: No data extracted! Check that samples exist in YACHT database.", 
              file=sys.stderr)
        sys.exit(1)
    
    print(f"\nSuccessfully extracted data for {len(results_df):,} samples")
    print(f"Samples with no data: {len(sample_ids) - len(results_df):,}")
    
    # Merge with mbases from original samples
    results_df = results_df.merge(
        samples_df[['acc', 'mbases']], 
        left_on='sample_id', 
        right_on='acc', 
        how='left'
    ).drop('acc', axis=1)
    
    # Calculate normalized metrics
    results_df['hashes_per_mb'] = results_df['num_hashes'] / results_df['mbases']
    results_df['diversity_per_mb'] = results_df['alpha_diversity'] / results_df['mbases']
    
    print(f"\nData summary:")
    print(f"  Hashes range: {results_df['num_hashes'].min():,.0f} - {results_df['num_hashes'].max():,.0f}")
    print(f"  Alpha diversity range: {results_df['alpha_diversity'].min():.0f} - {results_df['alpha_diversity'].max():.0f}")
    print(f"  Hashes per Mb range: {results_df['hashes_per_mb'].min():.2f} - {results_df['hashes_per_mb'].max():.2f}")
    print(f"  Diversity per Mb range: {results_df['diversity_per_mb'].min():.4f} - {results_df['diversity_per_mb'].max():.4f}")
    
    # Save raw data
    data_file = config.output_dir / "data" / "hash_diversity_data.csv"
    results_df.to_csv(data_file, index=False)
    print(f"\nSaved raw data to: {data_file}")
    
    return results_df


def perform_correlation_analysis(df: pd.DataFrame, config: Config) -> Dict:
    """
    Perform correlation analysis between hashes and alpha diversity.
    
    Args:
        df: DataFrame with hash and diversity data
        config: Configuration object
        
    Returns:
        Dictionary of analysis results
    """
    print("\n" + "="*70)
    print("STEP 3: Correlation analysis")
    print("="*70)
    
    # Remove any rows with missing data
    df_clean = df.dropna(subset=['hashes_per_mb', 'diversity_per_mb'])
    n_samples = len(df_clean)
    
    print(f"Analyzing {n_samples:,} samples with complete data")
    
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
    
    slope = model.coef_[0]
    intercept = model.intercept_
    
    results = {
        'n_samples': n_samples,
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p,
        'r2': r2,
        'slope': slope,
        'intercept': intercept,
    }
    
    print(f"\nCorrelation Results:")
    print(f"  Pearson r:  {pearson_r:.4f} (p={pearson_p:.2e})")
    print(f"  Spearman ρ: {spearman_r:.4f} (p={spearman_p:.2e})")
    print(f"  R²:         {r2:.4f}")
    print(f"  Regression: y = {slope:.6f}x + {intercept:.6f}")
    
    # Interpret results
    print(f"\nInterpretation:")
    if pearson_p < 0.001:
        sig_str = "highly significant (p < 0.001)"
    elif pearson_p < 0.01:
        sig_str = "significant (p < 0.01)"
    elif pearson_p < 0.05:
        sig_str = "significant (p < 0.05)"
    else:
        sig_str = "not significant (p >= 0.05)"
    
    if abs(pearson_r) > 0.7:
        strength = "strong"
    elif abs(pearson_r) > 0.4:
        strength = "moderate"
    elif abs(pearson_r) > 0.2:
        strength = "weak"
    else:
        strength = "very weak"
    
    direction = "positive" if pearson_r > 0 else "negative"
    
    print(f"  The correlation is {strength} and {direction}, and {sig_str}")
    print(f"  {r2*100:.1f}% of variance in diversity per Mb is explained by hashes per Mb")
    
    return results


def create_visualizations(df: pd.DataFrame, results: Dict, config: Config):
    """
    Create publication-quality visualizations.
    
    Args:
        df: DataFrame with hash and diversity data
        results: Correlation analysis results
        config: Configuration object
    """
    print("\n" + "="*70)
    print("STEP 4: Creating visualizations")
    print("="*70)
    
    df_clean = df.dropna(subset=['hashes_per_mb', 'diversity_per_mb'])
    
    # Main scatter plot with regression line
    fig, ax = plt.subplots(figsize=config.figure_size)
    
    # Scatter plot with transparency for dense regions
    ax.scatter(
        df_clean['hashes_per_mb'], 
        df_clean['diversity_per_mb'],
        alpha=0.3,
        s=20,
        c='steelblue',
        edgecolors='none',
        rasterized=True  # For better PDF rendering with many points
    )
    
    # Add regression line
    X = df_clean['hashes_per_mb'].values
    y_pred = results['slope'] * X + results['intercept']
    ax.plot(X, y_pred, 'r-', linewidth=2, label='Linear fit', zorder=10)
    
    # Labels and title
    ax.set_xlabel('Distinct hashes per megabase', fontsize=14, fontweight='bold')
    ax.set_ylabel('Alpha diversity per megabase', fontsize=14, fontweight='bold')
    ax.set_title('Correlation between Hash Density and Alpha Diversity\n' + 
                 f'in WGS Metagenomic Samples (coverage ≥ {config.coverage})',
                 fontsize=16, fontweight='bold', pad=20)
    
    # Add statistics box
    stats_text = (
        f'n = {results["n_samples"]:,}\n'
        f'Pearson r = {results["pearson_r"]:.4f}\n'
        f'p-value = {results["pearson_p"]:.2e}\n'
        f'R² = {results["r2"]:.4f}\n'
        f'Spearman ρ = {results["spearman_r"]:.4f}'
    )
    
    # Position stats box in upper left
    ax.text(
        0.05, 0.95, stats_text,
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    )
    
    # Add legend for regression line
    ax.legend(loc='lower right', fontsize=12)
    
    # Grid
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    plot_file = config.output_dir / "plots" / "hash_diversity_correlation.png"
    plt.savefig(plot_file, dpi=config.dpi, bbox_inches='tight')
    print(f"Saved main plot to: {plot_file}")
    plt.close()
    
    # Create additional diagnostic plots
    create_diagnostic_plots(df_clean, config)


def create_diagnostic_plots(df: pd.DataFrame, config: Config):
    """
    Create additional diagnostic plots.
    """
    # 1. Hexbin plot for dense data
    fig, ax = plt.subplots(figsize=config.figure_size)
    hexbin = ax.hexbin(
        df['hashes_per_mb'],
        df['diversity_per_mb'],
        gridsize=50,
        cmap='YlOrRd',
        mincnt=1,
        bins='log'
    )
    ax.set_xlabel('Distinct hashes per megabase', fontsize=14, fontweight='bold')
    ax.set_ylabel('Alpha diversity per megabase', fontsize=14, fontweight='bold')
    ax.set_title('Hash-Diversity Correlation (Hexbin Density)', fontsize=16, fontweight='bold')
    plt.colorbar(hexbin, ax=ax, label='log10(count)')
    plt.tight_layout()
    plt.savefig(config.output_dir / "plots" / "hash_diversity_hexbin.png", 
                dpi=config.dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved hexbin plot")
    
    # 2. Distribution plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Hashes per Mb distribution
    axes[0, 0].hist(df['hashes_per_mb'], bins=50, edgecolor='black', alpha=0.7)
    axes[0, 0].set_xlabel('Hashes per Mb')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Distribution of Hashes per Mb')
    axes[0, 0].axvline(df['hashes_per_mb'].median(), color='red', linestyle='--', 
                       label=f'Median: {df["hashes_per_mb"].median():.2f}')
    axes[0, 0].legend()
    
    # Diversity per Mb distribution
    axes[0, 1].hist(df['diversity_per_mb'], bins=50, edgecolor='black', alpha=0.7)
    axes[0, 1].set_xlabel('Diversity per Mb')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Distribution of Diversity per Mb')
    axes[0, 1].axvline(df['diversity_per_mb'].median(), color='red', linestyle='--',
                       label=f'Median: {df["diversity_per_mb"].median():.4f}')
    axes[0, 1].legend()
    
    # Raw hashes distribution
    axes[1, 0].hist(df['num_hashes'], bins=50, edgecolor='black', alpha=0.7)
    axes[1, 0].set_xlabel('Total distinct hashes')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Distribution of Total Hashes')
    axes[1, 0].ticklabel_format(style='scientific', axis='x', scilimits=(0,0))
    
    # Raw diversity distribution
    axes[1, 1].hist(df['alpha_diversity'], bins=50, edgecolor='black', alpha=0.7)
    axes[1, 1].set_xlabel('Alpha diversity (species count)')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].set_title('Distribution of Alpha Diversity')
    
    plt.tight_layout()
    plt.savefig(config.output_dir / "plots" / "distributions.png", 
                dpi=config.dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved distribution plots")
    
    # 3. Residual plot
    fig, ax = plt.subplots(figsize=config.figure_size)
    X = df['hashes_per_mb'].values
    y = df['diversity_per_mb'].values
    
    # Calculate residuals
    from sklearn.linear_model import LinearRegression
    model = LinearRegression()
    model.fit(X.reshape(-1, 1), y)
    y_pred = model.predict(X.reshape(-1, 1))
    residuals = y - y_pred
    
    ax.scatter(y_pred, residuals, alpha=0.3, s=20, c='steelblue', edgecolors='none')
    ax.axhline(y=0, color='red', linestyle='--', linewidth=2)
    ax.set_xlabel('Fitted values', fontsize=14, fontweight='bold')
    ax.set_ylabel('Residuals', fontsize=14, fontweight='bold')
    ax.set_title('Residual Plot', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(config.output_dir / "plots" / "residuals.png", 
                dpi=config.dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved residual plot")


def generate_report(df: pd.DataFrame, results: Dict, config: Config):
    """
    Generate a comprehensive text report of the analysis.
    """
    print("\n" + "="*70)
    print("STEP 5: Generating report")
    print("="*70)
    
    report_file = config.output_dir / "reports" / "analysis_report.txt"
    
    with open(report_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("HASH-DIVERSITY CORRELATION ANALYSIS REPORT\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Coverage Threshold: {config.coverage}\n")
        f.write(f"Minimum Mbases: {config.min_mbases}\n\n")
        
        f.write("-"*70 + "\n")
        f.write("DATA SUMMARY\n")
        f.write("-"*70 + "\n")
        f.write(f"Total samples analyzed: {len(df):,}\n")
        f.write(f"Samples with complete data: {results['n_samples']:,}\n\n")
        
        f.write("Raw Metrics:\n")
        f.write(f"  Distinct hashes:\n")
        f.write(f"    Range: {df['num_hashes'].min():,.0f} - {df['num_hashes'].max():,.0f}\n")
        f.write(f"    Mean: {df['num_hashes'].mean():,.0f}\n")
        f.write(f"    Median: {df['num_hashes'].median():,.0f}\n\n")
        
        f.write(f"  Alpha diversity (species count):\n")
        f.write(f"    Range: {df['alpha_diversity'].min():.0f} - {df['alpha_diversity'].max():.0f}\n")
        f.write(f"    Mean: {df['alpha_diversity'].mean():.2f}\n")
        f.write(f"    Median: {df['alpha_diversity'].median():.0f}\n\n")
        
        f.write(f"  Sequencing depth (Mbases):\n")
        f.write(f"    Range: {df['mbases'].min():.1f} - {df['mbases'].max():.1f}\n")
        f.write(f"    Mean: {df['mbases'].mean():.1f}\n")
        f.write(f"    Median: {df['mbases'].median():.1f}\n\n")
        
        f.write("Normalized Metrics (per Mb):\n")
        f.write(f"  Hashes per Mb:\n")
        f.write(f"    Range: {df['hashes_per_mb'].min():.2f} - {df['hashes_per_mb'].max():.2f}\n")
        f.write(f"    Mean: {df['hashes_per_mb'].mean():.2f}\n")
        f.write(f"    Median: {df['hashes_per_mb'].median():.2f}\n\n")
        
        f.write(f"  Diversity per Mb:\n")
        f.write(f"    Range: {df['diversity_per_mb'].min():.6f} - {df['diversity_per_mb'].max():.6f}\n")
        f.write(f"    Mean: {df['diversity_per_mb'].mean():.6f}\n")
        f.write(f"    Median: {df['diversity_per_mb'].median():.6f}\n\n")
        
        f.write("-"*70 + "\n")
        f.write("CORRELATION ANALYSIS\n")
        f.write("-"*70 + "\n")
        f.write(f"Pearson correlation coefficient: r = {results['pearson_r']:.6f}\n")
        f.write(f"  p-value: {results['pearson_p']:.4e}\n")
        f.write(f"  Significance: {'***' if results['pearson_p'] < 0.001 else '**' if results['pearson_p'] < 0.01 else '*' if results['pearson_p'] < 0.05 else 'ns'}\n\n")
        
        f.write(f"Spearman correlation coefficient: ρ = {results['spearman_r']:.6f}\n")
        f.write(f"  p-value: {results['spearman_p']:.4e}\n")
        f.write(f"  Significance: {'***' if results['spearman_p'] < 0.001 else '**' if results['spearman_p'] < 0.01 else '*' if results['spearman_p'] < 0.05 else 'ns'}\n\n")
        
        f.write(f"Linear Regression:\n")
        f.write(f"  R² (coefficient of determination): {results['r2']:.6f}\n")
        f.write(f"  Equation: diversity_per_mb = {results['slope']:.8f} × hashes_per_mb + {results['intercept']:.8f}\n")
        f.write(f"  Variance explained: {results['r2']*100:.2f}%\n\n")
        
        f.write("-"*70 + "\n")
        f.write("INTERPRETATION\n")
        f.write("-"*70 + "\n")
        
        # Determine strength
        r = abs(results['pearson_r'])
        if r > 0.7:
            strength = "strong"
        elif r > 0.4:
            strength = "moderate"
        elif r > 0.2:
            strength = "weak"
        else:
            strength = "very weak"
        
        direction = "positive" if results['pearson_r'] > 0 else "negative"
        
        if results['pearson_p'] < 0.001:
            significance = "highly statistically significant (p < 0.001)"
        elif results['pearson_p'] < 0.01:
            significance = "statistically significant (p < 0.01)"
        elif results['pearson_p'] < 0.05:
            significance = "statistically significant (p < 0.05)"
        else:
            significance = "not statistically significant (p ≥ 0.05)"
        
        f.write(f"The hypothesis that distinct hashes per basepair correlates with alpha diversity\n")
        f.write(f"shows a {strength} {direction} correlation that is {significance}.\n\n")
        
        f.write(f"The R² value of {results['r2']:.4f} indicates that {results['r2']*100:.2f}% of the variance\n")
        f.write(f"in normalized alpha diversity can be explained by the normalized hash count.\n\n")
        
        if results['pearson_r'] > 0:
            f.write("The positive correlation suggests that samples with more distinct k-mer hashes\n")
            f.write("per megabase tend to have higher taxonomic diversity per megabase. This supports\n")
            f.write("the hypothesis that hash-based metrics capture meaningful information about\n")
            f.write("community composition and diversity.\n\n")
        else:
            f.write("The negative correlation is unexpected and warrants further investigation.\n")
            f.write("This could indicate technical artifacts or biases in the data.\n\n")
        
        f.write("-"*70 + "\n")
        f.write("RECOMMENDATIONS\n")
        f.write("-"*70 + "\n")
        
        if results['r2'] < 0.3:
            f.write("• The moderate R² suggests other factors also influence alpha diversity.\n")
            f.write("  Consider investigating:\n")
            f.write("  - Sample type/environment effects\n")
            f.write("  - Sequencing platform biases\n")
            f.write("  - Coverage threshold sensitivity\n\n")
        
        f.write("• Examine the residual plot to identify potential outliers or systematic biases\n")
        f.write("• Consider stratifying analysis by environment type or sequencing platform\n")
        f.write("• Test sensitivity to coverage threshold parameter\n")
        f.write("• Investigate samples with high leverage (extreme hash counts)\n\n")
        
        f.write("-"*70 + "\n")
        f.write("OUTPUT FILES\n")
        f.write("-"*70 + "\n")
        f.write(f"Data:\n")
        f.write(f"  {config.output_dir / 'data' / 'selected_samples.csv'}\n")
        f.write(f"  {config.output_dir / 'data' / 'hash_diversity_data.csv'}\n\n")
        f.write(f"Plots:\n")
        f.write(f"  {config.output_dir / 'plots' / 'hash_diversity_correlation.png'}\n")
        f.write(f"  {config.output_dir / 'plots' / 'hash_diversity_hexbin.png'}\n")
        f.write(f"  {config.output_dir / 'plots' / 'distributions.png'}\n")
        f.write(f"  {config.output_dir / 'plots' / 'residuals.png'}\n\n")
        f.write(f"Reports:\n")
        f.write(f"  {config.output_dir / 'reports' / 'analysis_report.txt'}\n")
        f.write(f"  {config.output_dir / 'reports' / 'statistics_summary.csv'}\n\n")
        
        f.write("="*70 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*70 + "\n")
    
    print(f"Saved report to: {report_file}")
    
    # Also save statistics as CSV for easy import
    stats_df = pd.DataFrame([results])
    stats_file = config.output_dir / "reports" / "statistics_summary.csv"
    stats_df.to_csv(stats_file, index=False)
    print(f"Saved statistics to: {stats_file}")


def main():
    """Main analysis pipeline."""
    parser = argparse.ArgumentParser(
        description="Test correlation between distinct hashes and alpha diversity",
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
        '--coverage', '-c',
        type=float,
        default=0.0625,
        help='Minimum coverage threshold for YACHT algorithm'
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
    parser.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='DPI for saved figures'
    )
    
    args = parser.parse_args()
    
    # Create configuration
    config = Config(
        output_dir=args.output_dir,
        n_samples=args.n_samples,
        coverage=args.coverage,
        n_jobs=args.n_jobs,
        min_mbases=args.min_mbases,
        random_seed=args.random_seed,
        dpi=args.dpi
    )
    
    print("\n" + "="*70)
    print("HASH-DIVERSITY CORRELATION ANALYSIS")
    print("="*70)
    print(f"Configuration:")
    print(f"  Output directory: {config.output_dir}")
    print(f"  Number of samples: {'all' if config.n_samples is None else f'{config.n_samples:,}'}")
    print(f"  Coverage threshold: {config.coverage}")
    print(f"  Parallel workers: {config.n_jobs}")
    print(f"  Minimum mbases: {config.min_mbases}")
    print(f"  Random seed: {config.random_seed}")
    print("="*70)
    
    start_time = time.time()
    
    # Step 1: Get WGS samples
    samples_df = get_wgs_samples(config)
    
    # Step 2: Extract hash and diversity data
    data_df = extract_hash_and_diversity_data(samples_df, config)
    
    # Step 3: Perform correlation analysis
    results = perform_correlation_analysis(data_df, config)
    
    # Step 4: Create visualizations
    create_visualizations(data_df, results, config)
    
    # Step 5: Generate report
    generate_report(data_df, results, config)
    
    elapsed = time.time() - start_time
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"Total time: {elapsed/60:.1f} minutes")
    print(f"\nKey findings:")
    print(f"  Pearson r = {results['pearson_r']:.4f} (p = {results['pearson_p']:.2e})")
    print(f"  R² = {results['r2']:.4f}")
    print(f"  Samples analyzed: {results['n_samples']:,}")
    print(f"\nOutputs saved to: {config.output_dir}/")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()
