#!/usr/bin/env python3
"""
Functional Diversity Downstream Analysis Helper Script
=======================================================
Demonstrates how to use the functional diversity parquet files for custom analysis,
filtering, and visualization with metadata.

Features:
- Apply quality filters (min hashes, min diversity, min mbases)
- Join with metadata database
- Generate filtered correlation plots for ALL diversity metrics
- Create correlation plots colored by categorical metadata variables:
  * center_name, instrument, platform, organism
  * librarylayout, libraryselection, librarysource
  * mbytes (binned by order of magnitude)
  * mbases (binned by order of magnitude)
  * avgspotlen (binned by read length)
  * releasedate (binned by year)
  * biome (parsed from attributes JSON)

Diversity Metrics Analyzed:
- observed_richness (number of KOs)
- shannon_index
- simpson_index / gini_simpson
- hill_2 (effective number of species)
- berger_parker (dominance)
- pielou_evenness

Usage:
    python3 functional_analyze_parquet_data.py --input results/data/functional_hash_diversity_data.parquet --output filtered_analysis --join-metadata
"""

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import duckdb
import numpy as np
import time
import logging
from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Set matplotlib to non-interactive backend for parallel processing
import matplotlib
matplotlib.use('Agg')

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.3)

# Import optional libraries
try:
    import dcor
    HAS_DCOR = True
except ImportError:
    HAS_DCOR = False

try:
    from minepy import MINE
    HAS_MINE = True
except ImportError:
    HAS_MINE = False

# Diversity metrics available in the data
DIVERSITY_METRICS = [
    ('observed_richness_per_mb', 'Observed Richness per Mb (KO count)'),
    ('shannon_index', 'Shannon Index'),
    ('gini_simpson', 'Gini-Simpson Index'),
    ('hill_2', 'Hill Number Order 2'),
    ('berger_parker', 'Berger-Parker Index (Dominance)'),
    ('pielou_evenness', "Pielou's Evenness"),
]


def log_step(step_name: str, start: bool = True):
    """Log the start or end of a processing step with timestamp."""
    if start:
        logger.info(f"{'='*60}")
        logger.info(f"STARTING: {step_name}")
        logger.info(f"{'='*60}")
    else:
        logger.info(f"COMPLETED: {step_name}")
        logger.info(f"{'-'*60}")


def load_parquet_data(parquet_file: Path) -> pd.DataFrame:
    """Load the parquet file."""
    log_step("Loading parquet data")
    start = time.time()
    
    logger.info(f"Loading data from: {parquet_file}")
    df = pd.read_parquet(parquet_file)
    
    elapsed = time.time() - start
    logger.info(f"Loaded {len(df):,} samples in {elapsed:.1f}s")
    logger.info(f"Columns: {df.columns.tolist()}")
    
    log_step(f"Loading parquet data (took {elapsed:.1f}s)", start=False)
    return df


def apply_filters(df: pd.DataFrame, min_hashes: int = 1000, min_diversity: int = 1, 
                  min_mbases: int = 0) -> pd.DataFrame:
    """Apply quality filters to the data."""
    log_step("Applying filters")
    start = time.time()
    
    logger.info(f"Filter settings:")
    logger.info(f"  Minimum total hashes: {min_hashes:,}")
    logger.info(f"  Minimum observed richness (KOs): {min_diversity}")
    if min_mbases > 0:
        logger.info(f"  Minimum mbases: {min_mbases:,}")
    
    original_count = len(df)
    
    # Build filter conditions
    filter_mask = (df['total_distinct_hashes'] >= min_hashes)
    
    # Use observed_richness or alpha_diversity
    if 'observed_richness' in df.columns:
        filter_mask = filter_mask & (df['observed_richness'] >= min_diversity)
    elif 'alpha_diversity' in df.columns:
        filter_mask = filter_mask & (df['alpha_diversity'] >= min_diversity)
    
    # Add mbases filter if column exists and threshold > 0
    if min_mbases > 0 and 'mbases' in df.columns:
        filter_mask = filter_mask & (df['mbases'] >= min_mbases)
    elif min_mbases > 0:
        logger.warning("mbases column not found, skipping mbases filter")
    
    df_filtered = df[filter_mask].copy()
    
    filtered_count = len(df_filtered)
    elapsed = time.time() - start
    
    logger.info(f"Filtering results:")
    logger.info(f"  Original samples: {original_count:,}")
    logger.info(f"  After filtering: {filtered_count:,}")
    logger.info(f"  Removed: {original_count - filtered_count:,} ({100*(original_count-filtered_count)/original_count:.1f}%)")
    
    log_step(f"Applying filters (took {elapsed:.1f}s)", start=False)
    return df_filtered


def join_metadata(df: pd.DataFrame, metadata_db: str) -> pd.DataFrame:
    """Join with metadata database to get additional sample information."""
    log_step("Joining with metadata")
    start = time.time()
    
    logger.info(f"Metadata database: {metadata_db}")
    
    conn = duckdb.connect(metadata_db, read_only=True, config={'threads': 1})
    
    # Get metadata for samples
    accessions = df['accession'].tolist()
    accessions_str = "','".join(accessions)
    
    query = f"""
    SELECT 
        acc,
        biosample,
        organism,
        instrument,
        platform,
        sra_study,
        bioproject,
        center_name,
        librarylayout,
        libraryselection,
        librarysource,
        releasedate,
        mbytes,
        mbases,
        avgspotlen,
        biome,
        attributes
    FROM metadata_geo_joined 
    WHERE acc IN ('{accessions_str}')
    """
    
    logger.info(f"Querying metadata for {len(accessions):,} samples...")
    metadata_df = conn.execute(query).fetchdf()
    conn.close()
    
    # Merge with original data
    logger.info("Merging with original data...")
    df_merged = df.merge(metadata_df, left_on='accession', right_on='acc', how='left')
    df_merged = df_merged.drop('acc', axis=1)
    
    elapsed = time.time() - start
    logger.info(f"Successfully joined metadata for {len(df_merged):,} samples")
    logger.info(f"New columns: {[c for c in df_merged.columns if c not in df.columns]}")
    
    log_step(f"Joining with metadata (took {elapsed:.1f}s)", start=False)
    return df_merged


def calculate_metrics(x: np.ndarray, y: np.ndarray, max_samples_expensive: int = 50000) -> dict:
    """
    Calculate all correlation metrics.
    
    Args:
        x, y: Arrays to correlate
        max_samples_expensive: Max samples for expensive calculations (dcor, MIC)
    """
    n = len(x)
    
    # Fast correlations on full data
    pearson_r, pearson_p = pearsonr(x, y)
    spearman_r, spearman_p = spearmanr(x, y)
    kendall_tau, kendall_p = kendalltau(x, y)
    
    # Subsample for expensive correlations if needed
    if n > max_samples_expensive:
        np.random.seed(42)
        indices = np.random.choice(n, max_samples_expensive, replace=False)
        x_sub = x[indices]
        y_sub = y[indices]
    else:
        x_sub = x
        y_sub = y
    
    if HAS_DCOR:
        try:
            distance_corr = dcor.distance_correlation(x_sub, y_sub, method='AVL')
        except:
            distance_corr = np.nan
    else:
        distance_corr = np.nan
    
    if HAS_MINE and len(x_sub) >= 50:
        try:
            mine = MINE(alpha=0.6, c=15)
            mine.compute_score(x_sub, y_sub)
            mic_score = mine.mic()
        except:
            mic_score = np.nan
    else:
        mic_score = np.nan
    
    X = x.reshape(-1, 1)
    model = LinearRegression()
    model.fit(X, y)
    y_pred = model.predict(X)
    r2 = r2_score(y, y_pred)
    slope = model.coef_[0]
    intercept = model.intercept_
    
    return {
        'n': n,
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p,
        'kendall_tau': kendall_tau,
        'kendall_p': kendall_p,
        'distance_corr': distance_corr,
        'mic': mic_score,
        'r2': r2,
        'slope': slope,
        'intercept': intercept
    }


def plot_filtered_correlation(df: pd.DataFrame, output_dir: Path, label: str = "filtered",
                              max_samples_expensive: int = 50000):
    """Create correlation plots for all diversity metrics."""
    log_step(f"Creating correlation plots ({label})")
    step_start = time.time()
    
    all_stats = {}
    
    for metric_col, metric_name in DIVERSITY_METRICS:
        if metric_col not in df.columns:
            logger.warning(f"Skipping {metric_name}: column not found")
            continue
        
        df_clean = df.dropna(subset=['hashes_per_mb', metric_col])
        
        if len(df_clean) < 10:
            logger.warning(f"Skipping {metric_name}: insufficient data (n={len(df_clean)})")
            continue
        
        metric_start = time.time()
        logger.info(f"Processing: {metric_name} (n={len(df_clean):,})")
        
        x = df_clean['hashes_per_mb'].values
        y = df_clean[metric_col].values
        
        metrics = calculate_metrics(x, y, max_samples_expensive)
        all_stats[metric_col] = metrics
        
        # Create plot
        fig, ax = plt.subplots(figsize=(12, 9))
        
        ax.scatter(x, y, alpha=0.3, s=20, c='steelblue', edgecolors='none', rasterized=True)
        
        # Regression line
        X = x.reshape(-1, 1)
        model = LinearRegression()
        model.fit(X, y)
        y_pred = model.predict(X)
        sort_idx = np.argsort(x)
        ax.plot(x[sort_idx], y_pred[sort_idx], 'r-', linewidth=2, label='Linear fit', zorder=10)
        
        # Labels
        ax.set_xlabel('Functional k-mer hash density per megabase\n(11-mer AA FracMinHash, scale=1000)', 
                      fontsize=14, fontweight='bold')
        ax.set_ylabel(metric_name, fontsize=14, fontweight='bold')
        ax.set_title(f'Correlation: Hash Density vs {metric_name}\n({label})\nn = {len(df_clean):,}',
                     fontsize=16, fontweight='bold', pad=20)
        
        # Stats box
        stats_text = f'n = {metrics["n"]:,}\n\n'
        stats_text += 'Linear:\n'
        stats_text += f'  Pearson r = {metrics["pearson_r"]:.4f}\n'
        stats_text += f'  R² = {metrics["r2"]:.4f}\n'
        stats_text += f'  p = {metrics["pearson_p"]:.2e}\n\n'
        stats_text += 'Monotonic:\n'
        stats_text += f'  Spearman ρ = {metrics["spearman_r"]:.4f}\n'
        stats_text += f'  Kendall τ = {metrics["kendall_tau"]:.4f}\n\n'
        stats_text += 'Non-linear:\n'
        if not np.isnan(metrics['distance_corr']):
            stats_text += f'  Distance corr = {metrics["distance_corr"]:.4f}\n'
        if not np.isnan(metrics['mic']):
            stats_text += f'  MIC = {metrics["mic"]:.4f}'
        
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                fontsize=11, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                family='monospace')
        
        ax.legend(loc='lower right', fontsize=12)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        safe_name = metric_col.replace('/', '_').replace(' ', '_')
        plot_file = output_dir / f"correlation_{label}_{safe_name}.png"
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        metric_elapsed = time.time() - metric_start
        logger.info(f"  Saved: {plot_file.name} (took {metric_elapsed:.1f}s)")
        logger.info(f"  Spearman ρ={metrics['spearman_r']:.4f}, Pearson r={metrics['pearson_r']:.4f}, R²={metrics['r2']:.4f}")
    
    # Create summary multi-panel figure
    logger.info("Creating summary panel...")
    create_summary_panel(df, output_dir, label, all_stats)
    
    step_elapsed = time.time() - step_start
    log_step(f"Creating correlation plots (took {step_elapsed:.1f}s)", start=False)
    
    # Return stats for primary metric (observed_richness_per_mb) for compatibility
    primary_stats = all_stats.get('observed_richness_per_mb', 
                                   all_stats.get(list(all_stats.keys())[0] if all_stats else 'n/a', {}))
    return primary_stats


def create_summary_panel(df: pd.DataFrame, output_dir: Path, label: str, all_stats: dict):
    """Create a summary multi-panel figure for all metrics."""
    
    n_metrics = len(all_stats)
    if n_metrics == 0:
        return
    
    n_cols = min(3, n_metrics)
    n_rows = int(np.ceil(n_metrics / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_metrics == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    for idx, ((metric_col, metrics), (_, metric_name)) in enumerate(zip(all_stats.items(), 
                                                                        [(k, v) for k, v in DIVERSITY_METRICS if k in all_stats])):
        ax = axes[idx]
        df_clean = df.dropna(subset=['hashes_per_mb', metric_col])
        
        if len(df_clean) < 10:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
            continue
        
        x = df_clean['hashes_per_mb'].values
        y = df_clean[metric_col].values
        
        ax.scatter(x, y, alpha=0.3, s=10, c='steelblue', edgecolors='none', rasterized=True)
        
        # Regression line
        X = x.reshape(-1, 1)
        model = LinearRegression()
        model.fit(X, y)
        y_pred = model.predict(X)
        sort_idx = np.argsort(x)
        ax.plot(x[sort_idx], y_pred[sort_idx], 'r-', linewidth=2, alpha=0.8)
        
        # Stats
        stats_text = f"ρ={metrics['spearman_r']:.3f}\nr={metrics['pearson_r']:.3f}\nR²={metrics['r2']:.3f}"
        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                family='monospace')
        
        ax.set_xlabel('Hashes per Mb', fontsize=10)
        short_name = metric_name.replace(' per Mb', '').replace(' (KO count)', '')
        ax.set_ylabel(short_name, fontsize=10)
        ax.set_title(short_name, fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
    
    # Hide unused subplots
    for idx in range(len(all_stats), len(axes)):
        axes[idx].set_visible(False)
    
    fig.suptitle(f'Functional Hash-Diversity Correlations ({label})', fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    plt.savefig(output_dir / f"summary_all_metrics_{label}.png", dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved summary figure: summary_all_metrics_{label}.png")


def parse_biome(attributes_json):
    """Extract biome information from attributes JSON."""
    import json
    try:
        if pd.isna(attributes_json):
            return 'unknown'
        attrs = json.loads(attributes_json) if isinstance(attributes_json, str) else attributes_json
        for attr in attrs:
            tag = attr.get('tag', '').lower()
            if 'env_biome' in tag or 'biome' in tag:
                return attr.get('value', 'unknown')
            if 'environment' in tag and 'environmental' not in tag:
                return attr.get('value', 'unknown')
        return 'unknown'
    except:
        return 'unknown'


def bin_continuous_variable(series: pd.Series, var_name: str) -> pd.Series:
    """Bin continuous variables intelligently."""
    if var_name in ['mbytes', 'mbases']:
        bins = [0, 10, 100, 1000, 10000, 100000, float('inf')]
        if var_name == 'mbytes':
            labels = ['<10 MB', '10-100 MB', '100 MB-1 GB', '1-10 GB', '10-100 GB', '>100 GB']
        else:
            labels = ['<10 Mb', '10-100 Mb', '100-1000 Mb', '1-10 Gb', '10-100 Gb', '>100 Gb']
        return pd.cut(series, bins=bins, labels=labels, include_lowest=True)
    
    elif var_name == 'avgspotlen':
        bins = [0, 50, 100, 150, 250, 500, float('inf')]
        labels = ['<50 bp', '50-100 bp', '100-150 bp', '150-250 bp', '250-500 bp', '>500 bp']
        return pd.cut(series, bins=bins, labels=labels, include_lowest=True)
    
    elif var_name == 'releasedate':
        try:
            years = pd.to_datetime(series, errors='coerce').dt.year
            return years.apply(lambda y: f"{int(y)}" if pd.notna(y) else 'unknown')
        except:
            return pd.Series(['unknown'] * len(series))
    
    else:
        return pd.qcut(series, q=4, labels=['Q1', 'Q2', 'Q3', 'Q4'], duplicates='drop')


def consolidate_rare_categories(series: pd.Series, top_n: int = 10, other_label: str = 'other') -> pd.Series:
    """Group rare categories into 'other' category."""
    value_counts = series.value_counts()
    top_categories = value_counts.head(top_n).index.tolist()
    return series.apply(lambda x: x if x in top_categories else other_label)


def calculate_category_statistics(df: pd.DataFrame, category_col: str, 
                                  diversity_col: str, min_samples: int = 30,
                                  max_samples_expensive: int = 50000) -> pd.DataFrame:
    """
    Calculate correlation statistics for each value within a categorical variable.
    
    Args:
        df: DataFrame with data and category column
        category_col: Name of the categorical column
        diversity_col: Name of the diversity metric column
        min_samples: Minimum number of samples required for statistical calculation
        max_samples_expensive: Max samples for expensive correlations
        
    Returns:
        DataFrame with statistics for each category value
    """
    if category_col not in df.columns:
        return pd.DataFrame()
    
    df_clean = df[df[category_col].notna()].copy()
    
    if len(df_clean) == 0:
        return pd.DataFrame()
    
    results = []
    
    for category_value in df_clean[category_col].unique():
        category_data = df_clean[df_clean[category_col] == category_value]
        
        if len(category_data) < min_samples:
            continue
        
        category_data = category_data[
            category_data['hashes_per_mb'].notna() & 
            category_data[diversity_col].notna()
        ]
        
        if len(category_data) < min_samples:
            continue
        
        try:
            x = category_data['hashes_per_mb'].values
            y = category_data[diversity_col].values
            
            metrics = calculate_metrics(x, y, max_samples_expensive)
            
            # Add additional statistics
            results.append({
                'variable': category_col,
                'category': str(category_value),
                'diversity_metric': diversity_col,
                'n_samples': len(category_data),
                'pearson_r': metrics['pearson_r'],
                'pearson_p': metrics['pearson_p'],
                'r_squared': metrics['r2'],
                'spearman_r': metrics['spearman_r'],
                'spearman_p': metrics['spearman_p'],
                'kendall_tau': metrics['kendall_tau'],
                'kendall_p': metrics['kendall_p'],
                'distance_corr': metrics['distance_corr'],
                'mic': metrics['mic'],
                'slope': metrics['slope'],
                'intercept': metrics['intercept'],
                'mean_hashes_per_mb': x.mean(),
                'mean_diversity': y.mean()
            })
            
        except Exception:
            continue
    
    if not results:
        return pd.DataFrame()
    
    return pd.DataFrame(results)


def apply_multiple_testing_correction(stats_df: pd.DataFrame, 
                                      alpha: float = 0.05) -> pd.DataFrame:
    """
    Apply multiple testing correction to p-values.
    
    Args:
        stats_df: DataFrame with statistics including p-values
        alpha: Significance level
        
    Returns:
        DataFrame with corrected p-values and significance flags
    """
    from statsmodels.stats.multitest import multipletests
    
    if len(stats_df) == 0:
        return stats_df
    
    # Bonferroni correction
    stats_df['pearson_p_bonferroni'] = stats_df['pearson_p'] * len(stats_df)
    stats_df['pearson_p_bonferroni'] = stats_df['pearson_p_bonferroni'].clip(upper=1.0)
    
    # Benjamini-Hochberg FDR correction
    reject, pvals_corrected, _, _ = multipletests(
        stats_df['pearson_p'].values, 
        alpha=alpha, 
        method='fdr_bh'
    )
    
    stats_df['pearson_p_fdr'] = pvals_corrected
    stats_df['significant_fdr'] = reject
    stats_df['significant_bonferroni'] = stats_df['pearson_p_bonferroni'] < alpha
    
    return stats_df


def calculate_all_category_statistics(df: pd.DataFrame, 
                                      categorical_vars: list,
                                      diversity_col: str,
                                      min_samples: int = 30,
                                      max_samples_expensive: int = 50000) -> pd.DataFrame:
    """
    Calculate statistics for all categorical variables and combine.
    
    Args:
        df: DataFrame with data and metadata
        categorical_vars: List of categorical variable names
        diversity_col: Name of the diversity metric column
        min_samples: Minimum samples per category
        max_samples_expensive: Max samples for expensive correlations
        
    Returns:
        Combined DataFrame with all statistics
    """
    all_stats = []
    
    logger.info(f"Calculating statistics for each category ({diversity_col})...")
    
    for var in categorical_vars:
        if var not in df.columns:
            continue
        
        var_stats = calculate_category_statistics(df, var, diversity_col, min_samples, max_samples_expensive)
        
        if not var_stats.empty:
            all_stats.append(var_stats)
            logger.info(f"  {var}: {len(var_stats)} categories with ≥{min_samples} samples")
    
    if not all_stats:
        return pd.DataFrame()
    
    combined_stats = pd.concat(all_stats, ignore_index=True)
    
    # Apply multiple testing correction
    logger.info(f"Applying multiple testing correction (n={len(combined_stats)} comparisons)...")
    combined_stats = apply_multiple_testing_correction(combined_stats)
    
    # Sort by R-squared
    combined_stats = combined_stats.sort_values('r_squared', ascending=False)
    
    return combined_stats


def create_multipanel_figure(df: pd.DataFrame, category_col: str, 
                            stats_df: pd.DataFrame, output_dir: Path,
                            diversity_col: str = 'observed_richness_per_mb',
                            max_panels: int = 12):
    """
    Create multi-panel figure showing each category separately.
    
    Args:
        df: DataFrame with data
        category_col: Categorical column name
        stats_df: Statistics for this variable's categories
        output_dir: Output directory
        diversity_col: Which diversity metric to plot
        max_panels: Maximum number of panels to show
    """
    if category_col not in df.columns:
        return
    
    # Get categories from stats, sorted by FREQUENCY (n_samples), not R²
    var_stats = stats_df[stats_df['variable'] == category_col].copy()
    var_stats = var_stats.sort_values('n_samples', ascending=False)
    categories = var_stats['category'].tolist()
    
    if not categories:
        return
    
    # Limit to max_panels
    categories = categories[:max_panels]
    n_categories = len(categories)
    
    if n_categories == 0:
        return
    
    logger.info(f"Creating multi-panel figure for: {category_col} ({n_categories} panels)")
    
    # Calculate grid dimensions
    n_cols = min(3, n_categories)
    n_rows = int(np.ceil(n_categories / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
    
    # Flatten axes for easier iteration
    if n_categories == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if n_rows > 1 or n_cols > 1 else [axes]
    
    for idx, category in enumerate(categories):
        ax = axes[idx]
        
        # Get data for this category
        cat_data = df[df[category_col] == category].copy()
        cat_data = cat_data[cat_data['hashes_per_mb'].notna() & 
                           cat_data[diversity_col].notna()]
        
        if len(cat_data) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', 
                   transform=ax.transAxes)
            ax.set_title(str(category))
            continue
        
        # Scatter plot
        ax.scatter(cat_data['hashes_per_mb'], cat_data[diversity_col],
                  alpha=0.4, s=15, c='steelblue', edgecolors='none', rasterized=True)
        
        # Regression line
        X = cat_data['hashes_per_mb'].values.reshape(-1, 1)
        y = cat_data[diversity_col].values
        model = LinearRegression()
        model.fit(X, y)
        y_pred = model.predict(X)
        
        # Sort for clean line
        sort_idx = np.argsort(X.flatten())
        ax.plot(X[sort_idx], y_pred[sort_idx], 'r-', linewidth=2, alpha=0.8)
        
        # Get statistics for this category
        cat_stats = stats_df[(stats_df['variable'] == category_col) & 
                            (stats_df['category'] == str(category))]
        
        if not cat_stats.empty:
            r = cat_stats.iloc[0]['pearson_r']
            r2 = cat_stats.iloc[0]['r_squared']
            spearman = cat_stats.iloc[0]['spearman_r']
            kendall = cat_stats.iloc[0]['kendall_tau']
            dcorr = cat_stats.iloc[0]['distance_corr']
            mic = cat_stats.iloc[0]['mic']
            p = cat_stats.iloc[0]['pearson_p_fdr']
            n = cat_stats.iloc[0]['n_samples']
            
            # Determine significance
            if p < 0.001:
                sig = '***'
            elif p < 0.01:
                sig = '**'
            elif p < 0.05:
                sig = '*'
            else:
                sig = 'ns'
            
            # Add stats box with all metrics
            stats_text = f'n={n:,}{sig}\n'
            stats_text += f'r={r:.3f}\n'
            stats_text += f'ρ={spearman:.3f}\n'
            stats_text += f'τ={kendall:.3f}\n'
            if not pd.isna(dcorr):
                stats_text += f'dCor={dcorr:.3f}\n'
            if not pd.isna(mic):
                stats_text += f'MIC={mic:.3f}\n'
            stats_text += f'R²={r2:.3f}'
            
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                   fontsize=8, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                   family='monospace')
        
        # Labels
        ax.set_xlabel('Hashes per Mb', fontsize=9)
        
        # Get nice y-label
        metric_names = dict(DIVERSITY_METRICS)
        y_label = metric_names.get(diversity_col, diversity_col).replace(' per Mb', '').replace(' (KO count)', '')
        ax.set_ylabel(y_label, fontsize=9)
        
        # Title with category name (truncate if too long)
        title = str(category)
        if len(title) > 30:
            title = title[:27] + '...'
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.grid(True, alpha=0.3)
    
    # Hide unused subplots
    for idx in range(n_categories, len(axes)):
        axes[idx].set_visible(False)
    
    # Overall title
    clean_name = category_col.replace('_', ' ').title()
    metric_short = diversity_col.replace('_per_mb', '').replace('_', ' ').title()
    fig.suptitle(f'Functional Hash-{metric_short} Correlation by {clean_name}\n(Most Frequent Categories)', 
                fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    
    # Save
    safe_name = category_col.replace('/', '_').replace(' ', '_')
    safe_metric = diversity_col.replace('/', '_').replace(' ', '_')
    plot_file = output_dir / f"multipanel_{safe_name}_{safe_metric}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.debug(f"  Saved: {plot_file.name}")


def generate_statistical_report(stats_df: pd.DataFrame, output_dir: Path, 
                                filters_applied: dict):
    """
    Generate comprehensive markdown report of statistical findings.
    
    Args:
        stats_df: DataFrame with all category statistics
        output_dir: Output directory
        filters_applied: Dictionary of filters that were applied
    """
    report_file = output_dir / "statistical_report.md"
    
    with open(report_file, 'w') as f:
        f.write("# Functional Diversity Correlation Analysis Report\n\n")
        f.write("## Analysis Parameters\n\n")
        f.write(f"- **Input file**: {filters_applied.get('input_file', 'N/A')}\n")
        f.write(f"- **Minimum hashes**: {filters_applied.get('min_hashes', 'N/A'):,}\n")
        f.write(f"- **Minimum diversity**: {filters_applied.get('min_diversity', 'N/A')}\n")
        f.write(f"- **Minimum mbases**: {filters_applied.get('min_mbases', 'N/A')}\n")
        f.write(f"- **Samples after filtering**: {filters_applied.get('n_samples_after_filtering', 'N/A'):,}\n\n")
        
        f.write("## Summary Statistics\n\n")
        f.write(f"- **Total categories analyzed**: {len(stats_df)}\n")
        f.write(f"- **Significant after FDR correction**: {stats_df['significant_fdr'].sum()}\n")
        f.write(f"- **Significant after Bonferroni correction**: {stats_df['significant_bonferroni'].sum()}\n\n")
        
        f.write("## Top Categories by R²\n\n")
        f.write("| Variable | Category | n | r | ρ | R² | FDR p |\n")
        f.write("|----------|----------|---|---|---|----|---------|\n")
        
        for _, row in stats_df.head(20).iterrows():
            sig = '***' if row['pearson_p_fdr'] < 0.001 else ('**' if row['pearson_p_fdr'] < 0.01 else ('*' if row['pearson_p_fdr'] < 0.05 else ''))
            f.write(f"| {row['variable']} | {row['category'][:30]} | {row['n_samples']:,} | "
                   f"{row['pearson_r']:.3f} | {row['spearman_r']:.3f} | {row['r_squared']:.3f}{sig} | "
                   f"{row['pearson_p_fdr']:.2e} |\n")
        
        f.write("\n## Interpretation\n\n")
        f.write("- **r**: Pearson correlation coefficient (linear relationship)\n")
        f.write("- **ρ (rho)**: Spearman rank correlation (monotonic relationship)\n")
        f.write("- **R²**: Coefficient of determination (variance explained)\n")
        f.write("- **FDR p**: False Discovery Rate corrected p-value\n")
        f.write("- Significance: *** p<0.001, ** p<0.01, * p<0.05\n\n")
        
        f.write("## Output Files\n\n")
        f.write("- `category_statistics.csv` - Full statistics for all categories\n")
        f.write("- `correlation_by_*.png` - Overlay plots colored by category\n")
        f.write("- `multipanel_*.png` - Multi-panel figures showing most frequent categories\n")
    
    logger.info(f"Saved statistical report to: {report_file}")


def plot_correlation_by_category(df: pd.DataFrame, category_col: str, diversity_col: str,
                                 output_dir: Path, max_categories: int = 10):
    """Create correlation plot colored by a categorical variable."""
    
    if category_col not in df.columns or diversity_col not in df.columns:
        logger.debug(f"Skipping {category_col}/{diversity_col}: column not found")
        return None
    
    df_plot = df[df[category_col].notna()].copy()
    
    if len(df_plot) == 0:
        logger.debug(f"Skipping {category_col}/{diversity_col}: no data available")
        return None
    
    category_counts = df_plot[category_col].value_counts()
    
    if len(category_counts) > max_categories:
        df_plot[category_col] = consolidate_rare_categories(df_plot[category_col], top_n=max_categories)
        category_counts = df_plot[category_col].value_counts()
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    colors = sns.color_palette("tab20", n_colors=min(20, len(category_counts)))
    
    for idx, (category, color) in enumerate(zip(category_counts.index, colors)):
        category_data = df_plot[df_plot[category_col] == category]
        
        if len(category_data) == 0:
            continue
        
        label = f"{category} (n={len(category_data):,})"
        
        ax.scatter(
            category_data['hashes_per_mb'],
            category_data[diversity_col],
            alpha=0.4,
            s=15,
            color=color,
            label=label,
            rasterized=True
        )
    
    ax.set_xlabel('Functional k-mer hash density per megabase\n(11-mer AA FracMinHash, scale=1000)', 
                  fontsize=12, fontweight='bold')
    
    # Get nice name for diversity metric
    metric_names = dict(DIVERSITY_METRICS)
    y_label = metric_names.get(diversity_col, diversity_col)
    ax.set_ylabel(y_label, fontsize=12, fontweight='bold')
    
    title_name = category_col.replace('_', ' ').title()
    ax.set_title(f'Functional Hash-Diversity Correlation by {title_name}',
                 fontsize=14, fontweight='bold', pad=15)
    
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    safe_cat = category_col.replace('/', '_').replace(' ', '_')
    safe_div = diversity_col.replace('/', '_').replace(' ', '_')
    plot_file = output_dir / f"correlation_by_{safe_cat}_{safe_div}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return plot_file.name


def _plot_category_task(args: Tuple) -> Optional[str]:
    """Worker function for parallel categorical plotting."""
    df_dict, category_col, diversity_col, output_dir_str, max_categories = args
    
    # Reconstruct DataFrame and Path in worker process
    df = pd.DataFrame(df_dict)
    output_dir = Path(output_dir_str)
    
    try:
        result = plot_correlation_by_category(df, category_col, diversity_col, output_dir, max_categories)
        return result
    except Exception as e:
        return None


def generate_all_categorical_plots(df: pd.DataFrame, output_dir: Path, filters_applied: dict,
                                   n_jobs: int = 8, max_samples_expensive: int = 50000):
    """Generate correlation plots for all categorical metadata variables and diversity metrics."""
    
    log_step("Generating categorical analysis")
    step_start = time.time()
    
    cat_plot_dir = output_dir / "categorical_plots"
    cat_plot_dir.mkdir(exist_ok=True)
    
    categorical_vars = [
        'center_name',
        'instrument', 
        'librarylayout',
        'libraryselection',
        'librarysource',
        'platform',
        'organism'
    ]
    
    # Parse biome if needed
    if 'biome' in df.columns:
        logger.info("Filling missing biome values from attributes...")
        if 'attributes' in df.columns:
            missing_biome = df['biome'].isna()
            df.loc[missing_biome, 'biome'] = df.loc[missing_biome, 'attributes'].apply(parse_biome)
            logger.info(f"  Biome filled for {missing_biome.sum():,} samples")
        categorical_vars.append('biome')
    elif 'attributes' in df.columns:
        logger.info("Parsing biome from attributes...")
        df['biome'] = df['attributes'].apply(parse_biome)
        categorical_vars.append('biome')
    
    # Bin continuous variables
    continuous_to_bin = {
        'mbytes': 'mbytes',
        'mbases': 'mbases', 
        'avgspotlen': 'avgspotlen',
        'releasedate': 'releasedate'
    }
    
    for var_name, col_name in continuous_to_bin.items():
        if col_name in df.columns:
            logger.info(f"Binning {var_name}...")
            df[f'{var_name}_binned'] = bin_continuous_variable(df[col_name], var_name)
            categorical_vars.append(f'{var_name}_binned')
    
    # Metrics to analyze
    metrics_to_analyze = ['observed_richness_per_mb', 'shannon_index']
    
    # Calculate statistics for all categories
    logger.info("-" * 60)
    logger.info("Calculating category statistics...")
    logger.info("-" * 60)
    
    all_stats_combined = []
    for metric in metrics_to_analyze:
        if metric not in df.columns:
            continue
        all_stats = calculate_all_category_statistics(
            df, categorical_vars, metric, 
            min_samples=30, max_samples_expensive=max_samples_expensive
        )
        if not all_stats.empty:
            all_stats_combined.append(all_stats)
    
    if all_stats_combined:
        combined_stats = pd.concat(all_stats_combined, ignore_index=True)
        
        # Save statistics to CSV
        stats_file = cat_plot_dir / "category_statistics.csv"
        combined_stats.to_csv(stats_file, index=False)
        logger.info(f"Saved category statistics to: {stats_file}")
        logger.info(f"  Total categories analyzed: {len(combined_stats)}")
        logger.info(f"  Significant after FDR correction: {combined_stats['significant_fdr'].sum()}")
    else:
        combined_stats = pd.DataFrame()
        logger.warning("No categories had sufficient samples for analysis")
    
    # Generate overlay plots in parallel
    logger.info("-" * 60)
    logger.info("Generating overlay plots...")
    logger.info("-" * 60)
    
    plot_tasks = []
    for metric in metrics_to_analyze:
        if metric not in df.columns:
            continue
        for var in categorical_vars:
            if var in df.columns:
                plot_tasks.append((var, metric))
    
    logger.info(f"Generating {len(plot_tasks)} overlay plots using {n_jobs} workers...")
    
    # Convert DataFrame to dict for pickling (needed for multiprocessing)
    needed_cols = set(['hashes_per_mb'] + metrics_to_analyze + categorical_vars)
    needed_cols = [c for c in needed_cols if c in df.columns]
    df_subset = df[needed_cols].copy()
    
    # Execute plots in parallel
    completed = 0
    failed = 0
    
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = {}
        for var, metric in plot_tasks:
            future = executor.submit(
                _plot_category_task,
                (df_subset.to_dict('list'), var, metric, str(cat_plot_dir), 10)
            )
            futures[future] = (var, metric)
        
        for future in as_completed(futures):
            var, metric = futures[future]
            try:
                result = future.result()
                if result:
                    completed += 1
                else:
                    failed += 1
            except Exception as e:
                failed += 1
                logger.warning(f"  Failed {var}/{metric}: {e}")
    
    logger.info(f"Overlay plots: {completed} completed, {failed} skipped/failed")
    
    # Generate multi-panel figures
    if not combined_stats.empty:
        logger.info("-" * 60)
        logger.info("Generating multi-panel figures...")
        logger.info("-" * 60)
        
        for metric in metrics_to_analyze:
            if metric not in df.columns:
                continue
            metric_stats = combined_stats[combined_stats['diversity_metric'] == metric]
            if metric_stats.empty:
                continue
            for var in categorical_vars:
                if var in df.columns:
                    create_multipanel_figure(df, var, metric_stats, cat_plot_dir, 
                                           diversity_col=metric, max_panels=12)
        
        # Generate statistical report
        logger.info("-" * 60)
        logger.info("Generating statistical report...")
        logger.info("-" * 60)
        generate_statistical_report(combined_stats, output_dir, filters_applied)
    
    step_elapsed = time.time() - step_start
    logger.info(f"Results saved to: {cat_plot_dir}/")
    
    log_step(f"Categorical analysis (took {step_elapsed:.1f}s)", start=False)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze functional diversity parquet data with custom filters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', '-i', required=True, help='Input parquet file')
    parser.add_argument('--output', '-o', required=True, help='Output directory')
    parser.add_argument('--min-hashes', type=int, default=1000, 
                       help='Minimum total distinct hashes')
    parser.add_argument('--min-diversity', type=int, default=1,
                       help='Minimum observed richness (KO count)')
    parser.add_argument('--min-mbases', type=int, default=0,
                       help='Minimum sequencing depth in megabases (0 = no filter)')
    parser.add_argument('--metadata-db', 
                       default='/scratch/shared_data_new/Logan_yacht_data/metadata/aws_sra_metadata/metadata_geo_joined_5M.duckdb',
                       help='Path to metadata database')
    parser.add_argument('--join-metadata', action='store_true',
                       help='Join with metadata database')
    parser.add_argument('--n-jobs', '-j', type=int, default=8,
                       help='Number of parallel workers for plotting')
    parser.add_argument('--max-corr-samples', type=int, default=50000,
                       help='Max samples for expensive correlations (dcor, MIC)')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    total_start = time.time()
    
    logger.info("="*70)
    logger.info("FUNCTIONAL DIVERSITY DOWNSTREAM ANALYSIS")
    logger.info("="*70)
    logger.info(f"Configuration:")
    logger.info(f"  Input: {args.input}")
    logger.info(f"  Output: {args.output}")
    logger.info(f"  Min hashes: {args.min_hashes:,}")
    logger.info(f"  Min diversity: {args.min_diversity}")
    logger.info(f"  Min mbases: {args.min_mbases}")
    logger.info(f"  Join metadata: {args.join_metadata}")
    logger.info(f"  Parallel workers: {args.n_jobs}")
    logger.info(f"  Max samples for dcor/MIC: {args.max_corr_samples:,}")
    logger.info("="*70)
    
    # Load data
    df = load_parquet_data(Path(args.input))
    
    # Apply filters
    df_filtered = apply_filters(df, args.min_hashes, args.min_diversity, args.min_mbases)
    
    # Join metadata if requested
    if args.join_metadata:
        df_filtered = join_metadata(df_filtered, args.metadata_db)
    
    # Create plots
    label_parts = [f"min{args.min_hashes}hashes"]
    if args.min_mbases > 0:
        label_parts.append(f"min{args.min_mbases}mbases")
    label = "_".join(label_parts)
    stats = plot_filtered_correlation(df_filtered, output_dir, label=label, 
                                      max_samples_expensive=args.max_corr_samples)
    
    if args.join_metadata:
        filters_applied = {
            'min_hashes': args.min_hashes,
            'min_diversity': args.min_diversity,
            'min_mbases': args.min_mbases,
            'n_samples_after_filtering': len(df_filtered),
            'input_file': args.input
        }
        generate_all_categorical_plots(df_filtered, output_dir, filters_applied,
                                       n_jobs=args.n_jobs, 
                                       max_samples_expensive=args.max_corr_samples)
    
    # Save filtered data
    log_step("Saving outputs")
    save_start = time.time()
    
    output_file = output_dir / "filtered_data.parquet"
    df_filtered.to_parquet(output_file, index=False)
    logger.info(f"Saved filtered data to: {output_file}")
    
    # Save summary
    summary_file = output_dir / "summary.txt"
    with open(summary_file, 'w') as f:
        f.write("FUNCTIONAL DIVERSITY FILTERED DATA SUMMARY\n")
        f.write("="*70 + "\n\n")
        f.write(f"Input file: {args.input}\n")
        f.write(f"Filters:\n")
        f.write(f"  Minimum hashes: {args.min_hashes:,}\n")
        f.write(f"  Minimum diversity (KOs): {args.min_diversity}\n")
        if args.min_mbases > 0:
            f.write(f"  Minimum mbases: {args.min_mbases:,}\n")
        f.write(f"\nResults:\n")
        f.write(f"  Samples after filtering: {stats.get('n', len(df_filtered)):,}\n")
        if stats:
            f.write(f"  Spearman ρ (richness): {stats.get('spearman_r', 'N/A')}\n")
            f.write(f"  Pearson r (richness): {stats.get('pearson_r', 'N/A')}\n")
            f.write(f"  R² (richness): {stats.get('r2', 'N/A')}\n")
    
    logger.info(f"Saved summary to: {summary_file}")
    
    save_elapsed = time.time() - save_start
    log_step(f"Saving outputs (took {save_elapsed:.1f}s)", start=False)
    
    total_elapsed = time.time() - total_start
    logger.info("="*70)
    logger.info("ANALYSIS COMPLETE")
    logger.info("="*70)
    logger.info(f"Total time: {total_elapsed/60:.1f} minutes")
    logger.info(f"Outputs saved to: {output_dir}/")


if __name__ == "__main__":
    main()
