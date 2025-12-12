#!/usr/bin/env python3
"""
Downstream Analysis Helper Script
==================================
Demonstrates how to use the parquet files for custom analysis, filtering,
and visualization with metadata.

Features:
- Apply quality filters (min hashes, min diversity)
- Join with metadata database
- Generate filtered correlation plots
- Create correlation plots colored by ALL categorical metadata variables:
  * center_name, instrument, platform, organism
  * librarylayout, libraryselection, librarysource
  * mbytes (binned by order of magnitude)
  * mbases (binned by order of magnitude)
  * avgspotlen (binned by read length)
  * releasedate (binned by year)
  * biome (parsed from attributes JSON)

Usage:
    python3 analyze_parquet_data.py --input results/data/hash_diversity_data.parquet --output filtered_analysis --join-metadata
"""

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import duckdb
import numpy as np
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.3)

def load_parquet_data(parquet_file: Path) -> pd.DataFrame:
    """Load the parquet file."""
    print(f"Loading data from: {parquet_file}")
    df = pd.read_parquet(parquet_file)
    print(f"Loaded {len(df):,} samples")
    print(f"Columns: {df.columns.tolist()}")
    return df

def apply_filters(df: pd.DataFrame, min_hashes: int = 1000, min_diversity: int = 1) -> pd.DataFrame:
    """Apply quality filters to the data."""
    print(f"\nApplying filters:")
    print(f"  Minimum total hashes: {min_hashes:,}")
    print(f"  Minimum alpha diversity: {min_diversity}")
    
    original_count = len(df)
    df_filtered = df[
        (df['total_distinct_hashes'] >= min_hashes) &
        (df['alpha_diversity'] >= min_diversity)
    ].copy()
    
    filtered_count = len(df_filtered)
    print(f"\nFiltering results:")
    print(f"  Original samples: {original_count:,}")
    print(f"  After filtering: {filtered_count:,}")
    print(f"  Removed: {original_count - filtered_count:,} ({100*(original_count-filtered_count)/original_count:.1f}%)")
    
    return df_filtered

def join_metadata(df: pd.DataFrame, metadata_db: str) -> pd.DataFrame:
    """
    Join with metadata database to get additional sample information.
    
    Args:
        df: DataFrame with sample data
        metadata_db: Path to metadata DuckDB database
    """
    print(f"\nJoining with metadata from: {metadata_db}")
    
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
    
    metadata_df = conn.execute(query).fetchdf()
    conn.close()
    
    # Merge with original data
    df_merged = df.merge(metadata_df, left_on='accession', right_on='acc', how='left')
    df_merged = df_merged.drop('acc', axis=1)
    
    print(f"Successfully joined metadata for {len(df_merged):,} samples")
    print(f"New columns: {[c for c in df_merged.columns if c not in df.columns]}")
    
    return df_merged

def plot_filtered_correlation(df: pd.DataFrame, output_dir: Path, label: str = "filtered"):
    """Create correlation plot for filtered data."""
    from scipy.stats import pearsonr, spearmanr, kendalltau
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score
    
    # Import optional libraries
    try:
        import dcor
        has_dcor = True
    except ImportError:
        has_dcor = False
    
    try:
        from minepy import MINE
        has_mine = True
    except ImportError:
        has_mine = False
    
    print(f"\nCreating correlation plot for {label} data...")
    
    # Calculate all correlation metrics
    x = df['hashes_per_mb'].values
    y = df['diversity_per_mb'].values
    
    pearson_r, pearson_p = pearsonr(x, y)
    spearman_r, spearman_p = spearmanr(x, y)
    kendall_tau, kendall_p = kendalltau(x, y)
    
    if has_dcor:
        try:
            distance_corr = dcor.distance_correlation(x, y)
        except:
            distance_corr = np.nan
    else:
        distance_corr = np.nan
    
    if has_mine and len(df) >= 50:
        try:
            mine = MINE(alpha=0.6, c=15)
            mine.compute_score(x, y)
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
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 9))
    
    ax.scatter(
        df['hashes_per_mb'], 
        df['diversity_per_mb'],
        alpha=0.3,
        s=20,
        c='steelblue',
        edgecolors='none',
        rasterized=True
    )
    
    # Regression line
    ax.plot(X, y_pred, 'r-', linewidth=2, label='Linear fit', zorder=10)
    
    # Labels
    ax.set_xlabel('k-mer hash density per megabase\n(31-mer FracMinHash, scale=1000)', 
                  fontsize=14, fontweight='bold')
    ax.set_ylabel('Alpha diversity per megabase', fontsize=14, fontweight='bold')
    ax.set_title(f'Correlation after filtering ({label})\nn = {len(df):,}',
                 fontsize=16, fontweight='bold', pad=20)
    
    # Extended stats box with all metrics
    stats_text = f'n = {len(df):,}\n\n'
    stats_text += 'Linear:\n'
    stats_text += f'  Pearson r = {pearson_r:.4f}\n'
    stats_text += f'  R² = {r2:.4f}\n'
    stats_text += f'  p = {pearson_p:.2e}\n\n'
    stats_text += 'Monotonic:\n'
    stats_text += f'  Spearman ρ = {spearman_r:.4f}\n'
    stats_text += f'  Kendall τ = {kendall_tau:.4f}\n\n'
    stats_text += 'Non-linear:\n'
    if not np.isnan(distance_corr):
        stats_text += f'  Distance corr = {distance_corr:.4f}\n'
    if not np.isnan(mic_score):
        stats_text += f'  MIC = {mic_score:.4f}'
    
    ax.text(
        0.05, 0.95, stats_text,
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
        family='monospace'
    )
    
    ax.legend(loc='lower right', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_file = output_dir / f"correlation_{label}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved plot to: {plot_file}")
    print(f"  Pearson r = {pearson_r:.4f} (linear)")
    print(f"  Spearman ρ = {spearman_r:.4f} (monotonic)")
    print(f"  Kendall τ = {kendall_tau:.4f} (monotonic)")
    if not np.isnan(distance_corr):
        print(f"  Distance corr = {distance_corr:.4f} (non-linear)")
    if not np.isnan(mic_score):
        print(f"  MIC = {mic_score:.4f} (non-linear)")
    print(f"  R² = {r2:.4f}")
    
    return {
        'n': len(df), 
        'pearson_r': pearson_r, 
        'spearman_r': spearman_r,
        'kendall_tau': kendall_tau,
        'distance_corr': distance_corr,
        'mic': mic_score,
        'r2': r2,
        'p': pearson_p
    }

def calculate_category_statistics(df: pd.DataFrame, category_col: str, 
                                  min_samples: int = 30) -> pd.DataFrame:
    """
    Calculate correlation statistics for each value within a categorical variable.
    
    Args:
        df: DataFrame with data and category column
        category_col: Name of the categorical column
        min_samples: Minimum number of samples required for statistical calculation
        
    Returns:
        DataFrame with statistics for each category value
    """
    from scipy.stats import pearsonr, spearmanr, kendalltau
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score
    import numpy as np
    
    # Import optional libraries with error handling
    try:
        import dcor
        has_dcor = True
    except ImportError:
        has_dcor = False
        print("  Warning: dcor not installed - distance correlation will be skipped")
    
    try:
        from minepy import MINE
        has_mine = True
    except ImportError:
        has_mine = False
        print("  Warning: minepy not installed - MIC will be skipped")
    
    if category_col not in df.columns:
        return pd.DataFrame()
    
    # Remove rows with missing category values
    df_clean = df[df[category_col].notna()].copy()
    
    if len(df_clean) == 0:
        return pd.DataFrame()
    
    results = []
    
    for category_value in df_clean[category_col].unique():
        category_data = df_clean[df_clean[category_col] == category_value]
        
        # Skip if too few samples
        if len(category_data) < min_samples:
            continue
        
        # Remove any remaining NaN values in the variables we're correlating
        category_data = category_data[
            category_data['hashes_per_mb'].notna() & 
            category_data['diversity_per_mb'].notna()
        ]
        
        if len(category_data) < min_samples:
            continue
        
        try:
            x = category_data['hashes_per_mb'].values
            y = category_data['diversity_per_mb'].values
            
            # Pearson correlation (linear)
            pearson_r, pearson_p = pearsonr(x, y)
            
            # Spearman correlation (monotonic, rank-based)
            spearman_r, spearman_p = spearmanr(x, y)
            
            # Kendall's tau (monotonic, rank-based, more robust to outliers)
            kendall_tau, kendall_p = kendalltau(x, y)
            
            # Distance correlation (captures non-linear dependencies)
            if has_dcor:
                try:
                    distance_corr = dcor.distance_correlation(x, y)
                except:
                    distance_corr = np.nan
            else:
                distance_corr = np.nan
            
            # Maximal Information Coefficient (non-linear associations)
            if has_mine and len(category_data) >= 50:  # MIC needs more samples
                try:
                    mine = MINE(alpha=0.6, c=15)
                    mine.compute_score(x, y)
                    mic_score = mine.mic()
                except:
                    mic_score = np.nan
            else:
                mic_score = np.nan
            
            # Linear regression for R²
            X = x.reshape(-1, 1)
            model = LinearRegression()
            model.fit(X, y)
            y_pred = model.predict(X)
            r2 = r2_score(y, y_pred)
            
            # Additional statistics
            slope = model.coef_[0]
            intercept = model.intercept_
            mean_hashes = x.mean()
            mean_diversity = y.mean()
            
            results.append({
                'variable': category_col,
                'category': str(category_value),
                'n_samples': len(category_data),
                # Linear correlation
                'pearson_r': pearson_r,
                'pearson_p': pearson_p,
                'r_squared': r2,
                # Rank/monotonic correlations
                'spearman_r': spearman_r,
                'spearman_p': spearman_p,
                'kendall_tau': kendall_tau,
                'kendall_p': kendall_p,
                # Non-linear associations
                'distance_corr': distance_corr,
                'mic': mic_score,
                # Regression parameters
                'slope': slope,
                'intercept': intercept,
                'mean_hashes_per_mb': mean_hashes,
                'mean_diversity_per_mb': mean_diversity
            })
            
        except Exception as e:
            # Skip categories that cause errors
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
                                      min_samples: int = 30) -> pd.DataFrame:
    """
    Calculate statistics for all categorical variables and combine.
    
    Args:
        df: DataFrame with data and metadata
        categorical_vars: List of categorical variable names
        min_samples: Minimum samples per category
        
    Returns:
        Combined DataFrame with all statistics
    """
    all_stats = []
    
    print("\nCalculating statistics for each category...")
    
    for var in categorical_vars:
        if var not in df.columns:
            continue
        
        print(f"  Processing: {var}")
        var_stats = calculate_category_statistics(df, var, min_samples)
        
        if not var_stats.empty:
            all_stats.append(var_stats)
            print(f"    Found {len(var_stats)} categories with ≥{min_samples} samples")
    
    if not all_stats:
        return pd.DataFrame()
    
    combined_stats = pd.concat(all_stats, ignore_index=True)
    
    # Apply multiple testing correction
    print(f"\nApplying multiple testing correction...")
    print(f"  Total comparisons: {len(combined_stats)}")
    combined_stats = apply_multiple_testing_correction(combined_stats)
    
    # Sort by R-squared
    combined_stats = combined_stats.sort_values('r_squared', ascending=False)
    
    return combined_stats


def create_multipanel_figure(df: pd.DataFrame, category_col: str, 
                            stats_df: pd.DataFrame, output_dir: Path,
                            max_panels: int = 12):
    """
    Create multi-panel figure showing each category separately.
    
    Args:
        df: DataFrame with data
        category_col: Categorical column name
        stats_df: Statistics for this variable's categories
        output_dir: Output directory
        max_panels: Maximum number of panels to show
    """
    from scipy.stats import pearsonr
    from sklearn.linear_model import LinearRegression
    import numpy as np
    
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
    
    print(f"\nCreating multi-panel figure for: {category_col}")
    print(f"  Panels: {n_categories} (showing most frequent categories)")
    
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
                           cat_data['diversity_per_mb'].notna()]
        
        if len(cat_data) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', 
                   transform=ax.transAxes)
            ax.set_title(str(category))
            continue
        
        # Scatter plot
        ax.scatter(cat_data['hashes_per_mb'], cat_data['diversity_per_mb'],
                  alpha=0.4, s=15, c='steelblue', edgecolors='none', rasterized=True)
        
        # Regression line
        X = cat_data['hashes_per_mb'].values.reshape(-1, 1)
        y = cat_data['diversity_per_mb'].values
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
        ax.set_ylabel('Diversity per Mb', fontsize=9)
        
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
    fig.suptitle(f'Hash-Diversity Correlation by {clean_name} (Most Frequent Categories)', 
                fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    
    # Save
    safe_name = category_col.replace('/', '_').replace(' ', '_')
    plot_file = output_dir / f"multipanel_{safe_name}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {plot_file.name}")


def generate_statistical_report(stats_df: pd.DataFrame, output_dir: Path, 
                                filters_applied: dict):
    """
    Generate comprehensive markdown report of statistical findings.
    
    Args:
        stats_df: DataFrame with all category statistics
        output_dir: Output directory
        filters_applied: Dictionary of filters used in analysis
    """
    report_file = output_dir / "statistical_report.md"
    
    with open(report_file, 'w') as f:
        f.write("# Statistical Analysis Report: Hash-Diversity Correlation by Category\n\n")
        f.write("## Analysis Overview\n\n")
        
        f.write(f"**Date:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("**Filters Applied:**\n")
        for key, val in filters_applied.items():
            f.write(f"- {key}: {val}\n")
        f.write("\n")
        
        f.write(f"**Total Categories Analyzed:** {len(stats_df)}\n")
        f.write(f"**Variables Analyzed:** {stats_df['variable'].nunique()}\n")
        f.write(f"**Total Samples Across All Categories:** {stats_df['n_samples'].sum():,}\n\n")
        
        # Multiple testing correction info
        f.write("## Multiple Testing Correction\n\n")
        f.write(f"Given {len(stats_df)} statistical comparisons, multiple testing correction is essential.\n\n")
        f.write("**Methods Applied:**\n")
        f.write("- **Bonferroni correction:** Controls family-wise error rate (FWER)\n")
        f.write("- **Benjamini-Hochberg FDR:** Controls false discovery rate (recommended for exploratory analysis)\n\n")
        f.write(f"**Significance Threshold:** α = 0.05\n\n")
        
        n_sig_bonf = stats_df['significant_bonferroni'].sum()
        n_sig_fdr = stats_df['significant_fdr'].sum()
        
        f.write(f"**Significant Results:**\n")
        f.write(f"- Bonferroni-corrected: {n_sig_bonf}/{len(stats_df)} ({100*n_sig_bonf/len(stats_df):.1f}%)\n")
        f.write(f"- FDR-corrected: {n_sig_fdr}/{len(stats_df)} ({100*n_sig_fdr/len(stats_df):.1f}%)\n\n")
        
        # Highest correlations
        f.write("## Highest Correlations (Top 10)\n\n")
        f.write("Categories with the strongest hash-diversity relationships:\n\n")
        
        # Sort by Spearman (monotonic) since user cares about "more hashes → more diversity"
        top_10_spearman = stats_df.nlargest(10, 'spearman_r')
        f.write("### By Spearman ρ (Monotonic Relationship)\n\n")
        f.write("| Rank | Variable | Category | n | Spearman ρ | Pearson r | Kendall τ | dCor | MIC | R² | p (FDR) |\n")
        f.write("|------|----------|----------|---|------------|-----------|-----------|------|-----|----|---------|\n")
        
        for idx, row in enumerate(top_10_spearman.itertuples(), 1):
            sig = '***' if row.pearson_p_fdr < 0.001 else '**' if row.pearson_p_fdr < 0.01 else '*' if row.pearson_p_fdr < 0.05 else 'ns'
            dcor_str = f"{row.distance_corr:.4f}" if not pd.isna(row.distance_corr) else "N/A"
            mic_str = f"{row.mic:.4f}" if not pd.isna(row.mic) else "N/A"
            f.write(f"| {idx} | {row.variable} | {row.category} | {row.n_samples:,} | {row.spearman_r:.4f} | {row.pearson_r:.4f} | {row.kendall_tau:.4f} | {dcor_str} | {mic_str} | {row.r_squared:.4f} | {row.pearson_p_fdr:.2e}{sig} |\n")
        
        f.write("\n")
        
        # Also show top by distance correlation if available
        if not stats_df['distance_corr'].isna().all():
            top_10_dcor = stats_df.nlargest(10, 'distance_corr')
            f.write("### By Distance Correlation (Non-linear Dependencies)\n\n")
            f.write("| Rank | Variable | Category | n | dCor | Spearman ρ | Pearson r | MIC | R² |\n")
            f.write("|------|----------|----------|---|------|------------|-----------|-----|----|\n")
            
            for idx, row in enumerate(top_10_dcor.itertuples(), 1):
                mic_str = f"{row.mic:.4f}" if not pd.isna(row.mic) else "N/A"
                f.write(f"| {idx} | {row.variable} | {row.category} | {row.n_samples:,} | {row.distance_corr:.4f} | {row.spearman_r:.4f} | {row.pearson_r:.4f} | {mic_str} | {row.r_squared:.4f} |\n")
            
            f.write("\n")
        
        # Lowest correlations
        f.write("## Lowest Correlations (Bottom 10)\n\n")
        f.write("Categories with the weakest hash-diversity relationships:\n\n")
        
        bottom_10 = stats_df.nsmallest(10, 'spearman_r')
        f.write("| Rank | Variable | Category | n | Spearman ρ | Pearson r | Kendall τ | R² | p (FDR) |\n")
        f.write("|------|----------|----------|---|------------|-----------|-----------|----|---------|\n")
        
        for idx, row in enumerate(bottom_10.itertuples(), 1):
            sig = '***' if row.pearson_p_fdr < 0.001 else '**' if row.pearson_p_fdr < 0.01 else '*' if row.pearson_p_fdr < 0.05 else 'ns'
            f.write(f"| {idx} | {row.variable} | {row.category} | {row.n_samples:,} | {row.spearman_r:.4f} | {row.pearson_r:.4f} | {row.kendall_tau:.4f} | {row.r_squared:.4f} | {row.pearson_p_fdr:.2e}{sig} |\n")
        
        f.write("\n")
        
        # Analysis by variable
        f.write("## Analysis by Metadata Variable\n\n")
        
        for var in stats_df['variable'].unique():
            var_stats = stats_df[stats_df['variable'] == var].sort_values('spearman_r', ascending=False)
            
            f.write(f"### {var.replace('_', ' ').title()}\n\n")
            f.write(f"**Categories analyzed:** {len(var_stats)}\n")
            f.write(f"**Mean Spearman ρ:** {var_stats['spearman_r'].mean():.4f}\n")
            f.write(f"**Mean Pearson r:** {var_stats['pearson_r'].mean():.4f}\n")
            f.write(f"**Mean R²:** {var_stats['r_squared'].mean():.4f}\n")
            f.write(f"**Spearman ρ range:** {var_stats['spearman_r'].min():.4f} - {var_stats['spearman_r'].max():.4f}\n")
            f.write(f"**Significant (FDR < 0.05):** {var_stats['significant_fdr'].sum()}/{len(var_stats)}\n\n")
            
            # Top categories for this variable
            f.write("**Top categories (by Spearman ρ):**\n\n")
            f.write("| Category | n | Spearman ρ | Pearson r | Kendall τ | R² | p (FDR) |\n")
            f.write("|----------|---|------------|-----------|-----------|----|---------|\n")
            
            for row in var_stats.head(5).itertuples():
                sig = '***' if row.pearson_p_fdr < 0.001 else '**' if row.pearson_p_fdr < 0.01 else '*' if row.pearson_p_fdr < 0.05 else 'ns'
                f.write(f"| {row.category} | {row.n_samples:,} | {row.spearman_r:.4f} | {row.pearson_r:.4f} | {row.kendall_tau:.4f} | {row.r_squared:.4f} | {row.pearson_p_fdr:.2e}{sig} |\n")
            
            f.write("\n")
        
        # Key findings
        f.write("## Key Findings\n\n")
        
        # Overall correlation strength
        mean_spearman = stats_df['spearman_r'].mean()
        mean_pearson = stats_df['pearson_r'].mean()
        mean_r2 = stats_df['r_squared'].mean()
        f.write(f"1. **Overall Mean Correlations:**\n")
        f.write(f"   - Spearman ρ (monotonic): {mean_spearman:.4f}\n")
        f.write(f"   - Pearson r (linear): {mean_pearson:.4f}\n")
        f.write(f"   - R² (linear fit): {mean_r2:.4f}\n\n")
        
        # Comparison of metrics
        if not stats_df['distance_corr'].isna().all():
            mean_dcor = stats_df['distance_corr'].mean()
            f.write(f"   - Distance correlation (non-linear): {mean_dcor:.4f}\n\n")
        
        # Check if monotonic relationship is stronger than linear
        if mean_spearman > mean_pearson + 0.05:
            f.write(f"2. **Monotonic vs Linear:** Spearman ρ exceeds Pearson r, suggesting the relationship ")
            f.write(f"is monotonic but not strictly linear. More hashes consistently predict more diversity, ")
            f.write(f"but the relationship may be non-linear.\n\n")
        elif abs(mean_spearman - mean_pearson) < 0.05:
            f.write(f"2. **Monotonic vs Linear:** Spearman ρ and Pearson r are similar, suggesting a ")
            f.write(f"strong linear relationship between hash density and diversity.\n\n")
        else:
            f.write(f"2. **Monotonic vs Linear:** Pearson r exceeds Spearman ρ, which is unusual and ")
            f.write(f"may warrant further investigation.\n\n")
        
        # Variability across categories
        spearman_std = stats_df['spearman_r'].std()
        f.write(f"3. **Variability:** Standard deviation of Spearman ρ = {spearman_std:.4f}, indicating ")
        if spearman_std > 0.1:
            f.write("substantial heterogeneity across categories.\n\n")
        else:
            f.write("relatively consistent correlations across categories.\n\n")
        
        # Best performing variable
        best_var = stats_df.groupby('variable')['spearman_r'].mean().idxmax()
        best_var_spearman = stats_df.groupby('variable')['spearman_r'].mean().max()
        f.write(f"4. **Best-Performing Variable:** {best_var} (mean Spearman ρ = {best_var_spearman:.4f})\n\n")
        
        # Identify categories with weak correlations
        weak_threshold = 0.3
        weak_cats = stats_df[stats_df['spearman_r'] < weak_threshold]
        f.write(f"5. **Weak Correlations:** {len(weak_cats)} categories have Spearman ρ < {weak_threshold}, ")
        f.write("suggesting potential quality issues or distinct biological patterns.\n\n")
        
        # Identify categories with strong correlations
        strong_threshold = 0.6
        strong_cats = stats_df[stats_df['spearman_r'] > strong_threshold]
        f.write(f"6. **Strong Correlations:** {len(strong_cats)} categories have Spearman ρ > {strong_threshold}, ")
        f.write("indicating robust monotonic hash-diversity relationships.\n\n")
        
        # Statistical significance
        f.write(f"7. **Statistical Significance:** After FDR correction, {n_sig_fdr} categories ")
        f.write(f"remain significant (p < 0.05), representing {100*n_sig_fdr/len(stats_df):.1f}% of comparisons.\n\n")
        
        # Recommendations
        f.write("## Recommendations\n\n")
        
        f.write("1. **Focus on Monotonic Relationships:** Since you're interested in \"more hashes → more diversity\", ")
        f.write("prioritize Spearman ρ and Kendall τ over Pearson r and R². These capture the monotonic ")
        f.write("relationship without assuming linearity.\n\n")
        
        f.write("2. **Sample Selection:** Consider filtering out categories with Spearman ρ < 0.3 for cleaner analysis\n\n")
        
        if len(weak_cats) > 0:
            f.write("3. **Investigate Weak Correlations:** The following categories warrant further investigation:\n")
            for row in weak_cats.head(5).itertuples():
                f.write(f"   - {row.variable}: {row.category} (Spearman ρ = {row.spearman_r:.4f}, n = {row.n_samples:,})\n")
            f.write("\n")
        
        f.write("4. **Non-linear Patterns:** If distance correlation substantially exceeds Pearson correlation, ")
        f.write("consider non-linear modeling approaches.\n\n")
        
        f.write("5. **Technical Artifacts:** Categories with weak correlations may indicate:\n")
        f.write("   - Sequencing/processing issues\n")
        f.write("   - Contamination\n")
        f.write("   - Distinct biological patterns requiring separate modeling\n\n")
        
        # Methods section
        f.write("## Statistical Methods\n\n")
        f.write("**Correlation Metrics:**\n\n")
        f.write("*Linear relationships:*\n")
        f.write("- **Pearson correlation coefficient (r):** Measures linear correlation\n")
        f.write("- **R² (coefficient of determination):** Proportion of variance explained by linear regression\n\n")
        
        f.write("*Monotonic relationships (recommended for \"more X → more Y\" questions):*\n")
        f.write("- **Spearman's rank correlation (ρ):** Measures monotonic relationships without assuming linearity\n")
        f.write("- **Kendall's tau (τ):** Another rank-based correlation, more robust to outliers\n\n")
        
        f.write("*Non-linear dependencies:*\n")
        f.write("- **Distance correlation:** Captures both linear and non-linear dependencies (0 = independent, 1 = dependent)\n")
        f.write("- **Maximal Information Coefficient (MIC):** Detects various functional relationships including non-linear patterns\n\n")
        
        f.write("**Why Multiple Metrics?**\n\n")
        f.write("Different metrics capture different aspects of relationships:\n")
        f.write("- If Spearman ≈ Pearson: Relationship is approximately linear\n")
        f.write("- If Spearman > Pearson: Relationship is monotonic but non-linear\n")
        f.write("- If Distance corr > Pearson: Non-linear dependencies present\n")
        f.write("- MIC detects complex functional relationships\n\n")
        
        f.write("**Multiple Testing Correction:**\n")
        f.write("- Bonferroni: p_corrected = p_raw × n_comparisons\n")
        f.write("- Benjamini-Hochberg FDR: Controls expected proportion of false discoveries\n\n")
        
        f.write("**Minimum Sample Size:** 30 samples per category (to ensure statistical power)\n\n")
        
        f.write("**Significance Levels:**\n")
        f.write("- *** : p < 0.001 (highly significant)\n")
        f.write("- **  : p < 0.01 (significant)\n")
        f.write("- *   : p < 0.05 (significant)\n")
        f.write("- ns  : p ≥ 0.05 (not significant)\n\n")
        
        # Output files
        f.write("## Output Files\n\n")
        f.write("- `category_statistics.csv` - Complete statistics for all categories\n")
        f.write("- `multipanel_*.png` - Multi-panel figures showing the most frequent categories (not highest R²)\n")
        f.write("- `categorical_plots/correlation_by_*.png` - Individual overlay plots showing top 10 by frequency\n")
        f.write("- `statistical_report.md` - This report\n\n")
        
        f.write("---\n\n")
        f.write("*Report generated automatically by analyze_parquet_data.py*\n")
    
    print(f"\nGenerated statistical report: {report_file}")
    
    # Also create a plain text version
    txt_file = output_dir / "statistical_report.txt"
    with open(report_file, 'r') as md_f:
        content = md_f.read()
        # Remove markdown formatting for plain text
        content = content.replace('#', '').replace('**', '').replace('|', ' ')
    with open(txt_file, 'w') as txt_f:
        txt_f.write(content)
    print(f"Generated plain text report: {txt_file}")
    """Extract biome information from attributes JSON."""
    import json
    try:
        if pd.isna(attributes_json):
            return 'unknown'
        attrs = json.loads(attributes_json) if isinstance(attributes_json, str) else attributes_json
        # Look for common biome-related fields
        for attr in attrs:
            tag = attr.get('tag', '').lower()
            if 'env_biome' in tag or 'biome' in tag:
                return attr.get('value', 'unknown')
            if 'environment' in tag and 'environmental' not in tag:
                return attr.get('value', 'unknown')
        return 'unknown'
    except:
        return 'unknown'


def parse_biome(attributes_json):
    """Extract biome information from attributes JSON."""
    import json
    try:
        if pd.isna(attributes_json):
            return 'unknown'
        attrs = json.loads(attributes_json) if isinstance(attributes_json, str) else attributes_json
        # Look for common biome-related fields
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
    """
    Bin continuous variables intelligently.
    
    Args:
        series: Pandas series to bin
        var_name: Variable name for appropriate binning strategy
    """
    if var_name in ['mbytes', 'mbases']:
        # Bin by orders of magnitude
        bins = [0, 10, 100, 1000, 10000, 100000, float('inf')]
        if var_name == 'mbytes':
            labels = ['<10 MB', '10-100 MB', '100 MB-1 GB', '1-10 GB', '10-100 GB', '>100 GB']
        else:
            labels = ['<10 Mb', '10-100 Mb', '100-1000 Mb', '1-10 Gb', '10-100 Gb', '>100 Gb']
        return pd.cut(series, bins=bins, labels=labels, include_lowest=True)
    
    elif var_name == 'avgspotlen':
        # Bin by read length categories
        bins = [0, 50, 100, 150, 250, 500, float('inf')]
        labels = ['<50 bp', '50-100 bp', '100-150 bp', '150-250 bp', '250-500 bp', '>500 bp']
        return pd.cut(series, bins=bins, labels=labels, include_lowest=True)
    
    elif var_name == 'releasedate':
        # Extract year from date
        try:
            years = pd.to_datetime(series, errors='coerce').dt.year
            return years.apply(lambda y: f"{int(y)}" if pd.notna(y) else 'unknown')
        except:
            return pd.Series(['unknown'] * len(series))
    
    else:
        # Default: quartile binning
        return pd.qcut(series, q=4, labels=['Q1', 'Q2', 'Q3', 'Q4'], duplicates='drop')


def consolidate_rare_categories(series: pd.Series, top_n: int = 10, other_label: str = 'other') -> pd.Series:
    """
    Group rare categories into 'other' category.
    
    Args:
        series: Pandas series with categorical data
        top_n: Number of top categories to keep
        other_label: Label for consolidated rare categories
    """
    value_counts = series.value_counts()
    top_categories = value_counts.head(top_n).index.tolist()
    
    return series.apply(lambda x: x if x in top_categories else other_label)


def plot_correlation_by_category(df: pd.DataFrame, category_col: str, output_dir: Path, 
                                 max_categories: int = 10):
    """
    Create correlation plot colored by a categorical variable.
    
    Args:
        df: DataFrame with data and category column
        category_col: Name of the categorical column
        output_dir: Directory to save plots
        max_categories: Maximum number of categories to show
    """
    print(f"\nCreating plot colored by: {category_col}")
    
    if category_col not in df.columns:
        print(f"  Column '{category_col}' not found - skipping")
        return
    
    # Remove rows with missing category values
    df_plot = df[df[category_col].notna()].copy()
    
    if len(df_plot) == 0:
        print(f"  No data available for '{category_col}' - skipping")
        return
    
    # Get category counts
    category_counts = df_plot[category_col].value_counts()
    
    if len(category_counts) == 0:
        print(f"  No categories found in '{category_col}' - skipping")
        return
    
    # If too many categories, keep top N and group rest as 'other'
    if len(category_counts) > max_categories:
        print(f"  {len(category_counts)} categories found, keeping top {max_categories}")
        df_plot[category_col] = consolidate_rare_categories(df_plot[category_col], 
                                                            top_n=max_categories)
        category_counts = df_plot[category_col].value_counts()
    
    # Create plot
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Use a colormap with enough distinct colors
    colors = sns.color_palette("tab20", n_colors=min(20, len(category_counts)))
    
    for idx, (category, color) in enumerate(zip(category_counts.index, colors)):
        category_data = df_plot[df_plot[category_col] == category]
        
        # Skip if no data
        if len(category_data) == 0:
            continue
        
        label = f"{category} (n={len(category_data):,})"
        
        ax.scatter(
            category_data['hashes_per_mb'],
            category_data['diversity_per_mb'],
            alpha=0.4,
            s=15,
            color=color,
            label=label,
            rasterized=True
        )
    
    ax.set_xlabel('k-mer hash density per megabase\n(31-mer FracMinHash, scale=1000)', 
                  fontsize=12, fontweight='bold')
    ax.set_ylabel('Alpha diversity per megabase', fontsize=12, fontweight='bold')
    
    # Clean up column name for title
    title_name = category_col.replace('_', ' ').title()
    ax.set_title(f'Hash-Diversity Correlation by {title_name}',
                 fontsize=14, fontweight='bold', pad=15)
    
    # Legend outside plot area
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Safe filename
    safe_name = category_col.replace('/', '_').replace(' ', '_')
    plot_file = output_dir / f"correlation_by_{safe_name}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {plot_file.name}")
    print(f"  Categories plotted: {len(category_counts)}")


def generate_all_categorical_plots(df: pd.DataFrame, output_dir: Path, 
                                  filters_applied: dict):
    """
    Generate correlation plots and statistics for all categorical metadata variables.
    
    Args:
        df: DataFrame with metadata joined
        output_dir: Directory to save plots
        filters_applied: Dictionary of filters used
    """
    print("\n" + "="*70)
    print("GENERATING CATEGORICAL ANALYSIS")
    print("="*70)
    
    # Create subdirectory for categorical plots
    cat_plot_dir = output_dir / "categorical_plots"
    cat_plot_dir.mkdir(exist_ok=True)
    
    # List of categorical variables to plot
    categorical_vars = [
        'center_name',
        'instrument', 
        'librarylayout',
        'libraryselection',
        'librarysource',
        'platform',
        'organism'
    ]
    
    # Variables that need special handling
    # Parse biome from attributes only if database biome is NULL
    if 'biome' in df.columns:
        # Database has biome column - only parse from attributes for NULL values
        print("\nFilling missing biome values from attributes...")
        if 'attributes' in df.columns:
            missing_biome = df['biome'].isna()
            df.loc[missing_biome, 'biome'] = df.loc[missing_biome, 'attributes'].apply(parse_biome)
            print(f"  Biome filled for {missing_biome.sum():,} samples")
        categorical_vars.append('biome')
    elif 'attributes' in df.columns:
        # Database doesn't have biome - parse from attributes
        print("\nParsing biome from attributes...")
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
            print(f"Binning {var_name}...")
            df[f'{var_name}_binned'] = bin_continuous_variable(df[col_name], var_name)
            categorical_vars.append(f'{var_name}_binned')
    
    # Calculate statistics for all categories
    all_stats = calculate_all_category_statistics(df, categorical_vars, min_samples=30)
    
    if all_stats.empty:
        print("\nWARNING: No categories had sufficient samples for analysis")
        return
    
    # Save statistics to CSV
    stats_file = output_dir / "category_statistics.csv"
    all_stats.to_csv(stats_file, index=False)
    print(f"\nSaved category statistics to: {stats_file}")
    print(f"  Total categories analyzed: {len(all_stats)}")
    print(f"  Significant after FDR correction: {all_stats['significant_fdr'].sum()}")
    
    # Generate overlay plots for each categorical variable
    print("\n" + "-"*70)
    print("Generating overlay plots...")
    print("-"*70)
    for var in categorical_vars:
        plot_correlation_by_category(df, var, cat_plot_dir)
    
    # Generate multi-panel figures for each categorical variable
    print("\n" + "-"*70)
    print("Generating multi-panel figures...")
    print("-"*70)
    for var in categorical_vars:
        create_multipanel_figure(df, var, all_stats, cat_plot_dir, max_panels=12)
    
    # Generate comprehensive statistical report
    print("\n" + "-"*70)
    print("Generating statistical report...")
    print("-"*70)
    generate_statistical_report(all_stats, output_dir, filters_applied)
    
    print("\n" + "="*70)
    print(f"CATEGORICAL ANALYSIS COMPLETE")
    print(f"Results saved to: {output_dir}/")
    print("="*70)

def main():
    parser = argparse.ArgumentParser(
        description="Analyze parquet data with custom filters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', '-i', required=True, help='Input parquet file')
    parser.add_argument('--output', '-o', required=True, help='Output directory')
    parser.add_argument('--min-hashes', type=int, default=1000, 
                       help='Minimum total distinct hashes')
    parser.add_argument('--min-diversity', type=int, default=1,
                       help='Minimum alpha diversity')
    parser.add_argument('--metadata-db', 
                       default='/scratch/shared_data_new/Logan_yacht_data/metadata/aws_sra_metadata/metadata_geo_joined.duckdb',
                       help='Path to metadata database')
    parser.add_argument('--join-metadata', action='store_true',
                       help='Join with metadata database')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*70)
    print("DOWNSTREAM ANALYSIS")
    print("="*70)
    
    # Load data
    df = load_parquet_data(Path(args.input))
    
    # Apply filters
    df_filtered = apply_filters(df, args.min_hashes, args.min_diversity)
    
    # Join metadata if requested
    if args.join_metadata:
        df_filtered = join_metadata(df_filtered, args.metadata_db)
    
    # Create plots
    stats = plot_filtered_correlation(df_filtered, output_dir, 
                                     label=f"min{args.min_hashes}hashes")
    
    if args.join_metadata:
        # Prepare filter info for report
        filters_applied = {
            'min_hashes': args.min_hashes,
            'min_diversity': args.min_diversity,
            'min_mbases': 100,
            'n_samples_after_filtering': len(df_filtered),
            'input_file': args.input
        }
        generate_all_categorical_plots(df_filtered, output_dir, filters_applied)
    
    # Save filtered data
    output_file = output_dir / "filtered_data.parquet"
    df_filtered.to_parquet(output_file, index=False)
    print(f"\nSaved filtered data to: {output_file}")
    
    # Save summary
    summary_file = output_dir / "summary.txt"
    with open(summary_file, 'w') as f:
        f.write("FILTERED DATA SUMMARY\n")
        f.write("="*70 + "\n\n")
        f.write(f"Input file: {args.input}\n")
        f.write(f"Filters:\n")
        f.write(f"  Minimum hashes: {args.min_hashes:,}\n")
        f.write(f"  Minimum diversity: {args.min_diversity}\n\n")
        f.write(f"Results:\n")
        f.write(f"  Samples after filtering: {stats['n']:,}\n")
        f.write(f"  Pearson r: {stats['r']:.4f}\n")
        f.write(f"  p-value: {stats['p']:.2e}\n")
        f.write(f"  R²: {stats['r2']:.4f}\n")
    
    print(f"Saved summary to: {summary_file}")
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)

if __name__ == "__main__":
    main()
