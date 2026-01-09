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
from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

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


def load_parquet_data(parquet_file: Path) -> pd.DataFrame:
    """Load the parquet file."""
    print(f"Loading data from: {parquet_file}")
    df = pd.read_parquet(parquet_file)
    print(f"Loaded {len(df):,} samples")
    print(f"Columns: {df.columns.tolist()}")
    return df


def apply_filters(df: pd.DataFrame, min_hashes: int = 1000, min_diversity: int = 1, 
                  min_mbases: int = 0) -> pd.DataFrame:
    """Apply quality filters to the data."""
    print(f"\nApplying filters:")
    print(f"  Minimum total hashes: {min_hashes:,}")
    print(f"  Minimum observed richness (KOs): {min_diversity}")
    if min_mbases > 0:
        print(f"  Minimum mbases: {min_mbases:,}")
    
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
        print(f"  WARNING: mbases column not found, skipping mbases filter")
    
    df_filtered = df[filter_mask].copy()
    
    filtered_count = len(df_filtered)
    print(f"\nFiltering results:")
    print(f"  Original samples: {original_count:,}")
    print(f"  After filtering: {filtered_count:,}")
    print(f"  Removed: {original_count - filtered_count:,} ({100*(original_count-filtered_count)/original_count:.1f}%)")
    
    return df_filtered


def join_metadata(df: pd.DataFrame, metadata_db: str) -> pd.DataFrame:
    """Join with metadata database to get additional sample information."""
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


def calculate_metrics(x: np.ndarray, y: np.ndarray) -> dict:
    """Calculate all correlation metrics."""
    n = len(x)
    
    pearson_r, pearson_p = pearsonr(x, y)
    spearman_r, spearman_p = spearmanr(x, y)
    kendall_tau, kendall_p = kendalltau(x, y)
    
    if HAS_DCOR:
        try:
            distance_corr = dcor.distance_correlation(x, y, method='AVL')
        except:
            distance_corr = np.nan
    else:
        distance_corr = np.nan
    
    if HAS_MINE and n >= 50:
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


def plot_filtered_correlation(df: pd.DataFrame, output_dir: Path, label: str = "filtered"):
    """Create correlation plots for all diversity metrics."""
    
    print(f"\nCreating correlation plots for {label} data...")
    
    all_stats = {}
    
    for metric_col, metric_name in DIVERSITY_METRICS:
        if metric_col not in df.columns:
            print(f"  Skipping {metric_name}: column not found")
            continue
        
        df_clean = df.dropna(subset=['hashes_per_mb', metric_col])
        
        if len(df_clean) < 10:
            print(f"  Skipping {metric_name}: insufficient data (n={len(df_clean)})")
            continue
        
        print(f"\n  Processing: {metric_name}")
        
        x = df_clean['hashes_per_mb'].values
        y = df_clean[metric_col].values
        
        metrics = calculate_metrics(x, y)
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
        
        print(f"    Saved: {plot_file.name}")
        print(f"    Spearman ρ = {metrics['spearman_r']:.4f}")
        print(f"    Pearson r = {metrics['pearson_r']:.4f}")
        print(f"    R² = {metrics['r2']:.4f}")
    
    # Create summary multi-panel figure
    create_summary_panel(df, output_dir, label, all_stats)
    
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
    print(f"\nSaved summary figure: summary_all_metrics_{label}.png")


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
                                  diversity_col: str, min_samples: int = 30) -> pd.DataFrame:
    """Calculate correlation statistics for each value within a categorical variable."""
    
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
            
            metrics = calculate_metrics(x, y)
            
            results.append({
                'variable': category_col,
                'category': str(category_value),
                'diversity_metric': diversity_col,
                'n_samples': len(category_data),
                **metrics
            })
            
        except Exception:
            continue
    
    if not results:
        return pd.DataFrame()
    
    return pd.DataFrame(results)


def plot_correlation_by_category(df: pd.DataFrame, category_col: str, diversity_col: str,
                                 output_dir: Path, max_categories: int = 10):
    """Create correlation plot colored by a categorical variable."""
    
    print(f"\nCreating plot colored by: {category_col} for {diversity_col}")
    
    if category_col not in df.columns or diversity_col not in df.columns:
        print(f"  Column not found - skipping")
        return
    
    df_plot = df[df[category_col].notna()].copy()
    
    if len(df_plot) == 0:
        print(f"  No data available - skipping")
        return
    
    category_counts = df_plot[category_col].value_counts()
    
    if len(category_counts) > max_categories:
        print(f"  {len(category_counts)} categories found, keeping top {max_categories}")
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
    
    print(f"  Saved: {plot_file.name}")


def generate_all_categorical_plots(df: pd.DataFrame, output_dir: Path, filters_applied: dict):
    """Generate correlation plots for all categorical metadata variables and diversity metrics."""
    
    print("\n" + "="*70)
    print("GENERATING CATEGORICAL ANALYSIS")
    print("="*70)
    
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
        print("\nFilling missing biome values from attributes...")
        if 'attributes' in df.columns:
            missing_biome = df['biome'].isna()
            df.loc[missing_biome, 'biome'] = df.loc[missing_biome, 'attributes'].apply(parse_biome)
            print(f"  Biome filled for {missing_biome.sum():,} samples")
        categorical_vars.append('biome')
    elif 'attributes' in df.columns:
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
    
    # Generate plots for each diversity metric and categorical variable combination
    # Focus on the main diversity metric: observed_richness_per_mb
    primary_metric = 'observed_richness_per_mb'
    if primary_metric in df.columns:
        print(f"\nGenerating categorical plots for: {primary_metric}")
        for var in categorical_vars:
            plot_correlation_by_category(df, var, primary_metric, cat_plot_dir)
    
    # Also generate for Shannon index as it's commonly used
    shannon_metric = 'shannon_index'
    if shannon_metric in df.columns:
        print(f"\nGenerating categorical plots for: {shannon_metric}")
        for var in categorical_vars:
            plot_correlation_by_category(df, var, shannon_metric, cat_plot_dir)
    
    print("\n" + "="*70)
    print(f"CATEGORICAL ANALYSIS COMPLETE")
    print(f"Results saved to: {cat_plot_dir}/")
    print("="*70)


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
    
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*70)
    print("FUNCTIONAL DIVERSITY DOWNSTREAM ANALYSIS")
    print("="*70)
    
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
    stats = plot_filtered_correlation(df_filtered, output_dir, label=label)
    
    if args.join_metadata:
        filters_applied = {
            'min_hashes': args.min_hashes,
            'min_diversity': args.min_diversity,
            'min_mbases': args.min_mbases,
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
    
    print(f"Saved summary to: {summary_file}")
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)


if __name__ == "__main__":
    main()
