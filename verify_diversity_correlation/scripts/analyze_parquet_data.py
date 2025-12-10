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
    from scipy.stats import pearsonr, spearmanr
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score
    
    print(f"\nCreating correlation plot for {label} data...")
    
    # Calculate statistics
    pearson_r, pearson_p = pearsonr(df['hashes_per_mb'], df['diversity_per_mb'])
    spearman_r, spearman_p = spearmanr(df['hashes_per_mb'], df['diversity_per_mb'])
    
    X = df['hashes_per_mb'].values.reshape(-1, 1)
    y = df['diversity_per_mb'].values
    model = LinearRegression()
    model.fit(X, y)
    y_pred = model.predict(X)
    r2 = r2_score(y, y_pred)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
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
    
    # Stats box
    stats_text = (
        f'n = {len(df):,}\n'
        f'Pearson r = {pearson_r:.4f}\n'
        f'p-value = {pearson_p:.2e}\n'
        f'R² = {r2:.4f}\n'
        f'Spearman ρ = {spearman_r:.4f}'
    )
    
    ax.text(
        0.05, 0.95, stats_text,
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    )
    
    ax.legend(loc='lower right', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_file = output_dir / f"correlation_{label}.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved plot to: {plot_file}")
    print(f"  Pearson r = {pearson_r:.4f}")
    print(f"  R² = {r2:.4f}")
    
    return {'n': len(df), 'r': pearson_r, 'p': pearson_p, 'r2': r2}

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


def generate_all_categorical_plots(df: pd.DataFrame, output_dir: Path):
    """
    Generate correlation plots for all categorical metadata variables.
    
    Args:
        df: DataFrame with metadata joined
        output_dir: Directory to save plots
    """
    print("\n" + "="*70)
    print("GENERATING CATEGORICAL PLOTS")
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
    # Parse biome from attributes
    if 'attributes' in df.columns:
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
            print(f"\nBinning {var_name}...")
            df[f'{var_name}_binned'] = bin_continuous_variable(df[col_name], var_name)
            categorical_vars.append(f'{var_name}_binned')
    
    # Generate plots for each categorical variable
    for var in categorical_vars:
        plot_correlation_by_category(df, var, cat_plot_dir)
    
    print("\n" + "="*70)
    print(f"CATEGORICAL PLOTS COMPLETE")
    print(f"Saved to: {cat_plot_dir}/")
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
        generate_all_categorical_plots(df_filtered, output_dir)
    
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
