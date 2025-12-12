#!/usr/bin/env python3
"""
Custom Correlation Plot Generator
==================================
Generate correlation plots with flexible filtering and categorical coloring.

Features:
- Load pre-computed parquet files
- Filter by arbitrary metadata (equality, inequality, comparisons)
- Color points by categorical variables
- Compute all correlation metrics
- Save publication-quality plots

Filter Operators:
- = (equals)
- != (not equals)
- > (greater than)
- >= (greater than or equal)
- < (less than)
- <= (less than or equal)

Examples:
    # Basic plot for soil metagenome samples from JGI
    python3 plot_custom_correlation.py \\
        --input hash_diversity_data.parquet \\
        --output soil_jgi_analysis \\
        --filter "organism=soil metagenome" \\
        --filter "center_name=JGI"
    
    # Filter by file size (large samples only)
    python3 plot_custom_correlation.py \\
        --input hash_diversity_data.parquet \\
        --output large_samples \\
        --filter "mbytes>1024" \\
        --color-by platform
    
    # Filter by sequencing depth range
    python3 plot_custom_correlation.py \\
        --input hash_diversity_data.parquet \\
        --output medium_depth \\
        --filter "mbases>=100" \\
        --filter "mbases<=1000" \\
        --color-by organism
    
    # Multiple criteria with comparison operators
    python3 plot_custom_correlation.py \\
        --input filtered_data.parquet \\
        --output high_quality \\
        --filter "platform=ILLUMINA" \\
        --filter "mbases>500" \\
        --filter "alpha_diversity>10" \\
        --color-by organism \\
        --top-n 8

Usage:
    python3 plot_custom_correlation.py --input FILE --output DIR [OPTIONS]
"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr, kendalltau
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

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

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.3)


def parse_filter(filter_str: str) -> tuple:
    """
    Parse filter string with various operators.
    
    Supports:
        column=value       (equality)
        column!=value      (inequality)
        column>value       (greater than)
        column>=value      (greater than or equal)
        column<value       (less than)
        column<=value      (less than or equal)
    
    Returns:
        (column, operator, value)
    """
    # Check operators in order of precedence (longer operators first)
    if '>=' in filter_str:
        column, value = filter_str.split('>=', 1)
        return column.strip(), '>=', value.strip()
    elif '<=' in filter_str:
        column, value = filter_str.split('<=', 1)
        return column.strip(), '<=', value.strip()
    elif '!=' in filter_str:
        column, value = filter_str.split('!=', 1)
        return column.strip(), '!=', value.strip()
    elif '>' in filter_str:
        column, value = filter_str.split('>', 1)
        return column.strip(), '>', value.strip()
    elif '<' in filter_str:
        column, value = filter_str.split('<', 1)
        return column.strip(), '<', value.strip()
    elif '=' in filter_str:
        column, value = filter_str.split('=', 1)
        return column.strip(), '=', value.strip()
    else:
        raise ValueError(f"Invalid filter format: {filter_str}. Use operators: =, !=, >, >=, <, <=")


def apply_filters(df: pd.DataFrame, filters: list) -> pd.DataFrame:
    """
    Apply multiple filters to dataframe.
    
    Args:
        df: Input dataframe
        filters: List of (column, operator, value) tuples
        
    Returns:
        Filtered dataframe
    """
    df_filtered = df.copy()
    
    print("\nApplying filters:")
    for column, operator, value in filters:
        if column not in df_filtered.columns:
            print(f"  WARNING: Column '{column}' not found in data - skipping")
            continue
        
        before_count = len(df_filtered)
        
        # Try to convert value to numeric for comparison operators
        try:
            numeric_value = float(value)
            is_numeric = True
        except ValueError:
            is_numeric = False
            numeric_value = None
        
        # Apply filter based on operator
        if operator == '=':
            if is_numeric:
                df_filtered = df_filtered[df_filtered[column] == numeric_value]
            else:
                df_filtered = df_filtered[df_filtered[column] == value]
        
        elif operator == '!=':
            if is_numeric:
                df_filtered = df_filtered[df_filtered[column] != numeric_value]
            else:
                df_filtered = df_filtered[df_filtered[column] != value]
        
        elif operator == '>':
            if not is_numeric:
                print(f"  WARNING: Operator '>' requires numeric value, got '{value}' - skipping")
                continue
            df_filtered = df_filtered[df_filtered[column] > numeric_value]
        
        elif operator == '>=':
            if not is_numeric:
                print(f"  WARNING: Operator '>=' requires numeric value, got '{value}' - skipping")
                continue
            df_filtered = df_filtered[df_filtered[column] >= numeric_value]
        
        elif operator == '<':
            if not is_numeric:
                print(f"  WARNING: Operator '<' requires numeric value, got '{value}' - skipping")
                continue
            df_filtered = df_filtered[df_filtered[column] < numeric_value]
        
        elif operator == '<=':
            if not is_numeric:
                print(f"  WARNING: Operator '<=' requires numeric value, got '{value}' - skipping")
                continue
            df_filtered = df_filtered[df_filtered[column] <= numeric_value]
        
        after_count = len(df_filtered)
        print(f"  {column} {operator} {value}: {before_count:,} → {after_count:,} samples")
    
    return df_filtered


def consolidate_categories(df: pd.DataFrame, column: str, top_n: int = 10) -> pd.DataFrame:
    """
    Consolidate rare categories into 'other'.
    
    Args:
        df: Input dataframe
        column: Column to consolidate
        top_n: Number of top categories to keep
        
    Returns:
        Dataframe with consolidated categories
    """
    if column not in df.columns:
        return df
    
    value_counts = df[column].value_counts()
    top_categories = value_counts.head(top_n).index.tolist()
    
    df_copy = df.copy()
    df_copy[f'{column}_consolidated'] = df_copy[column].apply(
        lambda x: x if x in top_categories else 'other'
    )
    
    return df_copy


def calculate_metrics(x: np.ndarray, y: np.ndarray) -> dict:
    """
    Calculate all correlation metrics.
    
    Args:
        x: Independent variable (hashes per mb)
        y: Dependent variable (diversity per mb)
        
    Returns:
        Dictionary of metrics
    """
    n = len(x)
    
    # Linear correlations
    pearson_r, pearson_p = pearsonr(x, y)
    
    # Monotonic correlations
    spearman_r, spearman_p = spearmanr(x, y)
    kendall_tau, kendall_p = kendalltau(x, y)
    
    # Distance correlation
    if HAS_DCOR:
        try:
            distance_corr = dcor.distance_correlation(x, y)
        except:
            distance_corr = np.nan
    else:
        distance_corr = np.nan
    
    # MIC
    if HAS_MINE and n >= 50:
        try:
            mine = MINE(alpha=0.6, c=15)
            mine.compute_score(x, y)
            mic_score = mine.mic()
        except:
            mic_score = np.nan
    else:
        mic_score = np.nan
    
    # Linear regression
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


def create_plot(df: pd.DataFrame, metrics: dict, output_dir: Path, 
                base_name: str = 'custom', color_by: str = None, 
                filters_applied: list = None):
    """
    Create scatter plot with optional categorical coloring.
    
    Args:
        df: Dataframe with data
        metrics: Correlation metrics
        output_dir: Output directory
        color_by: Optional column for coloring
        filters_applied: List of filters applied
    """
    fig, ax = plt.subplots(figsize=(14, 10))
    
    x = df['hashes_per_mb'].values
    y = df['diversity_per_mb'].values
    
    if color_by and f'{color_by}_consolidated' in df.columns:
        # Colored by category
        # Sort categories by frequency (most frequent first)
        cat_col = f'{color_by}_consolidated'
        value_counts = df[cat_col].value_counts()
        categories = value_counts.index.tolist()  # Already sorted by frequency
        
        # Use maximally distinguishable colors
        # Combine multiple palettes to ensure distinct colors (no similar blues, etc.)
        distinct_colors = [
            '#e41a1c',  # red
            '#ff7f00',  # orange
            '#4daf4a',  # green
            '#984ea3',  # purple
            '#f781bf',  # pink
            '#a65628',  # brown
            '#ffff33',  # yellow
            '#377eb8',  # blue (only one blue)
            '#999999',  # gray
            '#e41a1c',  # cycling back with different saturation
            '#66c2a5',  # teal (only one teal-ish color)
            '#fc8d62',  # light orange
            '#8da0cb',  # light purple
            '#e78ac3',  # light pink
            '#a6d854',  # lime
            '#ffd92f',  # gold
            '#e5c494',  # tan
            '#b3b3b3',  # light gray
            '#8dd3c7',  # cyan
            '#bebada',  # lavender
        ]
        
        # Plot in order of frequency (most frequent first)
        # This ensures less frequent points are plotted on top
        for idx, category in enumerate(categories):
            mask = df[cat_col] == category
            cat_data = df[mask]
            
            label = f"{category} (n={len(cat_data):,})"
            
            ax.scatter(
                cat_data['hashes_per_mb'],
                cat_data['diversity_per_mb'],
                alpha=0.6,
                s=8,  # Smaller points to reduce overlap
                color=distinct_colors[idx % len(distinct_colors)],
                label=label,
                edgecolors='none',
                rasterized=True
            )
        
        # Add legend (categories already in frequency order)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
        
    else:
        # Single color
        ax.scatter(
            x, y,
            alpha=0.5,
            s=8,  # Smaller points to reduce overlap
            c='steelblue',
            edgecolors='none',
            rasterized=True
        )
    
    # Add regression line
    X_sorted = np.sort(x)
    y_pred = metrics['slope'] * X_sorted + metrics['intercept']
    ax.plot(X_sorted, y_pred, 'r-', linewidth=2.5, label='Linear fit', zorder=100)
    
    # Labels
    ax.set_xlabel('k-mer hash density per megabase\n(31-mer FracMinHash, scale=1000)',
                  fontsize=14, fontweight='bold')
    ax.set_ylabel('Alpha diversity per megabase',
                  fontsize=14, fontweight='bold')
    
    # Title
    title = 'Hash-Diversity Correlation'
    if filters_applied:
        # Format each filter with its actual operator
        filter_strs = []
        for c, op, v in filters_applied[:3]:
            filter_strs.append(f'{c}{op}{v}')
        title += '\n' + ' & '.join(filter_strs)
        if len(filters_applied) > 3:
            title += f' (+{len(filters_applied)-3} more)'
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    
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
        stats_text += f'  dCor = {metrics["distance_corr"]:.4f}\n'
    if not np.isnan(metrics['mic']):
        stats_text += f'  MIC = {metrics["mic"]:.4f}'
    
    # Position based on whether there's a legend
    if color_by:
        # Put stats box in upper left if legend is on right
        ax.text(0.02, 0.98, stats_text,
                transform=ax.transAxes,
                fontsize=11,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9),
                family='monospace')
    else:
        # Put stats box in upper left
        ax.text(0.05, 0.95, stats_text,
                transform=ax.transAxes,
                fontsize=11,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                family='monospace')
    
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save
    plot_file = output_dir / f"{base_name}_correlation.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nSaved plot to: {plot_file}")


def save_metrics(metrics: dict, output_dir: Path, base_name: str = 'custom',
                 filters_applied: list = None):
    """Save metrics to CSV and text files."""
    
    # CSV
    csv_file = output_dir / f"{base_name}_metrics.csv"
    metrics_df = pd.DataFrame([metrics])
    metrics_df.to_csv(csv_file, index=False)
    print(f"Saved metrics to: {csv_file}")
    
    # Text report
    txt_file = output_dir / f"{base_name}_report.txt"
    with open(txt_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("CUSTOM CORRELATION ANALYSIS REPORT\n")
        f.write("="*70 + "\n\n")
        
        if filters_applied:
            f.write("FILTERS APPLIED:\n")
            f.write("-"*70 + "\n")
            for column, operator, value in filters_applied:
                f.write(f"  {column} {operator} {value}\n")
            f.write("\n")
        
        f.write("SAMPLE SIZE:\n")
        f.write("-"*70 + "\n")
        f.write(f"  n = {metrics['n']:,} samples\n\n")
        
        f.write("CORRELATION METRICS:\n")
        f.write("-"*70 + "\n\n")
        
        f.write("Linear Correlations:\n")
        f.write(f"  Pearson r: {metrics['pearson_r']:.6f}\n")
        f.write(f"    p-value: {metrics['pearson_p']:.4e}\n")
        f.write(f"    Significance: {'***' if metrics['pearson_p'] < 0.001 else '**' if metrics['pearson_p'] < 0.01 else '*' if metrics['pearson_p'] < 0.05 else 'ns'}\n")
        f.write(f"  R²: {metrics['r2']:.6f}\n")
        f.write(f"    Variance explained: {metrics['r2']*100:.2f}%\n")
        f.write(f"  Regression: y = {metrics['slope']:.8f}x + {metrics['intercept']:.8f}\n\n")
        
        f.write("Monotonic Correlations:\n")
        f.write(f"  Spearman ρ: {metrics['spearman_r']:.6f}\n")
        f.write(f"    p-value: {metrics['spearman_p']:.4e}\n")
        f.write(f"    Significance: {'***' if metrics['spearman_p'] < 0.001 else '**' if metrics['spearman_p'] < 0.01 else '*' if metrics['spearman_p'] < 0.05 else 'ns'}\n")
        f.write(f"  Kendall τ: {metrics['kendall_tau']:.6f}\n")
        f.write(f"    p-value: {metrics['kendall_p']:.4e}\n\n")
        
        if not np.isnan(metrics['distance_corr']) or not np.isnan(metrics['mic']):
            f.write("Non-linear Dependencies:\n")
            if not np.isnan(metrics['distance_corr']):
                f.write(f"  Distance correlation: {metrics['distance_corr']:.6f}\n")
            if not np.isnan(metrics['mic']):
                f.write(f"  MIC: {metrics['mic']:.6f}\n")
            f.write("\n")
        
        f.write("-"*70 + "\n")
        f.write("INTERPRETATION:\n")
        f.write("-"*70 + "\n\n")
        
        spearman_r = abs(metrics['spearman_r'])
        if spearman_r > 0.7:
            strength = "strong"
        elif spearman_r > 0.4:
            strength = "moderate"
        elif spearman_r > 0.2:
            strength = "weak"
        else:
            strength = "very weak"
        
        direction = "positive" if metrics['spearman_r'] > 0 else "negative"
        
        f.write(f"The monotonic correlation is {strength} and {direction}.\n")
        f.write(f"More hashes {'consistently predict' if spearman_r > 0.5 else 'modestly predict'} ")
        f.write(f"{'higher' if metrics['spearman_r'] > 0 else 'lower'} diversity.\n\n")
        
        if abs(metrics['spearman_r'] - metrics['pearson_r']) > 0.1:
            if metrics['spearman_r'] > metrics['pearson_r']:
                f.write(f"Note: Spearman ρ ({metrics['spearman_r']:.3f}) > Pearson r ({metrics['pearson_r']:.3f})\n")
                f.write(f"This suggests a monotonic but non-linear relationship (e.g., diminishing returns).\n")
    
    print(f"Saved report to: {txt_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate custom correlation plots with filtering and coloring",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input parquet file (hash_diversity_data.parquet or filtered_data.parquet)'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output directory'
    )
    parser.add_argument(
        '--base-name', '-b',
        default='custom',
        help='Base name for output files (default: custom)'
    )
    parser.add_argument(
        '--filter', '-f',
        action='append',
        default=[],
        help='Filter condition using operators: =, !=, >, >=, <, <= (e.g., "mbytes>1024" or "organism=soil metagenome"). Can be used multiple times'
    )
    parser.add_argument(
        '--color-by', '-c',
        default=None,
        help='Categorical column to color points by (e.g., organism, platform, releasedate_binned)'
    )
    parser.add_argument(
        '--top-n', '-n',
        type=int,
        default=10,
        help='Number of top categories to show when coloring (rest grouped as "other")'
    )
    parser.add_argument(
        '--min-samples',
        type=int,
        default=10,
        help='Minimum number of samples required for analysis'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*70)
    print("CUSTOM CORRELATION ANALYSIS")
    print("="*70)
    
    # Load data
    print(f"\nLoading data from: {args.input}")
    df = pd.read_parquet(args.input)
    print(f"Loaded {len(df):,} samples")
    print(f"Columns available: {', '.join(df.columns.tolist())}")
    
    # Parse and apply filters
    filters = []
    if args.filter:
        for filter_str in args.filter:
            filters.append(parse_filter(filter_str))
        
        df = apply_filters(df, filters)
        
        if len(df) == 0:
            print("\nERROR: No samples remaining after filtering!")
            return
    
    # Check for required columns
    required_cols = ['hashes_per_mb', 'diversity_per_mb']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"\nERROR: Required columns missing: {missing_cols}")
        print("This script requires pre-computed normalized metrics.")
        return
    
    # Remove NaN values
    df_clean = df.dropna(subset=['hashes_per_mb', 'diversity_per_mb'])
    
    if len(df_clean) < args.min_samples:
        print(f"\nERROR: Only {len(df_clean)} samples with complete data (minimum: {args.min_samples})")
        return
    
    print(f"\nSamples with complete data: {len(df_clean):,}")
    
    # Handle coloring
    if args.color_by:
        if args.color_by not in df_clean.columns:
            print(f"\nWARNING: Column '{args.color_by}' not found - generating plot without coloring")
            args.color_by = None
        else:
            print(f"\nConsolidating '{args.color_by}' to top {args.top_n} categories...")
            df_clean = consolidate_categories(df_clean, args.color_by, args.top_n)
            
            # Print category distribution
            cat_col = f'{args.color_by}_consolidated'
            value_counts = df_clean[cat_col].value_counts()
            print(f"Categories ({len(value_counts)}):")
            for cat, count in value_counts.items():
                print(f"  {cat}: {count:,} samples ({100*count/len(df_clean):.1f}%)")
    
    # Calculate metrics
    print("\nCalculating correlation metrics...")
    x = df_clean['hashes_per_mb'].values
    y = df_clean['diversity_per_mb'].values
    metrics = calculate_metrics(x, y)
    
    print("\nMetrics:")
    print(f"  Spearman ρ = {metrics['spearman_r']:.4f} (p={metrics['spearman_p']:.2e})")
    print(f"  Pearson r = {metrics['pearson_r']:.4f} (p={metrics['pearson_p']:.2e})")
    print(f"  Kendall τ = {metrics['kendall_tau']:.4f}")
    print(f"  R² = {metrics['r2']:.4f}")
    if not np.isnan(metrics['distance_corr']):
        print(f"  Distance corr = {metrics['distance_corr']:.4f}")
    if not np.isnan(metrics['mic']):
        print(f"  MIC = {metrics['mic']:.4f}")
    
    # Create plot
    print("\nGenerating plot...")
    create_plot(df_clean, metrics, output_dir, args.base_name, args.color_by, filters)
    
    # Save metrics
    save_metrics(metrics, output_dir, args.base_name, filters)
    
    # Save filtered data
    data_file = output_dir / f"{args.base_name}_data.parquet"
    df_clean.to_parquet(data_file, index=False)
    print(f"Saved filtered data to: {data_file}")
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nOutput files in: {output_dir}/")
    print(f"  - {args.base_name}_correlation.png (plot)")
    print(f"  - {args.base_name}_metrics.csv (metrics)")
    print(f"  - {args.base_name}_report.txt (detailed report)")
    print(f"  - {args.base_name}_data.parquet (filtered data)")


if __name__ == "__main__":
    main()
