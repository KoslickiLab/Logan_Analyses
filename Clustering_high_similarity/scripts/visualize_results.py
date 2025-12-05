#!/usr/bin/env python3
"""
Visualization script for Jaccard matrix analysis results.
Creates detailed plots from the CSV outputs.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import sys
import os

sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.3)


def load_results(output_dir="/mnt/user-data/outputs"):
    """Load analysis results from CSV files."""
    results = {}
    
    edge_file = os.path.join(output_dir, "intermediate_similarity_edges.csv")
    comp_file = os.path.join(output_dir, "component_membership.csv")
    
    if os.path.exists(edge_file):
        print(f"Loading edges from {edge_file}...")
        results['edges'] = pd.read_csv(edge_file)
        print(f"  Loaded {len(results['edges']):,} edges")
    
    if os.path.exists(comp_file):
        print(f"Loading components from {comp_file}...")
        results['components'] = pd.read_csv(comp_file)
        print(f"  Loaded {len(results['components']):,} sample-component pairs")
    
    return results


def plot_similarity_distribution(edges_df, output_dir="/mnt/user-data/outputs"):
    """Plot the distribution of Jaccard similarities."""
    print("\nCreating similarity distribution plots...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    similarities = edges_df['jaccard_similarity'].values
    
    # 1. Full histogram
    axes[0, 0].hist(similarities, bins=100, edgecolor='black', alpha=0.7, color='steelblue')
    axes[0, 0].set_xlabel('Jaccard Similarity')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title(f'Distribution of Jaccard Similarities (n={len(similarities):,})')
    axes[0, 0].axvline(np.median(similarities), color='red', linestyle='--', 
                       label=f'Median: {np.median(similarities):.4f}')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. Log scale
    axes[0, 1].hist(similarities, bins=100, edgecolor='black', alpha=0.7, color='forestgreen')
    axes[0, 1].set_xlabel('Jaccard Similarity')
    axes[0, 1].set_ylabel('Frequency (log scale)')
    axes[0, 1].set_title('Distribution (Log Scale)')
    axes[0, 1].set_yscale('log')
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. Cumulative distribution
    sorted_sim = np.sort(similarities)
    cumulative = np.arange(1, len(sorted_sim) + 1) / len(sorted_sim)
    axes[1, 0].plot(sorted_sim, cumulative, linewidth=2, color='darkorange')
    axes[1, 0].set_xlabel('Jaccard Similarity')
    axes[1, 0].set_ylabel('Cumulative Probability')
    axes[1, 0].set_title('Cumulative Distribution Function')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].axhline(0.5, color='gray', linestyle='--', alpha=0.5)
    axes[1, 0].axvline(np.median(similarities), color='red', linestyle='--', alpha=0.5)
    
    # 4. Boxplot with quantiles
    axes[1, 1].boxplot([similarities], vert=False, widths=0.6)
    axes[1, 1].set_xlabel('Jaccard Similarity')
    axes[1, 1].set_title('Box Plot')
    axes[1, 1].set_yticks([])
    axes[1, 1].grid(True, alpha=0.3, axis='x')
    
    # Add statistics text
    stats_text = f"""Statistics:
    n = {len(similarities):,}
    min = {np.min(similarities):.4f}
    Q1 = {np.percentile(similarities, 25):.4f}
    median = {np.median(similarities):.4f}
    Q3 = {np.percentile(similarities, 75):.4f}
    max = {np.max(similarities):.4f}
    mean = {np.mean(similarities):.4f}
    std = {np.std(similarities):.4f}"""
    
    axes[1, 1].text(0.02, 0.98, stats_text, transform=axes[1, 1].transAxes,
                    verticalalignment='top', fontsize=9,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, "similarity_distribution_detailed.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved to: {output_file}")
    plt.close()


def plot_component_sizes(components_df, output_dir="/mnt/user-data/outputs"):
    """Plot component size distribution."""
    print("\nCreating component size distribution plots...")
    
    # Get component sizes
    component_sizes = components_df.groupby('component_id').size().values
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # 1. Histogram
    axes[0].hist(component_sizes, bins=50, edgecolor='black', alpha=0.7, color='purple')
    axes[0].set_xlabel('Component Size (number of samples)')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title(f'Component Size Distribution (n={len(component_sizes):,} components)')
    axes[0].grid(True, alpha=0.3)
    
    # 2. Log-log plot (for power law)
    size_counts = Counter(component_sizes)
    sizes = sorted(size_counts.keys())
    counts = [size_counts[s] for s in sizes]
    
    axes[1].scatter(sizes, counts, alpha=0.6, s=50, color='darkred')
    axes[1].set_xlabel('Component Size')
    axes[1].set_ylabel('Number of Components')
    axes[1].set_title('Size vs Frequency (Log-Log Scale)')
    axes[1].set_xscale('log')
    axes[1].set_yscale('log')
    axes[1].grid(True, alpha=0.3, which='both')
    
    # 3. Top 20 components
    top_components = components_df.groupby('component_id').size().nlargest(20)
    
    axes[2].barh(range(len(top_components)), top_components.values, color='teal')
    axes[2].set_xlabel('Number of Samples')
    axes[2].set_ylabel('Component Rank')
    axes[2].set_title('Top 20 Largest Components')
    axes[2].invert_yaxis()
    axes[2].grid(True, alpha=0.3, axis='x')
    
    # Add value labels
    for i, v in enumerate(top_components.values):
        axes[2].text(v, i, f' {v:,}', va='center', fontsize=9)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, "component_size_distribution.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved to: {output_file}")
    plt.close()
    
    # Print statistics
    print(f"\nComponent size statistics:")
    print(f"  Total components: {len(component_sizes):,}")
    print(f"  Min size: {np.min(component_sizes)}")
    print(f"  Max size: {np.max(component_sizes)}")
    print(f"  Mean size: {np.mean(component_sizes):.2f}")
    print(f"  Median size: {np.median(component_sizes):.0f}")
    print(f"\n  Singletons (size=1): {np.sum(component_sizes == 1):,}")
    print(f"  Pairs (size=2): {np.sum(component_sizes == 2):,}")
    print(f"  Large (size>10): {np.sum(component_sizes > 10):,}")
    print(f"  Very large (size>100): {np.sum(component_sizes > 100):,}")


def plot_degree_distribution(edges_df, output_dir="/mnt/user-data/outputs"):
    """Plot degree distribution (how many connections each sample has)."""
    print("\nCreating degree distribution plots...")
    
    # Calculate degrees
    sample_i_counts = edges_df['sample_i'].value_counts()
    sample_j_counts = edges_df['sample_j'].value_counts()
    
    all_samples = pd.concat([sample_i_counts, sample_j_counts]).groupby(level=0).sum()
    degrees = all_samples.values
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # 1. Histogram
    axes[0].hist(degrees, bins=50, edgecolor='black', alpha=0.7, color='indianred')
    axes[0].set_xlabel('Degree (number of connections)')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title(f'Degree Distribution (n={len(degrees):,} samples with connections)')
    axes[0].grid(True, alpha=0.3)
    axes[0].axvline(np.median(degrees), color='blue', linestyle='--', 
                    label=f'Median: {np.median(degrees):.0f}')
    axes[0].legend()
    
    # 2. Log-log scale
    degree_counts = Counter(degrees)
    deg_values = sorted(degree_counts.keys())
    freq_values = [degree_counts[d] for d in deg_values]
    
    axes[1].scatter(deg_values, freq_values, alpha=0.6, s=50, color='navy')
    axes[1].set_xlabel('Degree')
    axes[1].set_ylabel('Frequency')
    axes[1].set_title('Degree Distribution (Log-Log Scale)')
    axes[1].set_xscale('log')
    axes[1].set_yscale('log')
    axes[1].grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, "degree_distribution.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved to: {output_file}")
    plt.close()
    
    print(f"\nDegree statistics:")
    print(f"  Min: {np.min(degrees)}")
    print(f"  Max: {np.max(degrees)}")
    print(f"  Mean: {np.mean(degrees):.2f}")
    print(f"  Median: {np.median(degrees):.0f}")


def plot_similarity_vs_component(edges_df, components_df, output_dir="/mnt/user-data/outputs"):
    """Plot how similarity relates to component structure."""
    print("\nCreating similarity vs component plots...")
    
    # Merge edges with component info
    edges_with_comp = edges_df.merge(
        components_df[['sample_id', 'component_size']],
        left_on='sample_i', right_on='sample_id', how='left'
    ).rename(columns={'component_size': 'comp_size_i'})
    
    # Bin by component size
    size_bins = [0, 2, 5, 10, 20, 50, 100, 1000, 10000]
    bin_labels = ['1-2', '3-5', '6-10', '11-20', '21-50', '51-100', '101-1000', '>1000']
    
    edges_with_comp['size_bin'] = pd.cut(
        edges_with_comp['comp_size_i'], 
        bins=size_bins + [float('inf')],
        labels=bin_labels + ['>10000']
    )
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # 1. Violin plot
    valid_data = edges_with_comp.dropna(subset=['size_bin', 'jaccard_similarity'])
    if len(valid_data) > 0:
        # Only include bins with data
        bins_with_data = valid_data['size_bin'].value_counts().index[:8]
        plot_data = valid_data[valid_data['size_bin'].isin(bins_with_data)]
        
        parts = axes[0].violinplot(
            [plot_data[plot_data['size_bin'] == bin]['jaccard_similarity'].values 
             for bin in bins_with_data],
            positions=range(len(bins_with_data)),
            showmeans=True,
            showmedians=True
        )
        
        axes[0].set_xticks(range(len(bins_with_data)))
        axes[0].set_xticklabels(bins_with_data, rotation=45)
        axes[0].set_xlabel('Component Size')
        axes[0].set_ylabel('Jaccard Similarity')
        axes[0].set_title('Similarity Distribution by Component Size')
        axes[0].grid(True, alpha=0.3, axis='y')
    
    # 2. Scatter plot with trend
    sample_sizes = components_df.groupby('component_id')['component_size'].first()
    mean_similarities = edges_df.groupby('sample_i')['jaccard_similarity'].mean()
    
    # Match samples
    samples_with_both = set(sample_sizes.index) & set(mean_similarities.index)
    sizes = [sample_sizes[s] for s in samples_with_both]
    means = [mean_similarities[s] for s in samples_with_both]
    
    if len(sizes) > 0:
        axes[1].scatter(sizes, means, alpha=0.3, s=20)
        axes[1].set_xlabel('Component Size')
        axes[1].set_ylabel('Mean Jaccard Similarity')
        axes[1].set_title('Component Size vs Mean Similarity')
        axes[1].set_xscale('log')
        axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, "similarity_vs_components.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved to: {output_file}")
    plt.close()


def create_summary_report(results, output_dir="/mnt/user-data/outputs"):
    """Create a text summary report."""
    print("\nCreating summary report...")
    
    report_file = os.path.join(output_dir, "analysis_summary.txt")
    
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("JACCARD MATRIX ANALYSIS SUMMARY\n")
        f.write("="*80 + "\n\n")
        
        if 'edges' in results:
            edges_df = results['edges']
            f.write(f"EDGE ANALYSIS\n")
            f.write(f"-"*80 + "\n")
            f.write(f"Total edges (sample pairs with intermediate similarity): {len(edges_df):,}\n")
            f.write(f"\nSimilarity statistics:\n")
            f.write(f"  Min:    {edges_df['jaccard_similarity'].min():.6f}\n")
            f.write(f"  Q1:     {edges_df['jaccard_similarity'].quantile(0.25):.6f}\n")
            f.write(f"  Median: {edges_df['jaccard_similarity'].median():.6f}\n")
            f.write(f"  Q3:     {edges_df['jaccard_similarity'].quantile(0.75):.6f}\n")
            f.write(f"  Max:    {edges_df['jaccard_similarity'].max():.6f}\n")
            f.write(f"  Mean:   {edges_df['jaccard_similarity'].mean():.6f}\n")
            f.write(f"  Std:    {edges_df['jaccard_similarity'].std():.6f}\n")
            
            # Similarity ranges
            f.write(f"\nSimilarity ranges:\n")
            ranges = [
                (0.9, 1.0, "Very high (0.9-1.0)"),
                (0.8, 0.9, "High (0.8-0.9)"),
                (0.5, 0.8, "Medium (0.5-0.8)"),
                (0.0, 0.5, "Low (0.0-0.5)")
            ]
            for low, high, label in ranges:
                count = ((edges_df['jaccard_similarity'] >= low) & 
                        (edges_df['jaccard_similarity'] < high)).sum()
                pct = 100 * count / len(edges_df)
                f.write(f"  {label}: {count:,} ({pct:.2f}%)\n")
        
        if 'components' in results:
            components_df = results['components']
            f.write(f"\n\nCOMPONENT ANALYSIS\n")
            f.write(f"-"*80 + "\n")
            f.write(f"Total samples with connections: {len(components_df):,}\n")
            
            n_components = components_df['component_id'].nunique()
            f.write(f"Total components: {n_components:,}\n")
            
            component_sizes = components_df.groupby('component_id').size()
            f.write(f"\nComponent size statistics:\n")
            f.write(f"  Min:    {component_sizes.min()}\n")
            f.write(f"  Q1:     {component_sizes.quantile(0.25):.0f}\n")
            f.write(f"  Median: {component_sizes.median():.0f}\n")
            f.write(f"  Q3:     {component_sizes.quantile(0.75):.0f}\n")
            f.write(f"  Max:    {component_sizes.max()}\n")
            f.write(f"  Mean:   {component_sizes.mean():.2f}\n")
            
            f.write(f"\nTop 10 largest components:\n")
            top10 = component_sizes.nlargest(10)
            for i, (comp_id, size) in enumerate(top10.items(), 1):
                f.write(f"  {i:2d}. Component {comp_id}: {size:,} samples\n")
        
        f.write(f"\n" + "="*80 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*80 + "\n")
    
    print(f"Summary report saved to: {report_file}")


def main(output_dir="/mnt/user-data/outputs"):
    """Main visualization pipeline."""
    print("="*80)
    print("JACCARD MATRIX VISUALIZATION")
    print("="*80)
    
    # Load results
    results = load_results(output_dir)
    
    if not results:
        print("\nNo results found! Run the analysis script first.")
        print("Expected files in:", output_dir)
        return
    
    # Create visualizations
    if 'edges' in results:
        plot_similarity_distribution(results['edges'], output_dir)
        plot_degree_distribution(results['edges'], output_dir)
    
    if 'components' in results:
        plot_component_sizes(results['components'], output_dir)
    
    if 'edges' in results and 'components' in results:
        plot_similarity_vs_component(results['edges'], results['components'], output_dir)
    
    # Create summary report
    create_summary_report(results, output_dir)
    
    print("\n" + "="*80)
    print("VISUALIZATION COMPLETE")
    print("="*80)
    print(f"\nAll outputs saved to: {output_dir}/")
    print("\nGenerated files:")
    print("  - similarity_distribution_detailed.png")
    print("  - component_size_distribution.png")
    print("  - degree_distribution.png")
    print("  - similarity_vs_components.png")
    print("  - analysis_summary.txt")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        output_dir = sys.argv[1]
    else:
        output_dir = "/mnt/user-data/outputs"
    
    main(output_dir)
