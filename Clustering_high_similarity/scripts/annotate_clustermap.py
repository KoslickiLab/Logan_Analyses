#!/usr/bin/env python3
"""
Create annotated clustermaps from pre-computed linkage matrices with metadata decorations.

This script takes a pre-computed linkage matrix, similarity matrix, and ordered sample IDs,
then creates a clustermap with metadata annotations as horizontal color bars (col_colors).
"""

import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import duckdb
from pathlib import Path
import sys
import warnings
from scipy.sparse import load_npz, issparse
from scipy.cluster.hierarchy import dendrogram, linkage


def parse_args():
    parser = argparse.ArgumentParser(
        description='Create annotated clustermap with metadata column colors',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with one metadata column
  %(prog)s -l linkage.npy -s similarity.npz -i sample_ids.txt -o output.png \\
    -m country

  # Subsample to 5000 samples for visualization
  %(prog)s -l linkage.npy -s similarity.npz -i sample_ids.txt -o output.png \\
    -m country:10 -m organism --subsample 5000

  # Multiple metadata columns
  %(prog)s -l linkage.npy -s similarity.npz -i sample_ids.txt -o output.png \\
    -m country -m organism -m instrument

  # Limit to top 10 categories for a column
  %(prog)s -l linkage.npy -s similarity.npz -i sample_ids.txt -o output.png \\
    -m country:10

  # Mix of limited and unlimited columns  
  %(prog)s -l linkage.npy -s similarity.npz -i sample_ids.txt -o output.png \\
    -m country:5 -m organism -m instrument:3

  # Use bucketed scalar variables
  %(prog)s -l linkage.npy -s similarity.npz -i sample_ids.txt -o output.png \\
    -m mbytes -m avgspotlen

Metadata column format:
  column_name       - Include all unique values
  column_name:N     - Include only top N most frequent values (others become "Other")

Special handling:
  - mbytes: Automatically bucketed by order of magnitude
  - avgspotlen: Automatically bucketed into size ranges
  - Columns with only one unique value are automatically skipped
        """
    )
    
    parser.add_argument('-l', '--linkage', required=True, type=Path,
                        help='Path to linkage matrix (.npy file)')
    parser.add_argument('-s', '--similarity', required=True, type=Path,
                        help='Path to similarity matrix (.csv, .tsv, .npy, or .npz file)')
    parser.add_argument('-i', '--sample-ids', required=True, type=Path,
                        help='Path to ordered sample IDs file (one accession per line)')
    parser.add_argument('-d', '--database', type=Path,
                        default=Path('/scratch/shared_data_new/Logan_yacht_data/metadata/aws_sra_metadata/metadata_geo_joined.duckdb'),
                        help='Path to DuckDB metadata database (default: %(default)s)')
    parser.add_argument('-m', '--metadata', action='append', required=True, dest='metadata_cols',
                        help='Metadata column to include. Can specify multiple times. '
                             'Format: column_name or column_name:N to show only top N categories')
    parser.add_argument('-o', '--output', required=True, type=Path,
                        help='Output figure path (e.g., output.png, output.pdf)')
    parser.add_argument('--subsample', type=int, default=None,
                        help='Randomly subsample N samples for visualization (recommended for >10k samples)')
    parser.add_argument('--figsize', type=float, nargs=2, default=[20, 15],
                        help='Figure size in inches (width height). Default: 20 15')
    parser.add_argument('--cmap', default='viridis',
                        help='Colormap for the heatmap. Default: viridis')
    parser.add_argument('--dpi', type=int, default=300,
                        help='DPI for output figure. Default: 300')
    parser.add_argument('--vmin', type=float, default=None,
                        help='Minimum value for colormap scaling')
    parser.add_argument('--vmax', type=float, default=None,
                        help='Maximum value for colormap scaling')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for subsampling. Default: 42')
    
    return parser.parse_args()


def bucket_mbytes(value):
    """Bucket mbytes by order of magnitude."""
    if pd.isna(value) or value == 0:
        return '0 MB'
    
    try:
        value = float(value)
        if value == 0:
            return '0 MB'
        
        magnitude = int(np.floor(np.log10(value)))
        
        if magnitude < 0:
            return '<1 MB'
        elif magnitude == 0:
            return '1-10 MB'
        elif magnitude == 1:
            return '10-100 MB'
        elif magnitude == 2:
            return '100-1K MB'
        elif magnitude == 3:
            return '1-10 GB'
        elif magnitude == 4:
            return '10-100 GB'
        else:
            return '>100 GB'
    except (ValueError, TypeError):
        return 'NA'


def bucket_avgspotlen(value):
    """Bucket avgspotlen into predefined ranges."""
    if pd.isna(value):
        return 'NA'
    
    try:
        value = float(value)
        if value <= 100:
            return 'short [0-100]'
        elif value <= 250:
            return 'medium [101-250]'
        elif value <= 1000:
            return 'longer [251-1000]'
        else:
            return 'longest [1001+]'
    except (ValueError, TypeError):
        return 'NA'


def load_sample_ids(path):
    """Load sample IDs from file."""
    print(f"Loading sample IDs from {path}...")
    with open(path, 'r') as f:
        sample_ids = [line.strip() for line in f if line.strip()]
    print(f"  Loaded {len(sample_ids)} sample IDs")
    return sample_ids


def load_linkage(path):
    """Load linkage matrix."""
    print(f"Loading linkage matrix from {path}...")
    linkage = np.load(path)
    print(f"  Linkage matrix shape: {linkage.shape}")
    return linkage


def load_similarity_matrix(path, sample_ids):
    """Load similarity matrix from CSV, TSV, NPY, or sparse NPZ."""
    print(f"Loading similarity matrix from {path}...")
    
    if path.suffix == '.npz':
        # Load sparse matrix
        sparse_matrix = load_npz(path)
        print(f"  Loaded sparse matrix: shape={sparse_matrix.shape}, nnz={sparse_matrix.nnz}")
        print(f"  Sparsity: {1 - sparse_matrix.nnz / (sparse_matrix.shape[0] ** 2):.4f}")
        
        # Convert to dense for clustermap (only if reasonably sized)
        if sparse_matrix.shape[0] > 50000:
            warnings.warn(f"Matrix is very large ({sparse_matrix.shape[0]}x{sparse_matrix.shape[0]}). "
                         f"Consider using --subsample for better performance.")
        
        # Assume order matches sample_ids
        return sparse_matrix
    elif path.suffix == '.npy':
        matrix = np.load(path)
        print(f"  Matrix shape: {matrix.shape}")
        # Assume order matches sample_ids - DO NOT reorder
        return matrix
    elif path.suffix in ['.csv', '.tsv']:
        sep = '\t' if path.suffix == '.tsv' else ','
        df = pd.read_csv(path, sep=sep, index_col=0)
        print(f"  Matrix shape: {df.shape}")
        # Assume order matches sample_ids - DO NOT reorder
        return df.values
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")



def fetch_metadata(db_path, sample_ids, columns):
    """Fetch metadata for given sample IDs and columns from DuckDB."""
    print(f"Fetching metadata from {db_path}...")
    print(f"  Requesting columns: {', '.join(columns)}")
    
    # Connect in read-only mode
    con = duckdb.connect(str(db_path), read_only=True)
    
    # Build query - select only needed columns
    columns_str = ', '.join(['acc'] + columns)
    placeholders = ', '.join([f"'{sid}'" for sid in sample_ids])
    query = f"SELECT {columns_str} FROM metadata_geo_joined WHERE acc IN ({placeholders})"
    
    try:
        # Execute query
        df = con.execute(query).fetchdf()
        print(f"  Retrieved metadata for {len(df)} samples")
    finally:
        con.close()
    
    # Set acc as index
    df = df.set_index('acc')
    
    # Reindex to match sample_ids order (fills missing with NaN)
    df = df.reindex(sample_ids)
    
    # Report missing samples
    missing = df.index[df.isnull().all(axis=1)]
    if len(missing) > 0:
        print(f"  WARNING: {len(missing)} samples not found in database")
    
    return df


def process_metadata_column(series, column_name, top_n=None):
    """
    Process a metadata column for visualization.
    
    Returns:
        tuple: (processed_series, is_valid, color_map, category_counts)
        - is_valid is False if column should be skipped (e.g., only one unique value)
        - category_counts is a Series with value counts for reporting
    """
    # Handle special bucketing cases
    if column_name == 'mbytes':
        series = series.apply(bucket_mbytes)
    elif column_name == 'avgspotlen':
        series = series.apply(bucket_avgspotlen)
    
    # Convert to string and handle missing values
    series = series.astype(str)
    series = series.replace(['nan', 'None', '', 'NULL', 'null'], 'NA')
    
    # Get value counts before filtering
    value_counts = series.value_counts()
    
    # Get unique values (excluding NA for uniqueness check)
    unique_vals = series[series != 'NA'].unique()
    
    # Skip if only one unique value (excluding NA)
    if len(unique_vals) <= 1:
        return series, False, None, value_counts
    
    # Handle top_n filtering
    if top_n is not None and len(value_counts) > top_n:
        # Keep top N categories, make others "Other"
        # Don't count NA toward the top N
        non_na_counts = value_counts[value_counts.index != 'NA']
        top_categories = non_na_counts.head(top_n).index.tolist()
        
        def categorize(x):
            if x == 'NA':
                return 'NA'
            elif x in top_categories:
                return x
            else:
                return 'Other'
        
        series = series.apply(categorize)
        # Recalculate value counts
        value_counts = series.value_counts()
    
    # Get final unique values for color mapping
    unique_vals = series.unique()
    n_colors = len(unique_vals)
    
    # Create color palette (excluding NA which will be white)
    non_na_vals = [v for v in unique_vals if v != 'NA']
    n_colors_needed = len(non_na_vals)
    
    if n_colors_needed <= 10:
        colors = sns.color_palette('tab10', n_colors_needed)
    elif n_colors_needed <= 20:
        colors = sns.color_palette('tab20', n_colors_needed)
    else:
        colors = sns.color_palette('husl', n_colors_needed)
    
    # Sort non-NA values for consistent coloring
    sorted_vals = sorted(non_na_vals)
    
    # Create color map with white for NA
    color_map = dict(zip(sorted_vals, colors))
    if 'NA' in unique_vals:
        color_map['NA'] = (1.0, 1.0, 1.0)  # White for missing values
    
    return series, True, color_map, value_counts


def get_dendrogram_order(linkage_matrix, n_samples):
    """
    Get sample ordering from linkage matrix using dendrogram.
    
    Args:
        linkage_matrix: Linkage matrix from scipy
        n_samples: Number of samples
        
    Returns:
        Array of indices representing the dendrogram leaf order
    """
    # Create dendrogram to get leaf ordering
    dendro = dendrogram(linkage_matrix, no_plot=True)
    return np.array(dendro['leaves'])


def create_clustermap(similarity_matrix, linkage_matrix, metadata_df, metadata_specs, 
                     sample_ids, output_path, figsize=(20, 15), cmap='viridis', 
                     dpi=300, vmin=None, vmax=None, subsample=None, seed=42):
    """
    Create annotated clustermap using pre-computed linkage matrix.
    
    Args:
        similarity_matrix: Similarity matrix (sparse or dense, N x N)
        linkage_matrix: Pre-computed linkage matrix from scipy
        metadata_df: DataFrame with metadata (index = sample IDs)
        metadata_specs: list of tuples (column_name, top_n)
        sample_ids: List of sample IDs in original order
        output_path: Path to save figure
        figsize: Figure size tuple
        cmap: Colormap name
        dpi: DPI for output
        vmin, vmax: Color scale limits
        subsample: Number of samples to randomly select (None = use all)
        seed: Random seed for subsampling
    """
    print("\nProcessing metadata columns...")
    
    n_samples = len(sample_ids)
    print(f"Starting with {n_samples} samples")
    
    # Track whether we actually subsample
    did_subsample = False
    
    # CRITICAL: Get ordering from linkage matrix FIRST (before subsampling)
    # The linkage matrix refers to the full dataset
    print("Computing dendrogram ordering from linkage matrix...")
    leaf_order = get_dendrogram_order(linkage_matrix, n_samples)
    
    # Reorder everything according to dendrogram BEFORE subsampling
    print("Reordering data according to dendrogram...")
    if issparse(similarity_matrix):
        # For sparse, reorder then subsample then convert to dense
        similarity_matrix = similarity_matrix[leaf_order, :][:, leaf_order]
    else:
        similarity_matrix = similarity_matrix[np.ix_(leaf_order, leaf_order)]
    
    sample_ids = [sample_ids[i] for i in leaf_order]
    metadata_df = metadata_df.iloc[leaf_order].reset_index(drop=True)
    
    # NOW subsample if requested (from the already-ordered data)
    if subsample is not None and subsample < n_samples:
        did_subsample = True
        print(f"Subsampling {subsample} from {n_samples} samples (preserving dendrogram order)...")
        np.random.seed(seed)
        
        # Randomly select indices from the ordered data
        # Keep them sorted to maintain relative order
        indices = np.sort(np.random.choice(n_samples, subsample, replace=False))
        
        # Subsample
        sample_ids = [sample_ids[i] for i in indices]
        metadata_df = metadata_df.iloc[indices].reset_index(drop=True)
        
        if issparse(similarity_matrix):
            similarity_matrix = similarity_matrix[indices, :][:, indices].toarray()
        else:
            similarity_matrix = similarity_matrix[np.ix_(indices, indices)]
        
        # Compute NEW linkage matrix for the subsampled data
        # This allows us to show dendrograms that reflect the subsampled clustering
        print("  Computing linkage for subsampled data...")
        from scipy.spatial.distance import squareform
        distance_matrix = 1 - similarity_matrix
        condensed_dist = squareform(distance_matrix, checks=False)
        subsampled_linkage = linkage(condensed_dist, method='average')
        
        n_samples = subsample
    elif issparse(similarity_matrix):
        # No subsampling but need to convert sparse to dense
        print("Converting sparse matrix to dense...")
        similarity_matrix = similarity_matrix.toarray()
        subsampled_linkage = None  # Use original linkage
    else:
        subsampled_linkage = None  # Use original linkage
    
    print(f"Working with {n_samples} samples")
    
    # Data is now ordered and subsampled - use as-is (no further reordering)
    ordered_matrix = similarity_matrix
    ordered_ids = sample_ids
    ordered_metadata = metadata_df
    
    # Process metadata columns
    col_colors_dict = {}
    color_maps = {}
    all_value_counts = {}
    
    for col_name, top_n in metadata_specs:
        if col_name not in ordered_metadata.columns:
            warnings.warn(f"Column '{col_name}' not found in metadata. Skipping.")
            continue
        
        series, is_valid, color_map, value_counts = process_metadata_column(
            ordered_metadata[col_name].copy(), col_name, top_n
        )
        
        if not is_valid:
            warnings.warn(f"Column '{col_name}' has only one unique value. Skipping.")
            continue
        
        # Map values to colors
        col_colors_dict[col_name] = series.map(color_map)
        color_maps[col_name] = color_map
        all_value_counts[col_name] = value_counts
        
        print(f"  {col_name}: {len(color_map)} categories" + 
              (f" (limited to top {top_n})" if top_n else ""))
    
    if not col_colors_dict:
        print("ERROR: No valid metadata columns to display!", file=sys.stderr)
        sys.exit(1)
    
    # Create col_colors as list of lists (one list per metadata track)
    # Each track is a list of RGB tuples, one per sample
    # Seaborn expects format: [[color1_track1, color2_track1, ...], [color1_track2, ...]]
    col_colors_list = []
    for col_name in col_colors_dict:
        track_colors = col_colors_dict[col_name].tolist()
        col_colors_list.append(track_colors)
    
    print(f"\nCreating clustermap with {len(col_colors_dict)} metadata annotations...")
    print(f"  col_colors: {len(col_colors_list)} tracks, {len(col_colors_list[0])} samples each")
    
    # For very large figures, tight_layout can fail
    # Temporarily patch it to handle errors gracefully
    import matplotlib.figure
    original_tight_layout = matplotlib.figure.Figure.tight_layout
    
    def safe_tight_layout(self, *args, **kwargs):
        """Wrapper for tight_layout that catches and ignores layout errors."""
        try:
            return original_tight_layout(self, *args, **kwargs)
        except (ValueError, RuntimeError) as e:
            # tight_layout failed (often with very large figures or empty axes)
            # Just skip it and proceed
            warnings.warn(f"tight_layout failed ({e}), proceeding without it")
            pass
    
    # Temporarily replace tight_layout
    matplotlib.figure.Figure.tight_layout = safe_tight_layout
    
    try:
        # Determine which linkage to use for dendrograms
        if did_subsample and subsampled_linkage is not None:
            # Use linkage computed from subsampled data
            print("  Showing dendrogram computed from subsampled data")
            dendrogram_linkage = subsampled_linkage
            # Let seaborn cluster to show dendrograms (will reorder based on subsampled linkage)
            do_cluster = True
        else:
            # Use original linkage
            print("  Showing dendrogram from original linkage matrix")
            dendrogram_linkage = linkage_matrix
            # Let seaborn cluster to show dendrograms
            do_cluster = True
        
        # Create clustermap with dendrograms and external colorbar
        g = sns.clustermap(
            ordered_matrix,
            row_linkage=dendrogram_linkage,
            col_linkage=dendrogram_linkage,
            row_cluster=do_cluster,  # Allow clustering to draw dendrograms
            col_cluster=do_cluster,  # Allow clustering to draw dendrograms
            col_colors=col_colors_list,
            cmap=cmap,
            figsize=figsize,
            xticklabels=False,
            yticklabels=False,
            vmin=vmin,
            vmax=vmax,
            cbar_pos=(0.02, 0.8, 0.03, 0.15),  # (left, bottom, width, height) - position on left side
            cbar_kws={'label': 'Similarity'},
            dendrogram_ratio=0.15,  # Larger dendrograms for visibility
            colors_ratio=0.03,
        )
    finally:
        # Restore original tight_layout
        matplotlib.figure.Figure.tight_layout = original_tight_layout
    
    # Create custom legends for each metadata column
    # Position them vertically on the right side
    n_cols = len(color_maps)
    legend_x = 1.02
    legend_y_spacing = 1.0 / max(n_cols, 1)
    
    for idx, (col_name, color_map) in enumerate(color_maps.items()):
        # Sort categories for legend (NA last if present)
        if 'NA' in color_map:
            sorted_categories = sorted([k for k in color_map.keys() if k != 'NA']) + ['NA']
        else:
            sorted_categories = sorted(color_map.keys())
        
        # Create legend handles
        handles = [mpatches.Patch(facecolor=color_map[cat], edgecolor='black', linewidth=0.5)
                  for cat in sorted_categories]
        
        # Calculate legend position (top to bottom)
        legend_y = 0.98 - (idx * legend_y_spacing)
        
        # Add count information to labels
        labels = []
        for cat in sorted_categories:
            count = all_value_counts[col_name].get(cat, 0)
            labels.append(f'{cat} (n={count})')
        
        # Add legend
        legend = g.fig.legend(
            handles, labels,
            title=col_name,
            loc='upper left',
            bbox_to_anchor=(legend_x, legend_y),
            frameon=True,
            fontsize=7,
            title_fontsize=8,
            ncol=1
        )
        legend.get_frame().set_linewidth(0.5)
    
    # Save figure (bbox_inches='tight' handles layout automatically)
    print(f"\nSaving figure to {output_path}...")
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Clustermap saved successfully!")
    
    # Print detailed summary
    print("\n" + "="*60)
    print("METADATA SUMMARY")
    print("="*60)
    for col_name, top_n in metadata_specs:
        if col_name in color_maps:
            counts = all_value_counts[col_name]
            n_unique = len(color_maps[col_name])
            total = counts.sum()
            
            print(f"\n{col_name}:")
            print(f"  Total categories: {n_unique}" + 
                  (f" (limited to top {top_n})" if top_n else ""))
            print(f"  Total samples: {total}")
            print(f"  Top categories:")
            for cat, count in counts.head(10).items():
                pct = 100 * count / total
                print(f"    {cat:30s} {count:6d} ({pct:5.1f}%)")
            if len(counts) > 10:
                print(f"    ... and {len(counts) - 10} more categories")
    print("="*60)


def parse_metadata_spec(spec):
    """Parse metadata specification like 'column_name' or 'column_name:10'."""
    parts = spec.split(':', 1)
    col_name = parts[0].strip()
    top_n = int(parts[1]) if len(parts) > 1 else None
    return col_name, top_n


def main():
    args = parse_args()
    
    print("="*60)
    print("CLUSTERMAP METADATA ANNOTATION")
    print("="*60)
    
    # Load data
    sample_ids = load_sample_ids(args.sample_ids)
    linkage_matrix = load_linkage(args.linkage)
    similarity_matrix = load_similarity_matrix(args.similarity, sample_ids)
    
    # Validate dimensions
    n_samples = len(sample_ids)
    expected_linkage_shape = (n_samples - 1, 4)
    
    if linkage_matrix.shape != expected_linkage_shape:
        print(f"WARNING: Expected linkage shape {expected_linkage_shape}, got {linkage_matrix.shape}")
        print(f"  This may indicate a mismatch between the linkage matrix and sample IDs.")
        # Try to proceed anyway - the linkage might be from a subset
    
    # Parse metadata specifications
    metadata_specs = [parse_metadata_spec(spec) for spec in args.metadata_cols]
    all_columns = list(set(spec[0] for spec in metadata_specs))  # Unique columns
    
    # Fetch metadata
    metadata_df = fetch_metadata(args.database, sample_ids, all_columns)
    
    # Warn about subsampling
    if args.subsample and args.subsample < n_samples:
        print(f"\nNote: Will subsample to {args.subsample} samples for visualization")
        print(f"  Original dataset: {n_samples} samples")
        print(f"  Random seed: {args.seed}")
    
    # Create clustermap
    create_clustermap(
        similarity_matrix,
        linkage_matrix,
        metadata_df,
        metadata_specs,
        sample_ids,
        args.output,
        figsize=tuple(args.figsize),
        cmap=args.cmap,
        dpi=args.dpi,
        vmin=args.vmin,
        vmax=args.vmax,
        subsample=args.subsample,
        seed=args.seed
    )
    
    print("\nDone!")


if __name__ == '__main__':
    main()
