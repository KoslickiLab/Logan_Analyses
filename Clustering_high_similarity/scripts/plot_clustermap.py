#!/usr/bin/env python3
"""
Create clustermap visualization from saved clustering artifacts.

This script loads pre-computed clustering results and generates visualizations
without recomputing the expensive clustering step.
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import load_npz
from scipy.cluster.hierarchy import linkage as recompute_linkage
from scipy.spatial.distance import squareform


def setup_logging(verbose: bool = False) -> logging.Logger:
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Suppress verbose matplotlib font messages
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)
    
    return logging.getLogger(__name__)


def load_artifacts(artifacts_dir: Path, logger: logging.Logger):
    """
    Load saved clustering artifacts.
    
    Args:
        artifacts_dir: Directory containing clustering artifacts
        logger: Logger instance
        
    Returns:
        Tuple of (distance_matrix, sample_ids, linkage_matrix)
    """
    logger.info(f"Loading artifacts from {artifacts_dir}")
    
    # Check required files exist
    required_files = [
        'expanded_similarity_matrix.npz',
        'expanded_sample_ids.txt',
        'expanded_linkage_matrix.npy'
    ]
    
    for filename in required_files:
        filepath = artifacts_dir / filename
        if not filepath.exists():
            logger.error(f"Required file not found: {filepath}")
            sys.exit(1)
    
    # Load sparse similarity matrix and convert to distance
    similarity_matrix = load_npz(artifacts_dir / 'expanded_similarity_matrix.npz')
    distance_matrix = 1 - similarity_matrix.toarray()
    logger.info(f"Loaded distance matrix: {distance_matrix.shape}")
    
    # Load ordered sample IDs
    with open(artifacts_dir / 'expanded_sample_ids.txt', 'r') as f:
        sample_ids = [line.strip() for line in f]
    logger.info(f"Loaded {len(sample_ids)} sample IDs")
    
    # Load linkage matrix
    linkage_matrix = np.load(artifacts_dir / 'expanded_linkage_matrix.npy')
    logger.info(f"Loaded linkage matrix: {linkage_matrix.shape}")
    
    return distance_matrix, sample_ids, linkage_matrix


def create_clustermap(distance_matrix, sample_ids, linkage_matrix,
                     output_file, max_samples, dpi, figsize, logger):
    """
    Create and save clustermap visualization.
    
    Args:
        distance_matrix: Distance matrix (1 - Jaccard)
        sample_ids: Ordered list of sample IDs
        linkage_matrix: Precomputed linkage matrix
        output_file: Path for output plot
        max_samples: Maximum number of samples to plot
        dpi: Output DPI
        figsize: Figure size tuple (width, height)
        logger: Logger instance
    """
    n_samples = distance_matrix.shape[0]
    original_n_samples = n_samples
    
    # Downsample if necessary
    if n_samples > max_samples:
        logger.info(f"Downsampling from {n_samples} to {max_samples} samples")
        step = n_samples // max_samples
        indices = np.arange(0, n_samples, step)[:max_samples]
        
        distance_matrix = distance_matrix[np.ix_(indices, indices)]
        sample_ids = [sample_ids[i] for i in indices]
        n_samples = len(sample_ids)
        
        # Recompute linkage for downsampled data
        logger.info("Recomputing linkage for downsampled data")
        condensed = squareform(distance_matrix, checks=False)
        linkage_matrix = recompute_linkage(condensed, method='average')
    
    # Scale figure size for large heatmaps to avoid rasterization issues
    # Rule: need at least 0.5 pixels per sample for proper rendering
    min_pixels_per_sample = 0.5
    required_pixels = n_samples * min_pixels_per_sample
    current_pixels = min(figsize[0] * dpi, figsize[1] * dpi)
    
    if required_pixels > current_pixels:
        scale_factor = required_pixels / current_pixels
        figsize = (int(figsize[0] * scale_factor), int(figsize[1] * scale_factor))
        logger.info(f"Scaled figure size to {figsize[0]}×{figsize[1]} inches for {n_samples} samples")
    
    logger.info(f"Creating clustermap with {n_samples} samples")
    
    # Set matplotlib rasterization threshold for large heatmaps
    # This ensures proper rendering even with dense data
    import matplotlib
    matplotlib.rcParams['agg.path.chunksize'] = 10000
    
    # Create clustermap with precomputed linkage
    g = sns.clustermap(
        distance_matrix,
        row_linkage=linkage_matrix,
        col_linkage=linkage_matrix,
        cmap='viridis',
        xticklabels=False,
        yticklabels=False,
        cbar_kws={'label': 'Distance (1 - Jaccard)'},
        figsize=figsize
    )
    
    title = f'Hierarchical Clustering of {original_n_samples:,} Samples'
    if n_samples < original_n_samples:
        title += f'\n(showing {n_samples:,} samples)'
    title += '\n(Distance = 1 - Jaccard)'
    
    g.fig.suptitle(title, y=0.98, fontsize=14)
    
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved clustermap to {output_file} (DPI={dpi})")


def main():
    parser = argparse.ArgumentParser(
        description='Create clustermap visualization from saved clustering artifacts',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--artifacts-dir',
        type=Path,
        required=True,
        help='Directory containing clustering artifacts'
    )
    
    parser.add_argument(
        '--output-file',
        type=Path,
        required=True,
        help='Output file for the clustermap (e.g., clustermap.png)'
    )
    
    parser.add_argument(
        '--max-samples',
        type=int,
        default=5000,
        help='Maximum samples to plot (for performance)'
    )
    
    parser.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='Output DPI (dots per inch)'
    )
    
    parser.add_argument(
        '--figsize',
        nargs=2,
        type=int,
        default=[16, 14],
        help='Figure size: width height (in inches)'
    )
    
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Setup
    logger = setup_logging(args.verbose)
    
    logger.info("="*80)
    logger.info("Clustermap Visualization from Saved Artifacts")
    logger.info("="*80)
    
    # Load artifacts
    distance_matrix, sample_ids, linkage_matrix = load_artifacts(
        args.artifacts_dir, logger
    )
    
    # Create clustermap
    create_clustermap(
        distance_matrix,
        sample_ids,
        linkage_matrix,
        args.output_file,
        args.max_samples,
        args.dpi,
        tuple(args.figsize),
        logger
    )
    
    logger.info("="*80)
    logger.info("Visualization completed successfully")
    logger.info("="*80)


if __name__ == "__main__":
    main()
