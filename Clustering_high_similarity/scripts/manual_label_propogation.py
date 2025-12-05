#!/usr/bin/env python3
"""
Community-based hierarchical clustering for large-scale Jaccard similarity data.

This script processes CSV files containing Jaccard=1 similarities, identifies communities,
selects representatives, computes pairwise similarities, performs hierarchical clustering,
and creates visualizations.
"""

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
from scipy.sparse import csr_matrix, save_npz


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


def load_csv_files(input_dir: Path, limit: str, logger: logging.Logger) -> Dict[str, Set[str]]:
    """
    Load CSV files and build adjacency information for Jaccard=1 connections.

    Args:
        input_dir: Directory containing CSV files
        limit: Number of files to process ('all' or integer)
        logger: Logger instance

    Returns:
        Dictionary mapping sample IDs to sets of connected IDs (Jaccard=1)
    """
    csv_files = sorted(input_dir.glob("*_jaccard_1.csv"))

    if not csv_files:
        logger.error(f"No *_jaccard_1.csv files found in {input_dir}")
        sys.exit(1)

    if limit != 'all':
        try:
            n_files = int(limit)
            csv_files = csv_files[:n_files]
            logger.info(f"Processing first {n_files} files")
        except ValueError:
            logger.error(f"Invalid limit value: {limit}. Use 'all' or an integer.")
            sys.exit(1)
    else:
        logger.info(f"Processing all {len(csv_files)} files")

    # Build adjacency: for each ID, track all IDs with Jaccard=1
    adjacency: Dict[str, Set[str]] = {}

    for i, csv_file in enumerate(csv_files, 1):
        if i % 10000 == 0:
            logger.info(f"Processed {i}/{len(csv_files)} files")

        # Extract source ID from filename
        source_id = csv_file.stem.replace("_jaccard_1", "")

        try:
            df = pd.read_csv(csv_file)
            if df.empty:
                continue

            # Get IDs with Jaccard=1
            connected_ids = set(df['ID'].astype(str).tolist())

            # Add to adjacency
            if source_id not in adjacency:
                adjacency[source_id] = set()
            adjacency[source_id].update(connected_ids)

            # Add bidirectional connections
            for conn_id in connected_ids:
                if conn_id not in adjacency:
                    adjacency[conn_id] = set()
                adjacency[conn_id].add(source_id)

        except Exception as e:
            logger.warning(f"Error processing {csv_file}: {e}")
            continue

    logger.info(f"Loaded {len(adjacency)} unique sample IDs")
    return adjacency


def find_communities(adjacency: Dict[str, Set[str]], logger: logging.Logger) -> Dict[int, List[str]]:
    """
    Find connected components (communities) in the adjacency graph using iterative DFS.

    Args:
        adjacency: Dictionary mapping IDs to sets of connected IDs
        logger: Logger instance

    Returns:
        Dictionary mapping community ID to list of member IDs
    """
    visited = set()
    communities = {}
    community_id = 0

    # Find all connected components using iterative DFS (avoids recursion limit)
    for start_node in adjacency:
        if start_node not in visited:
            current_community = []
            stack = [start_node]

            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    current_community.append(node)

                    # Add unvisited neighbors to stack
                    for neighbor in adjacency.get(node, set()):
                        if neighbor not in visited:
                            stack.append(neighbor)

            communities[community_id] = sorted(current_community)
            community_id += 1

    # Include isolated nodes (those only in CSV files but not in adjacency keys)
    all_nodes = set(adjacency.keys())
    for node_set in adjacency.values():
        all_nodes.update(node_set)

    isolated = all_nodes - visited
    for node in isolated:
        communities[community_id] = [node]
        community_id += 1

    logger.info(f"Identified {len(communities)} communities")

    # Log community size statistics
    sizes = [len(comm) for comm in communities.values()]
    logger.info(f"Community size: min={min(sizes)}, max={max(sizes)}, "
                f"mean={np.mean(sizes):.2f}, median={np.median(sizes):.0f}")

    return communities


def save_community_outputs(communities: Dict[int, List[str]],
                           output_csv: Path,
                           output_json: Path,
                           logger: logging.Logger) -> None:
    """
    Save community assignments to CSV and JSON files.

    Args:
        communities: Dictionary mapping community ID to member IDs
        output_csv: Path for CSV output
        output_json: Path for JSON output
        logger: Logger instance
    """
    # Create CSV
    records = []
    for comm_id, members in communities.items():
        for member in members:
            records.append({'sample_id': member, 'community': comm_id})

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    logger.info(f"Saved community CSV to {output_csv}")

    # Create JSON
    json_data = {f"community_{k}": v for k, v in communities.items()}
    with open(output_json, 'w') as f:
        json.dump(json_data, f, indent=2)
    logger.info(f"Saved community JSON to {output_json}")


def select_representatives(communities: Dict[int, List[str]],
                           output_file: Path,
                           logger: logging.Logger) -> Tuple[List[str], Dict[str, int]]:
    """
    Select one representative per community and save to file.

    Args:
        communities: Dictionary mapping community ID to member IDs
        output_file: Path for representatives text file
        logger: Logger instance

    Returns:
        Tuple of (list of representative IDs, mapping from representative to community ID)
    """
    representatives = []
    rep_to_community = {}

    for comm_id, members in communities.items():
        rep = members[0]  # Arbitrary selection (first member)
        representatives.append(rep)
        rep_to_community[rep] = comm_id

    # Save to file
    with open(output_file, 'w') as f:
        for rep in representatives:
            f.write(f"{rep}\n")

    logger.info(f"Selected {len(representatives)} representatives")
    logger.info(f"Saved representatives to {output_file}")

    return representatives, rep_to_community


def compute_representative_similarities(representatives_file: Path,
                                        output_matrix: Path,
                                        query_pc_mat_path: str,
                                        matrix_path: str,
                                        db_path: str,
                                        logger: logging.Logger) -> None:
    """
    Call C++ executable to compute pairwise similarities between representatives.

    Args:
        representatives_file: Path to representatives text file
        output_matrix: Path for output similarity matrix
        query_pc_mat_path: Path to query_pc_mat executable
        matrix_path: Path to matrix directory
        db_path: Path to database directory
        logger: Logger instance
    """
    cmd = [
        query_pc_mat_path,
        "--matrix", matrix_path,
        "--db", db_path,
        "--row_file", str(representatives_file),
        "--col_file", str(representatives_file),
        "--show_all",
        "--write_to_file", str(output_matrix)
    ]

    logger.info(f"Computing pairwise similarities for {len(open(representatives_file).readlines())} representatives")
    logger.info(f"Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if result.stdout:
            logger.debug(f"stdout: {result.stdout}")
        if result.stderr:
            logger.debug(f"stderr: {result.stderr}")
        logger.info(f"Saved representative similarity matrix to {output_matrix}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running query_pc_mat: {e}")
        logger.error(f"stdout: {e.stdout}")
        logger.error(f"stderr: {e.stderr}")
        sys.exit(1)


def perform_hierarchical_clustering(similarity_matrix: np.ndarray,
                                    logger: logging.Logger) -> np.ndarray:
    """
    Perform hierarchical clustering with average linkage.

    Args:
        similarity_matrix: Square similarity matrix
        logger: Logger instance

    Returns:
        Linkage matrix from hierarchical clustering
    """
    # Convert similarity to distance
    distance_matrix = 1 - similarity_matrix

    # Ensure diagonal is zero
    np.fill_diagonal(distance_matrix, 0)

    # Convert to condensed form for linkage
    condensed_dist = squareform(distance_matrix, checks=False)

    logger.info("Performing hierarchical clustering (average linkage)")
    linkage_matrix = linkage(condensed_dist, method='average')

    return linkage_matrix


def expand_matrix(similarity_matrix: np.ndarray,
                  representatives: List[str],
                  communities: Dict[int, List[str]],
                  rep_to_community: Dict[str, int],
                  linkage_matrix: np.ndarray,
                  output_expanded: Path,
                  logger: logging.Logger) -> Tuple[np.ndarray, List[str]]:
    """
    Expand similarity matrix to include all community members.

    Args:
        similarity_matrix: Representative similarity matrix
        representatives: List of representative IDs (in matrix order)
        communities: Dictionary mapping community ID to member IDs
        rep_to_community: Mapping from representative to community ID
        linkage_matrix: Linkage matrix from hierarchical clustering
        output_expanded: Path for expanded matrix output
        logger: Logger instance

    Returns:
        Tuple of (expanded similarity matrix, ordered list of all IDs)
    """
    # Get optimal leaf ordering from clustering
    from scipy.cluster.hierarchy import optimal_leaf_ordering

    logger.info("Computing optimal leaf ordering from hierarchical clustering")
    linkage_ordered = optimal_leaf_ordering(linkage_matrix, squareform(1 - similarity_matrix))

    # Get the dendrogram to extract leaf order
    with plt.ioff():
        dend = dendrogram(linkage_ordered, no_plot=True)
        rep_order = dend['leaves']

    # Reorder representatives according to clustering
    ordered_representatives = [representatives[i] for i in rep_order]

    # Build expanded ID list maintaining clustering order
    expanded_ids = []
    rep_to_idx = {}  # Track where each representative's block starts

    for idx, rep in enumerate(ordered_representatives):
        comm_id = rep_to_community[rep]
        members = communities[comm_id]
        rep_to_idx[rep] = len(expanded_ids)
        expanded_ids.extend(members)

    n_expanded = len(expanded_ids)
    logger.info(f"Expanding matrix from {len(representatives)}x{len(representatives)} "
                f"to {n_expanded}x{n_expanded}")

    # Create expanded matrix (using dense for now, will convert to sparse)
    expanded_matrix = np.zeros((n_expanded, n_expanded), dtype=np.float32)

    # Fill in values
    current_row = 0
    for i, rep_i in enumerate(ordered_representatives):
        comm_i = rep_to_community[rep_i]
        members_i = communities[comm_i]
        n_members_i = len(members_i)

        current_col = 0
        for j, rep_j in enumerate(ordered_representatives):
            comm_j = rep_to_community[rep_j]
            members_j = communities[comm_j]
            n_members_j = len(members_j)

            # Get similarity between representatives
            sim_value = similarity_matrix[rep_order[i], rep_order[j]]

            # Broadcast to all members of both communities
            expanded_matrix[current_row:current_row + n_members_i,
            current_col:current_col + n_members_j] = sim_value

            current_col += n_members_j

        current_row += n_members_i

    # Convert to sparse and save
    sparse_matrix = csr_matrix(expanded_matrix)
    save_npz(output_expanded, sparse_matrix)
    logger.info(f"Saved expanded sparse matrix to {output_expanded}")
    logger.info(f"Matrix sparsity: {1 - sparse_matrix.nnz / (n_expanded ** 2):.4f}")

    return expanded_matrix, expanded_ids


def create_clustermap(distance_matrix: np.ndarray,
                      sample_ids: List[str],
                      linkage_matrix: np.ndarray,
                      output_plot: Path,
                      max_samples: int,
                      logger: logging.Logger) -> None:
    """
    Create and save clustermap visualization.

    Args:
        distance_matrix: Distance matrix (1 - Jaccard)
        sample_ids: Ordered list of sample IDs
        linkage_matrix: Precomputed linkage matrix
        output_plot: Path for output plot
        max_samples: Maximum number of samples to plot (for performance)
        logger: Logger instance
    """
    n_samples = distance_matrix.shape[0]

    # Downsample if necessary
    if n_samples > max_samples:
        logger.info(f"Downsampling from {n_samples} to {max_samples} samples for visualization")
        step = n_samples // max_samples
        indices = np.arange(0, n_samples, step)[:max_samples]
        distance_matrix = distance_matrix[np.ix_(indices, indices)]
        sample_ids = [sample_ids[i] for i in indices]
        # Need to recompute linkage for downsampled data
        from scipy.spatial.distance import squareform
        condensed = squareform(distance_matrix, checks=False)
        linkage_matrix = linkage(condensed, method='average')

    logger.info(f"Creating clustermap with {len(sample_ids)} samples")

    # Create clustermap with precomputed linkage
    g = sns.clustermap(
        distance_matrix,
        row_linkage=linkage_matrix,
        col_linkage=linkage_matrix,
        cmap='viridis',
        xticklabels=False,
        yticklabels=False,
        cbar_kws={'label': 'Distance (1 - Jaccard)'},
        figsize=(16, 14)
    )

    g.fig.suptitle(f'Hierarchical Clustering of {n_samples} Samples\n(Distance = 1 - Jaccard)',
                   y=0.98, fontsize=14)

    plt.savefig(output_plot, dpi=1500, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved clustermap to {output_plot}")


def main():
    parser = argparse.ArgumentParser(
        description='Community-based hierarchical clustering for Jaccard similarity data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--input-dir',
        type=Path,
        required=True,
        help='Directory containing *_jaccard_1.csv files'
    )

    parser.add_argument(
        '--output-dir',
        type=Path,
        required=True,
        help='Directory for output files'
    )

    parser.add_argument(
        '--limit',
        default='all',
        help='Number of CSV files to process (default: all)'
    )

    parser.add_argument(
        '--query-pc-mat',
        default='/scratch/dmk333_new/Logan/Logan_Analyses/metagenome_vector_sketches/build/query_pc_mat',
        help='Path to query_pc_mat executable'
    )

    parser.add_argument(
        '--matrix-dir',
        default='/scratch/mgs_project/matrix_unzipped/',
        help='Path to matrix directory for query_pc_mat'
    )

    parser.add_argument(
        '--db-dir',
        default='/scratch/mgs_project/db/',
        help='Path to database directory for query_pc_mat'
    )

    parser.add_argument(
        '--max-plot-samples',
        type=int,
        default=5000,
        help='Maximum samples to include in clustermap (for performance)'
    )

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    # Setup
    logger = setup_logging(args.verbose)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 80)
    logger.info("Community-based Hierarchical Clustering")
    logger.info("=" * 80)

    # Define output paths
    communities_csv = args.output_dir / "communities.csv"
    communities_json = args.output_dir / "communities.json"
    representatives_file = args.output_dir / "representatives.txt"
    rep_similarity_matrix = args.output_dir / "pw_of_representatives_jaccard_1_no_diag.npy"
    expanded_matrix_file = args.output_dir / "expanded_similarity_matrix.npz"
    clustermap_plot = args.output_dir / "clustermap.png"

    # Step 1: Load CSV files
    logger.info("Step 1: Loading CSV files")
    adjacency = load_csv_files(args.input_dir, args.limit, logger)

    # Step 2: Find communities
    logger.info("Step 2: Identifying communities")
    communities = find_communities(adjacency, logger)

    # Step 3: Save community outputs
    logger.info("Step 3: Saving community assignments")
    save_community_outputs(communities, communities_csv, communities_json, logger)

    # Step 4: Select representatives
    logger.info("Step 4: Selecting community representatives")
    representatives, rep_to_community = select_representatives(
        communities, representatives_file, logger
    )

    # Step 5: Compute representative similarities
    logger.info("Step 5: Computing pairwise similarities between representatives")
    compute_representative_similarities(
        representatives_file,
        rep_similarity_matrix,
        args.query_pc_mat,
        args.matrix_dir,
        args.db_dir,
        logger
    )

    # Step 6: Load similarity matrix and perform clustering
    logger.info("Step 6: Loading similarity matrix and performing hierarchical clustering")
    similarity_matrix = np.load(rep_similarity_matrix)
    linkage_matrix = perform_hierarchical_clustering(similarity_matrix, logger)

    # Step 7: Expand matrix
    logger.info("Step 7: Expanding matrix to include all community members")
    expanded_matrix, expanded_ids = expand_matrix(
        similarity_matrix,
        representatives,
        communities,
        rep_to_community,
        linkage_matrix,
        expanded_matrix_file,
        logger
    )

    # Save ordered sample IDs for future plotting
    expanded_ids_file = args.output_dir / "expanded_sample_ids.txt"
    with open(expanded_ids_file, 'w') as f:
        for sample_id in expanded_ids:
            f.write(f"{sample_id}\n")
    logger.info(f"Saved ordered sample IDs to {expanded_ids_file}")

    # Step 8: Create clustermap
    logger.info("Step 8: Creating clustermap visualization")
    distance_matrix = 1 - expanded_matrix

    # Recompute linkage on the full expanded matrix for visualization
    from scipy.spatial.distance import squareform
    condensed_dist = squareform(distance_matrix, checks=False)
    expanded_linkage = linkage(condensed_dist, method='average')

    # Save linkage matrix for future plotting
    expanded_linkage_file = args.output_dir / "expanded_linkage_matrix.npy"
    np.save(expanded_linkage_file, expanded_linkage)
    logger.info(f"Saved linkage matrix to {expanded_linkage_file}")

    create_clustermap(
        distance_matrix,
        expanded_ids,
        expanded_linkage,
        clustermap_plot,
        args.max_plot_samples,
        logger
    )

    logger.info("=" * 80)
    logger.info("Pipeline completed successfully")
    logger.info(f"Output files in: {args.output_dir}")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()