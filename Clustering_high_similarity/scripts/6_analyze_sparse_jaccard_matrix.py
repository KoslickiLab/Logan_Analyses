#!/usr/bin/env python3
"""
Optimized analysis for SPARSE Jaccard matrix.
Works with scipy.sparse CSR/CSC format matrices.
"""

import numpy as np
import scipy.sparse as sp
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import time
import os

# Optional advanced clustering
try:
    import networkx as nx
    from networkx.algorithms import community

    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    print("Note: networkx not available (install for community detection)")

try:
    from sklearn.cluster import AgglomerativeClustering, DBSCAN
    from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
    from scipy.spatial.distance import squareform

    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    print("Note: sklearn not available (install for hierarchical clustering)")

sns.set_style("whitegrid")


def load_sparse_matrix(filepath):
    """Load sparse matrix file."""
    print("=" * 80)
    print("LOADING SPARSE MATRIX")
    print("=" * 80)
    print(f"\nFile: {filepath}")
    print(f"Size: {os.path.getsize(filepath) / (1024 ** 3):.2f} GB")

    print("\nLoading sparse matrix...")
    start_time = time.time()

    matrix = sp.load_npz(filepath)

    elapsed = time.time() - start_time

    print(f"\n✓ Loaded in {elapsed:.1f} seconds")
    print(f"  Shape: {matrix.shape}")
    print(f"  Format: {type(matrix).__name__}")
    print(f"  Dtype: {matrix.dtype}")
    print(f"  Non-zero elements: {matrix.nnz:,}")

    sparsity = 100 * (1 - matrix.nnz / (matrix.shape[0] * matrix.shape[1]))
    print(f"  Sparsity: {sparsity:.4f}%")

    # Verify it's a square matrix
    if matrix.shape[0] != matrix.shape[1]:
        print(f"  WARNING: Matrix is not square!")
    else:
        print(f"  ✓ Square matrix: {matrix.shape[0]:,} × {matrix.shape[1]:,}")

    # Convert to CSR for efficient row operations
    if not isinstance(matrix, sp.csr_matrix):
        print(f"  Converting to CSR format...")
        matrix = matrix.tocsr()

    return matrix


def analyze_value_distribution(matrix):
    """Analyze the distribution of values in the sparse matrix."""
    print("\n" + "=" * 80)
    print("VALUE DISTRIBUTION ANALYSIS")
    print("=" * 80)

    print("\nAnalyzing sparse matrix values...")
    start_time = time.time()

    # For sparse matrices, we analyze the explicit values
    data_values = matrix.data

    total_elements = matrix.shape[0] * matrix.shape[1]
    n_explicit = len(data_values)
    n_implicit_zeros = total_elements - n_explicit

    print(f"Total elements: {total_elements:,}")
    print(f"Explicit non-zeros: {n_explicit:,}")
    print(f"Implicit zeros: {n_implicit_zeros:,}")

    # Count special values in explicit entries
    n_ones = np.sum(np.abs(data_values - 1.0) < 1e-10)
    n_explicit_zeros = np.sum(np.abs(data_values) < 1e-10)
    n_other = n_explicit - n_ones - n_explicit_zeros

    # Total zeros includes implicit
    n_total_zeros = n_implicit_zeros + n_explicit_zeros

    elapsed = time.time() - start_time
    print(f"Computed in {elapsed:.1f} seconds")

    print(f"\nValue counts:")
    print(f"  Zeros (0.0):  {n_total_zeros:15,} ({100 * n_total_zeros / total_elements:6.3f}%)")
    print(f"  Ones (1.0):   {n_ones:15,} ({100 * n_ones / total_elements:6.3f}%)")
    print(f"  Other values: {n_other:15,} ({100 * n_other / total_elements:6.3f}%)")

    print(f"\n{'*' * 80}")
    print(f"*** NON-ZERO, NON-ONE ENTRIES: {n_other:,} ***")
    print(f"{'*' * 80}")

    if n_other == 0:
        print("\nNo intermediate values found!")
        print("All duplicated samples only connect to their exact duplicates.")
        return False, None

    # Extract and analyze non-zero, non-one values
    print(f"\nExtracting {n_other:,} intermediate values...")
    start_time = time.time()

    mask = (np.abs(data_values - 1.0) >= 1e-10) & (np.abs(data_values) >= 1e-10)
    other_values = data_values[mask]

    elapsed = time.time() - start_time
    print(f"Extracted in {elapsed:.1f} seconds")

    print(f"\nStatistics for intermediate values:")
    print(f"  Count:  {len(other_values):,}")
    print(f"  Min:    {np.min(other_values):.8f}")
    print(f"  Max:    {np.max(other_values):.8f}")
    print(f"  Mean:   {np.mean(other_values):.8f}")
    print(f"  Median: {np.median(other_values):.8f}")
    print(f"  Std:    {np.std(other_values):.8f}")

    # Percentiles
    percentiles = [1, 5, 10, 25, 50, 75, 90, 95, 99]
    print(f"\nPercentiles:")
    for p in percentiles:
        val = np.percentile(other_values, p)
        print(f"  {p:3d}%: {val:.8f}")

    # Value range distribution
    print(f"\nValue range distribution:")
    ranges = [
        (0.0, 0.1, "0.0 - 0.1"),
        (0.1, 0.2, "0.1 - 0.2"),
        (0.2, 0.3, "0.2 - 0.3"),
        (0.3, 0.4, "0.3 - 0.4"),
        (0.4, 0.5, "0.4 - 0.5"),
        (0.5, 0.6, "0.5 - 0.6"),
        (0.6, 0.7, "0.6 - 0.7"),
        (0.7, 0.8, "0.7 - 0.8"),
        (0.8, 0.9, "0.8 - 0.9"),
        (0.9, 1.0, "0.9 - 1.0"),
    ]

    for low, high, label in ranges:
        count = np.sum((other_values >= low) & (other_values < high))
        if count > 0:
            pct = 100 * count / len(other_values)
            print(f"  {label}: {count:12,} ({pct:6.2f}%)")

    # Plot distribution
    plot_value_distribution(other_values)

    return True, other_values


def plot_value_distribution(values, output_dir="/scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/pw_analysis_jacc_1"):
    """Create detailed plots of value distribution."""
    print("\nCreating value distribution plots...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Histogram - linear scale
    axes[0, 0].hist(values, bins=100, edgecolor='black', alpha=0.7, color='steelblue')
    axes[0, 0].set_xlabel('Jaccard Similarity')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title(f'Distribution of Intermediate Similarities (n={len(values):,})')
    axes[0, 0].axvline(np.median(values), color='red', linestyle='--',
                       label=f'Median: {np.median(values):.4f}')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    # 2. Histogram - log scale
    axes[0, 1].hist(values, bins=100, edgecolor='black', alpha=0.7, color='forestgreen')
    axes[0, 1].set_xlabel('Jaccard Similarity')
    axes[0, 1].set_ylabel('Frequency (log scale)')
    axes[0, 1].set_title('Distribution (Log Y-axis)')
    axes[0, 1].set_yscale('log')
    axes[0, 1].grid(True, alpha=0.3)

    # 3. Cumulative distribution
    sorted_vals = np.sort(values)
    cumulative = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
    axes[1, 0].plot(sorted_vals, cumulative, linewidth=2, color='darkorange')
    axes[1, 0].set_xlabel('Jaccard Similarity')
    axes[1, 0].set_ylabel('Cumulative Probability')
    axes[1, 0].set_title('Cumulative Distribution Function')
    axes[1, 0].grid(True, alpha=0.3)

    # 4. Box plot
    axes[1, 1].boxplot([values], vert=False, widths=0.6)
    axes[1, 1].set_xlabel('Jaccard Similarity')
    axes[1, 1].set_title('Box Plot with Quartiles')
    axes[1, 1].set_yticks([])
    axes[1, 1].grid(True, alpha=0.3, axis='x')

    # Add stats
    stats_text = f"""n = {len(values):,}
min = {np.min(values):.6f}
Q1  = {np.percentile(values, 25):.6f}
med = {np.median(values):.6f}
Q3  = {np.percentile(values, 75):.6f}
max = {np.max(values):.6f}"""
    axes[1, 1].text(0.02, 0.98, stats_text, transform=axes[1, 1].transAxes,
                    verticalalignment='top', fontsize=10, family='monospace',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    output_file = os.path.join(output_dir, "value_distribution.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()


def build_edge_list_sparse(matrix, threshold=0.0):
    """Build edge list from sparse matrix efficiently."""
    print("\n" + "=" * 80)
    print("BUILDING EDGE LIST FROM SPARSE MATRIX")
    print("=" * 80)

    print(f"\nExtracting edges with similarity > {threshold}")
    start_time = time.time()

    # Convert to COO format for easy iteration
    matrix_coo = matrix.tocoo()

    edges = []
    n_processed = 0

    print(f"Processing {matrix_coo.nnz:,} non-zero entries...")

    for i, j, val in zip(matrix_coo.row, matrix_coo.col, matrix_coo.data):
        # Only upper triangle to avoid duplicates
        if i < j and val > threshold:
            edges.append((i, j, val))

        n_processed += 1
        if n_processed % 1000000 == 0:
            elapsed = time.time() - start_time
            progress = n_processed / matrix_coo.nnz
            eta = elapsed / progress - elapsed if progress > 0 else 0
            print(f"  Progress: {n_processed:,}/{matrix_coo.nnz:,} ({100 * progress:.1f}%) - "
                  f"ETA: {eta / 60:.1f} min - Edges: {len(edges):,}", end='\r')

    elapsed = time.time() - start_time
    print(f"\n\n✓ Found {len(edges):,} edges in {elapsed:.1f} seconds")

    if len(edges) == 0:
        print("No intermediate similarity edges found!")
        return None

    # Statistics
    weights = np.array([e[2] for e in edges])
    print(f"\nEdge weight statistics:")
    print(f"  Min:    {np.min(weights):.8f}")
    print(f"  Max:    {np.max(weights):.8f}")
    print(f"  Mean:   {np.mean(weights):.8f}")
    print(f"  Median: {np.median(weights):.8f}")

    # Save to CSV
    save_edges_csv(edges)

    return edges


def save_edges_csv(edges, output_dir="/scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/pw_analysis_jacc_1"):
    """Save edge list to CSV."""
    output_file = os.path.join(output_dir, "intermediate_similarity_edges.csv")
    print(f"\nSaving edges to {output_file}...")

    with open(output_file, 'w') as f:
        f.write("sample_i,sample_j,jaccard_similarity\n")
        for i, j, v in edges:
            f.write(f"{i},{j},{v:.8f}\n")

    print(f"✓ Saved {len(edges):,} edges")


def analyze_connectivity(edges, n_samples):
    """Analyze connectivity patterns."""
    print("\n" + "=" * 80)
    print("CONNECTIVITY ANALYSIS")
    print("=" * 80)

    # Build adjacency list
    print("\nBuilding adjacency list...")
    adj_list = defaultdict(list)
    for i, j, v in edges:
        adj_list[i].append((j, v))
        adj_list[j].append((i, v))

    print(f"Samples with connections: {len(adj_list):,} / {n_samples:,} ({100 * len(adj_list) / n_samples:.2f}%)")

    # Degree distribution
    degrees = [len(neighbors) for neighbors in adj_list.values()]

    print(f"\nDegree distribution:")
    print(f"  Min:    {min(degrees)}")
    print(f"  Max:    {max(degrees)}")
    print(f"  Mean:   {np.mean(degrees):.2f}")
    print(f"  Median: {np.median(degrees):.0f}")

    degree_counts = Counter(degrees)
    print(f"\nTop degree values:")
    for degree, count in degree_counts.most_common(20):
        print(f"  Degree {degree:5d}: {count:8,} samples")

    return adj_list


def find_connected_components(adj_list):
    """Find connected components using iterative DFS (avoids recursion limit)."""
    print("\n" + "=" * 80)
    print("CONNECTED COMPONENTS")
    print("=" * 80)

    print("\nFinding connected components...")
    start_time = time.time()

    nodes = set(adj_list.keys())
    visited = set()
    components = []

    for start_node in nodes:
        if start_node in visited:
            continue

        # Iterative DFS using a stack
        component = set()
        stack = [start_node]

        while stack:
            node = stack.pop()
            if node in visited:
                continue

            visited.add(node)
            component.add(node)

            # Add unvisited neighbors to stack
            for neighbor, _ in adj_list.get(node, []):
                if neighbor not in visited:
                    stack.append(neighbor)

        components.append(component)

    elapsed = time.time() - start_time
    print(f"✓ Found {len(components):,} components in {elapsed:.1f} seconds")

    # Statistics
    sizes = [len(c) for c in components]
    print(f"\nComponent size statistics:")
    print(f"  Min:    {min(sizes)}")
    print(f"  Max:    {max(sizes)}")
    print(f"  Mean:   {np.mean(sizes):.2f}")
    print(f"  Median: {np.median(sizes):.0f}")

    size_counts = Counter(sizes)
    print(f"\nSize distribution:")
    for size in sorted(size_counts.keys())[:25]:
        count = size_counts[size]
        print(f"  Size {size:5d}: {count:8,} components")
    if len(size_counts) > 25:
        print(f"  ... ({len(size_counts) - 25} more unique sizes)")

    print(f"\nLargest 10 components:")
    largest = sorted(components, key=len, reverse=True)[:10]
    for i, comp in enumerate(largest, 1):
        print(f"  {i:2d}. {len(comp):8,} samples")

    # Save component membership
    save_components_csv(components)

    return components


def save_components_csv(components, output_dir="/scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/pw_analysis_jacc_1"):
    """Save component membership to CSV."""
    output_file = os.path.join(output_dir, "component_membership.csv")
    print(f"\nSaving component membership to {output_file}...")

    with open(output_file, 'w') as f:
        f.write("sample_id,component_id,component_size\n")
        for comp_id, comp in enumerate(components):
            for node in sorted(comp):
                f.write(f"{node},{comp_id},{len(comp)}\n")

    total_samples = sum(len(c) for c in components)
    print(f"✓ Saved {total_samples:,} sample-component pairs")


def network_analysis(edges):
    """Perform network-based community detection."""
    if not NETWORKX_AVAILABLE:
        print("\nSkipping network analysis (networkx not available)")
        return

    print("\n" + "=" * 80)
    print("NETWORK COMMUNITY DETECTION")
    print("=" * 80)

    print("\nBuilding NetworkX graph...")
    G = nx.Graph()
    G.add_weighted_edges_from(edges)

    print(f"Graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

    # Analyze largest component
    if not nx.is_connected(G):
        print("\nGraph is disconnected, analyzing largest component...")
        largest_cc = max(nx.connected_components(G), key=len)
        G_largest = G.subgraph(largest_cc).copy()
    else:
        G_largest = G

    print(f"Analyzing component: {G_largest.number_of_nodes():,} nodes")

    # Louvain community detection
    try:
        print("\nRunning Louvain community detection...")
        start_time = time.time()
        communities = community.louvain_communities(G_largest, resolution=1.0, seed=42)
        elapsed = time.time() - start_time

        print(f"✓ Found {len(communities)} communities in {elapsed:.1f} seconds")

        sizes = [len(c) for c in communities]
        print(f"\nCommunity sizes:")
        print(f"  Min:    {min(sizes)}")
        print(f"  Max:    {max(sizes)}")
        print(f"  Mean:   {np.mean(sizes):.2f}")
        print(f"  Median: {np.median(sizes):.0f}")

        # Modularity
        mod = community.modularity(G_largest, communities)
        print(f"\nModularity: {mod:.6f}")

    except Exception as e:
        print(f"Louvain failed: {e}")


def hierarchical_clustering_sample(matrix, edges, max_samples=10000):
    """Perform hierarchical clustering on a sample of highly connected nodes."""
    if not SKLEARN_AVAILABLE:
        print("\nSkipping hierarchical clustering (sklearn/scipy not available)")
        return

    print("\n" + "=" * 80)
    print("HIERARCHICAL CLUSTERING (SAMPLE)")
    print("=" * 80)

    # Find most connected samples
    print(f"\nSelecting up to {max_samples} most connected samples...")
    sample_degrees = defaultdict(int)
    for i, j, _ in edges:
        sample_degrees[i] += 1
        sample_degrees[j] += 1

    # Top samples by degree
    top_samples = sorted(sample_degrees.keys(),
                         key=lambda x: sample_degrees[x],
                         reverse=True)[:max_samples]

    print(f"Selected {len(top_samples)} samples")

    # Extract submatrix - handle sparse matrix
    print("Extracting submatrix...")
    sub_matrix = matrix[np.ix_(top_samples, top_samples)]

    # Convert to dense for hierarchical clustering
    if sp.issparse(sub_matrix):
        print("Converting to dense...")
        sub_matrix = sub_matrix.toarray()

    # Convert similarity to distance
    distance_matrix = 1 - sub_matrix
    np.fill_diagonal(distance_matrix, 0)

    # Hierarchical clustering
    print("Performing hierarchical clustering...")
    start_time = time.time()

    condensed_dist = squareform(distance_matrix, checks=False)
    linkage_matrix = linkage(condensed_dist, method='average')

    elapsed = time.time() - start_time
    print(f"✓ Completed in {elapsed:.1f} seconds")

    # Cut at various thresholds
    print("\nClusters at different distance thresholds:")
    for threshold in [0.01, 0.05, 0.1, 0.2, 0.5]:
        clusters = fcluster(linkage_matrix, threshold, criterion='distance')
        n_clusters = len(np.unique(clusters))
        print(f"  Distance {threshold:.2f}: {n_clusters:5d} clusters")


def main(filepath):
    """Main analysis pipeline."""
    print("=" * 80)
    print("SPARSE JACCARD MATRIX ANALYSIS")
    print("=" * 80)
    print(f"\nInput: {filepath}\n")

    # Ensure output directory exists
    os.makedirs("/scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/pw_analysis_jacc_1", exist_ok=True)

    # Load sparse matrix
    matrix = load_sparse_matrix(filepath)

    # Analyze values
    has_intermediate, other_values = analyze_value_distribution(matrix)

    if not has_intermediate:
        print("\n" + "=" * 80)
        print("CONCLUSION: No intermediate similarities")
        print("=" * 80)
        return

    # Build edge list
    edges = build_edge_list_sparse(matrix, threshold=0.0)

    if edges is None:
        return

    # Connectivity analysis
    adj_list = analyze_connectivity(edges, matrix.shape[0])

    # Find connected components
    components = find_connected_components(adj_list)

    # Network analysis
    network_analysis(edges)

    # Hierarchical clustering on sample
    hierarchical_clustering_sample(matrix, edges, max_samples=10000)

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"\nKey Results:")
    print(f"  - Matrix: {matrix.shape[0]:,} × {matrix.shape[1]:,}")
    print(f"  - Non-zero, non-one entries: {len(other_values):,}")
    print(f"  - Intermediate similarity edges: {len(edges):,}")
    print(f"  - Samples with connections: {len(adj_list):,}")
    print(f"  - Connected components: {len(components):,}")
    print(f"  - Largest component: {max(len(c) for c in components):,} samples")
    print(f"\nOutputs saved to: /scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/pw_analysis_jacc_1")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python analyze_sparse_jaccard_matrix.py <path_to_sparse_matrix.npz>")
        sys.exit(1)

    filepath = sys.argv[1]
    main(filepath)