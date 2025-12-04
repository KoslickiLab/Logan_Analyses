#!/usr/bin/env python3
"""
Optimized analysis for SPARSE Jaccard matrix with visualizations.
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
import json
import pickle

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


def plot_value_distribution(values, output_dir="/mnt/user-data/outputs"):
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

    print(f"\nExtracting edges with similarity > {threshold}...")
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
        print("No similarity edges found!")
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


def save_edges_csv(edges, output_dir="/mnt/user-data/outputs"):
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
    save_components_ragged(components)

    return components


def save_components_csv(components, output_dir="/mnt/user-data/outputs"):
    """Save component membership to CSV (long format)."""
    output_file = os.path.join(output_dir, "component_membership.csv")
    print(f"\nSaving component membership to {output_file}...")

    with open(output_file, 'w') as f:
        f.write("sample_id,component_id,component_size\n")
        for comp_id, comp in enumerate(components):
            for node in sorted(comp):
                f.write(f"{node},{comp_id},{len(comp)}\n")

    total_samples = sum(len(c) for c in components)
    print(f"✓ Saved {total_samples:,} sample-component pairs")


def save_components_ragged(components, output_dir="/mnt/user-data/outputs"):
    """Save components in ragged CSV format (components as columns)."""
    output_file = os.path.join(output_dir, "components_ragged.csv")
    print(f"\nSaving ragged component format to {output_file}...")

    # Sort components by size (largest first)
    sorted_components = sorted(components, key=len, reverse=True)

    # Convert to lists and pad
    max_size = max(len(c) for c in sorted_components)

    with open(output_file, 'w') as f:
        # Header
        header = ','.join([f'component_{i}' for i in range(len(sorted_components))])
        f.write(header + '\n')

        # Data rows
        for row_idx in range(max_size):
            row = []
            for comp in sorted_components:
                comp_list = sorted(list(comp))
                if row_idx < len(comp_list):
                    row.append(str(comp_list[row_idx]))
                else:
                    row.append('')  # Empty for shorter components
            f.write(','.join(row) + '\n')

    print(f"✓ Saved {len(sorted_components)} components in ragged format")

    # Also save as JSON for easier parsing
    json_file = os.path.join(output_dir, "components.json")
    components_dict = {
        f"component_{i}": sorted(list(comp))
        for i, comp in enumerate(sorted_components)
    }
    with open(json_file, 'w') as f:
        json.dump(components_dict, f, indent=2)
    print(f"✓ Saved components to {json_file}")


def visualize_components(edges, components, output_dir="/mnt/user-data/outputs", max_components=10,
                         max_nodes_per_plot=100):
    """Visualize connected components with edge weights."""
    if not NETWORKX_AVAILABLE:
        print("\nSkipping component visualization (networkx not available)")
        return

    print("\n" + "=" * 80)
    print("VISUALIZING CONNECTED COMPONENTS")
    print("=" * 80)

    # Build graph
    G = nx.Graph()
    G.add_weighted_edges_from(edges)

    # Sort components by size
    sorted_components = sorted(components, key=len, reverse=True)

    # Visualize largest components
    n_to_plot = min(max_components, len(sorted_components))
    print(f"\nCreating visualizations for {n_to_plot} largest components...")

    for idx, component in enumerate(sorted_components[:n_to_plot]):
        comp_size = len(component)

        # Skip if too large
        if comp_size > max_nodes_per_plot:
            print(f"  Component {idx} (size {comp_size}): Too large, skipping visualization")
            continue

        print(f"  Visualizing component {idx} (size {comp_size})...")

        # Extract subgraph
        subgraph = G.subgraph(component).copy()

        # Get edge weights
        edges_in_comp = list(subgraph.edges(data=True))
        weights = np.array([e[2]['weight'] for e in edges_in_comp])

        # Normalize weights for visualization
        if len(weights) > 0:
            min_weight = weights.min()
            max_weight = weights.max()

            # Edge widths (1-5 range)
            if max_weight > min_weight:
                edge_widths = 1 + 4 * (weights - min_weight) / (max_weight - min_weight)
            else:
                edge_widths = np.full(len(weights), 2.5)

            # Edge alphas (0.3-1.0 range)
            if max_weight > min_weight:
                edge_alphas = 0.3 + 0.7 * (weights - min_weight) / (max_weight - min_weight)
            else:
                edge_alphas = np.full(len(weights), 0.65)
        else:
            edge_widths = []
            edge_alphas = []

        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))

        # Layout
        if comp_size <= 20:
            pos = nx.spring_layout(subgraph, k=2, iterations=50, seed=42)
        else:
            pos = nx.spring_layout(subgraph, k=1, iterations=30, seed=42)

        # Draw edges with varying width and alpha
        for (u, v, data), width, alpha in zip(edges_in_comp, edge_widths, edge_alphas):
            nx.draw_networkx_edges(
                subgraph, pos,
                edgelist=[(u, v)],
                width=width,
                alpha=alpha,
                edge_color='steelblue',
                ax=ax
            )

        # Draw nodes
        nx.draw_networkx_nodes(
            subgraph, pos,
            node_size=300,
            node_color='lightcoral',
            edgecolors='black',
            linewidths=1.5,
            ax=ax
        )

        # Draw labels
        if comp_size <= 50:
            nx.draw_networkx_labels(
                subgraph, pos,
                font_size=8,
                font_weight='bold',
                ax=ax
            )

        # Title and info
        if len(weights) > 0:
            title = f'Component {idx} ({comp_size} nodes, {len(edges_in_comp)} edges)\n'
            title += f'Jaccard range: [{min_weight:.3f}, {max_weight:.3f}]'
        else:
            title = f'Component {idx} ({comp_size} nodes, isolated)'

        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.axis('off')

        # Add legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color='steelblue', linewidth=1, alpha=0.3, label='Low similarity'),
            Line2D([0], [0], color='steelblue', linewidth=3, alpha=0.65, label='Medium similarity'),
            Line2D([0], [0], color='steelblue', linewidth=5, alpha=1.0, label='High similarity'),
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

        plt.tight_layout()
        output_file = os.path.join(output_dir, f"component_{idx}_visualization.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"    Saved: {output_file}")

    print(
        f"\n✓ Created {min(n_to_plot, sum(1 for c in sorted_components[:n_to_plot] if len(c) <= max_nodes_per_plot))} component visualizations")


def network_analysis(edges):
    """Perform network-based community detection."""
    if not NETWORKX_AVAILABLE:
        print("\nSkipping network analysis (networkx not available)")
        return None

    print("\n" + "=" * 80)
    print("NETWORK COMMUNITY DETECTION")
    print("=" * 80)

    print("\nBuilding NetworkX graph...")
    G = nx.Graph()
    G.add_weighted_edges_from(edges)

    print(f"Graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

    # Save graph for future use
    save_graph(G)

    # Analyze largest component
    if not nx.is_connected(G):
        print("\nGraph is disconnected, analyzing largest component...")
        largest_cc = max(nx.connected_components(G), key=len)
        G_largest = G.subgraph(largest_cc).copy()
    else:
        G_largest = G

    print(f"Analyzing component: {G_largest.number_of_nodes():,} nodes")

    # Louvain community detection
    communities = None
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

        # Save communities
        save_communities(communities)

    except Exception as e:
        print(f"Louvain failed: {e}")

    return communities


def save_communities(communities, output_dir="/mnt/user-data/outputs"):
    """Save community assignments."""
    # Long format CSV
    output_file = os.path.join(output_dir, "community_membership.csv")
    print(f"\nSaving community membership to {output_file}...")

    with open(output_file, 'w') as f:
        f.write("sample_id,community_id,community_size\n")
        for comm_id, comm in enumerate(communities):
            for node in sorted(comm):
                f.write(f"{node},{comm_id},{len(comm)}\n")

    total_samples = sum(len(c) for c in communities)
    print(f"✓ Saved {total_samples:,} sample-community pairs")

    # Ragged format
    output_file_ragged = os.path.join(output_dir, "communities_ragged.csv")
    print(f"Saving ragged community format to {output_file_ragged}...")

    sorted_communities = sorted(communities, key=len, reverse=True)
    max_size = max(len(c) for c in sorted_communities)

    with open(output_file_ragged, 'w') as f:
        header = ','.join([f'community_{i}' for i in range(len(sorted_communities))])
        f.write(header + '\n')

        for row_idx in range(max_size):
            row = []
            for comm in sorted_communities:
                comm_list = sorted(list(comm))
                if row_idx < len(comm_list):
                    row.append(str(comm_list[row_idx]))
                else:
                    row.append('')
            f.write(','.join(row) + '\n')

    print(f"✓ Saved {len(sorted_communities)} communities in ragged format")

    # JSON format
    json_file = os.path.join(output_dir, "communities.json")
    communities_dict = {
        f"community_{i}": sorted(list(comm))
        for i, comm in enumerate(sorted_communities)
    }
    with open(json_file, 'w') as f:
        json.dump(communities_dict, f, indent=2)
    print(f"✓ Saved communities to {json_file}")


def save_graph(G, output_dir="/mnt/user-data/outputs"):
    """Save NetworkX graph for future analysis."""
    output_file = os.path.join(output_dir, "graph.pickle")
    print(f"\nSaving graph to {output_file}...")

    with open(output_file, 'wb') as f:
        pickle.dump(G, f)

    print(f"✓ Saved graph with {G.number_of_nodes():,} nodes and {G.number_of_edges():,} edges")


def hierarchical_clustering_sample(matrix, edges, max_samples=10000, output_dir="/mnt/user-data/outputs"):
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

    # Save sample indices
    sample_file = os.path.join(output_dir, "hierarchical_clustering_samples.txt")
    with open(sample_file, 'w') as f:
        for sample in top_samples:
            f.write(f"{sample}\n")
    print(f"Saved sample indices to {sample_file}")

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

    # Save linkage matrix
    linkage_file = os.path.join(output_dir, "linkage_matrix.npy")
    np.save(linkage_file, linkage_matrix)
    print(f"Saved linkage matrix to {linkage_file}")

    # Cut at various thresholds
    print("\nClusters at different distance thresholds:")
    threshold_results = {}
    for threshold in [0.01, 0.05, 0.1, 0.2, 0.5]:
        clusters = fcluster(linkage_matrix, threshold, criterion='distance')
        n_clusters = len(np.unique(clusters))
        print(f"  Distance {threshold:.2f}: {n_clusters:5d} clusters")
        threshold_results[threshold] = clusters

    # Save clustering results
    clusters_file = os.path.join(output_dir, "hierarchical_clusters.csv")
    with open(clusters_file, 'w') as f:
        f.write("sample_id," + ",".join([f"threshold_{t}" for t in [0.01, 0.05, 0.1, 0.2, 0.5]]) + "\n")
        for i, sample_id in enumerate(top_samples):
            row = [str(sample_id)]
            for threshold in [0.01, 0.05, 0.1, 0.2, 0.5]:
                row.append(str(threshold_results[threshold][i]))
            f.write(",".join(row) + "\n")
    print(f"Saved clustering results to {clusters_file}")

    # Plot dendrogram
    print("\nCreating dendrogram visualization...")
    plot_dendrogram(linkage_matrix, top_samples, output_dir)

    return linkage_matrix


def plot_dendrogram(linkage_matrix, sample_ids, output_dir="/mnt/user-data/outputs", max_labels=50):
    """Plot hierarchical clustering dendrogram."""
    fig, ax = plt.subplots(figsize=(16, 10))

    # If too many samples, don't show labels
    if len(sample_ids) > max_labels:
        dendrogram(
            linkage_matrix,
            ax=ax,
            no_labels=True,
            color_threshold=0.2,
            above_threshold_color='gray'
        )
        ax.set_xlabel('Sample Index (labels omitted for clarity)', fontsize=12)
    else:
        dendrogram(
            linkage_matrix,
            labels=sample_ids,
            ax=ax,
            leaf_font_size=8,
            color_threshold=0.2,
            above_threshold_color='gray'
        )
        ax.set_xlabel('Sample ID', fontsize=12)
        plt.xticks(rotation=90)

    ax.set_ylabel('Distance (1 - Jaccard)', fontsize=12)
    ax.set_title(f'Hierarchical Clustering Dendrogram ({len(sample_ids)} samples)',
                 fontsize=14, fontweight='bold')
    ax.axhline(y=0.2, color='red', linestyle='--', linewidth=1, alpha=0.5, label='Distance threshold = 0.2')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    output_file = os.path.join(output_dir, "dendrogram.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved dendrogram to {output_file}")
    plt.close()


def main(filepath):
    """Main analysis pipeline."""
    print("=" * 80)
    print("SPARSE JACCARD MATRIX ANALYSIS WITH VISUALIZATIONS")
    print("=" * 80)
    print(f"\nInput: {filepath}\n")

    # Ensure output directory exists
    output_dir = "/mnt/user-data/outputs"
    os.makedirs(output_dir, exist_ok=True)

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

    # Visualize components
    visualize_components(edges, components, output_dir, max_components=20, max_nodes_per_plot=100)

    # Network analysis and community detection
    communities = network_analysis(edges)

    # Hierarchical clustering on sample
    hierarchical_clustering_sample(matrix, edges, max_samples=10000, output_dir=output_dir)

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
    if communities:
        print(f"  - Communities detected: {len(communities)}")

    print(f"\nOutputs saved to: {output_dir}/")
    print("\nGenerated files:")
    print("  - intermediate_similarity_edges.csv")
    print("  - component_membership.csv")
    print("  - components_ragged.csv")
    print("  - components.json")
    print("  - community_membership.csv (if communities found)")
    print("  - communities_ragged.csv (if communities found)")
    print("  - communities.json (if communities found)")
    print("  - component_*_visualization.png (for visualizable components)")
    print("  - dendrogram.png")
    print("  - hierarchical_clusters.csv")
    print("  - linkage_matrix.npy")
    print("  - graph.pickle (NetworkX graph object)")
    print("  - value_distribution.png")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python analyze_sparse_jaccard_matrix.py <path_to_sparse_matrix.npz>")
        sys.exit(1)

    filepath = sys.argv[1]
    main(filepath)