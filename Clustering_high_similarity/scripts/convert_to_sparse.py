#!/usr/bin/env python3
"""
Convert dense Jaccard matrix to sparse format.
"""

import numpy as np
import scipy.sparse as sp
import sys

if len(sys.argv) != 3:
    print("Usage: python convert_to_sparse.py <input.npz> <output_sparse.npz>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

print("Loading dense matrix...")
data = np.load(input_file, allow_pickle=False)

if isinstance(data, np.ndarray):
    matrix = data
else:
    matrix = data[list(data.keys())[0]]

print(f"Matrix shape: {matrix.shape}")
print(f"Matrix size in memory: {matrix.nbytes / (1024**3):.2f} GB")

# Count zeros
n_zeros = np.sum(np.abs(matrix) < 1e-10)
sparsity = 100 * n_zeros / matrix.size
print(f"Sparsity: {sparsity:.2f}%")

print("\nConverting to sparse format (CSR)...")
sparse_matrix = sp.csr_matrix(matrix)

print(f"Sparse matrix non-zero elements: {sparse_matrix.nnz:,}")
print(f"Compression ratio: {matrix.nbytes / (sparse_matrix.data.nbytes + sparse_matrix.indices.nbytes + sparse_matrix.indptr.nbytes):.1f}x")

print(f"\nSaving to {output_file}...")
sp.save_npz(output_file, sparse_matrix)

import os
output_size = os.path.getsize(output_file) / (1024**3)
print(f"✓ Saved! Output file size: {output_size:.2f} GB")
print(f"Space saved: {(os.path.getsize(input_file) - os.path.getsize(output_file)) / (1024**3):.2f} GB")