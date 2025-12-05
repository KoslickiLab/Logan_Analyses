#!/bin/bash

# For each sample, get all other samples with jaccard = 1 (besides itself)
# nohup ./2_parallel_query_jaccard_1.sh > 2_parallel_query_jaccard_1.log 2>&1

# Get some stats about how big are the groupings of such duplicate samples 
#nohup ./2.1.analysis.matches_per_csv.sh > 2.1.analysis.matches_per_csv.log 2>&1

# Get all the unique accessions in all of these
#nohup ./3_collect_accessions_with_other_jaccard_1.sh > 3_collect_accessions_with_other_jaccard_1.log 2>&1

# Get all pw Jaccards between all of these samples
# nohup ./4_pairwise_of_acc_with_other_jaccard_1.sh > 4_pairwise_of_acc_with_other_jaccard_1.log 2>&1

# Convert this to sparse
# nohup ./5_convert_to_sparse.sh > 5_convert_to_sparse.log 2>&1

# Cluster, group, and visualize these.
nohup python  analyze_sparse_jaccard_matrix.py ../data/pw_of_acc_with_other_jaccard_1.npz ../data/pw_jaccard_clustering_of_duplicates > 6_analyze_sparse_jaccard_matrix.log 2>&1
