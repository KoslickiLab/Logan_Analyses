#!/bin/bash
cd /scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/output
/scratch/dmk333_new/Logan/Logan_Analyses/metagenome_vector_sketches/build/./query_pc_mat --matrix /scratch/mgs_project/matrix_unzipped/ --db /scratch/mgs_project/db/ --query_file /scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/all_acc_100000.txt --write_to_file test.csv --no_self --min_jaccard 1
