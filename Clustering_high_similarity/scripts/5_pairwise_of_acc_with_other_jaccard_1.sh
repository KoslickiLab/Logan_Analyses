#!/bin/bash
ACC_WITH_JACC_1=/scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/accessions_mbases_geq_10_with_other_Jaccard_1.txt
/scratch/dmk333_new/Logan/Logan_Analyses/metagenome_vector_sketches/build/./query_pc_mat --matrix /scratch/mgs_project/matrix_unzipped/ --db /scratch/mgs_project/db/ --row_file ${ACC_WITH_JACC_1} --col_file ${ACC_WITH_JACC_1} --write_to_file /scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/pw_of_acc_with_other_jaccard_1.npz --show_all
