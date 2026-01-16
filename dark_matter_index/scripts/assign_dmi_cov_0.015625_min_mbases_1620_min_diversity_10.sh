#!/bin/bash
#python batch_compute_dmi.py \
#	--reference /scratch/dmk333_new/Logan/Logan_Analyses/dark_matter_index/data_orig/reference_hashes_k31.bin \
#	--input /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.015625/analysis/filtered_analysis_min_mbases_1620_min_diversity_10/filtered_data.parquet \
#	--output /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.015625/analysis/filtered_analysis_min_mbases_1620_min_diversity_10/filtered_data_with_dmi.parquet \
#	--database /scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db \
#	--ksize 31 \
#	--workers 1 \
#	--chunk-size 500 \
#	--checkpoint /scratch/dmk333_new/Logan/Logan_Analyses/dark_matter_index/data/checkpoint.parquet \
#	--checkpoint-interval 5

# Takes Waaaaay tooo looong
#python batch_compute_dmi_optimized.py \
#        --reference /scratch/dmk333_new/Logan/Logan_Analyses/dark_matter_index/data_orig/reference_hashes_k31.bin \
#        --input /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.015625/analysis/filtered_analysis_min_mbases_1620_min_diversity_10/filtered_data.parquet \
#        --output /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.015625/analysis/filtered_analysis_min_mbases_1620_min_diversity_10/filtered_data_with_dmi.parquet \
#        --database /scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db \
#        --ksize 31 \
#	--batch-size 500000
