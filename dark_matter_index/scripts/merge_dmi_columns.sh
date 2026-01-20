#!/bin/bash

if false; then
python merge_dmi_columns.py \
	--parquet-with-dmi /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.0/data/hash_diversity_data_with_dmi.parquet \
    --parquet-without-dmi /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.0/analysis/filtered_analysis_min_hash_10000_min_diversity_10/filtered_data.parquet \
    --output /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.0/analysis/filtered_analysis_min_hash_10000_min_diversity_10/filtered_data_with_dmi_merged.parquet
fi
