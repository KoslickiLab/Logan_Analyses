#!/bin/bash
python compute_dmi_native.py \
	--database ../data/samples_with_reference_hashes_cov_0.015625_min_mbases_1620_min_diversity_10.db \
	--input /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.015625/analysis/filtered_analysis_min_mbases_1620_min_diversity_10/filtered_data.parquet \
       	--output /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.015625/analysis/filtered_analysis_min_mbases_1620_min_diversity_10/filtered_data_with_dmi.parquet \
	--chunked \
	--chunk-size 10000 \
	--threads 200 \
	--memory 3000GB
