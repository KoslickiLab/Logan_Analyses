#!/bin/bash
python compute_dmi_native.py \
	--database /scratch/dmk333_new/Logan/Logan_Analyses/dark_matter_index/data/samples_with_reference_hashes_cov_0.0_min_mbases_100.db \
	--input /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.0/data/hash_diversity_data.parquet \
       	--output /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.0/data/hash_diversity_data_with_dmi.parquet \
	--chunked \
	--chunk-size 10000 \
	--threads 200 \
	--memory 3000GB
