#!/bin/bash
#python create_dmi_database.py \
#	--source /scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db \
#	--samples /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.015625/analysis/filtered_analysis_min_mbases_1620_min_diversity_10/filtered_data.parquet \
#	--reference /scratch/dmk333_new/Logan/Logan_Analyses/dark_matter_index/data_orig/reference_hashes_k31.bin \
#	--output /scratch/dmk333_new/Logan/Logan_Analyses/dark_matter_index/data/samples_with_reference_hashes_cov_0.015625_min_mbases_1620_min_diversity_10.db \
#	--ksize 31

# crashed at the indexing step, trying again
python create_dmi_database.py \
        --source /scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db \
        --samples /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.015625/analysis/filtered_analysis_min_mbases_1620_min_diversity_10/filtered_data.parquet \
        --reference /scratch/dmk333_new/Logan/Logan_Analyses/dark_matter_index/data_orig/reference_hashes_k31.bin \
        --output /scratch/dmk333_new/Logan/Logan_Analyses/dark_matter_index/data/samples_with_reference_hashes_cov_0.015625_min_mbases_1620_min_diversity_10.db \
        --ksize 31 \
	--skip-reference \
	--skip-samples
