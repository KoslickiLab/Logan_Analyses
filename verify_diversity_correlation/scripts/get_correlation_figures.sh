#!/bin/bash

parallel --jobs 8 'python3 analyze_parquet_data.py --input /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_{}/data/hash_diversity_data.parquet --output /scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_{}/analysis/filtered_analysis_min_hash_1000_min_diversity_100 --min-hashes 1000 --min-diversity 100 --join-metadata' ::: 0.0 1.0 0.5 0.25 0.125 0.0625 0.03125 0.015625
