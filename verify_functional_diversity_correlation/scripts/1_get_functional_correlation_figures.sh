#!/bin/bash
# Generate filtered correlation figures for functional diversity analysis

# Run downstream analysis with various filter settings
# These mirror the filters used for taxonomic analysis

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="/scratch/dmk333_new/Logan/Logan_Analyses/verify_functional_diversity_correlation/data"

# Standard filter: min 1620 mbases, min 10 KOs (functional diversity)
echo "Running filtered analysis with min_mbases=1620, min_diversity=10..."
python3 "${SCRIPT_DIR}/functional_analyze_parquet_data.py" \
    --input "${BASE_DIR}/functional_hash_diversity_results/data/functional_hash_diversity_data.parquet" \
    --output "${BASE_DIR}/functional_hash_diversity_results/analysis/filtered_analysis_min_mbases_1620_min_diversity_10" \
    --min-mbases 1620 \
    --min-diversity 10 \
    --join-metadata

echo ""
echo "Running filtered analysis with min_hashes=10000, min_diversity=10..."
python3 "${SCRIPT_DIR}/functional_analyze_parquet_data.py" \
    --input "${BASE_DIR}/functional_hash_diversity_results/data/functional_hash_diversity_data.parquet" \
    --output "${BASE_DIR}/functional_hash_diversity_results/analysis/filtered_analysis_min_hash_10000_min_diversity_10" \
    --min-hashes 10000 \
    --min-diversity 10 \
    --join-metadata

echo ""
echo "Running filtered analysis with min_hashes=1000, min_diversity=100..."
python3 "${SCRIPT_DIR}/functional_analyze_parquet_data.py" \
    --input "${BASE_DIR}/functional_hash_diversity_results/data/functional_hash_diversity_data.parquet" \
    --output "${BASE_DIR}/functional_hash_diversity_results/analysis/filtered_analysis_min_hash_1000_min_diversity_100" \
    --min-hashes 1000 \
    --min-diversity 100 \
    --join-metadata

echo ""
echo "All analyses complete!"
echo "Results saved to: ${BASE_DIR}/functional_hash_diversity_results/analysis/"
