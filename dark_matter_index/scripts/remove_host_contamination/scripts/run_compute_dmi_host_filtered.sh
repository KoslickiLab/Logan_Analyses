#!/bin/bash
# Compute DMI with host contamination removal for human-associated metagenomes
#
# This script:
# 1. Filters samples to only human-associated metagenome types
# 2. Excludes human genome hashes from the DMI calculation
# 3. Outputs a parquet file compatible with existing analysis scripts
#
# Prerequisites:
#   1. Run extract_human_hashes.py to create human_hashes_k31.db
#   2. Have the main DMI database (from create_dmi_database.py)

set -e

# Configuration - Update these paths as needed
DMI_DATABASE="/scratch/dmk333_new/Logan/Logan_Analyses/dark_matter_index/data/samples_with_reference_hashes_cov_0.015625_min_mbases_1620_min_diversity_10.db"
HUMAN_HASHES="/scratch/dmk333_new/Logan/Logan_Analyses/dark_matter_index/data/human_hashes_k31.db"
INPUT_PARQUET="/scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.015625/analysis/filtered_analysis_min_mbases_1620_min_diversity_10/filtered_data.parquet"
OUTPUT_PARQUET="/scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.015625/analysis/filtered_analysis_min_mbases_1620_min_diversity_10/filtered_data_with_dmi_host_removed.parquet"

# Processing parameters
THREADS=200
MEMORY="3000GB"
CHUNK_SIZE=5000

# Organism filter - add or remove organisms as needed
# Common human-associated metagenome types:
#   - "human gut metagenome"
#   - "human feces metagenome" 
#   - "human skin metagenome"
#   - "human oral metagenome"
#   - "human saliva metagenome"
#   - "human metagenome"
#   - "human nasopharyngeal metagenome"
#   - "human vaginal metagenome"

ORGANISMS=(
    "human gut metagenome"
    "human metagenome"
    "human oral metagenome"
    "human skin metagenome"
    "human vaginal metagenome"
    "human feces metagenome"
    "human lung metagenome"
    "human nasopharyngeal metagenome"
    "human saliva metagenome"
    "human reproductive system metagenome"
    "human blood metagenome"
    "human urinary tract metagenome"
    "human milk metagenome"
    "human sputum metagenome"
    "human eye metagenome"
    "human viral metagenome"
    "human bile metagenome"
    "human tracheal metagenome"
    "uncultured human oral bacterium"
    "human brain metagenome"
    "Homo sapiens"
)

echo "============================================================"
echo "Computing DMI with Host Contamination Removal"
echo "============================================================"
echo "DMI Database:   ${DMI_DATABASE}"
echo "Human Hashes:   ${HUMAN_HASHES}"
echo "Input:          ${INPUT_PARQUET}"
echo "Output:         ${OUTPUT_PARQUET}"
echo "Organisms:      ${ORGANISMS[*]}"
echo "============================================================"

# Build organism arguments
ORG_ARGS=""
for org in "${ORGANISMS[@]}"; do
    ORG_ARGS="${ORG_ARGS} \"${org}\""
done

# Run the DMI computation
python compute_dmi_host_filtered.py \
    --database "${DMI_DATABASE}" \
    --host-hashes "${HUMAN_HASHES}" \
    --input "${INPUT_PARQUET}" \
    --output "${OUTPUT_PARQUET}" \
    --organisms ${ORGANISMS[@]} \
    --chunked \
    --chunk-size ${CHUNK_SIZE} \
    --threads ${THREADS} \
    --memory "${MEMORY}"

echo "============================================================"
echo "Done! Output written to: ${OUTPUT_PARQUET}"
echo "============================================================"
