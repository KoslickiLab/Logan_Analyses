#!/bin/bash
# Extract human genome hashes from sourmash signatures to DuckDB
#
# Prerequisites:
#   - sourmash >= 4.8.0
#   - Human genome sketches file (human_genome_sketches_all.sig.zip)
#
# The output database will be used by compute_dmi_host_filtered.py

set -e

# Configuration
INPUT_SIG="/scratch/dmk333_new/Logan/Logan_Analyses/dark_matter_index/scripts/remove_host_contamination/data/human_genome_sketches_all.sig.zip"
OUTPUT_DB="/scratch/dmk333_new/Logan/Logan_Analyses/dark_matter_index/data/human_hashes_k31.db"
KSIZE=31
SCALED=1000
THREADS=64
MEMORY="256GB"

echo "============================================================"
echo "Extracting human genome hashes"
echo "============================================================"
echo "Input:  ${INPUT_SIG}"
echo "Output: ${OUTPUT_DB}"
echo "ksize:  ${KSIZE}"
echo "============================================================"

python extract_human_hashes.py \
    --input "${INPUT_SIG}" \
    --output "${OUTPUT_DB}" \
    --ksize ${KSIZE} \
    --scaled ${SCALED} \
    --threads ${THREADS} \
    --memory "${MEMORY}"

echo "============================================================"
echo "Done! Human hash database created at: ${OUTPUT_DB}"
echo "============================================================"
