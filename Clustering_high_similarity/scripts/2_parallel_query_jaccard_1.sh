#!/bin/bash

set -euo pipefail

# Configuration
INPUT_FILE="/scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/accessions_mbases_geq_10_with_other_Jaccard_1.txt"
MATRIX_DIR="/scratch/mgs_project/matrix_unzipped/"
DB_DIR="/scratch/mgs_project/db/"
QUERY_PC_MAT="/scratch/dmk333_new/Logan/Logan_Analyses/metagenome_vector_sketches/build/query_pc_mat"
OUTPUT_SUFFIX="jaccard_1.csv"
NUM_CHUNKS=128
NUM_JOBS=128

# Hack to get around bug in query_pc_mat not respecting output paths
cd /scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/output

# Create temporary directory for chunks
TEMP_DIR="./chunks_tmp_$$"
mkdir -p "$TEMP_DIR"

# Log file
LOGFILE="./parallel_processing_$$.log"

echo "Starting parallel processing at $(date)"
echo "Input file: $INPUT_FILE"
echo "Number of chunks: $NUM_CHUNKS"
echo "Parallel jobs: $NUM_JOBS"
echo "Temporary directory: $TEMP_DIR"
echo "Log file: $LOGFILE"
echo ""

# Cleanup function
cleanup() {
    local exit_code=$?
    echo ""
    echo "Cleaning up temporary files..."
    rm -rf "$TEMP_DIR"
    if [ $exit_code -eq 0 ]; then
        echo "Processing completed successfully at $(date)"
    else
        echo "Processing failed with exit code $exit_code at $(date)"
        echo "Check $LOGFILE for details"
    fi
}

# Set trap to cleanup on exit
trap cleanup EXIT

# Split the input file into chunks
echo "Splitting $INPUT_FILE into $NUM_CHUNKS chunks..."
total_lines=$(wc -l < "$INPUT_FILE")
lines_per_chunk=$(( (total_lines + NUM_CHUNKS - 1) / NUM_CHUNKS ))
split -l "$lines_per_chunk" -d -a 4 "$INPUT_FILE" "$TEMP_DIR/chunk_"
echo "Created $(ls -1 "$TEMP_DIR"/chunk_* | wc -l) chunk files"
echo ""

# Function to process a single chunk
process_chunk() {
    local chunk_file=$1
    "$QUERY_PC_MAT" \
        --matrix "$MATRIX_DIR" \
        --db "$DB_DIR" \
        --query_file "$chunk_file" \
        --write_to_file "$OUTPUT_SUFFIX" \
        --no_self \
        --min_jaccard 1 \
	--show_all
}

export -f process_chunk
export QUERY_PC_MAT MATRIX_DIR DB_DIR OUTPUT_SUFFIX

# Run parallel processing
echo "Starting parallel processing of chunks..."
parallel \
    --jobs "$NUM_JOBS" \
    --halt now,fail=1 \
    --joblog "$LOGFILE" \
    --progress \
    process_chunk ::: "$TEMP_DIR"/chunk_*

echo ""
echo "All chunks processed successfully!"
