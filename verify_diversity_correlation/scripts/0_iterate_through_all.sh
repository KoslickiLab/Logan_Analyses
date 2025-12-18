#!/bin/bash
# Run hash-diversity correlation analysis for all coverage thresholds
# Coverage values: 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 0

# Configuration
BASE_OUTPUT_DIR="/scratch/dmk333_new/Logan/Logan_Analyses/verify_diversity_correlation/data"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
N_JOBS=200  # Leave some cores free for system
LOG_FILE="${BASE_OUTPUT_DIR}/run_all_coverages.log"

# Create base output directory if needed
mkdir -p "$BASE_OUTPUT_DIR"

# Start logging
exec > >(tee -a "$LOG_FILE")
exec 2>&1

# Coverage values (1, 1/2, 1/2^2, ..., 1/2^6, 0)
COVERAGES=(1.0 0.5 0.25 0.125 0.0625 0.03125 0.015625 0.0)

# Track successes and failures
declare -a SUCCESSFUL=()
declare -a FAILED=()

echo "======================================================================="
echo "HASH-DIVERSITY CORRELATION ANALYSIS - ALL COVERAGE THRESHOLDS"
echo "======================================================================="
echo "Starting: $(date)"
echo "Coverage values: ${COVERAGES[@]}"
echo "Number of coverage thresholds: ${#COVERAGES[@]}"
echo "Parallel jobs: $N_JOBS"
echo "Using ALL available samples (no subsampling)"
echo ""
echo "ESTIMATED RUNTIME:"
echo "  Per coverage: ~4-8 hours (depends on sample count and database speed)"
echo "  Total: ~32-64 hours for all 8 coverage values"
echo ""
echo "Output directories will be created in:"
echo "  $BASE_OUTPUT_DIR/"
echo "======================================================================="
echo ""
echo "Starting analysis..."
echo ""

# Loop through each coverage value
for COV in "${COVERAGES[@]}"; do
    OUTPUT_DIR="${BASE_OUTPUT_DIR}/hash_diversity_results_full_cov_${COV}"
    
    echo ""
    echo "-----------------------------------------------------------------------"
    echo "Processing coverage = $COV"
    echo "Output directory: $OUTPUT_DIR"
    echo "Started: $(date)"
    echo "-----------------------------------------------------------------------"
    
    # Run the analysis
    if python3 "${SCRIPT_DIR}/hash_diversity_correlation.py" \
        --output-dir "$OUTPUT_DIR" \
        --coverage "$COV" \
        --n-jobs "$N_JOBS" \
        --min-mbases 100; then
        
        SUCCESSFUL+=("$COV")
        echo ""
        echo "✓ SUCCESS: Completed coverage = $COV at $(date)"
        echo "-----------------------------------------------------------------------"
    else
        FAILED+=("$COV")
        echo ""
        echo "✗ FAILED: Coverage = $COV failed at $(date)"
        echo "-----------------------------------------------------------------------"
        echo "Continuing with remaining coverage values..."
    fi
done

echo ""
echo "======================================================================="
echo "ALL ANALYSES COMPLETE"
echo "======================================================================="
echo "Finished: $(date)"
echo ""
echo "Summary:"
echo "  Successful: ${#SUCCESSFUL[@]}/${#COVERAGES[@]}"
echo "  Failed: ${#FAILED[@]}/${#COVERAGES[@]}"
echo ""

if [ ${#SUCCESSFUL[@]} -gt 0 ]; then
    echo "✓ Successful coverage values:"
    for COV in "${SUCCESSFUL[@]}"; do
        echo "    $COV → ${BASE_OUTPUT_DIR}/hash_diversity_results_full_cov_${COV}/"
    done
    echo ""
fi

if [ ${#FAILED[@]} -gt 0 ]; then
    echo "✗ Failed coverage values:"
    for COV in "${FAILED[@]}"; do
        echo "    $COV"
    done
    echo ""
fi

echo "Key outputs for each successful run:"
echo "  - plots/hash_diversity_correlation.png"
echo "  - reports/analysis_report.txt"
echo "  - data/hash_diversity_data.csv"
echo ""
echo "Full log saved to: $LOG_FILE"
echo "======================================================================="

# Exit with error if any failed
if [ ${#FAILED[@]} -gt 0 ]; then
    exit 1
fi
