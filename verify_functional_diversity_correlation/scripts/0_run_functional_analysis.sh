#!/bin/bash
# Run functional hash-diversity correlation analysis
# Uses FracMinHash amino acid sketches (k=11, scale=1000) from fmh-funprofiler

# Configuration
BASE_OUTPUT_DIR="/scratch/dmk333_new/Logan/Logan_Analyses/verify_functional_diversity_correlation/data"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
N_JOBS=200  # Leave some cores free for system
LOG_FILE="${BASE_OUTPUT_DIR}/run_functional_analysis.log"

# Create base output directory if needed
mkdir -p "$BASE_OUTPUT_DIR"

# Start logging
exec > >(tee -a "$LOG_FILE")
exec 2>&1

echo "======================================================================="
echo "FUNCTIONAL HASH-DIVERSITY CORRELATION ANALYSIS"
echo "======================================================================="
echo "Starting: $(date)"
echo "Parallel jobs: $N_JOBS"
echo "Using ALL available samples (no subsampling)"
echo ""
echo "ANALYSIS OVERVIEW:"
echo "  - Data source: fmh-funprofiler functional profiles"
echo "  - FracMinHash: k=11 (amino acid), scale=1000"
echo "  - Diversity metrics: Richness, Shannon, Simpson, Hill2, Berger-Parker, Pielou"
echo ""
echo "Output directory will be created in:"
echo "  $BASE_OUTPUT_DIR/"
echo "======================================================================="
echo ""
echo "Starting analysis..."
echo ""

OUTPUT_DIR="${BASE_OUTPUT_DIR}/functional_hash_diversity_results"

echo "-----------------------------------------------------------------------"
echo "Processing functional diversity analysis"
echo "Output directory: $OUTPUT_DIR"
echo "Started: $(date)"
echo "-----------------------------------------------------------------------"

# Run the analysis
if python3 "${SCRIPT_DIR}/functional_hash_diversity_correlation.py" \
    --output-dir "$OUTPUT_DIR" \
    --n-jobs "$N_JOBS" \
    --min-mbases 100; then
    
    echo ""
    echo "✓ SUCCESS: Completed functional analysis at $(date)"
    echo "-----------------------------------------------------------------------"
else
    echo ""
    echo "✗ FAILED: Functional analysis failed at $(date)"
    echo "-----------------------------------------------------------------------"
    exit 1
fi

echo ""
echo "======================================================================="
echo "ANALYSIS COMPLETE"
echo "======================================================================="
echo "Finished: $(date)"
echo ""
echo "Key outputs:"
echo "  - plots/summary_all_metrics.png"
echo "  - plots/hash_vs_*_correlation.png (one per metric)"
echo "  - reports/analysis_report.txt"
echo "  - data/functional_hash_diversity_data.parquet"
echo ""
echo "Full log saved to: $LOG_FILE"
echo "======================================================================="
