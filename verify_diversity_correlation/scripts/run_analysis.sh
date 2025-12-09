#!/bin/bash
# Wrapper script for hash-diversity correlation analysis with common presets

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default values
N_JOBS=128
OUTPUT_BASE="hash_diversity_results"

show_help() {
    cat << EOF
Hash-Diversity Correlation Analysis Wrapper
===========================================

This wrapper provides convenient presets for common analysis scenarios.

USAGE:
    $0 [PRESET] [OPTIONS]

PRESETS:
    quick       - Quick test with 1,000 samples (for testing pipeline)
    small       - Small analysis with 10,000 samples  
    medium      - Medium analysis with 50,000 samples (recommended for initial analysis)
    large       - Large analysis with 100,000 samples
    full        - Full analysis with all available samples (may take hours)
    sensitivity - Sensitivity analysis across coverage thresholds (10,000 samples)

OPTIONS:
    -j, --jobs N        Number of parallel jobs (default: 128)
    -o, --output DIR    Output directory (default: hash_diversity_results_PRESET)
    -c, --coverage VAL  Coverage threshold (default: 0.0625)
    -h, --help          Show this help message

EXAMPLES:
    # Quick test run
    $0 quick

    # Medium analysis with custom job count
    $0 medium --jobs 256

    # Large analysis with custom coverage
    $0 large --coverage 0.125 --output my_results

    # Sensitivity analysis
    $0 sensitivity

    # Custom run (not using presets)
    python3 hash_diversity_correlation.py \\
        --output-dir custom_results \\
        --n-samples 25000 \\
        --coverage 0.0625 \\
        --n-jobs 128

EOF
}

# Parse arguments
PRESET=""
COVERAGE=0.0625
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        quick|small|medium|large|full|sensitivity)
            PRESET=$1
            shift
            ;;
        -j|--jobs)
            N_JOBS=$2
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR=$2
            shift 2
            ;;
        -c|--coverage)
            COVERAGE=$2
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Check if preset was provided
if [ -z "$PRESET" ]; then
    echo "Error: No preset specified"
    show_help
    exit 1
fi

# Set default output directory if not specified
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="${OUTPUT_BASE}_${PRESET}_$(date +%Y%m%d_%H%M%S)"
fi

# Execute based on preset
echo "=============================================="
echo "Hash-Diversity Correlation Analysis"
echo "=============================================="
echo "Preset:         $PRESET"
echo "Output:         $OUTPUT_DIR"
echo "Parallel jobs:  $N_JOBS"
echo "Coverage:       $COVERAGE"
echo "=============================================="
echo ""

case $PRESET in
    quick)
        echo "Running QUICK test (1,000 samples)..."
        python3 "$SCRIPT_DIR/hash_diversity_correlation.py" \
            --output-dir "$OUTPUT_DIR" \
            --n-samples 1000 \
            --coverage "$COVERAGE" \
            --n-jobs "$N_JOBS"
        ;;
    
    small)
        echo "Running SMALL analysis (10,000 samples)..."
        python3 "$SCRIPT_DIR/hash_diversity_correlation.py" \
            --output-dir "$OUTPUT_DIR" \
            --n-samples 10000 \
            --coverage "$COVERAGE" \
            --n-jobs "$N_JOBS"
        ;;
    
    medium)
        echo "Running MEDIUM analysis (50,000 samples)..."
        python3 "$SCRIPT_DIR/hash_diversity_correlation.py" \
            --output-dir "$OUTPUT_DIR" \
            --n-samples 50000 \
            --coverage "$COVERAGE" \
            --n-jobs "$N_JOBS"
        ;;
    
    large)
        echo "Running LARGE analysis (100,000 samples)..."
        python3 "$SCRIPT_DIR/hash_diversity_correlation.py" \
            --output-dir "$OUTPUT_DIR" \
            --n-samples 100000 \
            --coverage "$COVERAGE" \
            --n-jobs "$N_JOBS"
        ;;
    
    full)
        echo "Running FULL analysis (all available samples)..."
        echo "Warning: This may take several hours!"
        read -p "Continue? (y/n) " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            python3 "$SCRIPT_DIR/hash_diversity_correlation.py" \
                --output-dir "$OUTPUT_DIR" \
                --coverage "$COVERAGE" \
                --n-jobs "$N_JOBS"
        else
            echo "Cancelled."
            exit 0
        fi
        ;;
    
    sensitivity)
        echo "Running SENSITIVITY analysis (10,000 samples, multiple coverages)..."
        python3 "$SCRIPT_DIR/hash_diversity_sensitivity.py" \
            --output-dir "$OUTPUT_DIR" \
            --n-samples 10000 \
            --n-jobs "$N_JOBS"
        ;;
esac

echo ""
echo "=============================================="
echo "Analysis complete!"
echo "Results saved to: $OUTPUT_DIR"
echo "=============================================="
echo ""
echo "Key outputs:"
echo "  Main plot:    $OUTPUT_DIR/plots/hash_diversity_correlation.png"
echo "  Report:       $OUTPUT_DIR/reports/analysis_report.txt"
echo "  Data:         $OUTPUT_DIR/data/hash_diversity_data.csv"
echo ""
