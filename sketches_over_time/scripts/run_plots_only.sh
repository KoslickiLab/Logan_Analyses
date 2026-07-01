#!/usr/bin/env bash
# =============================================================================
# run_plots_only.sh — Re-run all three plotting scripts (no DB rebuild)
# =============================================================================
#
# Usage
# -----
#   # Auto-detect the most recent 'full' intermediate DB:
#   bash run_plots_only.sh
#
#   # Point at a specific DB:
#   bash run_plots_only.sh --db /scratch/dmk333/Logan_Analyses/sketches_over_time/data/intermediate_full.duckdb
#
#   # Smoke-test variant:
#   bash run_plots_only.sh --db data/intermediate_test500.duckdb
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/../logs"
mkdir -p "$LOG_DIR"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"

# Collect any extra args (e.g. --db /path/to/file) to forward to each plot script
EXTRA=("$@")

echo "================================================================"
echo " Logan — plots only"
echo " Start  : $(date)"
echo " Extra  : ${EXTRA[*]:-'(none)'}"
echo "================================================================"

run_plot() {
    local name="$1"; local script="$2"; shift 2
    local log="$LOG_DIR/${TIMESTAMP}_${name}.log"
    echo ""
    echo ">>> $name"
    python "$SCRIPT_DIR/$script" "${EXTRA[@]}" "$@" 2>&1 | tee "$log"
    local rc=${PIPESTATUS[0]}
    [ $rc -eq 0 ] || { echo "!!! $name FAILED (exit $rc)"; exit $rc; }
    echo "<<< $name done"
}

run_plot "02_plot_cumulative" "02_plot_cumulative_hashes.py"
run_plot "03_plot_new_hashes" "03_plot_new_hashes_per_sample.py"
run_plot "04_plot_survival"   "04_plot_survival.py"

echo ""
echo "================================================================"
echo " Plots complete  $(date)"
echo " Output → $SCRIPT_DIR/../data/plots/"
echo "================================================================"
