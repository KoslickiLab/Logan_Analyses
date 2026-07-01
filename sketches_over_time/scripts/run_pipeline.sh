#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Logan FracMinHash Timeline Analysis
# =============================================================================
#
# Runs all four stages in order:
#   01  Build intermediate DuckDB  (heavy — hours on full dataset)
#   02  Plot 1: cumulative distinct hashes
#   03  Plot 2: new hashes per sample
#   04  Plot 3: survival curves + birth/death heatmap
#
# All arguments after the script name are forwarded to 01_build_intermediate_db.py.
#
# Usage
# -----
#   # Smoke test — 500 samples (~minutes):
#   bash run_pipeline.sh --test --num-test-samples 500
#
#   # Full run  (~hours, use nohup / tmux):
#   bash run_pipeline.sh
#
#   # ILLUMINA only:
#   bash run_pipeline.sh --metadata-filter "platform='ILLUMINA'"
#
#   # Resume after interruption:
#   bash run_pipeline.sh --resume
#
#   # Only build stage 1, then stop:
#   bash run_pipeline.sh --stage 1
#
# Environment
# -----------
#   Assumes the 'logan' conda environment is active.
#   If not, activate it first:  conda activate logan
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/../logs"
mkdir -p "$LOG_DIR"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
BUILD_ARGS=("$@")

# Extract --output-db from BUILD_ARGS so the plot scripts read the same DB
# that step 01 wrote.  When --output-db is absent the plot scripts auto-detect.
PLOT_ARGS=()
i=0
while [ $i -lt ${#BUILD_ARGS[@]} ]; do
    if [ "${BUILD_ARGS[$i]}" = "--output-db" ]; then
        i=$((i+1))
        PLOT_ARGS=(--db "${BUILD_ARGS[$i]}")
    fi
    i=$((i+1))
done

echo "================================================================"
echo " Logan FracMinHash Timeline Analysis"
echo " Start  : $(date)"
echo " Args   : ${BUILD_ARGS[*]}"
echo " Plot DB: ${PLOT_ARGS[*]:-(auto-detect)}"
echo " Logs   : $LOG_DIR"
echo "================================================================"
echo ""

run_step() {
    local name="$1"
    local script="$2"
    shift 2
    local log_file="$LOG_DIR/${TIMESTAMP}_${name}.log"
    echo ">>> STEP $name"
    echo "    Script : $script"
    echo "    Log    : $log_file"
    python "$SCRIPT_DIR/$script" "$@" 2>&1 | tee "$log_file"
    local rc=${PIPESTATUS[0]}
    if [ $rc -ne 0 ]; then
        echo "!!! STEP $name FAILED  (exit $rc)" >&2
        exit $rc
    fi
    echo "<<< STEP $name done"
    echo ""
}

# ── Step 1: build the intermediate database ────────────────────────────────────
run_step "01_build" "01_build_intermediate_db.py" "${BUILD_ARGS[@]}"

# ── Steps 2-4: plots  (pass --db only if a specific DB was chosen) ─────────────
# The plotting scripts auto-detect the most recent 'full' DB in data/.
# If you want to point at a specific file, pass --db <path> explicitly:
#   python 02_plot_cumulative_hashes.py --db data/intermediate_full.duckdb

run_step "02_plot_cumulative"    "02_plot_cumulative_hashes.py"       "${PLOT_ARGS[@]}"
run_step "03_plot_new_hashes"    "03_plot_new_hashes_per_sample.py"   "${PLOT_ARGS[@]}"
run_step "04_plot_survival"      "04_plot_survival.py"                 "${PLOT_ARGS[@]}"

echo "================================================================"
echo " Pipeline complete  $(date)"
echo " Plots  → $SCRIPT_DIR/../data/plots/"
echo " CSVs   → $SCRIPT_DIR/../data/csv/"
echo " Logs   → $LOG_DIR"
echo "================================================================"
