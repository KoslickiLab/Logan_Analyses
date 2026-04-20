#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
FMH_BIN="$REPO_ROOT/kmer-sketch/bin/fmh_collisions"
ACCESSIONS_CSV="$REPO_ROOT/data/sorted_WGS_accessions.csv"

N=0   # 0 = process all
T=4   # parallel jobs
OUT_DIR="$REPO_ROOT/results"
TMP_DIR="/tmp"

usage() {
    cat >&2 <<EOF
Usage: $0 [-N num_accessions] [-t threads] [-o output_dir] [-d tmp_dir]
  -N  Number of accessions to process (default: all)
  -t  Parallel jobs (default: 4)
  -o  Output directory (default: <repo>/results)
  -d  Directory for temporary files (default: /tmp)
EOF
    exit 1
}

while getopts ":N:t:o:d:h" opt; do
    case $opt in
        N) N="$OPTARG" ;;
        t) T="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        d) TMP_DIR="$OPTARG" ;;
        h) usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage ;;
    esac
done

mkdir -p "$TMP_DIR"

mkdir -p "$OUT_DIR"

# ---------------------------------------------------------------------------
# Per-accession worker — called by GNU parallel.
# Downloads the unitig FASTA, runs fmh_collisions, emits one TSV row to stdout.
# Errors are reported to stderr and the job returns 0 so parallel keeps going.
# ---------------------------------------------------------------------------
run_one() {
    local acc="$1" out_dir="$2" fmh_bin="$3" tmp_dir="$4"

    local tmpfa sum_tmp col_out
    tmpfa=$(mktemp "${tmp_dir}/${acc}.XXXXXX.fa")
    sum_tmp=$(mktemp "${tmp_dir}/${acc}_sum.XXXXXX.tsv")
    col_out="$out_dir/${acc}.collisions.tsv"
    trap "rm -f '$tmpfa' '$sum_tmp'" EXIT

    local url="https://s3.amazonaws.com/logan-pub/u/${acc}/${acc}.unitigs.fa.zst"

    if ! wget -q -O - "$url" | zstd -d -q -f - -o "$tmpfa"; then
        echo "ERROR: download/decompress failed for $acc" >&2
        return 0
    fi

    if ! "$fmh_bin" \
            --input "$tmpfa" --kmer 31 --scale 0.001 --seed 42 \
            --summary "$sum_tmp" --collisions "$col_out" 2>/dev/null; then
        echo "ERROR: fmh_collisions failed for $acc" >&2
        return 0
    fi

    # Parse the summary file and emit one tab-separated row.
    awk -v acc="$acc" '
        /^# total_kmer_positions_scanned=/          { split($0,a,"="); pos=a[2] }
        /^# total_distinct_canonical_kmers_selected=/{ split($0,a,"="); dist=a[2] }
        /^sketch_size\t/            { split($0,a,"\t"); sz=a[2] }
        /^n_collision_hashes\t/     { split($0,a,"\t"); nh=a[2] }
        /^n_kmers_in_collisions\t/  { split($0,a,"\t"); nk=a[2] }
        END { printf "%s\t%s\t%s\t%s\t%s\t%s\n", acc, sz, nh, nk, pos, dist }
    ' "$sum_tmp"
}

export -f run_one

# ---------------------------------------------------------------------------
# Build accession list
# ---------------------------------------------------------------------------
# The CSV is sorted largest→smallest, so -N takes from the bottom (smallest/fastest).
if (( N > 0 )); then
    mapfile -t ACCS < <(tail -n +2 "$ACCESSIONS_CSV" | tail -n "$N")
else
    mapfile -t ACCS < <(tail -n +2 "$ACCESSIONS_CSV")
fi

echo "Processing ${#ACCS[@]} accessions with $T parallel jobs → $OUT_DIR (tmp: $TMP_DIR)" >&2

# ---------------------------------------------------------------------------
# Run in parallel; capture all TSV rows (with header) into summary.tsv
# ---------------------------------------------------------------------------
{
    printf 'acc\tsketch_size\tn_collision_hashes\tn_kmers_in_collisions\ttotal_kmer_positions\ttotal_distinct_kmers_selected\n'
    printf '%s\n' "${ACCS[@]}" \
        | parallel -j "$T" --line-buffer run_one {} "$OUT_DIR" "$FMH_BIN" "$TMP_DIR"
} > "$OUT_DIR/summary.tsv"

echo "Done. Summary: $OUT_DIR/summary.tsv" >&2
