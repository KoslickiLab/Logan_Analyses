#!/bin/bash

# Calculate statistics on the number of rows (matches) in each CSV file
# Uses parallel processing for speed

INPUT_DIR="/scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/output"
OUTPUT_FILE="/scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/jaccard_1_match_statistics.txt"
TEMP_DIR="/scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/scripts/temp"
NUM_CORES=64

# Create temporary directory for worker outputs
mkdir -p "$TEMP_DIR"
rm -f "$TEMP_DIR"/counts_*.txt

echo "Counting rows in parallel with $NUM_CORES cores..."

# Count rows in parallel (excluding header)
find "$INPUT_DIR" -maxdepth 1 -name "*_jaccard_1.csv" -type f | \
    parallel -j "$NUM_CORES" --line-buffer \
    "tail -n +2 {} | wc -l" > /scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/scripts/temp/all_counts.txt

echo "Calculating statistics..."

# Use Python to calculate statistics
python3 << 'EOF'
import sys
import statistics

# Read all counts
with open('/scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/scripts/temp/all_counts.txt', 'r') as f:
    counts = [int(line.strip()) for line in f if line.strip()]

if not counts:
    print("No data found!")
    sys.exit(1)

# Calculate statistics
counts_sorted = sorted(counts)
n = len(counts)

stats = {
    'Total files': n,
    'Min': min(counts),
    'Max': max(counts),
    'Mean': statistics.mean(counts),
    'Median': statistics.median(counts),
    'Std Dev': statistics.stdev(counts) if n > 1 else 0,
    'Q1 (25th percentile)': counts_sorted[n // 4],
    'Q3 (75th percentile)': counts_sorted[3 * n // 4],
}

# Print and save results
output = []
output.append("=" * 50)
output.append("Statistics: Number of matches per CSV file")
output.append("=" * 50)
for key, value in stats.items():
    if isinstance(value, float):
        line = f"{key:.<30} {value:>15,.2f}"
    else:
        line = f"{key:.<30} {value:>15,}"
    output.append(line)
output.append("=" * 50)

result = '\n'.join(output)
print(result)

# Save to file
with open('$OUTPUT_FILE', 'w') as f:
    f.write(result + '\n')

EOF

# Cleanup
rm -rf "$TEMP_DIR"

echo ""
echo "Results saved to: $OUTPUT_FILE"
