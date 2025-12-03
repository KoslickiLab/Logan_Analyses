#!/bin/bash

# Extract all unique accessions from CSV files and their filenames
# Input: CSV files in format {ACCESSION}_jaccard_1.csv with ID,Jaccard columns
# Output: One unique accession per line

INPUT_DIR="/scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/output"
OUTPUT_FILE="/scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/accessions_mbases_geq_10_with_other_Jaccard_1.txt"

# Remove output file if it exists
rm -f "$OUTPUT_FILE"

# Process all *_jaccard_1.csv files using find to avoid ARG_MAX issues
find "$INPUT_DIR" -maxdepth 1 -name "*_jaccard_1.csv" -type f | while read -r csv_file; do
    # Extract accession from filename (remove path and _jaccard_1.csv suffix)
    filename=$(basename "$csv_file")
    accession="${filename%_jaccard_1.csv}"
    echo "$accession" >> "$OUTPUT_FILE"
    
    # Extract accessions from CSV body (skip header, get first column)
    tail -n +2 "$csv_file" | cut -d',' -f1 >> "$OUTPUT_FILE"
done

# Sort and remove duplicates
sort -u "$OUTPUT_FILE" > "${OUTPUT_FILE}.copy"
mv "${OUTPUT_FILE}.copy" "$OUTPUT_FILE"

echo "Done! Unique accessions written to: $OUTPUT_FILE"
echo "Total unique accessions: $(wc -l < "$OUTPUT_FILE")"
