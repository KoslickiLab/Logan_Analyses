#!/bin/bash

# Path to the DuckDB database
DB_PATH="/scratch/shared_data_new/Logan_yacht_data/metadata/aws_sra_metadata/metagenome_metadata.duckdb"

# Output file
OUTPUT_FILE="/scratch/dmk333_new/Logan/Logan_Analyses/Clustering_high_similarity/data/accessions_mbases_geq_10.txt"

# Query DuckDB and save results
duckdb -readonly "$DB_PATH" -noheader -list "SELECT acc FROM metadata WHERE mbases >= 10;" > "${OUTPUT_FILE}"

echo "Query completed. Results saved to $OUTPUT_FILE"
echo "Total accessions found: $(wc -l < "$OUTPUT_FILE")"
