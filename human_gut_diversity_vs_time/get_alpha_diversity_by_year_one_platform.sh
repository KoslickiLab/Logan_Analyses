#!/bin/bash
for year in {2012..2023}
do
  echo "Processing year $year..."
  start_time=$(date +%s)
  duckdb -readonly /scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db <<EOF
SET enable_progress_bar=true;

-- Attach the metadata database (read-only)
ATTACH '/scratch/shared_data_new/Logan_yacht_data/metadata/aws_sra_metadata/metagenome_metadata.duckdb'
    AS meta (READ_ONLY);

COPY (
    SELECT
        p.sample_id as accession,
        COUNT(DISTINCT CASE WHEN p.min_coverage = 1 THEN p.tax_id END)        as taxa_coverage_1,
        COUNT(DISTINCT CASE WHEN p.min_coverage = 0.5 THEN p.tax_id END)      as taxa_coverage_0_5,
        COUNT(DISTINCT CASE WHEN p.min_coverage = 0.25 THEN p.tax_id END)     as taxa_coverage_0_25,
        COUNT(DISTINCT CASE WHEN p.min_coverage = 0.125 THEN p.tax_id END)    as taxa_coverage_0_125,
        COUNT(DISTINCT CASE WHEN p.min_coverage = 0.0625 THEN p.tax_id END)   as taxa_coverage_0_0625,
        COUNT(DISTINCT CASE WHEN p.min_coverage = 0.03125 THEN p.tax_id END)  as taxa_coverage_0_03125,
        COUNT(DISTINCT CASE WHEN p.min_coverage = 0.015625 THEN p.tax_id END) as taxa_coverage_0_015625,
        COUNT(DISTINCT CASE WHEN p.min_coverage = 0 THEN p.tax_id END)        as taxa_coverage_0
    FROM taxa_profiles.profiles p
    INNER JOIN (
        -- Year-specific sample IDs, filtered to Illumina MiSeq in metadata
        SELECT s.sample_id
        FROM read_csv(
            '/scratch/shared_data_new/Logan_yacht_data/metadata/slices/human_gut_metagenomes_partitioned_by_year/human_gut_metagenomes_${year}.csv',
            header=false,
            columns={'sample_id': 'VARCHAR'}
        ) s
        JOIN meta.metadata m
          ON s.sample_id = m.acc      -- adjust this column name if needed
        WHERE m.instrument = 'Illumina MiSeq'
           AND m.libraryselection = 'RANDOM'
    ) s ON p.sample_id = s.sample_id
    WHERE p.min_coverage IN (1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0)
    GROUP BY p.sample_id
) TO '/scratch/dmk333_new/Logan/Logan_Analyses/human_gut_diversity_vs_time/results/Illumina-MiSeq/alpha_diversity_${year}_Illumina_MiSeq.csv' (HEADER);
EOF
  end_time=$(date +%s)
  elapsed=$((end_time - start_time))
  echo "Completed year $year in ${elapsed} seconds"
  echo "---"
done

