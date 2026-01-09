# Functional Hash-Diversity Correlation Analysis

This set of scripts analyzes the relationship between FracMinHash hash diversity and functional diversity (KEGG KO-based) in metagenomic samples. This is analogous to the existing taxonomic hash-diversity analysis but uses functional profiler data.

## Overview

The analysis uses functional profile data from fmh-funprofiler stored in the YACHT DuckDB:
- **Hash counts**: `functional_profile_data.gather_data` table
- **KO profiles**: `functional_profile.profiles` table
- **FracMinHash parameters**: k=11 (amino acid), scale=1000

## Diversity Metrics Computed

Unlike the taxonomic analysis which only computed species richness, the functional analysis includes multiple abundance-aware diversity metrics:

| Metric | Description | Formula |
|--------|-------------|---------|
| **Observed Richness** | Count of distinct KOs | S |
| **Shannon Index (H')** | Entropy-based diversity | -Σ(pᵢ × ln(pᵢ)) |
| **Simpson's Index (D)** | Probability two random individuals are same species | Σ(pᵢ²) |
| **Gini-Simpson** | Complement of Simpson | 1 - D |
| **Hill Number Order 2** | Effective number of species | 1/D |
| **Berger-Parker** | Dominance index (abundance of most common) | max(pᵢ) |
| **Pielou's Evenness (J')** | How evenly distributed abundances are | H'/ln(S) |

Where pᵢ is the proportional abundance of KO i.

## Scripts

### 0_run_functional_analysis.sh
Main driver script that runs the full analysis pipeline.

```bash
./0_run_functional_analysis.sh
```

Output: `/scratch/dmk333_new/Logan/Logan_Analyses/verify_functional_correlation/data/functional_hash_diversity_results/`

### functional_hash_diversity_correlation.py
Core analysis script. Extracts data, computes diversity metrics, performs correlation analysis, generates plots.

```bash
python3 functional_hash_diversity_correlation.py \
    --output-dir results \
    --n-jobs 200 \
    --min-mbases 100
```

Options:
- `--output-dir, -o`: Output directory (required)
- `--n-samples, -n`: Limit samples (default: all)
- `--n-jobs, -j`: Parallel workers (default: 64)
- `--min-mbases, -m`: Minimum sequencing depth (default: 100)

### 1_get_functional_correlation_figures.sh
Runs filtered downstream analysis with multiple filter combinations.

```bash
./1_get_functional_correlation_figures.sh
```

### functional_analyze_parquet_data.py
Downstream analysis with metadata joining, filtering, and categorical plots.

```bash
python3 functional_analyze_parquet_data.py \
    --input results/data/functional_hash_diversity_data.parquet \
    --output filtered_analysis \
    --min-mbases 1620 \
    --min-diversity 10 \
    --join-metadata
```

### 2_plot_functional_custom_correlation.py
Generate custom correlation plots with arbitrary filtering and categorical coloring.

```bash
# Plot Shannon index for soil samples, colored by platform
python3 2_plot_functional_custom_correlation.py \
    --input filtered_data.parquet \
    --output soil_shannon \
    --filter "organism=soil metagenome" \
    --diversity-metric shannon_index \
    --color-by platform

# Available diversity metrics:
#   observed_richness_per_mb, shannon_index, simpson_index,
#   gini_simpson, hill_2, berger_parker, pielou_evenness
```

## Output Files

### Data
- `functional_hash_diversity_data.csv` - Full results as CSV
- `functional_hash_diversity_data.parquet` - Optimized parquet for analysis

Parquet columns:
- `accession` - SRA accession
- `total_distinct_hashes` - Total FracMinHash hashes
- `hashes_per_mb` - Normalized hash density
- `mbases` - Sequencing depth
- `observed_richness`, `shannon_index`, `simpson_index`, etc. - Raw metrics
- `observed_richness_per_mb`, `shannon_index_per_mb`, etc. - Normalized metrics
- `alpha_diversity`, `diversity_per_mb` - Aliases for compatibility

### Plots
- `summary_all_metrics.png` - Multi-panel summary of all diversity metrics
- `hash_vs_[metric]_correlation.png` - Individual metric correlation plots
- `hash_diversity_hexbin.png` - Density plot for richness
- `distributions.png` - Distribution of all metrics

### Reports
- `analysis_report.txt` - Comprehensive text report
- `statistics_summary.csv` - Correlation statistics for all metrics

## Key Differences from Taxonomic Analysis

1. **Data source**: Uses `functional_profile.*` tables instead of `taxa_profiles.*`
2. **No coverage parameter**: Functional profiles don't have a coverage threshold like taxonomic profiles
3. **Multiple diversity metrics**: Includes abundance-aware metrics (Shannon, Simpson, Hill, etc.)
4. **k-mer size**: Uses k=11 amino acid (33 nucleotide equivalent) vs k=31 for taxonomic

## Database Queries Used

Hash count:
```sql
SELECT query_n_hashes 
FROM functional_profile_data.gather_data 
WHERE sample_id = '{accession}' LIMIT 1
```

KO abundances:
```sql
SELECT ko_id, abundance
FROM functional_profile.profiles 
WHERE sample_id = '{accession}'
```

## Dependencies

- Python 3.8+
- duckdb
- pandas
- numpy
- matplotlib
- seaborn
- scipy
- scikit-learn
- tqdm
- Optional: dcor (distance correlation), minepy (MIC)

## Example Workflow

```bash
# Step 0: Run full analysis
./0_run_functional_analysis.sh

# Step 1: Generate filtered figures
./1_get_functional_correlation_figures.sh

# Step 2: Custom analysis
python3 2_plot_functional_custom_correlation.py \
    --input /path/to/functional_hash_diversity_data.parquet \
    --output my_analysis \
    --filter "platform=ILLUMINA" \
    --filter "mbases>500" \
    --diversity-metric hill_2 \
    --color-by organism
```

## Interpreting Results

- **Observed Richness**: Direct count of functional annotations
- **Shannon Index**: Higher values = more diverse, accounts for evenness
- **Simpson/Hill-2**: "Effective" number of species; less sensitive to rare KOs
- **Berger-Parker**: Lower values = less dominance by single KO
- **Pielou Evenness**: 1.0 = perfectly even, 0 = highly uneven

For correlation interpretation:
- **Spearman ρ**: Use for monotonic relationships (recommended for "more hashes → more diversity")
- **Pearson r**: Use for linear relationships
- **Distance correlation**: Captures non-linear dependencies
- **MIC**: Detects complex functional relationships
