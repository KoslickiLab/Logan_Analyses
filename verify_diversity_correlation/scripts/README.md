# Hash-Diversity Correlation Analysis

A comprehensive toolkit for testing the hypothesis that distinct k-mer hashes per basepair correlate with alpha diversity in WGS metagenomic samples.

## Overview

This analysis pipeline:
1. Queries WGS metagenomic samples from NCBI's SRA metadata
2. Extracts k-mer hash counts and taxonomic data from YACHT database
3. Calculates normalized metrics (per million bases)
4. Performs correlation analysis with publication-quality visualizations
5. Supports sensitivity analysis across coverage thresholds
6. Utilizes parallel processing for high-performance computing

## Quick Start

### Recommended: Use the wrapper script

```bash
# Quick test (1,000 samples)
bash run_analysis.sh quick

# Medium analysis (50,000 samples, ~30 minutes)
bash run_analysis.sh medium

# Full analysis (all samples, may take hours)
bash run_analysis.sh full

# Sensitivity analysis
bash run_analysis.sh sensitivity
```

### Manual execution

```bash
# Basic analysis
python3 hash_diversity_correlation.py \
    --output-dir results \
    --n-samples 10000 \
    --coverage 0.0625 \
    --n-jobs 128

# Sensitivity analysis
python3 hash_diversity_sensitivity.py \
    --output-dir results_sensitivity \
    --n-samples 5000 \
    --coverages 0.0625,0.125,0.25,0.5,1.0 \
    --n-jobs 128
```

## Requirements

### Python Dependencies
```bash
pip install duckdb pandas numpy matplotlib seaborn scipy scikit-learn tqdm
```

Or use the provided requirements file:
```bash
pip install -r requirements.txt
```

### Data Requirements
- YACHT database: `/scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db`
- Metadata database: `/scratch/shared_data_new/Logan_yacht_data/metadata/aws_sra_metadata/metadata_geo_joined.duckdb`

Both databases must be accessible in read-only mode.

### Hardware Requirements
- Recommended: 256+ CPU cores, 4TB RAM
- Minimum: 64 CPU cores, 256GB RAM
- Storage: ~5GB free space for outputs

## Scripts

### 1. `hash_diversity_correlation.py`

Main analysis script that tests the correlation hypothesis.

**Key features:**
- Parallel processing for efficient data extraction
- Automatic normalization by sequencing depth
- Publication-quality visualizations
- Comprehensive statistical analysis
- Detailed reporting

**Command-line options:**
```
--output-dir, -o      Output directory (required)
--n-samples, -n       Number of samples to analyze (default: all)
--coverage, -c        Coverage threshold (default: 0.0625)
--n-jobs, -j          Number of parallel workers (default: 64)
--min-mbases, -m      Minimum megabases for inclusion (default: 100)
--random-seed         Random seed for reproducibility (default: 42)
--dpi                 DPI for figures (default: 300)
```

**Example:**
```bash
python3 hash_diversity_correlation.py \
    --output-dir results_test \
    --n-samples 10000 \
    --coverage 0.0625 \
    --n-jobs 128 \
    --min-mbases 100
```

### 2. `hash_diversity_sensitivity.py`

Sensitivity analysis script that tests robustness across coverage thresholds.

**Key features:**
- Tests multiple coverage thresholds in parallel
- Generates comparative visualizations
- Assesses robustness of correlation
- Identifies optimal parameter choices

**Command-line options:**
```
--output-dir, -o      Output directory (required)
--n-samples, -n       Number of samples (default: all)
--coverages           Comma-separated coverage values 
                      (default: 0.015625,0.03125,0.0625,0.125,0.25,0.5,1.0)
--n-jobs, -j          Number of parallel workers (default: 64)
--min-mbases, -m      Minimum megabases (default: 100)
--random-seed         Random seed (default: 42)
```

**Example:**
```bash
python3 hash_diversity_sensitivity.py \
    --output-dir results_sensitivity \
    --n-samples 5000 \
    --coverages 0.0625,0.125,0.25,0.5,1.0 \
    --n-jobs 128
```

### 3. `run_analysis.sh`

Convenient wrapper script with presets for common scenarios.

**Presets:**
- `quick` - 1,000 samples (for testing)
- `small` - 10,000 samples
- `medium` - 50,000 samples (recommended)
- `large` - 100,000 samples
- `full` - All available samples
- `sensitivity` - Sensitivity analysis with 10,000 samples

**Options:**
```
-j, --jobs N          Number of parallel jobs (default: 128)
-o, --output DIR      Output directory
-c, --coverage VAL    Coverage threshold (default: 0.0625)
-h, --help            Show help message
```

**Examples:**
```bash
# Quick test
bash run_analysis.sh quick

# Medium analysis with custom job count
bash run_analysis.sh medium --jobs 256

# Large analysis with custom output directory
bash run_analysis.sh large --output my_results --jobs 256
```

## Output Structure

```
output_directory/
├── data/
│   ├── selected_samples.csv           # List of samples used in analysis
│   ├── hash_diversity_data.csv        # Main data with all metrics
│   └── sensitivity_results.csv        # Sensitivity analysis results
├── plots/
│   ├── hash_diversity_correlation.png # Main scatter plot with regression
│   ├── hash_diversity_hexbin.png      # Density visualization
│   ├── distributions.png              # Distribution plots
│   ├── residuals.png                  # Residual analysis
│   ├── sensitivity_analysis.png       # Sensitivity results (if applicable)
│   ├── mean_values_vs_coverage.png    # Coverage trends (if applicable)
│   └── metrics_heatmap.png            # Summary heatmap (if applicable)
└── reports/
    ├── analysis_report.txt            # Comprehensive text report
    ├── statistics_summary.csv         # Statistics in CSV format
    └── sensitivity_report.txt         # Sensitivity analysis report
```

## Key Outputs Explained

### Main Correlation Plot
- **File:** `plots/hash_diversity_correlation.png`
- **Content:** Scatter plot of hashes per Mb vs diversity per Mb
- **Features:** 
  - Regression line
  - Statistics box (r, p-value, R², n)
  - Publication-quality formatting

### Analysis Report
- **File:** `reports/analysis_report.txt`
- **Content:**
  - Data summary statistics
  - Correlation coefficients (Pearson, Spearman)
  - Linear regression results
  - Interpretation of findings
  - Recommendations for further analysis

### Data Files
- **File:** `data/hash_diversity_data.csv`
- **Columns:**
  - `sample_id` - Sample accession
  - `num_hashes` - Total distinct k-mer hashes
  - `alpha_diversity` - Number of distinct taxa
  - `mbases` - Sequencing depth in megabases
  - `hashes_per_mb` - Normalized hash density
  - `diversity_per_mb` - Normalized diversity
  - `num_records` - Number of taxonomic hits

## Methodology

### Sample Selection
1. Query metadata database for WGS metagenomic samples
2. Apply filters:
   - Assay type: WGS
   - Library selection: RANDOM
   - Minimum sequencing depth: 100 Mbases
3. Optional: Random sampling for manageable dataset size

### Data Extraction
1. Query YACHT database for each sample at specified coverage threshold
2. Extract:
   - `num_total_kmers_in_sample_sketch` - distinct hash count
   - Unique `organism_id` / `tax_id` - for alpha diversity
3. Parallel processing across CPU cores for efficiency

### Normalization
To control for sequencing depth variation:
- **Hashes per Mb** = total_hashes / mbases
- **Diversity per Mb** = alpha_diversity / mbases

This approach is analogous to RNA-seq normalization (RPKM/TPM), accounting for the fact that deeper sequencing naturally yields more hashes and detects more species.

### Statistical Analysis
1. **Pearson correlation** - Tests linear relationship
2. **Spearman correlation** - Tests monotonic relationship
3. **Linear regression** - Models relationship and calculates R²
4. **Residual analysis** - Checks for systematic biases

### Sensitivity Analysis
Tests correlation across multiple coverage thresholds to assess:
- Robustness of findings
- Optimal parameter choice
- Trade-offs between sensitivity and specificity

## Interpretation Guide

### Correlation Strength (|r|)
- 0.7-1.0: Strong correlation
- 0.4-0.7: Moderate correlation  
- 0.2-0.4: Weak correlation
- 0.0-0.2: Very weak correlation

### Statistical Significance (p-value)
- p < 0.001: Highly significant (★★★)
- p < 0.01: Significant (★★)
- p < 0.05: Significant (★)
- p ≥ 0.05: Not significant (ns)

### R² (Coefficient of Determination)
Percentage of variance in diversity explained by hash count:
- R² = 0.7: 70% of variance explained (strong model)
- R² = 0.4: 40% of variance explained (moderate model)
- R² = 0.1: 10% of variance explained (weak model)

### Positive Correlation
If hashes per Mb correlate positively with diversity per Mb:
- Supports hypothesis that hash-based metrics capture diversity
- Suggests k-mer complexity reflects taxonomic richness
- Validates use of hash counts as proxy for diversity

### Negative Correlation
Would indicate unexpected relationship and warrant investigation of:
- Technical artifacts in sequencing or analysis
- Biases in sample types or environments
- Coverage threshold effects

## Performance Notes

### Expected Runtime
(with 128 cores, depends on database speed)

| Preset      | Samples  | Approximate Time |
|-------------|----------|------------------|
| Quick       | 1,000    | 2-5 minutes      |
| Small       | 10,000   | 10-20 minutes    |
| Medium      | 50,000   | 30-60 minutes    |
| Large       | 100,000  | 1-2 hours        |
| Full        | ~500,000 | 4-8 hours        |
| Sensitivity | varies   | 30-90 minutes    |

### Optimization Tips
1. **Increase parallelism:** Use more cores with `--n-jobs`
2. **Sample strategically:** Start with smaller n-samples for pilot analysis
3. **Database indexing:** If possible, index the YACHT database (currently not indexed)
4. **Batch size:** Script auto-adjusts, but can be tuned in code
5. **Memory:** Ensure sufficient RAM; script uses in-memory processing

## Troubleshooting

### "No data extracted" error
- Check that YACHT database is accessible at specified path
- Verify samples exist in database at specified coverage threshold
- Try different coverage threshold or reduce n-samples

### "Database locked" error
- Ensure databases are opened in read-only mode
- Check that no other process has exclusive lock
- Wait and retry

### Slow performance
- Increase `--n-jobs` to use more CPU cores
- Reduce `--n-samples` for faster testing
- Check database indexing status
- Monitor I/O with `iostat -x 1`

### Memory errors
- Reduce `--n-samples`
- Increase available RAM
- Check for memory leaks with `top` or `htop`

### Missing dependencies
```bash
pip install duckdb pandas numpy matplotlib seaborn scipy scikit-learn tqdm
```

## Citation

If you use this analysis in your research, please cite:
- YACHT algorithm paper (reference to be added)
- This analysis pipeline (GitHub URL)

## Advanced Usage

### Custom Coverage Thresholds
```bash
python3 hash_diversity_correlation.py \
    --output-dir results_custom \
    --n-samples 10000 \
    --coverage 0.125 \
    --n-jobs 128
```

### Stratified Analysis by Environment
Modify script to filter samples by environment type from metadata.

### Integration with Existing Pipelines
Import functions from scripts:
```python
from hash_diversity_correlation import (
    get_wgs_samples, 
    extract_hash_and_diversity_data,
    perform_correlation_analysis
)
```

### Batch Processing Multiple Coverage Values
```bash
for cov in 0.0625 0.125 0.25 0.5 1.0; do
    python3 hash_diversity_correlation.py \
        --output-dir "results_coverage_${cov}" \
        --n-samples 10000 \
        --coverage $cov \
        --n-jobs 128
done
```

### Custom Visualization
Use `data/hash_diversity_data.csv` for custom plotting in R, Python, or other tools.

## Contact & Support

For questions, issues, or contributions:
- Open an issue on GitHub
- Contact the research team
- Check documentation for updates

## License

[To be determined based on your project requirements]

---

**Last Updated:** December 2024
