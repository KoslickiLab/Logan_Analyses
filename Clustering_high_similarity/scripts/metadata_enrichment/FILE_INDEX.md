# Component Metadata Analysis Toolkit - File Index

**Created:** December 4, 2025  
**Purpose:** Analyze metadata patterns in clustered duplicated genomic samples

## Quick Start

```bash
# 1. Quick overview (30 seconds)
python quick_insights.py

# 2. Complete analysis (10-30 minutes)
python run_all_analyses.py

# Or use Make for convenience:
make quick    # Quick overview
make all      # Complete analysis
```

## Documentation Files

### START HERE
- **README.md** - Quick start guide and overview
- **WORKFLOW_SUMMARY.md** - Recommended workflow and use cases

### Reference
- **ANALYSIS_GUIDE.md** - Detailed interpretation guide
- **FILE_INDEX.md** - This file

## Analysis Scripts

### Main Pipeline
- **run_all_analyses.py** - Master script (runs all 5 analysis steps)
- **Makefile** - Convenient shortcuts

### Individual Analysis Steps
1. **quick_insights.py** - Rapid overview (30 sec)
2. **analyze_component_metadata.py** - Comprehensive analysis (5-10 min)
3. **statistical_enrichment_analysis.py** - Statistical testing (5-10 min)
4. **explore_components_interactive.py** - Outlier detection (2-5 min)
5. **example_queries.py** - Pre-built queries (1 min)

## What Each Script Produces

### quick_insights.py
**Output:** Console only  
**Purpose:** Fast orientation
- Component size distribution
- Most predictive fields
- Geographic/temporal diversity
- Initial red flags

### analyze_component_metadata.py
**Outputs:**
- `component_metadata_summary.csv` - Per-component statistics
- `merged_components_metadata.csv` - Full merged dataset

**Purpose:** Comprehensive pattern analysis
- Single-field homogeneity
- Multi-field combinations
- Outlier detection in homogeneous groups

### statistical_enrichment_analysis.py
**Outputs:**
- `statistical_enrichment_results.csv` - Chi-square test results
- `pairwise_field_interactions.csv` - Field combination effects

**Purpose:** Formal hypothesis testing
- Significant field associations
- Effect sizes (Cramér's V)
- Multiple testing correction (FDR)

### explore_components_interactive.py
**Outputs:**
- `geographic_outliers.csv` - Multi-country components
- `bioproject_mixtures.csv` - Cross-study duplications
- `temporal_patterns.csv` - Temporal clustering

**Purpose:** Find specific interesting patterns
- Geographic diversity
- Cross-bioproject duplications
- Temporal spread

### example_queries.py
**Output:** Console summaries (also importable for interactive use)

**Purpose:** Pre-built investigation queries
- Large single-center components
- Geographically diverse groups
- Cross-study duplications
- Temporally spread components

## Input Files Required

Your working directory must contain:
1. `components.json` - Component membership in JSON format
   - Format: `{"component_0": [0, 1, 2, ...], "component_1": [100, 101, ...]}`
2. `accessions_mbases_geq_10.txt` - List of accessions  
3. `metagenome_metadata_with_geo.duckdb` - Metadata database

## Output Files Generated

After running the complete analysis:

### Primary Results
- `component_metadata_summary.csv` - One row per component
- `statistical_enrichment_results.csv` - Significant associations
- `geographic_outliers.csv` - Multi-country components
- `bioproject_mixtures.csv` - Cross-study duplications
- `temporal_patterns.csv` - Time-based clustering

### Supporting Data
- `merged_components_metadata.csv` - Full dataset (can be large!)
- `pairwise_field_interactions.csv` - Field combination analysis

## Typical Workflow

### Phase 1: Reconnaissance (5 min)
```bash
make quick
```
Review console output for orientation.

### Phase 2: Full Analysis (10-30 min)
```bash
make all
```
Generates all output files.

### Phase 3: Review Results (15-30 min)
1. Open `statistical_enrichment_results.csv` - Check significant fields
2. Open `geographic_outliers.csv` - Find surprising patterns
3. Open `component_metadata_summary.csv` - Browse all components

### Phase 4: Deep Investigation (as needed)
```python
from example_queries import *
df = load_data()
investigate_component(df, component_id=12345)
```

## Key Questions Answered

### By this toolkit:
✓ Why do these samples cluster together?  
✓ Which metadata fields explain component structure?  
✓ Are there geographic outliers (same cluster, different countries)?  
✓ Are there cross-study duplications?  
✓ Which patterns are expected vs surprising?  

### Common findings:
- Technical batch effects (same center + same time)
- Data reuse (same samples in multiple studies)
- Reference samples (used across years)
- Mislabeling (geographic inconsistencies)
- Contamination (organism mismatches)

## Interpreting Results

### Expected Patterns (normal)
- Same center + same bioproject → Technical replicates
- Same instrument + same year → Batch effects
- Perfect homogeneity (100%) → Exact duplicates

### Surprising Patterns (investigate!)
- Geographic outliers → Mislabeling or contamination
- Cross-bioproject → Data reuse or shared controls
- Temporal spread → Reference samples
- Single outliers → High priority investigation

### Key Metrics
- **Entropy**: 0 = homogeneous, higher = diverse
- **Frequency**: >0.95 = nearly pure, <0.7 = diverse
- **Cramér's V**: >0.3 = strong association
- **p_value_fdr**: <0.05 = significant

## Advanced Usage

### Custom analyses
Edit field lists in scripts to add your own metadata fields.

### Interactive investigation
Import functions from `example_queries.py` for custom queries.

### Visualizations
Use matplotlib/seaborn with the merged dataset.

## Troubleshooting

**Memory issues**
- Load specific columns: `pd.read_csv('file.csv', usecols=['col1', 'col2'])`
- Filter to specific components first

**"Too sparse" warnings**
- Normal for very diverse fields
- Automatically skipped in statistical tests

**No significant results**
- Try lowering `min_component_size`
- Focus on descriptive analysis

## Dependencies

```bash
pip install pandas numpy scipy statsmodels matplotlib seaborn duckdb
```

## Support

- **Usage questions:** See README.md
- **Interpretation:** See ANALYSIS_GUIDE.md
- **Workflows:** See WORKFLOW_SUMMARY.md
- **Error messages:** Check file paths and formats

## Summary

This toolkit provides everything needed to:
1. Understand component structure
2. Identify metadata patterns
3. Find outliers and anomalies
4. Make data quality decisions

**Recommended starting point:** `make quick` then `make all`

---

**Author:** [Your name]  
**Date:** December 4, 2025  
**Version:** 1.0
