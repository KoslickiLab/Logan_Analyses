# Component Metadata Enrichment Analysis

**Analyze metadata patterns in clustered duplicated genomic samples**

This toolkit helps you understand **why samples with Jaccard similarity = 1.0 cluster together** by analyzing their metadata. Find both expected patterns (same center, same study) and surprising anomalies (geographic outliers, cross-study duplications).

## Quick Start

### Prerequisites

Your working directory should contain:
- `components.json` - Component membership (JSON format)
- `accessions_mbases_geq_10.txt` - List of accessions
- `metagenome_metadata_with_geo.duckdb` - Metadata database

Required Python packages:
```bash
pip install pandas numpy scipy statsmodels matplotlib seaborn duckdb
```

### Run Complete Analysis

```bash
# Run everything at once (5-30 minutes)
python run_all_analyses.py
```

This executes the full pipeline and generates all output files.

### Test Your Files

Before running the full analysis, verify your JSON file is valid:

```bash
# Quick test
python test_json_loading.py

# Detailed validation
python convert_format.py validate components.json
```

### Or Run Step-by-Step

```bash
# 1. Quick overview (fast, ~30 seconds)
python quick_insights.py

# 2. Comprehensive analysis (~5-10 minutes)
python analyze_component_metadata.py

# 3. Statistical testing (~5-10 minutes)  
python statistical_enrichment_analysis.py

# 4. Outlier detection (~2-5 minutes)
python explore_components_interactive.py

# 5. Example queries (~1 minute)
python example_queries.py
```

## What Each Script Does

### 1. `quick_insights.py`
**Purpose**: Rapid overview before detailed analysis  
**Runtime**: ~30 seconds  
**Outputs**: Console summary only  
**Use when**: You want a quick scan of the data

Shows:
- Component size distribution
- Most predictive metadata fields
- Geographic/temporal diversity
- Red flags (very large homogeneous groups)

### 2. `analyze_component_metadata.py`
**Purpose**: Comprehensive metadata pattern analysis  
**Runtime**: ~5-10 minutes  
**Outputs**: 
- `component_metadata_summary.csv` - Per-component statistics
- `merged_components_metadata.csv` - Full dataset

Shows:
- Single-field homogeneity (which fields are pure within components)
- Multi-field combinations (e.g., center + bioproject)
- Outliers in homogeneous components

### 3. `statistical_enrichment_analysis.py`
**Purpose**: Formal statistical hypothesis testing  
**Runtime**: ~5-10 minutes  
**Outputs**:
- `statistical_enrichment_results.csv` - Chi-square test results
- `pairwise_field_interactions.csv` - Field combination effects

Shows:
- Statistically significant field associations
- Effect sizes (Cramér's V)
- Best field combinations for predicting components

### 4. `explore_components_interactive.py`
**Purpose**: Find specific types of interesting patterns  
**Runtime**: ~2-5 minutes  
**Outputs**:
- `geographic_outliers.csv` - Multi-country components
- `bioproject_mixtures.csv` - Cross-study duplications
- `temporal_patterns.csv` - Year-based clustering

Shows:
- Components with samples from multiple countries
- Components spanning multiple bioprojects
- Temporal clustering patterns

### 5. `example_queries.py`
**Purpose**: Pre-built queries for common investigations  
**Runtime**: ~1 minute  
**Outputs**: Console summaries + can be imported for interactive use

Shows:
- Large single-center components (batch effects)
- Geographic diversity (surprising)
- Cross-study duplications (data reuse)
- Temporal spread (reference samples)

## Understanding Results

### Expected Patterns ✓

These are **normal** and indicate technical duplications:

1. **Same center + same bioproject** → Lab sequenced samples multiple times
2. **Same instrument + same year** → Batch effects from a sequencing run  
3. **Perfect homogeneity (100%)** → Exact technical replicates

### Surprising Patterns ⚠

These are **interesting** and warrant investigation:

1. **Geographic outliers** → Samples from USA + China in same component
   - Could be: mislabeling, contamination, or biological similarity
   
2. **Cross-bioproject** → Samples from different studies clustering
   - Could be: sample reuse, shared controls, data contamination
   
3. **Temporal spread** → Old (2012) + new (2023) samples together
   - Could be: reference strains used across years
   
4. **Single outlier** → 99 from China, 1 from USA
   - Most interesting! Priority for investigation

## Key Output Files

### `component_metadata_summary.csv`
One row per component with summary statistics:
- `size` - Number of samples
- `{field}_top` - Most common value for each field
- `{field}_top_freq` - Frequency of most common value
- `{field}_n_unique` - Number of unique values

**Use for**: Quickly scanning all components, sorting by interesting properties

### `statistical_enrichment_results.csv`  
Statistical test results for each metadata field:
- `p_value_fdr` - Corrected p-value (< 0.05 = significant)
- `cramers_v` - Effect size (> 0.3 = strong association)
- `pct_dominated_components` - % of components with >80% one value

**Use for**: Identifying which fields most strongly predict components

### `geographic_outliers.csv`
Components with samples from multiple countries:
- `dominant_country` - Most common country
- `dominant_pct` - Percentage from dominant country
- `outlier_countries` - Other countries represented

**Use for**: Finding geographic anomalies

### `merged_components_metadata.csv`
Complete dataset (warning: can be large!):
- All samples with all metadata fields
- Use for custom analyses

**Use for**: Interactive investigation, custom queries

## Interactive Investigation

For deeper investigation, use Python interactively:

```python
# Load the analysis functions
from example_queries import *

# Load data
df = load_data()

# Investigate a specific component
investigate_component(df, component_id=12345)

# Compare two components
compare_two_components(df, comp_id1=100, comp_id2=200)

# Find components with specific properties
large_centers = find_large_single_center_components(df, min_size=100)
geo_diverse = find_geographically_diverse_components(df, min_size=20)
cross_study = find_cross_study_components(df, min_size=10)
```

## Interpretation Guide

### Metrics Explained

**Entropy**: Measure of diversity within a component
- 0.0 = perfectly homogeneous (all same value)
- Higher = more diversity
- Use to find "pure" components

**Frequency**: Proportion with most common value
- >0.95 = nearly pure (look for outliers)
- 0.7-0.9 = dominated but mixed (interesting!)
- <0.7 = diverse

**Cramér's V**: Effect size for chi-square test
- 0.0-0.1 = negligible association
- 0.1-0.3 = small association  
- 0.3-0.5 = medium association
- >0.5 = strong association

**Component Purity**: For field combinations
- >0.8 = strong predictor
- Find best multi-field signatures

### Prioritizing Investigations

**High Priority:**
1. Geographic outliers in large components
2. Cross-bioproject duplications (20-100 samples)
3. Single outliers in pure components

**Medium Priority:**
4. Large single-center components (document for filtering)
5. Temporal anomalies (old + new samples)
6. Multi-field pattern exceptions

**Low Priority:**
7. Small components (<10 samples) - limited statistical power
8. Expected homogeneity - normal technical duplicates

## Advanced Usage

### Custom Field Analysis

Edit the field lists in any script to add your own metadata fields:

```python
# In analyze_component_metadata.py
single_value_fields = [
    'center_name',
    'bioproject',
    'YOUR_CUSTOM_FIELD',  # Add here
]
```

### Custom Queries

Create your own queries using pandas:

```python
import pandas as pd

df = pd.read_csv('merged_components_metadata.csv')

# Find components with specific property
interesting = df[
    (df['component_size'] >= 50) &
    (df['YOUR_CONDITION'])
].groupby('component_id')

# Analyze specific metadata
for comp_id, comp_df in interesting:
    print(f"Component {comp_id}:")
    print(comp_df['your_field'].value_counts())
```

### Visualization

Create custom plots:

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Component size distribution
plt.figure(figsize=(12, 6))
sizes = df.groupby('component_id')['component_size'].first()
plt.hist(sizes, bins=50, log=True)
plt.xlabel('Component Size')
plt.ylabel('Count (log scale)')
plt.title('Component Size Distribution')
plt.savefig('component_sizes.png')
```

## Troubleshooting

**"File not found" errors**
- Ensure all input files are in the current directory
- Check file names match exactly (case-sensitive)

**Memory errors with merged_components_metadata.csv**
- Load specific columns only: `pd.read_csv('file.csv', usecols=['col1', 'col2'])`
- Filter to specific components before loading full data

**"Too sparse" warnings in statistical tests**
- Normal for very diverse fields (e.g., sample names)
- These fields are automatically skipped

**No significant results after FDR correction**
- Try lowering `min_component_size` to include more data
- Or focus on descriptive analysis rather than statistical testing

## Files in This Package

### Analysis Scripts
- `run_all_analyses.py` - Master pipeline (run everything)
- `quick_insights.py` - Rapid overview  
- `analyze_component_metadata.py` - Comprehensive analysis
- `statistical_enrichment_analysis.py` - Statistical testing
- `explore_components_interactive.py` - Outlier detection
- `example_queries.py` - Pre-built investigation queries

### Documentation
- `README.md` - This file (quick start)
- `ANALYSIS_GUIDE.md` - Detailed interpretation guide

## Citation

If you use this analysis toolkit in your research, please cite:

[Your paper reference here]

## Support

For detailed interpretation guidance, see `ANALYSIS_GUIDE.md`.

For questions or issues, check:
1. Input file formats match expected structure
2. All required Python packages are installed
3. Sufficient memory available (8GB+ recommended for large datasets)

## License

[Your license here]
