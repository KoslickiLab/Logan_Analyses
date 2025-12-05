# Component Metadata Analysis Guide

## Overview

This analysis suite helps identify **what metadata patterns explain the clustering of duplicated samples** (samples with Jaccard similarity = 1.0 with at least one other sample). The goal is to find both expected patterns (e.g., samples from the same sequencing center) and surprising outliers (e.g., geographically distant samples clustering together).

## Analysis Philosophy

### Three-Tier Approach

1. **Homogeneity Analysis**: Find components where samples share common metadata values
   - High homogeneity → likely technical duplicates or batch effects
   - Examples: all from same center, same bioproject, same instrument

2. **Statistical Enrichment**: Formal hypothesis testing to identify which metadata fields 
   are significantly associated with component structure
   - Controls for multiple testing
   - Quantifies effect sizes

3. **Outlier Detection**: Find samples that break otherwise consistent patterns
   - Geographic outliers (different countries in same component)
   - Bioproject mixtures (samples from multiple studies)
   - Temporal anomalies (old and new samples together)

## Quick Start

### Step 1: Initial Analysis

Run the main metadata analysis:

```bash
python analyze_component_metadata.py
```

**What it does:**
- Loads and merges all data sources
- Analyzes component size distribution
- Tests single-field enrichment (which metadata fields are homogeneous within components)
- Tests multi-field combinations (e.g., center + bioproject)
- Detects outliers in highly homogeneous components
- Saves: `component_metadata_summary.csv`, `merged_components_metadata.csv`

**What to look for:**
- Fields with high "frequency" values (>0.8) = very homogeneous components
- Low entropy = predictable patterns
- Large components with high homogeneity = systematic duplications

### Step 2: Statistical Testing

Run formal statistical tests:

```bash
python statistical_enrichment_analysis.py
```

**What it does:**
- Chi-square tests for each metadata field
- Multiple testing correction (FDR)
- Effect size calculation (Cramér's V)
- Pairwise field interaction analysis
- Saves: `statistical_enrichment_results.csv`, `pairwise_field_interactions.csv`

**What to look for:**
- Significant fields (p_value_fdr < 0.05) with high Cramér's V = strong predictors
- High "pct_dominated_components" = field explains many components well
- Field pairs with high "mean_purity" = combinations that strongly predict components

### Step 3: Outlier Detection

Run targeted outlier analyses:

```bash
python explore_components_interactive.py
```

**What it does:**
- Geographic outliers (samples from different countries in same component)
- Bioproject mixtures (cross-study duplications)
- Temporal clustering patterns
- Detailed exploration of specific components
- Saves: `geographic_outliers.csv`, `bioproject_mixtures.csv`, `temporal_patterns.csv`

**What to look for:**
- Components with 70-99% from one country (the exceptions are interesting!)
- Multi-bioproject components (data reuse or contamination?)
- Tight temporal clustering (batch effects?) vs wide temporal spread (ongoing duplications?)

## Interpreting Results

### Expected Patterns (Not Surprising)

These are *expected* and indicate technical/administrative duplications:

1. **Same sequencing center + bioproject**: Likely the same lab sequencing samples multiple times
2. **Same instrument + same year**: Batch effects from a sequencing run
3. **Same bioproject + same country**: Coordinated study with duplicate uploads
4. **Perfect homogeneity (100%)**: Exact technical replicates

### Surprising Patterns (Interesting!)

These warrant deeper investigation:

1. **Geographic mixing**: Samples from USA + China in same component
   - Could indicate: sample mislabeling, contamination, or genuine biological similarity
   - Look at: organism type, isolation source, collection context

2. **Cross-bioproject duplications**: Samples from completely different studies
   - Could indicate: sample reuse, shared controls, or data contamination
   - Check: dates, centers, sample names for clues

3. **Temporal spread**: Old (2012) and new (2023) samples clustering together
   - Could indicate: reference samples used across years, or long-running studies
   - Check: if it's a standard strain or control

4. **Single outlier in homogeneous component**: 99 samples from China, 1 from USA
   - Most interesting! Could be mislabeled, contaminated, or biologically remarkable
   - Priority for manual investigation

### Key Metrics Explained

**Entropy**: 
- 0 = perfectly homogeneous (all samples same value)
- Higher = more diversity
- Use to rank how "pure" components are

**Frequency**:
- Proportion of samples with the most common value
- >0.8 = dominated by one value
- >0.95 = nearly pure

**Cramér's V**:
- Effect size for chi-square test
- 0.0-0.1 = negligible, 0.1-0.3 = small, 0.3-0.5 = medium, >0.5 = large
- Tells you how strongly a field predicts component membership

**Component Purity** (for field pairs):
- Mean fraction of samples in components sharing the same combination
- >0.8 = strong predictor
- Use to find best field combinations

## Investigation Workflow

### For a Specific Component

```python
# In Python/IPython:
import pandas as pd
from explore_components_interactive import *

df = load_merged_data()

# Explore a specific component
explore_component(df, component_id=12345)

# Compare two components
compare_components(df, comp_id1=100, comp_id2=200)
```

### Finding Interesting Components

```python
# Load summary
summary = pd.read_csv('component_metadata_summary.csv')

# Find large components with mixed geography
mixed_geo = summary[
    (summary['size'] >= 50) & 
    (summary['geo_loc_name_country_calc_n_unique'] > 3)
]

# Find single-bioproject components (pure studies)
pure_studies = summary[
    (summary['size'] >= 20) &
    (summary['bioproject_n_unique'] == 1)
]

# Find cross-bioproject duplications
cross_study = summary[
    (summary['size'] >= 20) &
    (summary['bioproject_n_unique'] > 2)
]
```

## Output Files Reference

### From analyze_component_metadata.py
- `component_metadata_summary.csv`: One row per component with summary statistics
- `merged_components_metadata.csv`: Full dataset (all samples with all metadata)

### From statistical_enrichment_analysis.py
- `statistical_enrichment_results.csv`: Chi-square test results for each metadata field
- `pairwise_field_interactions.csv`: Results for field combinations

### From explore_components_interactive.py
- `geographic_outliers.csv`: Components with samples from multiple countries
- `bioproject_mixtures.csv`: Components spanning multiple bioprojects
- `temporal_patterns.csv`: Temporal clustering analysis

## Common Questions

**Q: Why are some components so large (>100,000 samples)?**

A: Likely represents systematic duplication at scale - same samples uploaded many times,
   or a very common control/reference being reused. Check if it's dominated by one center
   or bioproject.

**Q: What's the difference between a "surprising" and "expected" pattern?**

A: Expected patterns have >95% homogeneity in predictable fields (center, bioproject).
   Surprising patterns show either: (1) high homogeneity in unexpected fields, or 
   (2) outliers breaking otherwise consistent patterns.

**Q: How do I prioritize which components to investigate?**

A: Priority order:
   1. Geographic outliers in large components
   2. Cross-bioproject duplications in medium-sized components (20-100 samples)
   3. Large components with temporal spread
   4. Single outliers in otherwise pure components

**Q: What if I want to analyze a specific metadata field not in the scripts?**

A: Edit the `single_value_fields` or `comparison_fields` lists in the scripts to add
   your field. The analysis will automatically include it.

## Advanced Usage

### Custom Enrichment Tests

```python
from statistical_enrichment_analysis import test_field_enrichment, load_data

df = load_data()

# Test a specific field with custom parameters
result = test_field_enrichment(
    df, 
    field='your_custom_field',
    min_component_size=50,  # Only large components
    max_categories=20       # Limit field values
)
```

### Component-Specific Analysis

```python
# Analyze metadata co-occurrence within a component
comp_df = df[df['component_id'] == 12345]

# Cross-tabulation
pd.crosstab(comp_df['center_name'], comp_df['bioproject'])

# Find unique sample characteristics
for field in ['organism', 'country', 'instrument']:
    print(f"\n{field}:")
    print(comp_df[field].value_counts())
```

## Troubleshooting

**"Too sparse (>80% zeros), skipping..."**
- This field has too many categories for reliable statistical testing
- This is expected for very diverse fields (like individual sample names)

**"No results to analyze"**
- No metadata fields passed the sparsity filter
- Try lowering min_component_size or focusing on a subset of components

**Memory issues with merged_components_metadata.csv**
- The full merged dataset can be large
- Load specific columns only: `pd.read_csv('file.csv', usecols=['col1', 'col2'])`
- Or filter to specific components first

## Tips for Biological Interpretation

1. **Check organism field**: Human gut microbiome duplications have different implications
   than soil or marine samples

2. **Check isolation source**: Duplications of reference materials (mock communities, 
   standards) are expected and not concerning

3. **Check sample names/accessions**: Sometimes the naming pattern gives clues
   (e.g., "replicate", "control", "standard")

4. **Cross-reference with BioProject descriptions**: The original study design
   might explain apparent duplications

5. **Consider technical vs biological**: Perfect duplicates (J=1.0) are almost always
   technical, not biological replicates

## Next Steps After Analysis

1. **Export candidate components for manual review**: Focus on surprising patterns
2. **Validate with original data**: Check if FASTQ files are actually identical
3. **Contact submitters**: For suspicious duplications, notify data submitters
4. **Document patterns**: Create a report of common duplication causes
5. **Develop filters**: Use learned patterns to flag duplicates in new data

## Contact & Support

For questions about:
- Analysis interpretation: Review this guide's "Interpreting Results" section
- Script errors: Check file paths and data formats match expected structure
- Custom analyses: Modify scripts or create new ones following existing patterns
