# Component Metadata Analysis - Complete Workflow Summary

## What You Have

A complete analysis toolkit to investigate metadata patterns in your clustered duplicated samples. This addresses your core question: **"Why do these samples cluster together?"**

## The Complete Toolkit

### Analysis Scripts (5)
1. **quick_insights.py** - Rapid 30-second overview
2. **analyze_component_metadata.py** - Comprehensive pattern analysis  
3. **statistical_enrichment_analysis.py** - Formal hypothesis testing
4. **explore_components_interactive.py** - Outlier detection
5. **example_queries.py** - Pre-built investigation queries

### Utilities
6. **run_all_analyses.py** - Master pipeline (runs all 5 scripts)
7. **Makefile** - Convenient shortcuts (e.g., `make quick`)

### Documentation
8. **README.md** - Quick start guide
9. **ANALYSIS_GUIDE.md** - Detailed interpretation
10. **WORKFLOW_SUMMARY.md** - This file

## Recommended Workflow

### Phase 1: Initial Reconnaissance (5 minutes)

```bash
# Get oriented with a quick overview
make quick
# or: python quick_insights.py
```

**What you'll learn:**
- How many components you have
- Typical component sizes
- Which metadata fields appear most homogeneous
- Initial red flags

**Decision point:** Does this match your expectations? Any surprises?

### Phase 2: Comprehensive Analysis (10-20 minutes)

```bash
# Run the full pipeline
make all
# or: python run_all_analyses.py
```

**What you'll get:**
- 7 CSV files with detailed results
- Console output with key findings
- Ready for investigation

**Outputs:**
- `component_metadata_summary.csv` - Component-level statistics
- `merged_components_metadata.csv` - Full dataset
- `statistical_enrichment_results.csv` - Significant associations
- `pairwise_field_interactions.csv` - Field combinations
- `geographic_outliers.csv` - Multi-country components
- `bioproject_mixtures.csv` - Cross-study duplications
- `temporal_patterns.csv` - Time-based patterns

### Phase 3: Review Results (15-30 minutes)

**Priority order:**

1. **Start with statistical results**
   ```bash
   # Open in your favorite spreadsheet program or:
   import pandas as pd
   stats = pd.read_csv('statistical_enrichment_results.csv')
   stats[stats['significant_fdr']]  # Significant fields only
   ```
   
   Look for:
   - Fields with high Cramér's V (>0.3) = strong predictors
   - High `pct_dominated_components` = explains many groups

2. **Check for surprises**
   ```bash
   # Geographic outliers
   geo = pd.read_csv('geographic_outliers.csv')
   
   # Cross-study duplications  
   cross = pd.read_csv('bioproject_mixtures.csv')
   ```
   
   Look for:
   - Large components with diverse geography
   - Components spanning many bioprojects

3. **Review component summary**
   ```bash
   summary = pd.read_csv('component_metadata_summary.csv')
   
   # Sort by size to see largest groups
   summary.sort_values('size', ascending=False)
   
   # Find highly homogeneous groups
   summary[summary['center_name_top_freq'] > 0.95]
   ```

### Phase 4: Deep Investigation (as needed)

For interesting components, investigate interactively:

```python
# In Python/IPython
from example_queries import *

df = load_data()

# Investigate specific components
investigate_component(df, component_id=12345)

# Compare different groups
compare_two_components(df, comp_id1=100, comp_id2=200)

# Run custom queries
df[(df['component_size'] >= 50) & 
   (df['center_name'] == 'YOUR_CENTER')]
```

### Phase 5: Document Findings

Create a findings report with:

1. **Expected patterns** (document but not surprising)
   - "Components 1-500: All from sequencing center X, single bioproject"
   
2. **Systematic issues** (important for data quality)
   - "Component 0: 145,770 samples, 95% from center Y - likely batch upload"
   
3. **Surprising patterns** (warrant investigation)
   - "Component 1234: 50 samples, 48 from China, 2 from USA"
   - "Component 5678: Spans 3 bioprojects from different years"

4. **Recommendations**
   - Which components to flag as technical duplicates
   - Which samples might be mislabeled
   - Data quality issues to report

## Common Use Cases

### Use Case 1: Finding Technical Batch Effects

**Goal:** Identify large groups of duplicates from same center/run

```python
# Use example_queries.py
large_centers = find_large_single_center_components(df, min_size=100)

# These can probably be filtered out as batch duplications
```

### Use Case 2: Detecting Data Contamination

**Goal:** Find samples that shouldn't cluster together

```python
# Geographic mixing
geo_mixed = find_geographically_diverse_components(df, min_size=20)

# Cross-study duplications
cross_study = find_cross_study_components(df, min_size=10)

# Mixed organisms (very suspicious!)
mixed_org = find_mixed_organism_components(df, min_size=10)
```

### Use Case 3: Identifying Reference Samples

**Goal:** Find samples reused across time/studies

```python
# Temporally spread
temporal = find_temporally_spread_components(df, min_size=20, min_year_range=5)

# These might be reference strains or mock communities
```

### Use Case 4: Investigating Specific Outliers

**Goal:** Understand why one sample doesn't fit its group

```python
# Investigate the component
comp_df = investigate_component(df, component_id=12345)

# Look at the outlier sample specifically
outlier = comp_df[comp_df['accession'] == 'SRR123456']

# Compare all metadata fields
for col in ['center_name', 'bioproject', 'country', 'organism']:
    print(f"{col}: {outlier[col].iloc[0]}")
```

## Tips for Success

### 1. Start Small
Don't try to understand everything at once. Focus on:
- Largest components first
- Most statistically significant fields
- Most obvious outliers

### 2. Look for Patterns
Common explanations for duplications:
- Same center + same month = batch upload
- Same bioproject + same country = coordinated study
- Same instrument + same year = sequencing run issue
- Different bioprojects + same samples = data reuse

### 3. Validate Findings
For suspicious patterns:
- Check original BioProject descriptions
- Look at sample metadata in SRA
- Compare FASTQ checksums if available
- Contact data submitters if needed

### 4. Document Everything
Keep notes on:
- Interesting components and why they're interesting
- Explanations for patterns
- Components to filter out
- Follow-up investigations needed

## Quick Reference Commands

```bash
# Complete analysis
make all

# Quick overview only
make quick

# Individual steps
make full      # Comprehensive analysis
make stats     # Statistical testing
make explore   # Outlier detection
make queries   # Example queries

# Clean up
make clean
```

## Interpreting Results - Decision Tree

```
For each component, ask:

1. Is it HOMOGENEOUS (>95% same center/bioproject)?
   YES → Likely technical duplication
         └─> Document and potentially filter
   NO  → Continue to #2

2. Does it span MULTIPLE BIOPROJECTS?
   YES → Data reuse or contamination
         └─> Priority investigation
   NO  → Continue to #3

3. Does it have GEOGRAPHIC OUTLIERS?
   YES → Mislabeling or biological similarity
         └─> High priority investigation
   NO  → Continue to #4

4. Does it span MANY YEARS?
   YES → Reference sample or ongoing duplication
         └─> Check if it's a known standard
   NO  → Continue to #5

5. Is it SMALL (<10 samples)?
   YES → Limited statistical power
         └─> Lower priority
   NO  → Medium-sized heterogeneous group
         └─> Manual investigation needed
```

## Next Steps After Analysis

1. **Create filtered dataset**: Remove obvious technical duplicates
2. **Flag suspicious samples**: Mark for manual review
3. **Report to databases**: Notify submitters of apparent issues
4. **Document patterns**: Use for future data QC
5. **Publish findings**: Contribute to data quality literature

## Getting Help

1. **For usage questions**: See README.md
2. **For interpretation**: See ANALYSIS_GUIDE.md  
3. **For error messages**: Check file paths and data formats
4. **For custom analyses**: Modify scripts or create new ones

## Summary

You now have everything needed to:
- ✓ Understand why duplicated samples cluster together
- ✓ Identify expected vs surprising patterns
- ✓ Find outliers and anomalies
- ✓ Document findings systematically
- ✓ Make data quality decisions

**Start with:** `make quick` to get oriented, then `make all` for complete analysis.

Good luck with your investigation!
