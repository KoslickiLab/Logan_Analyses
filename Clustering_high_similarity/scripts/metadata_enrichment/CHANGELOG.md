# Changelog - Bug Fixes and Improvements

## Version 1.1 - Bug Fixes and Command-Line Arguments

**Date:** December 4, 2025

### Summary

Fixed critical bugs, removed unused code, and added full command-line argument support for custom file names and output directories.

---

## 🐛 Bug Fixes

### 1. **CRITICAL: Missing location_metadata in quick_insights.py**

**Issue:** The script attempted to use geographic fields (`geo_loc_name_country_calc`, `country`, etc.) but never loaded the `location_metadata` table from the database.

**Impact:** Script would fail with KeyError when trying to access geographic fields.

**Fix:** Added location_metadata loading and merge:
```python
location_df = conn.execute("SELECT * FROM location_metadata").df()
merged = merged.merge(location_df, left_on='biosample', right_on='accession', 
                     how='left', suffixes=('', '_loc'))
```

**Files Modified:** `quick_insights.py`

---

## 🧹 Code Cleanup

### 2. **Removed unused `component_sizes` return value**

**Issue:** `analyze_component_sizes()` returned a value that was never used by the caller.

**Fix:** Removed the return statement and the assignment:
```python
# Before:
component_sizes = analyze_component_sizes(df)  # Never used

# After:
analyze_component_sizes(df)  # Just prints analysis
```

**Files Modified:** `analyze_component_metadata.py`

### 3. **Removed unused `date_fields` variable**

**Issue:** Variable was defined but never used anywhere in the code.

**Fix:** Deleted the line entirely:
```python
# Removed this line:
date_fields = ['releasedate', 'collection_date_sam']
```

**Files Modified:** `analyze_component_metadata.py`

---

## ✨ New Features

### 4. **Command-Line Arguments for All Scripts**

**Feature:** All scripts now accept command-line arguments for input files and output directories.

**Benefits:**
- Use custom file names without editing code
- Organize outputs into separate directories
- Run multiple analyses with different parameters
- Integrate into automated pipelines

**Scripts Updated:**
1. `quick_insights.py`
2. `analyze_component_metadata.py`
3. `statistical_enrichment_analysis.py`
4. `explore_components_interactive.py`
5. `example_queries.py`
6. `run_all_analyses.py` (master pipeline)

---

## 📖 Usage Examples

### Basic Usage (Defaults)

```bash
# Uses default files and current directory
python quick_insights.py
python analyze_component_metadata.py
python run_all_analyses.py
```

### Custom Input Files

```bash
# Use different component file
python analyze_component_metadata.py \
  --components components_subset.json \
  --accessions accessions_mbases_geq_10.txt \
  --database metagenome_metadata_with_geo.duckdb
```

### Custom Output Directory

```bash
# Save all outputs to a specific directory
python analyze_component_metadata.py \
  --output-dir results/run1/

# Master pipeline with custom output
python run_all_analyses.py \
  --components my_components.json \
  --output-dir results/experiment_2025_12_04/
```

### Scripts Reading Merged CSV

```bash
# Statistical analysis with custom input/output
python statistical_enrichment_analysis.py \
  --input results/run1/merged_components_metadata.csv \
  --output-dir results/run1/

# Interactive exploration
python explore_components_interactive.py \
  --input results/run1/merged_components_metadata.csv \
  --output-dir results/run1/

# Example queries
python example_queries.py \
  --input results/run1/merged_components_metadata.csv
```

### Help and Documentation

Every script now has built-in help:

```bash
python analyze_component_metadata.py --help
```

Output:
```
usage: analyze_component_metadata.py [-h] 
       [--components COMPONENTS]
       [--accessions ACCESSIONS]
       [--database DATABASE]
       [--output-dir OUTPUT_DIR]

Comprehensive metadata analysis for duplicated sample components

optional arguments:
  -h, --help            show this help message and exit
  --components COMPONENTS
                        Component membership file (JSON format) (default: components.json)
  --accessions ACCESSIONS
                        Accessions list file (default: accessions_mbases_geq_10.txt)
  --database DATABASE   Metadata database file (default: metagenome_metadata_with_geo.duckdb)
  --output-dir OUTPUT_DIR
                        Output directory for results (default: ./)
```

---

## 🔧 Technical Details

### Command-Line Argument Implementation

All scripts now use Python's `argparse` module with:
- Default values for backward compatibility
- Help text for each parameter
- ArgumentDefaultsHelpFormatter for clear documentation

### Output Directory Handling

- Automatically creates directories if they don't exist
- Uses `pathlib.Path` for cross-platform compatibility
- All output files are saved to the specified directory
- Preserves filename conventions

### Master Pipeline Integration

`run_all_analyses.py` now:
- Passes arguments to all sub-scripts
- Constructs correct file paths for intermediate files
- Checks for output files in the specified directory
- Supports end-to-end custom workflows

---

## 📊 Complete Argument Reference

### Arguments by Script

| Script | Input Arguments | Output Arguments |
|--------|----------------|------------------|
| `quick_insights.py` | `--components`, `--accessions`, `--database` | - |
| `analyze_component_metadata.py` | `--components`, `--accessions`, `--database` | `--output-dir` |
| `statistical_enrichment_analysis.py` | `--input` | `--output-dir` |
| `explore_components_interactive.py` | `--input` | `--output-dir` |
| `example_queries.py` | `--input` | - |
| `run_all_analyses.py` | `--components`, `--accessions`, `--database` | `--output-dir` |

### Default Values

| Argument | Default Value |
|----------|---------------|
| `--components` | `components.json` |
| `--accessions` | `accessions_mbases_geq_10.txt` |
| `--database` | `metagenome_metadata_with_geo.duckdb` |
| `--input` | `merged_components_metadata.csv` |
| `--output-dir` | `./` (current directory) |

---

## 🧪 Testing

### Verify Bug Fixes

```bash
# Test that geographic fields work in quick_insights
python quick_insights.py
# Should show geographic distribution without errors
```

### Test Custom Arguments

```bash
# Create test output directory
mkdir -p test_results

# Run with custom output
python run_all_analyses.py --output-dir test_results/

# Check outputs
ls -lh test_results/
```

### Test Different Input Files

```bash
# If you have alternative component files
python analyze_component_metadata.py \
  --components communities.json \
  --output-dir communities_analysis/
```

---

## 🔄 Migration Guide

### For Existing Workflows

**No changes required!** All scripts maintain backward compatibility with default arguments.

### To Use New Features

1. **Organize by experiment:**
   ```bash
   python run_all_analyses.py --output-dir results/experiment_001/
   python run_all_analyses.py --output-dir results/experiment_002/
   ```

2. **Use different input files:**
   ```bash
   python analyze_component_metadata.py \
     --components subset_components.json \
     --output-dir subset_analysis/
   ```

3. **Integrate into pipelines:**
   ```bash
   #!/bin/bash
   for subset in subset1 subset2 subset3; do
     python run_all_analyses.py \
       --components ${subset}_components.json \
       --output-dir results/${subset}/
   done
   ```

---

## 📝 Summary of Changes

### Files Modified
- ✅ `quick_insights.py` - Bug fix + arguments
- ✅ `analyze_component_metadata.py` - Cleanup + arguments
- ✅ `statistical_enrichment_analysis.py` - Arguments
- ✅ `explore_components_interactive.py` - Arguments
- ✅ `example_queries.py` - Arguments
- ✅ `run_all_analyses.py` - Arguments + propagation

### Lines of Code
- **Bugs Fixed:** 1 critical bug
- **Code Removed:** ~10 lines (unused variables)
- **Code Added:** ~150 lines (argparse + output handling)

### Backward Compatibility
- ✅ 100% backward compatible
- ✅ All existing commands work unchanged
- ✅ Default behavior preserved

---

## 🎯 Next Steps

You can now:

1. **Run analyses with custom files:**
   ```bash
   python run_all_analyses.py \
     --components your_components.json \
     --output-dir your_results/
   ```

2. **Organize results by experiment:**
   ```bash
   python run_all_analyses.py --output-dir exp_20251204_01/
   python run_all_analyses.py --output-dir exp_20251204_02/
   ```

3. **Process subsets separately:**
   ```bash
   python analyze_component_metadata.py \
     --components large_components.json \
     --output-dir analysis_large/
   ```

4. **Integrate into automated workflows:**
   All scripts are now pipeline-friendly with full argument support!

---

**Version:** 1.1  
**Date:** December 4, 2025  
**Status:** ✅ Complete and tested
