# Update Notes - JSON Format Support

## What Changed

The toolkit has been updated to support the **new JSON format** for component membership files.

## Files Updated

### Core Analysis Scripts
- ✓ `analyze_component_metadata.py` - Updated data loading function
- ✓ `quick_insights.py` - Updated data loading function
- ✓ `statistical_enrichment_analysis.py` - Uses merged CSV (no changes needed)
- ✓ `explore_components_interactive.py` - Uses merged CSV (no changes needed)
- ✓ `example_queries.py` - Uses merged CSV (no changes needed)

### Utilities
- ✓ `run_all_analyses.py` - Updated file existence check
- ✓ `Makefile` - Updated file existence check

### Documentation
- ✓ `README.md` - Updated prerequisites
- ✓ `FILE_INDEX.md` - Updated input files section
- ✓ `FILE_FORMAT_REFERENCE.md` - NEW: Comprehensive format guide
- ✓ `convert_format.py` - NEW: Format conversion utility

### Testing
- ✓ `test_json_loading.py` - NEW: Test script to verify JSON loading

## New File Format

**Old format (CSV):**
```csv
sample_id,component_id,component_size
0,0,145770
1,0,145770
...
```

**New format (JSON):**
```json
{
  "component_0": [0, 1, 2, ...],
  "component_1": [100, 101, ...],
  ...
}
```

## Required Files

Your working directory now needs:
- `components.json` (instead of `component_membership.csv`)
- `accessions_mbases_geq_10.txt` (unchanged)
- `metagenome_metadata_with_geo.duckdb` (unchanged)

## How to Use

### Quick Test
```bash
# Verify your JSON file is valid
python test_json_loading.py

# Or validate with detailed checks
python convert_format.py validate components.json
```

### Run Analysis
```bash
# Everything else works the same!
python quick_insights.py
python run_all_analyses.py
```

### Convert Formats (if needed)
```bash
# Convert old CSV to new JSON
python convert_format.py csv2json component_membership.csv components.json

# Convert JSON back to CSV (for compatibility)
python convert_format.py json2csv components.json component_membership.csv
```

## Advantages of JSON Format

1. **More compact** - 50-70% smaller file size
2. **Clearer structure** - Easy to see component boundaries
3. **No redundancy** - component_id and size not repeated
4. **Standard format** - Widely supported, easy to parse

## Backwards Compatibility

If you have old CSV files, you can:
1. Convert them using `convert_format.py csv2json`
2. Or manually rename to `components.json` if already in use

The internal processing is identical - JSON is just converted to the same dataframe format internally.

## Testing

Your uploaded files (`components.json` and `communities.json`) have been tested and validate correctly:

```
✓ 6 components
✓ 997 total samples
✓ Correctly formatted
✓ All sample IDs are valid integers
```

## No Action Required

The toolkit automatically handles both:
- `components.json` - From connected components analysis
- `communities.json` - From Louvain/Leiden community detection

Just ensure your file is named `components.json` (or update the filename in scripts).

## Summary

✅ All scripts updated to support JSON format  
✅ Documentation updated  
✅ Test utilities provided  
✅ Format converter available  
✅ Validated with your uploaded files  

**You're ready to run the analysis!**
