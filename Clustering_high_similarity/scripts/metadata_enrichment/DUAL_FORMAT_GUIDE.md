# Dual JSON Format Support

## Overview

The toolkit now supports **TWO JSON formats** for component membership:

1. **Index-based** (original): Integer indices that reference an accessions file
2. **String-based** (new): Accession strings directly embedded in JSON

**Format is automatically detected** - no configuration needed!

---

## Format 1: Index-Based (Original)

### Structure
```json
{
  "component_0": [0, 1, 2, 3, 4, 6],
  "component_1": [7, 8, 9]
}
```

### Requirements
- **JSON file**: `components.json` with integer arrays
- **Accessions file**: `accessions_mbases_geq_10.txt` (one accession per line)

### How it works
1. Load components.json
2. Detect: first value is `int` → index-based format
3. Load accessions file
4. Map: `component[0]` → `accessions[0]` → `'SRR2846706'`

### Example
**components.json:**
```json
{
  "component_0": [0, 1, 2]
}
```

**accessions_mbases_geq_10.txt:**
```
SRR2846706
SRR2846720
SRR2846725
```

**Result:** Component 0 contains: `['SRR2846706', 'SRR2846720', 'SRR2846725']`

### Usage
```bash
python run_all_analyses.py \
  --components components.json \
  --accessions accessions_mbases_geq_10.txt
```

---

## Format 2: String-Based (New)

### Structure
```json
{
  "component_0": ["DRR008435", "DRR008436", "SRR001234"],
  "component_1": ["ERR123456", "SRR987654"]
}
```

### Requirements
- **JSON file**: `components.json` with string arrays
- **Accessions file**: NOT NEEDED! ✨

### How it works
1. Load components.json
2. Detect: first value is `str` → string-based format
3. Skip accessions file (not needed)
4. Use strings directly: `'DRR008435'` → `'DRR008435'`

### Example
**components.json:**
```json
{
  "component_0": ["DRR008435", "DRR008436", "DRR008437"]
}
```

**Result:** Component 0 contains: `['DRR008435', 'DRR008436', 'DRR008437']`

### Usage
```bash
# Accessions file not needed!
python run_all_analyses.py \
  --components components.json
```

---

## Automatic Format Detection

The code **automatically detects** which format you're using:

```python
# Detection logic (internal)
first_value = next(iter(components_json.values()))[0]
is_index_format = isinstance(first_value, int)

if is_index_format:
    # Load accessions file and map integers
    accessions = load_accessions_file()
    mapped = [accessions[i] for i in sample_ids]
else:
    # Use strings directly
    mapped = sample_ids
```

**Detection method:**
- Checks the type of the **first value** in the first component
- `int` → Index-based format (load accessions file)
- `str` → String-based format (skip accessions file)

---

## Advantages of Each Format

### Index-Based Format

**Advantages:**
- ✅ Smaller JSON file size (integers < strings)
- ✅ Single source of truth for accessions
- ✅ Easy to update accessions without changing JSON

**Use when:**
- Working with very large datasets (millions of samples)
- Accession list might change/update
- Want to minimize JSON file size

### String-Based Format

**Advantages:**
- ✅ Self-contained (no separate file needed)
- ✅ More explicit (see accessions directly)
- ✅ Easier to share/transfer (one file)
- ✅ No risk of accessions file mismatch

**Use when:**
- Working with moderate datasets
- Sharing analysis with others
- Want explicit accession identifiers
- Prefer simplicity over file size

---

## Migration Guide

### From Index-Based to String-Based

If you have index-based JSON and want to convert:

```python
import json

# Load index-based format
with open('components.json', 'r') as f:
    components = json.load(f)

# Load accessions
with open('accessions_mbases_geq_10.txt', 'r') as f:
    accessions = [line.strip() for line in f]

# Convert to string-based format
string_based = {}
for comp_name, indices in components.items():
    string_based[comp_name] = [accessions[i] for i in indices]

# Save
with open('components_string.json', 'w') as f:
    json.dump(string_based, f, indent=2)
```

### From String-Based to Index-Based

If you want to go the other way (less common):

```python
import json

# Load string-based format
with open('components.json', 'r') as f:
    components = json.load(f)

# Collect all unique accessions
all_accessions = []
for accessions in components.values():
    all_accessions.extend(accessions)
all_accessions = list(dict.fromkeys(all_accessions))  # Deduplicate, preserve order

# Create index mapping
acc_to_idx = {acc: idx for idx, acc in enumerate(all_accessions)}

# Convert to index-based format
index_based = {}
for comp_name, accessions in components.items():
    index_based[comp_name] = [acc_to_idx[acc] for acc in accessions]

# Save
with open('components_index.json', 'w') as f:
    json.dump(index_based, f, indent=2)

with open('accessions_mbases_geq_10.txt', 'w') as f:
    f.write('\n'.join(all_accessions))
```

---

## Command-Line Usage

### Index-Based Format
```bash
# Provide both files
python run_all_analyses.py \
  --components components.json \
  --accessions accessions_mbases_geq_10.txt \
  --database metagenome_metadata_with_geo.duckdb
```

### String-Based Format
```bash
# Accessions file not needed (will show warning but continues)
python run_all_analyses.py \
  --components components.json \
  --database metagenome_metadata_with_geo.duckdb

# Or explicitly skip accessions
python run_all_analyses.py \
  --components components_string.json \
  --accessions "" \
  --database metagenome_metadata_with_geo.duckdb
```

---

## File Size Comparison

**Example with 100,000 samples:**

### Index-Based
- `components.json`: ~1.5 MB (integers)
- `accessions_mbases_geq_10.txt`: ~1.0 MB (strings)
- **Total**: 2.5 MB

### String-Based
- `components.json`: ~2.0 MB (strings embedded)
- **Total**: 2.0 MB

**Conclusion:** String-based can actually be smaller for moderate datasets due to no duplication!

---

## Troubleshooting

### Error: "list indices must be integers or slices, not str"

**Cause:** Index-based JSON format but accessions file not provided or wrong format

**Solution:** Provide correct accessions file:
```bash
python run_all_analyses.py \
  --components components.json \
  --accessions accessions_mbases_geq_10.txt
```

### Warning: "Accessions file not found"

**Cause:** Using string-based format without accessions file (this is fine!)

**Solution:** This is just a warning. If your JSON has strings, it will work fine. You can ignore this warning.

### Mixed Format Error

**Problem:** Some values are integers, some are strings in the same JSON

**Solution:** Choose one format consistently. Don't mix integers and strings in the same component.

---

## Best Practices

1. **Choose one format** - Don't mix formats in the same project
2. **Document your choice** - Note in README which format you're using
3. **Validate JSON** - Use `python test_dual_format.py` to verify format
4. **Keep backup** - When converting, keep original format as backup

---

## Testing

Test that both formats work correctly:

```bash
python test_dual_format.py
```

Expected output:
```
✓ PASS   - Index format detection
✓ PASS   - String format detection
✓ PASS   - Empty component handling

🎉 All format tests passed!
```

---

## Summary

| Feature | Index-Based | String-Based |
|---------|-------------|--------------|
| **Format** | `[0, 1, 2]` | `["SRR001", "SRR002"]` |
| **Accessions File** | Required | Not needed |
| **File Size** | Smaller JSON | Smaller total |
| **Simplicity** | Two files | One file |
| **Explicit** | Indirect | Direct |
| **Best For** | Very large datasets | Moderate datasets |

**Both formats are fully supported and automatically detected!**

Choose the format that best fits your workflow and dataset size.

---

**Version:** 1.3  
**Date:** December 4, 2025  
**Status:** ✅ Production Ready
