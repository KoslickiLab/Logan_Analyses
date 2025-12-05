# Component File Format Reference

## Current Format (JSON)

The toolkit now expects component membership in **JSON format**.

### File: `components.json`

**Structure:**
```json
{
  "component_0": [6, 7, 8, 9, 10, ...],
  "component_1": [644, 645, 646, ...],
  "component_2": [26, 27, 28, ...],
  ...
}
```

**Key points:**
- Each key is a component name in format `component_N` where N is an integer
- Each value is a list of sample IDs (integers)
- Sample IDs correspond to line numbers (0-indexed) in `accessions_mbases_geq_10.txt`
- Component size is automatically calculated from the list length

### Alternative: Communities File

If you have Louvain/Leiden community detection results, you can also use:

**File: `communities.json`**
```json
{
  "community_0": [6, 7, 12, 13, ...],
  "community_1": [113, 115, 118, ...],
  ...
}
```

Same format as components.json, just with "community" prefix instead of "component".

## Converting Between Formats

### From Old CSV Format to New JSON Format

If you have the old CSV format:
```csv
sample_id,component_id,component_size
0,0,145770
1,0,145770
2,0,145770
```

Convert to JSON:
```python
import pandas as pd
import json

# Read CSV
df = pd.read_csv('component_membership.csv')

# Group by component
components = {}
for comp_id in df['component_id'].unique():
    sample_ids = df[df['component_id'] == comp_id]['sample_id'].tolist()
    components[f'component_{comp_id}'] = sample_ids

# Save as JSON
with open('components.json', 'w') as f:
    json.dump(components, f, indent=2)
```

### From JSON Back to CSV (for compatibility)

If you need CSV format for some reason:
```python
import json
import pandas as pd

# Load JSON
with open('components.json', 'r') as f:
    components = json.load(f)

# Convert to dataframe
rows = []
for comp_name, sample_ids in components.items():
    comp_id = int(comp_name.split('_')[1])
    comp_size = len(sample_ids)
    for sample_id in sample_ids:
        rows.append({
            'sample_id': sample_id,
            'component_id': comp_id,
            'component_size': comp_size
        })

df = pd.DataFrame(rows)
df.to_csv('component_membership.csv', index=False)
```

## How the Toolkit Processes JSON

The analysis scripts automatically:

1. **Load JSON file**
   ```python
   with open('components.json', 'r') as f:
       components_json = json.load(f)
   ```

2. **Extract component ID from name**
   ```python
   component_id = int(component_name.split('_')[1])
   # "component_0" -> 0
   # "component_123" -> 123
   ```

3. **Calculate component size**
   ```python
   component_size = len(sample_ids)
   ```

4. **Map sample IDs to accessions**
   ```python
   accessions = [line.strip() for line in open('accessions_mbases_geq_10.txt')]
   accession = accessions[sample_id]
   ```

5. **Create internal dataframe**
   Same structure as old CSV format for compatibility with existing code.

## Advantages of JSON Format

1. **More compact** - No repetition of component_id and component_size
2. **Easier to inspect** - Clear component boundaries in the file
3. **Standard format** - JSON is widely supported
4. **Flexible** - Easy to add metadata at component level if needed

## Testing Your File

To verify your components.json file is formatted correctly:

```bash
python test_json_loading.py
```

This will:
- Load your JSON file
- Convert to internal format
- Show component statistics
- Verify the structure is correct

## Common Issues

### Issue: "KeyError: invalid literal for int()"
**Cause:** Component names not in format `component_N`  
**Fix:** Ensure all keys are like `component_0`, `component_1`, etc.

### Issue: "sample_id out of range"
**Cause:** Sample ID in JSON doesn't exist in accessions file  
**Fix:** Verify all sample IDs are valid indices (0 to len(accessions)-1)

### Issue: "JSON decode error"
**Cause:** Invalid JSON syntax  
**Fix:** Validate JSON with `python -m json.tool components.json`

## Example Files

See the uploaded example files:
- `components.json` - Component membership from connected components
- `communities.json` - Community detection results from Louvain

Both formats work with the toolkit!

## Summary

**Required file:** `components.json` (or `communities.json`)  
**Format:** JSON dictionary mapping component names to lists of sample IDs  
**Component names:** Must be `component_N` where N is an integer  
**Sample IDs:** Must be valid indices into `accessions_mbases_geq_10.txt`

All analysis scripts have been updated to handle this format automatically.
