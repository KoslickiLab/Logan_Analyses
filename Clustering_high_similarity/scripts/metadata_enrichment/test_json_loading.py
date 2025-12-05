#!/usr/bin/env python3
"""
Test script to verify JSON component loading works correctly
"""

import json
import pandas as pd

def test_json_loading():
    """Test loading the new JSON format"""
    
    print("Testing JSON component loading...")
    print("="*60)
    
    # Test with the uploaded file
    component_file = '/mnt/user-data/uploads/components.json'
    
    # Load JSON
    print(f"\n1. Loading {component_file}...")
    with open(component_file, 'r') as f:
        components_json = json.load(f)
    
    print(f"   Found {len(components_json)} components")
    
    # Show component names
    print("\n2. Component names:")
    for name in list(components_json.keys())[:10]:
        print(f"   - {name}")
    if len(components_json) > 10:
        print(f"   ... and {len(components_json) - 10} more")
    
    # Convert to dataframe
    print("\n3. Converting to dataframe format...")
    component_data = []
    for component_name, sample_ids in components_json.items():
        component_id = int(component_name.split('_')[1])
        component_size = len(sample_ids)
        
        for sample_id in sample_ids:
            component_data.append({
                'sample_id': sample_id,
                'component_id': component_id,
                'component_size': component_size
            })
    
    df = pd.DataFrame(component_data)
    
    print(f"   Created dataframe with {len(df):,} rows")
    print(f"   Columns: {list(df.columns)}")
    
    # Show statistics
    print("\n4. Component statistics:")
    print(f"   Total samples: {len(df):,}")
    print(f"   Unique components: {df['component_id'].nunique()}")
    print(f"   Component sizes:")
    
    size_summary = df.groupby('component_id')['component_size'].first()
    print(f"     Min: {size_summary.min()}")
    print(f"     Max: {size_summary.max():,}")
    print(f"     Mean: {size_summary.mean():.1f}")
    print(f"     Median: {size_summary.median():.0f}")
    
    # Show sample of data
    print("\n5. Sample of dataframe:")
    print(df.head(10))
    
    print("\n" + "="*60)
    print("✓ JSON loading test PASSED")
    print("="*60)
    
    return df


if __name__ == '__main__':
    test_json_loading()
