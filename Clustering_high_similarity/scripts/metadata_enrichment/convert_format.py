#!/usr/bin/env python3
"""
Component File Format Converter

Converts between CSV and JSON formats for component membership.
"""

import json
import pandas as pd
import argparse
from pathlib import Path


def csv_to_json(csv_file, json_file):
    """Convert CSV format to JSON format"""
    print(f"Converting {csv_file} to {json_file}...")
    
    # Read CSV
    df = pd.read_csv(csv_file)
    print(f"  Loaded {len(df):,} rows")
    
    # Group by component
    components = {}
    for comp_id in sorted(df['component_id'].unique()):
        sample_ids = df[df['component_id'] == comp_id]['sample_id'].tolist()
        components[f'component_{comp_id}'] = sample_ids
        print(f"  Component {comp_id}: {len(sample_ids)} samples")
    
    # Save as JSON
    with open(json_file, 'w') as f:
        json.dump(components, f, indent=2)
    
    print(f"✓ Saved {len(components)} components to {json_file}")


def json_to_csv(json_file, csv_file):
    """Convert JSON format to CSV format"""
    print(f"Converting {json_file} to {csv_file}...")
    
    # Load JSON
    with open(json_file, 'r') as f:
        components = json.load(f)
    
    print(f"  Loaded {len(components)} components")
    
    # Convert to dataframe
    rows = []
    for comp_name, sample_ids in components.items():
        # Extract component ID from name
        comp_id = int(comp_name.split('_')[1])
        comp_size = len(sample_ids)
        
        for sample_id in sample_ids:
            rows.append({
                'sample_id': sample_id,
                'component_id': comp_id,
                'component_size': comp_size
            })
        
        print(f"  Component {comp_id}: {comp_size} samples")
    
    df = pd.DataFrame(rows)
    df = df.sort_values(['component_id', 'sample_id'])
    df.to_csv(csv_file, index=False)
    
    print(f"✓ Saved {len(df):,} rows to {csv_file}")


def validate_json(json_file):
    """Validate JSON format"""
    print(f"Validating {json_file}...")
    
    errors = []
    warnings = []
    
    # Load JSON
    try:
        with open(json_file, 'r') as f:
            components = json.load(f)
    except json.JSONDecodeError as e:
        print(f"✗ Invalid JSON: {e}")
        return False
    
    # Check structure
    if not isinstance(components, dict):
        errors.append("Root element must be a dictionary")
    
    # Check each component
    for comp_name, sample_ids in components.items():
        # Check name format
        if not comp_name.startswith('component_') and not comp_name.startswith('community_'):
            warnings.append(f"Unusual component name: {comp_name}")
        
        try:
            comp_id = int(comp_name.split('_')[1])
        except (ValueError, IndexError):
            errors.append(f"Cannot extract integer ID from: {comp_name}")
        
        # Check sample IDs
        if not isinstance(sample_ids, list):
            errors.append(f"{comp_name}: value must be a list")
            continue
        
        if len(sample_ids) == 0:
            warnings.append(f"{comp_name}: empty component")
        
        for sample_id in sample_ids:
            if not isinstance(sample_id, int):
                errors.append(f"{comp_name}: sample_id {sample_id} is not an integer")
                break
            if sample_id < 0:
                errors.append(f"{comp_name}: negative sample_id {sample_id}")
                break
    
    # Report results
    if errors:
        print(f"\n✗ Validation FAILED with {len(errors)} errors:")
        for error in errors:
            print(f"  - {error}")
        return False
    
    if warnings:
        print(f"\n⚠ Validation passed with {len(warnings)} warnings:")
        for warning in warnings:
            print(f"  - {warning}")
    else:
        print("\n✓ Validation PASSED - no errors")
    
    # Show summary
    print(f"\nSummary:")
    print(f"  Total components: {len(components)}")
    print(f"  Total samples: {sum(len(v) for v in components.values()):,}")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description='Convert between CSV and JSON component formats',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert CSV to JSON
  python convert_format.py csv2json component_membership.csv components.json
  
  # Convert JSON to CSV
  python convert_format.py json2csv components.json component_membership.csv
  
  # Validate JSON file
  python convert_format.py validate components.json
        """
    )
    
    parser.add_argument('command', choices=['csv2json', 'json2csv', 'validate'],
                       help='Conversion command')
    parser.add_argument('input_file', help='Input file')
    parser.add_argument('output_file', nargs='?', help='Output file (not needed for validate)')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not Path(args.input_file).exists():
        print(f"Error: Input file '{args.input_file}' not found")
        return 1
    
    # Execute command
    try:
        if args.command == 'csv2json':
            if not args.output_file:
                print("Error: output_file required for csv2json")
                return 1
            csv_to_json(args.input_file, args.output_file)
        
        elif args.command == 'json2csv':
            if not args.output_file:
                print("Error: output_file required for json2csv")
                return 1
            json_to_csv(args.input_file, args.output_file)
        
        elif args.command == 'validate':
            if not validate_json(args.input_file):
                return 1
        
        return 0
    
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    exit(main())
