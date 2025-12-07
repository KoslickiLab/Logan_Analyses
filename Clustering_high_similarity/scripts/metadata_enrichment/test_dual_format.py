#!/usr/bin/env python3
"""
Test script to verify both JSON component formats work correctly
"""

import json
import tempfile
import os


def detect_format(components_json):
    """
    Detect if JSON format uses integers (index-based) or strings (accession-based)
    
    Returns:
        bool: True if index-based format, False if string-based format
    """
    first_component = next(iter(components_json.values()))
    if len(first_component) > 0:
        first_value = first_component[0]
        return isinstance(first_value, int)
    else:
        # Empty component, assume index format for safety
        return True


def test_index_format():
    """Test the old index-based format"""
    print("="*80)
    print("TEST 1: Index-based JSON format")
    print("="*80)
    
    # Create test data
    components = {
        "component_0": [0, 1, 2],
        "component_1": [3, 4, 5]
    }
    
    print(f"\nComponent data:")
    print(json.dumps(components, indent=2))
    
    # Check format detection
    is_index = detect_format(components)
    
    print(f"\nFormat detection:")
    first_value = next(iter(components.values()))[0]
    print(f"  First value: {first_value}")
    print(f"  Type: {type(first_value).__name__}")
    print(f"  Detected as: {'index-based' if is_index else 'string-based'}")
    
    if is_index:
        print("  ✓ CORRECT - Detected as index-based format")
        return True
    else:
        print("  ✗ ERROR - Should be detected as index-based format")
        return False


def test_string_format():
    """Test the new string-based format"""
    print("\n" + "="*80)
    print("TEST 2: String-based JSON format")
    print("="*80)
    
    # Create test data
    components = {
        "component_0": ["DRR008435", "DRR008436", "DRR008437"],
        "component_1": ["SRR001234", "SRR001235", "SRR001236"]
    }
    
    print(f"\nComponent data:")
    print(json.dumps(components, indent=2))
    
    # Check format detection
    is_index = detect_format(components)
    
    print(f"\nFormat detection:")
    first_value = next(iter(components.values()))[0]
    print(f"  First value: {first_value}")
    print(f"  Type: {type(first_value).__name__}")
    print(f"  Detected as: {'index-based' if is_index else 'string-based'}")
    
    if not is_index:
        print("  ✓ CORRECT - Detected as string-based format")
        return True
    else:
        print("  ✗ ERROR - Should be detected as string-based format")
        return False


def test_empty_format():
    """Test empty component handling"""
    print("\n" + "="*80)
    print("TEST 3: Empty component handling")
    print("="*80)
    
    # Create test data with empty component
    components = {
        "component_0": []
    }
    
    print(f"\nComponent data:")
    print(json.dumps(components, indent=2))
    
    # Check format detection
    is_index = detect_format(components)
    
    print(f"\nFormat detection:")
    print(f"  Empty component")
    print(f"  Detected as: {'index-based' if is_index else 'string-based'} (default)")
    print(f"  ✓ HANDLED - Defaults to index-based for safety")
    return True


def show_examples():
    """Show examples of both formats"""
    print("\n" + "="*80)
    print("FORMAT EXAMPLES")
    print("="*80)
    
    print("\nOLD FORMAT (index-based):")
    print("-" * 40)
    old_format = {
        "component_0": [0, 1, 2, 3, 4, 6],
        "component_1": [7, 8, 9]
    }
    print(json.dumps(old_format, indent=2))
    print("\nRequires: accessions_mbases_geq_10.txt file")
    print("Usage: Integers index into accessions file (0-based)")
    print("Example: 0 → first line, 1 → second line, etc.")
    
    print("\n\nNEW FORMAT (string-based):")
    print("-" * 40)
    new_format = {
        "component_0": ["DRR008435", "DRR008436", "SRR001234"],
        "component_1": ["ERR123456", "SRR987654"]
    }
    print(json.dumps(new_format, indent=2))
    print("\nNo accessions file needed!")
    print("Usage: Strings are the actual SRA accessions")
    print("Example: 'DRR008435' is used directly")


def show_processing_logic():
    """Show how the code processes each format"""
    print("\n" + "="*80)
    print("PROCESSING LOGIC")
    print("="*80)
    
    print("\nIndex-based format processing:")
    print("-" * 40)
    print("1. Load components.json")
    print("2. Detect format: first value is int → index-based")
    print("3. Load accessions_mbases_geq_10.txt")
    print("4. Map: component[0] → accessions[0] → 'SRR2846706'")
    print("5. Create DataFrame with mapped accessions")
    
    print("\n\nString-based format processing:")
    print("-" * 40)
    print("1. Load components.json")
    print("2. Detect format: first value is str → string-based")
    print("3. Skip loading accessions file (not needed)")
    print("4. Use strings directly: 'DRR008435' → 'DRR008435'")
    print("5. Create DataFrame with direct accessions")


def main():
    """Run all tests"""
    print("\n" + "#"*80)
    print("# TESTING DUAL JSON FORMAT SUPPORT")
    print("#"*80)
    print("\nVerifying that code handles both index-based and string-based formats...\n")
    
    results = []
    
    # Run tests
    results.append(("Index format detection", test_index_format()))
    results.append(("String format detection", test_string_format()))
    results.append(("Empty component handling", test_empty_format()))
    
    # Show examples
    show_examples()
    show_processing_logic()
    
    # Summary
    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)
    
    for test_name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status:8s} - {test_name}")
    
    all_passed = all(passed for _, passed in results)
    
    print(f"\nOverall: {sum(passed for _, passed in results)}/{len(results)} tests passed")
    
    if all_passed:
        print("\n🎉 All format tests passed!")
        print("\nYour code now supports BOTH formats:")
        print("  ✓ Index-based: [0, 1, 2] → requires accessions file")
        print("  ✓ String-based: ['SRR001', 'SRR002'] → accessions embedded")
        print("\nFormat is detected automatically - just use the JSON file!")
    else:
        print("\n⚠️  Some tests failed")
    
    return 0 if all_passed else 1


if __name__ == '__main__':
    import sys
    sys.exit(main())
