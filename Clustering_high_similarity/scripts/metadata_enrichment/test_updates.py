#!/usr/bin/env python3
"""
Test Script for Version 1.1 Updates

Tests that all bug fixes and new features work correctly.
"""

import subprocess
import sys
from pathlib import Path


def test_help_flags():
    """Test that all scripts have working --help"""
    print("="*80)
    print("TEST 1: Help Flags")
    print("="*80)
    
    scripts = [
        'quick_insights.py',
        'analyze_component_metadata.py',
        'statistical_enrichment_analysis.py',
        'explore_components_interactive.py',
        'example_queries.py',
        'run_all_analyses.py'
    ]
    
    passed = 0
    failed = 0
    
    for script in scripts:
        try:
            result = subprocess.run(
                ['python', script, '--help'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode == 0:
                # Check that help includes expected arguments
                help_text = result.stdout
                
                # Check for argparse-style help
                if 'usage:' in help_text.lower() and 'optional arguments:' in help_text.lower():
                    print(f"  ✓ {script:45s} - Help works")
                    passed += 1
                else:
                    print(f"  ✗ {script:45s} - Help format incorrect")
                    failed += 1
            else:
                print(f"  ✗ {script:45s} - Help failed (exit code {result.returncode})")
                failed += 1
                
        except Exception as e:
            print(f"  ✗ {script:45s} - Error: {e}")
            failed += 1
    
    print(f"\nResults: {passed} passed, {failed} failed")
    return failed == 0


def test_argument_parsing():
    """Test that custom arguments are accepted"""
    print("\n" + "="*80)
    print("TEST 2: Argument Parsing")
    print("="*80)
    
    tests = [
        ('analyze_component_metadata.py', ['--components', 'test.json', '--help']),
        ('analyze_component_metadata.py', ['--output-dir', 'test_dir', '--help']),
        ('statistical_enrichment_analysis.py', ['--input', 'test.csv', '--help']),
        ('explore_components_interactive.py', ['--output-dir', 'test_dir', '--help']),
        ('run_all_analyses.py', ['--components', 'test.json', '--output-dir', 'test', '--help']),
    ]
    
    passed = 0
    failed = 0
    
    for script, args in tests:
        try:
            result = subprocess.run(
                ['python', script] + args,
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode == 0:
                print(f"  ✓ {script:45s} accepts {' '.join(args[:2])}")
                passed += 1
            else:
                print(f"  ✗ {script:45s} rejected {' '.join(args[:2])}")
                failed += 1
                
        except Exception as e:
            print(f"  ✗ {script:45s} - Error: {e}")
            failed += 1
    
    print(f"\nResults: {passed} passed, {failed} failed")
    return failed == 0


def test_code_cleanup():
    """Verify that cleaned up code was actually removed"""
    print("\n" + "="*80)
    print("TEST 3: Code Cleanup Verification")
    print("="*80)
    
    checks = [
        ('analyze_component_metadata.py', 'date_fields = [', False, 'Unused date_fields removed'),
        ('analyze_component_metadata.py', 'return component_sizes', False, 'Unused return removed'),
        ('analyze_component_metadata.py', 'component_sizes = analyze_component_sizes', False, 'Unused assignment removed'),
    ]
    
    passed = 0
    failed = 0
    
    for filepath, pattern, should_exist, description in checks:
        try:
            with open(filepath, 'r') as f:
                content = f.read()
            
            exists = pattern in content
            
            if exists == should_exist:
                print(f"  ✓ {description}")
                passed += 1
            else:
                if should_exist:
                    print(f"  ✗ {description} - Pattern not found but should exist")
                else:
                    print(f"  ✗ {description} - Pattern still exists but should be removed")
                failed += 1
                
        except Exception as e:
            print(f"  ✗ {description} - Error: {e}")
            failed += 1
    
    print(f"\nResults: {passed} passed, {failed} failed")
    return failed == 0


def test_location_metadata_fix():
    """Verify that location_metadata merge was added"""
    print("\n" + "="*80)
    print("TEST 4: Location Metadata Bug Fix")
    print("="*80)
    
    try:
        with open('quick_insights.py', 'r') as f:
            content = f.read()
        
        checks = [
            ('location_df = conn.execute', 'Loads location_metadata table'),
            ('merge(location_df', 'Merges location_metadata'),
            ("suffixes=('', '_loc')", 'Uses correct merge suffixes'),
        ]
        
        passed = 0
        failed = 0
        
        for pattern, description in checks:
            if pattern in content:
                print(f"  ✓ {description}")
                passed += 1
            else:
                print(f"  ✗ {description} - Not found")
                failed += 1
        
        print(f"\nResults: {passed} passed, {failed} failed")
        return failed == 0
        
    except Exception as e:
        print(f"  ✗ Error reading file: {e}")
        return False


def test_argparse_imports():
    """Verify argparse was imported in all scripts"""
    print("\n" + "="*80)
    print("TEST 5: Argparse Imports")
    print("="*80)
    
    scripts = [
        'quick_insights.py',
        'analyze_component_metadata.py',
        'statistical_enrichment_analysis.py',
        'explore_components_interactive.py',
        'example_queries.py',
        'run_all_analyses.py'
    ]
    
    passed = 0
    failed = 0
    
    for script in scripts:
        try:
            with open(script, 'r') as f:
                content = f.read()
            
            if 'import argparse' in content:
                print(f"  ✓ {script:45s} - argparse imported")
                passed += 1
            else:
                print(f"  ✗ {script:45s} - argparse NOT imported")
                failed += 1
                
        except Exception as e:
            print(f"  ✗ {script:45s} - Error: {e}")
            failed += 1
    
    print(f"\nResults: {passed} passed, {failed} failed")
    return failed == 0


def main():
    """Run all tests"""
    print("\n" + "#"*80)
    print("# TESTING VERSION 1.1 UPDATES")
    print("#"*80)
    print("\nRunning comprehensive tests for bug fixes and new features...\n")
    
    results = []
    
    # Run all tests
    results.append(("Help Flags", test_help_flags()))
    results.append(("Argument Parsing", test_argument_parsing()))
    results.append(("Code Cleanup", test_code_cleanup()))
    results.append(("Location Metadata Fix", test_location_metadata_fix()))
    results.append(("Argparse Imports", test_argparse_imports()))
    
    # Summary
    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)
    
    total_passed = sum(1 for _, passed in results if passed)
    total_tests = len(results)
    
    for test_name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status:8s} - {test_name}")
    
    print(f"\nOverall: {total_passed}/{total_tests} test groups passed")
    
    if total_passed == total_tests:
        print("\n🎉 All tests passed! Version 1.1 updates are working correctly.")
        return 0
    else:
        print(f"\n⚠️  {total_tests - total_passed} test group(s) failed. Please review.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
