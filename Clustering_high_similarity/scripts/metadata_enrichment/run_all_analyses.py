#!/usr/bin/env python3
"""
Master Analysis Pipeline

Runs the complete metadata enrichment analysis workflow:
1. Quick insights summary
2. Comprehensive metadata analysis  
3. Statistical enrichment testing
4. Interactive exploration and outlier detection
5. Example queries

Run this to execute the full analysis suite.
"""

import sys
import subprocess
from pathlib import Path
import argparse


def check_files_exist(component_file='components.json',
                     accession_file='accessions_mbases_geq_10.txt',
                     db_file='metagenome_metadata_with_geo.duckdb'):
    """Check that required input files exist"""
    # Always required
    required_files = [component_file, db_file]
    
    # Accessions file is optional (only needed for index-based JSON format)
    # We'll check if it exists if it's provided, but won't fail if it's missing
    # (the actual code will handle the error gracefully if needed)
    
    missing = []
    for f in required_files:
        if not Path(f).exists():
            missing.append(f)
    
    if missing:
        print("ERROR: Missing required input files:")
        for f in missing:
            print(f"  - {f}")
        print("\nPlease ensure these files are in the current directory.")
        return False
    
    # Check accessions file but only warn if missing (might not be needed)
    if not Path(accession_file).exists():
        print(f"Note: Accessions file '{accession_file}' not found.")
        print("      This is OK if your components JSON uses string accessions directly.")
        print()
    
    return True


def run_script(script_name, description, extra_args=None):
    """Run a Python script and handle errors"""
    print("\n" + "="*80)
    print(f"RUNNING: {description}")
    print("="*80)
    print(f"Script: {script_name}\n")
    
    cmd = ['python', script_name]
    if extra_args:
        cmd.extend(extra_args)
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=False,
            text=True
        )
        print(f"\n✓ {description} completed successfully")
        return True
    
    except subprocess.CalledProcessError as e:
        print(f"\n✗ {description} failed with error code {e.returncode}")
        print("See error output above for details")
        return False
    
    except Exception as e:
        print(f"\n✗ Unexpected error running {description}:")
        print(str(e))
        return False


def main(component_file='components.json',
         accession_file='accessions_mbases_geq_10.txt',
         db_file='metagenome_metadata_with_geo.duckdb',
         output_dir='./'):
    """
    Run the complete analysis pipeline
    
    Args:
        component_file: Path to component membership JSON
        accession_file: Path to accessions list
        db_file: Path to metadata database
        output_dir: Directory for output files
    """
    print("="*80)
    print("COMPONENT METADATA ANALYSIS - MASTER PIPELINE")
    print("="*80)
    print("\nThis will run the complete analysis suite.")
    print("Estimated time: 5-30 minutes depending on data size")
    print()
    
    # Check files
    print("Checking for required input files...")
    if not check_files_exist(component_file, accession_file, db_file):
        sys.exit(1)
    print("✓ All required files found\n")
    
    # Build common arguments
    file_args = [
        '--components', component_file,
        '--accessions', accession_file,
        '--database', db_file
    ]
    
    output_args = ['--output-dir', output_dir]
    
    # For scripts that read from merged CSV
    merged_csv = str(Path(output_dir) / 'merged_components_metadata.csv')
    merged_args = ['--input', merged_csv]
    
    # Define analysis pipeline
    pipeline = [
        ('quick_insights.py', 
         'Quick Insights Summary',
         file_args),
        
        ('analyze_component_metadata.py', 
         'Comprehensive Metadata Analysis',
         file_args + output_args),
        
        ('statistical_enrichment_analysis.py', 
         'Statistical Enrichment Testing',
         merged_args + output_args),
        
        ('explore_components_interactive.py', 
         'Interactive Exploration & Outlier Detection',
         merged_args + output_args),
        
        ('example_queries.py', 
         'Example Queries',
         merged_args),
    ]
    
    # Run pipeline
    results = []
    
    for i, (script, description, extra_args) in enumerate(pipeline, 1):
        print(f"\n{'#'*80}")
        print(f"# STEP {i}/{len(pipeline)}")
        print(f"{'#'*80}")
        
        success = run_script(script, description, extra_args)
        results.append((description, success))
        
        if not success:
            print(f"\nWARNING: {description} did not complete successfully")
            response = input("Continue with remaining steps? (y/n): ")
            if response.lower() != 'y':
                print("Analysis pipeline stopped by user")
                break
    
    # Summary
    print("\n" + "="*80)
    print("PIPELINE SUMMARY")
    print("="*80)
    
    for description, success in results:
        status = "✓" if success else "✗"
        print(f"{status} {description}")
    
    # List output files
    print("\n" + "="*80)
    print("OUTPUT FILES GENERATED")
    print("="*80)
    
    output_path = Path(output_dir)
    
    output_files = [
        ('component_metadata_summary.csv', 'Summary statistics for each component'),
        ('merged_components_metadata.csv', 'Full merged dataset'),
        ('statistical_enrichment_results.csv', 'Statistical test results'),
        ('pairwise_field_interactions.csv', 'Field combination analysis'),
        ('geographic_outliers.csv', 'Components with geographic diversity'),
        ('geographic_distances.csv', 'Geographic distance analysis with lat/lon'),
        ('bioproject_mixtures.csv', 'Cross-study duplications'),
        ('temporal_patterns.csv', 'Temporal clustering patterns'),
    ]
    
    for filename, description in output_files:
        filepath = output_path / filename
        if filepath.exists():
            size = filepath.stat().st_size / (1024*1024)  # MB
            print(f"✓ {filename:40s} ({size:>6.1f} MB) - {description}")
        else:
            print(f"✗ {filename:40s} - Not generated")
    
    # Next steps
    print("\n" + "="*80)
    print("NEXT STEPS")
    print("="*80)
    print("""
1. Review ANALYSIS_GUIDE.md for interpretation guidance

2. Open the CSV files in a spreadsheet program or pandas:
   - Start with component_metadata_summary.csv for overview
   - Check geographic_outliers.csv for surprising patterns
   - Review statistical_enrichment_results.csv for significant fields

3. For interactive investigation, use Python:
   
   from example_queries import *
   df = load_data()
   investigate_component(df, component_id=YOUR_ID)

4. Create visualizations or custom analyses based on your findings

5. Document interesting patterns for your research
    """)
    
    print("="*80)
    print("Analysis pipeline complete!")
    print("="*80)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run the complete component metadata analysis pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--components', default='components.json',
                       help='Component membership file (JSON format)')
    parser.add_argument('--accessions', default='accessions_mbases_geq_10.txt',
                       help='Accessions list file (only needed for index-based JSON format)')
    parser.add_argument('--database', default='metagenome_metadata_with_geo.duckdb',
                       help='Metadata database file')
    parser.add_argument('--output-dir', default='./',
                       help='Output directory for all results')
    
    args = parser.parse_args()
    
    main(
        component_file=args.components,
        accession_file=args.accessions,
        db_file=args.database,
        output_dir=args.output_dir
    )
