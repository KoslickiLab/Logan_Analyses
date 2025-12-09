#!/usr/bin/env python3
"""
Setup Verification Script
=========================
Checks that all requirements are met before running the analysis.

Usage:
    python3 verify_setup.py
"""

import sys
from pathlib import Path

# Color codes for output
GREEN = '\033[92m'
RED = '\033[91m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
RESET = '\033[0m'

def print_status(message, status):
    """Print a status message with color."""
    if status == 'ok':
        print(f"{GREEN}✓{RESET} {message}")
        return True
    elif status == 'warning':
        print(f"{YELLOW}⚠{RESET} {message}")
        return True
    elif status == 'error':
        print(f"{RED}✗{RESET} {message}")
        return False
    else:
        print(f"{BLUE}ℹ{RESET} {message}")
        return True

def check_python_version():
    """Check Python version."""
    version = sys.version_info
    if version.major >= 3 and version.minor >= 8:
        return print_status(f"Python version: {version.major}.{version.minor}.{version.micro}", 'ok')
    else:
        return print_status(f"Python version: {version.major}.{version.minor}.{version.micro} (3.8+ required)", 'error')

def check_dependencies():
    """Check if all required Python packages are installed."""
    required = [
        'duckdb',
        'pandas',
        'numpy',
        'matplotlib',
        'seaborn',
        'scipy',
        'sklearn',
        'tqdm'
    ]
    
    all_ok = True
    for package in required:
        try:
            __import__(package)
            print_status(f"Package '{package}' installed", 'ok')
        except ImportError:
            print_status(f"Package '{package}' NOT installed", 'error')
            all_ok = False
    
    return all_ok

def check_databases():
    """Check if databases are accessible."""
    yacht_db = Path("/scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db")
    metadata_db = Path("/scratch/shared_data_new/Logan_yacht_data/metadata/aws_sra_metadata/metadata_geo_joined.duckdb")
    
    all_ok = True
    
    if yacht_db.exists():
        size = yacht_db.stat().st_size / (1024**3)  # GB
        print_status(f"YACHT database found ({size:.1f} GB)", 'ok')
    else:
        print_status(f"YACHT database NOT found at: {yacht_db}", 'error')
        all_ok = False
    
    if metadata_db.exists():
        size = metadata_db.stat().st_size / (1024**3)  # GB
        print_status(f"Metadata database found ({size:.1f} GB)", 'ok')
    else:
        print_status(f"Metadata database NOT found at: {metadata_db}", 'error')
        all_ok = False
    
    return all_ok

def check_database_access():
    """Check if we can actually query the databases."""
    try:
        import duckdb
        
        # Test metadata database
        metadata_db = "/scratch/shared_data_new/Logan_yacht_data/metadata/aws_sra_metadata/metadata_geo_joined.duckdb"
        conn = duckdb.connect(metadata_db, read_only=True)
        result = conn.execute("SELECT COUNT(*) FROM metadata_geo_joined WHERE assay_type='WGS'").fetchone()
        conn.close()
        print_status(f"Metadata database accessible ({result[0]:,} WGS samples)", 'ok')
        
        # Test YACHT database
        yacht_db = "/scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db"
        conn = duckdb.connect(yacht_db, read_only=True)
        result = conn.execute("SELECT COUNT(DISTINCT sample_id) FROM taxa_profiles.profiles").fetchone()
        conn.close()
        print_status(f"YACHT database accessible ({result[0]:,} samples)", 'ok')
        
        return True
    except Exception as e:
        print_status(f"Database access error: {e}", 'error')
        return False

def check_system_resources():
    """Check available system resources."""
    try:
        import os
        import psutil
        
        # CPU cores
        cpu_count = os.cpu_count()
        if cpu_count >= 64:
            print_status(f"CPU cores: {cpu_count} (excellent)", 'ok')
        elif cpu_count >= 32:
            print_status(f"CPU cores: {cpu_count} (good)", 'ok')
        else:
            print_status(f"CPU cores: {cpu_count} (may be slow, 64+ recommended)", 'warning')
        
        # RAM
        ram_gb = psutil.virtual_memory().total / (1024**3)
        if ram_gb >= 500:
            print_status(f"RAM: {ram_gb:.0f} GB (excellent)", 'ok')
        elif ram_gb >= 256:
            print_status(f"RAM: {ram_gb:.0f} GB (good)", 'ok')
        else:
            print_status(f"RAM: {ram_gb:.0f} GB (may be insufficient, 256+ GB recommended)", 'warning')
        
        # Disk space
        disk = psutil.disk_usage('.')
        free_gb = disk.free / (1024**3)
        if free_gb >= 50:
            print_status(f"Free disk space: {free_gb:.0f} GB (sufficient)", 'ok')
        elif free_gb >= 10:
            print_status(f"Free disk space: {free_gb:.0f} GB (minimal)", 'warning')
        else:
            print_status(f"Free disk space: {free_gb:.0f} GB (insufficient, 10+ GB required)", 'error')
            return False
        
        return True
    except ImportError:
        print_status("psutil not installed (system resource check skipped)", 'warning')
        return True
    except Exception as e:
        print_status(f"System resource check error: {e}", 'warning')
        return True

def check_scripts():
    """Check if analysis scripts are present."""
    scripts = [
        'hash_diversity_correlation.py',
        'hash_diversity_sensitivity.py',
        'run_analysis.sh'
    ]
    
    all_ok = True
    for script in scripts:
        if Path(script).exists():
            print_status(f"Script '{script}' found", 'ok')
        else:
            print_status(f"Script '{script}' NOT found", 'error')
            all_ok = False
    
    return all_ok

def main():
    """Run all verification checks."""
    print("\n" + "="*70)
    print("HASH-DIVERSITY CORRELATION ANALYSIS - SETUP VERIFICATION")
    print("="*70 + "\n")
    
    checks = [
        ("Python Version", check_python_version),
        ("Python Dependencies", check_dependencies),
        ("Analysis Scripts", check_scripts),
        ("Database Files", check_databases),
        ("Database Access", check_database_access),
        ("System Resources", check_system_resources),
    ]
    
    results = []
    for name, check_func in checks:
        print(f"\n{BLUE}Checking: {name}{RESET}")
        print("-" * 70)
        try:
            results.append(check_func())
        except Exception as e:
            print_status(f"Error during check: {e}", 'error')
            results.append(False)
    
    print("\n" + "="*70)
    print("VERIFICATION SUMMARY")
    print("="*70)
    
    if all(results):
        print(f"\n{GREEN}✓ All checks passed!{RESET}")
        print("You are ready to run the analysis.\n")
        print("Quick start:")
        print("  bash run_analysis.sh quick        # Test run with 1,000 samples")
        print("  bash run_analysis.sh medium       # Medium run with 50,000 samples")
        print("  bash run_analysis.sh sensitivity  # Sensitivity analysis")
        return 0
    else:
        print(f"\n{RED}✗ Some checks failed!{RESET}")
        print("Please resolve the issues above before running the analysis.\n")
        
        if not results[1]:  # Dependencies failed
            print("To install missing dependencies:")
            print("  pip install -r requirements.txt")
        
        return 1

if __name__ == "__main__":
    sys.exit(main())
