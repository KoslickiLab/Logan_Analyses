# Dark Matter Index (DMI) Calculator

A toolkit for computing the Dark Matter Index — a metric that quantifies the fraction of k-mer hashes in a metagenomic sample that do NOT appear in any reference database, indicating potentially novel or uncharacterized genetic material.

## Overview

The Dark Matter Index is defined as:

```
DMI = |sample_hashes - reference_hashes| / |sample_hashes|
```

- **DMI = 0.0**: All hashes are explained by reference databases
- **DMI = 1.0**: No hashes match any reference (100% "dark matter")

High DMI values indicate samples with significant uncharacterized genetic content — potential "microbial dark matter" not yet represented in public databases.

## Scripts Overview

| Script | Purpose | Audience |
|--------|---------|----------|
| `compute_dmi.py` | Compute DMI for individual .sig files | **End users** |
| `batch_compute_dmi.py` | Batch process SRA samples from DuckDB | Internal analysis |
| `build_reference_hashes.py` | Build the reference hash set | One-time setup |

## Reference Databases Included

The reference hash set combines unique k-mers (k=31, FracMinHash scale=1000) from:

1. **GTDB** — Genome Taxonomy Database
2. **GenBank WGS** — Whole Genome Shotgun sequences
3. **GenBank Genomes** — Complete genome assemblies
4. **AllTheBacteria** — Comprehensive bacterial genome collection
5. **SILVA** — Ribosomal RNA database
6. **GenBank TLS** — Targeted Locus Study sequences
7. **GenBank TSA** — Transcriptome Shotgun Assembly
8. **Logan Obelisks** — Novel RNA elements
9. **Logan Plasmids** — Plasmid sequences
10. **Serratus Viruses** — Viral sequences from Serratus project

Combined: **~20 billion unique hashes** (~160 GB binary file)

## Installation

### Requirements

```bash
pip install numpy pandas duckdb pyarrow
```

### File Structure

```
dark_matter_index/
├── README.md                    # This file
├── DB_info.json                 # Database configuration
├── compute_dmi.py               # User-facing DMI calculator
├── dmi_calculator.py            # Core DMI computation module
├── build_reference_hashes.py    # Script to build reference set
├── batch_compute_dmi.py         # Batch processing for SRA samples
└── demo_dmi.py                  # Demonstration script
```

## Quick Start

### 1. Build the Reference Hash Set (One-time)

For systems with 250GB+ RAM (recommended):

```bash
python build_reference_hashes.py \
    DB_info.json \
    /scratch/dmk333_new/dmi_reference/ \
    --ksize 31 \
    --method inmemory
```

For systems with less RAM, use chunked external merge sort:

```bash
python build_reference_hashes.py \
    DB_info.json \
    /scratch/dmk333_new/dmi_reference/ \
    --ksize 31 \
    --method chunked \
    --chunk-size 1.0
```

**Expected output:**
- `reference_hashes_k31.bin` (~160 GB)
- `reference_metadata.json` (statistics)

**Time:** 30-60 minutes with inmemory method on high-RAM system

**Note:** All source databases are connected in **read-only mode** to prevent accidental modifications. Individual databases are already deduplicated; only the union requires deduplication.

### 2. Compute DMI for Your Sample (User-Facing)

```bash
# Basic usage with a single signature file
python compute_dmi.py \
    --reference /path/to/reference_hashes_k31.bin \
    --input sample.sig \
    --ksize 31

# With gzipped signature
python compute_dmi.py -r reference_hashes_k31.bin -i sample.sig.gz -k 31

# Process multiple files and save results to JSON
python compute_dmi.py \
    -r reference_hashes_k31.bin \
    -i *.sig \
    -k 31 \
    -o results.json

# Quiet mode (just print DMI values)
python compute_dmi.py -r reference_hashes_k31.bin -i sample.sig -k 31 --quiet
```

**Example output:**
```
============================================================
Sample: SRR1234567
============================================================
Dark Matter Index (DMI): 0.3421
Dark Matter:             34.21%

Total hashes:            158,432
Mapped to reference:     104,234 (65.79%)
Unmapped (dark matter):  54,198 (34.21%)
============================================================
```

### 3. Batch Process SRA Samples (Internal Analysis)

```bash
python batch_compute_dmi.py \
    --reference /scratch/dmk333_new/dmi_reference/reference_hashes_k31.bin \
    --input /path/to/filtered_data.parquet \
    --output /path/to/output_with_dmi.parquet \
    --database /scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db \
    --ksize 31 \
    --workers 128 \
    --chunk-size 500
```

**Input parquet must have:** `accession` column  
**Output parquet adds:** `dmi`, `total_hashes_dmi`, `unmapped_hashes`, `mapped_hashes`

## Technical Details

### Data Structure

The reference hash set is stored as a **sorted array of 64-bit unsigned integers**:

```
[hash_1, hash_2, hash_3, ..., hash_20B]  (sorted ascending)
```

This enables:
- **Memory-mapped access**: No need to load 160GB into RAM
- **Binary search**: O(log n) membership test per hash
- **Vectorized operations**: NumPy's `searchsorted` for batch queries

### Membership Testing Algorithm

```python
# For each sample hash, binary search finds insertion point
indices = np.searchsorted(reference, sample_hashes)

# Check if the value at that index matches
matches = reference[indices] == sample_hashes

# Count non-matches
unmapped = len(sample_hashes) - np.sum(matches)
dmi = unmapped / len(sample_hashes)
```

### Performance Characteristics

| Operation | Complexity | Typical Time |
|-----------|------------|--------------|
| Load reference (mmap) | O(1) | < 1 second |
| Single hash lookup | O(log 20B) ≈ 34 ops | < 1 μs |
| Sample (100K hashes) | O(m log n) | 10-50 ms |
| Batch (1M samples, 128 workers) | Parallel | 1-2 hours |

### Memory Usage

| Mode | RAM Required |
|------|--------------|
| Memory-mapped (default) | ~100 MB base (page cache grows) |
| Full preload | ~160 GB |

### Recommended System Configuration

For batch processing millions of samples:

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| RAM | 32 GB | 256 GB+ |
| CPU Cores | 8 | 64-256 |
| Storage | SSD, 200GB free | NVMe SSD |
| Workers | cores/2 | cores * 0.5-0.75 |

**Note:** With 4TB RAM and 256 cores, use `--workers 128` for optimal throughput.

## Interpreting DMI Values

### Expected Ranges by Environment

Based on preliminary analysis:

| Environment | Typical DMI | Interpretation |
|-------------|-------------|----------------|
| Human gut | 0.05-0.20 | Well-characterized |
| Human feces | 0.10-0.25 | Mostly characterized |
| Soil | **0.40-0.80** | High dark matter! |
| Marine | 0.15-0.40 | Moderate dark matter |
| Wastewater | 0.20-0.40 | Mixed sources |

### Caveats

- **Low sequencing depth** can inflate DMI (rare k-mers less likely to match)
- **Contamination/artifacts** can increase DMI (non-biological sequences)
- **Novel taxa** genuinely have high DMI (the intended signal!)

Always cross-reference with:
- Sequencing depth / read count
- Quality control metrics
- Known sample characteristics

## Advanced Usage

### Build Methods

Three methods are available for building the reference hash set:

```bash
# In-memory (FASTEST, requires ~250GB RAM) - DEFAULT
python build_reference_hashes.py DB_info.json output/ --ksize 31 --method inmemory

# Chunked external merge sort (lower memory, slower)
python build_reference_hashes.py DB_info.json output/ --ksize 31 --method chunked --chunk-size 1.0

# DuckDB in-database (alternative approach)
python build_reference_hashes.py DB_info.json output/ --ksize 31 --method duckdb
```

### Resumable Batch Processing

```bash
# Start processing
python batch_compute_dmi.py ... --ksize 31 --output results.parquet

# If interrupted, resume:
python batch_compute_dmi.py ... --ksize 31 --resume-checkpoint checkpoint.parquet
```

### Custom Reference Sets

To use a subset of databases, modify `DB_info.json`:

```json
{
  "databases": [
    {
      "name": "GTDB",
      "set_id": "right",  // Include this database
      ...
    },
    {
      "name": "Some Other DB",
      "set_id": "exclude",  // Skip this database
      ...
    }
  ]
}
```

### Low-Memory Build

For systems with limited RAM, use the chunked method:

```bash
python build_reference_hashes.py DB_info.json output/ --ksize 31 --method chunked --chunk-size 0.5
```

## Output Schema

The output parquet file includes all original columns plus:

| Column | Type | Description |
|--------|------|-------------|
| `dmi` | float | Dark Matter Index (0.0-1.0) |
| `total_hashes_dmi` | int | Total hashes in sample |
| `unmapped_hashes` | int | Hashes not in reference |
| `mapped_hashes` | int | Hashes found in reference |

## Troubleshooting

### "Reference file not found"

Build the reference first:
```bash
python build_reference_hashes.py DB_info.json /path/to/output/
```

### "Memory error" during build

Reduce chunk size:
```bash
python build_reference_hashes.py DB_info.json output/ --chunk-size 0.25
```

### Slow performance

- Increase workers: `--workers 64`
- Use SSD for reference file
- Preload reference if RAM available (edit code: `preload=True`)

### DMI is NaN

Sample has no hashes. Check:
- Is the accession correct?
- Does the sample exist in the database?

## Citation

If you use this tool, please cite:

```
[Koslicki lab DMI tool - citation pending]
```

## License

[License information]

## Contact

David Koslicki Lab  
Penn State University
