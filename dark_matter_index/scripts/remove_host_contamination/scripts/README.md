# Script order
1. `get_T2T.sh`
2. `get_HPRC_Y1.sh`
3. `get_HPRC_Y2.sh`
4. `sketch_human_manysketch.sh`
5. `run_extract_human_hashes.sh`
6. `run_compute_dmi_host_filtered.sh`


# Host Contamination Removal for DMI Calculation

This directory contains tools to compute the Dark Matter Index (DMI) while excluding host-derived hashes from the calculation. This is essential for human-associated metagenome samples where host contamination would artificially inflate the mapped hash count.

## Overview

When computing DMI on human metagenomes (gut, skin, oral, etc.), a significant fraction of hashes may come from the host rather than the microbial community. If not removed, these host hashes get counted as "mapped" (since human genomes are in reference databases), which artificially **decreases** the DMI and underestimates the true microbial dark matter.

### The Modified DMI Calculation

```
Standard DMI:
    DMI = |sample_hashes - reference_hashes| / |sample_hashes|

Host-Filtered DMI:
    filtered_hashes = sample_hashes - host_hashes
    DMI = |filtered_hashes - reference_hashes| / |filtered_hashes|
```

## Quick Start

### Step 1: Extract Human Genome Hashes

First, extract hashes from your human genome sketches into a DuckDB database:

```bash
python extract_human_hashes.py \
    --input /path/to/human_genome_sketches_all.sig.zip \
    --output /path/to/human_hashes_k31.db \
    --ksize 31 \
    --scaled 1000
```

This only needs to be done once. The output database contains ~200-500 million unique human hashes (depending on how many assemblies you included).

### Step 2: Compute DMI with Host Removal

```bash
python compute_dmi_host_filtered.py \
    --database /path/to/dmi_database.db \
    --host-hashes /path/to/human_hashes_k31.db \
    --input /path/to/filtered_data.parquet \
    --output /path/to/filtered_data_with_dmi_host_removed.parquet \
    --organisms "human gut metagenome" "human feces metagenome"
```

## Scripts

| Script | Purpose |
|--------|---------|
| `extract_human_hashes.py` | Extract hashes from sourmash signatures to DuckDB |
| `compute_dmi_host_filtered.py` | Compute DMI with host hash removal and organism filtering |
| `run_extract_human_hashes.sh` | Example script for hash extraction |
| `run_compute_dmi_host_filtered.sh` | Example script for DMI computation |

## Output Format

The output parquet file is compatible with existing analysis scripts and includes all original columns plus:

| Column | Type | Description |
|--------|------|-------------|
| `dmi` | float | Dark Matter Index (computed on non-host hashes) |
| `total_hashes_dmi` | int | Non-host hashes used for DMI calculation |
| `total_hashes_original` | int | Original hash count before host removal |
| `host_hashes` | int | Number of host hashes removed |
| `host_fraction` | float | Fraction of hashes that were host-derived |
| `mapped_hashes` | int | Non-host hashes found in reference |
| `unmapped_hashes` | int | Non-host hashes not in reference (dark matter) |

## Organism Filtering

The script allows you to filter samples by organism type. This is useful for focusing on specific sample types:

```bash
# List available organisms in your data
python compute_dmi_host_filtered.py \
    --database db.db \
    --host-hashes human.db \
    --input data.parquet \
    --output out.parquet \
    --list-organisms

# Filter for specific organisms
python compute_dmi_host_filtered.py \
    ... \
    --organisms "human gut metagenome" "human feces metagenome" "human skin metagenome"
```

### Common Human-Associated Metagenome Types

- `human gut metagenome`
- `human feces metagenome`
- `human skin metagenome`
- `human oral metagenome`
- `human saliva metagenome`
- `human metagenome`
- `human nasopharyngeal metagenome`
- `human vaginal metagenome`
- `human milk metagenome`
- `human lung metagenome`

## Human Genome Sources

The human genome sketches were created from:

1. **HPRC (Human Pangenome Reference Consortium)**
   - Year 1 assemblies
   - Year 2 assemblies
   
2. **T2T (Telomere-to-Telomere) CHM13**
   - Complete human genome assembly

3. **NCBI Human Reference Genomes**
   - GRCh38/hg38
   - Additional assemblies

See `scripts/get_HPRC_Y1.sh`, `scripts/get_HPRC_Y2.sh`, and `scripts/get_T2T.sh` for download scripts.

## Performance

### Hash Extraction (`extract_human_hashes.py`)

| Metric | Value |
|--------|-------|
| Input | 503 sketches (~1.2B total hashes) |
| Output | ~200-500M unique hashes |
| Time | 10-30 minutes |
| Memory | 32-128 GB |

### DMI Computation (`compute_dmi_host_filtered.py`)

| Samples | Chunked | Time (estimate) |
|---------|---------|-----------------|
| 10,000 | No | 5-15 min |
| 100,000 | Yes | 30-90 min |
| 500,000 | Yes | 2-6 hours |

Memory usage depends on database size and chunk settings.

## Interpreting Results

### Expected Host Contamination Levels

| Sample Type | Typical Host Fraction |
|-------------|----------------------|
| Human gut | 1-5% |
| Human feces | 0.5-3% |
| Human skin | 5-20% |
| Human oral | 5-15% |
| Human saliva | 2-10% |

### Impact on DMI

Host removal typically **increases** DMI because:
1. Host hashes are removed from the denominator
2. Host hashes were likely mapped (counted against dark matter)

Example:
- Before: 1000 hashes total, 300 unmapped → DMI = 0.30
- Host removal: 50 host hashes removed → 950 non-host hashes
- Of those 950: 280 unmapped (some host hashes were also mapped)
- After: DMI = 280/950 = 0.295

The change may be small for low-contamination samples but significant for high-contamination samples.

## Troubleshooting

### "No samples match the organism filter"

Check available organisms with `--list-organisms`. Organism names are case-insensitive.

### Memory issues

- Use `--chunked` mode with smaller `--chunk-size`
- Reduce `--memory` setting
- Process fewer samples at once

### Slow performance

- Increase `--threads` (up to available cores)
- Use SSD storage for databases
- Pre-index the database (see `create_dmi_database.py`)

### High host fraction (>50%)

This may indicate:
- Sequencing from host tissue rather than microbial community
- Sample quality issues
- Incorrect sample annotation

## Citation

If you use this tool, please cite:

```
[Koslicki lab DMI tool - citation pending]
```

## Contact

David Koslicki Lab  
Penn State University
