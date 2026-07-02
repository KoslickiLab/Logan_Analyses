# FracMinHash Sketch Diversity Over Time

Analysis code accompanying the paper: **"[Title]"**.

This directory contains scripts that quantify how the diversity of FracMinHash
k-mer sketches in the NCBI Sequence Read Archive (SRA) has evolved from 2008
to 2023, using the [Logan](https://github.com/IndexThePlanet/Logan) large-scale
sketch database.

## Overview

Each publicly available SRA metagenome in the Logan dataset is represented as a
FracMinHash sketch (k = 31, scale 1/1000).  Each 64-bit hash value in a sketch
corresponds to a unique 31-mer.  The analysis tracks when each distinct hash
first appears ("birth year") and last appears ("death year") across the ~4.8
million metagenomes in the corpus.  From these records we compute cumulative
sequence diversity growth over fifteen years, per-sample novelty (the fraction
of k-mers not previously observed in any earlier-released run), and
Kaplan-Meier survival curves for individual hashes.

Two analysis strata are produced:
- **Full corpus** — all SRA metagenomes with a release date
- **WGS Illumina metagenomics** — filtered to `librarysource=METAGENOMIC`,
  `platform=ILLUMINA`, `assay_type=WGS`, `libraryselection IN (RANDOM, RANDOM PCR)`,
  and `mbases > 500`

## Repository layout

```
sketches_over_time/
|-- scripts/          analysis code (see below)
|-- data/
|   |-- *.duckdb      intermediate databases (produced by step 01)
|   |-- csv/          exported aggregate tables (gzip CSV)
|   \-- plots/        output figures and figure_captions.md
|-- logs/             timestamped run logs
\-- paper/            manuscript files
```

## Scripts

| Script | Purpose |
|---|---|
| `01_build_intermediate_db.py` | Build the intermediate DuckDB database (heavy computation) |
| `02_plot_cumulative_hashes.py` | Figure 1: cumulative distinct hashes over time |
| `03_plot_new_hashes_per_sample.py` | Figures 2a–2b: new hashes per sample over time |
| `04_plot_survival.py` | Figures 3–5: Kaplan-Meier survival and birth/death heatmap |
| `compute_application0_stats.py` | Reproduce every statistic cited in the paper |
| `run_pipeline.sh` | Orchestrate steps 01–04 end-to-end |
| `run_WGS.sh` | Example invocation for the WGS-filtered stratum |

### Step 01 — build intermediate database

`01_build_intermediate_db.py` performs three sequential stages over the full
Logan signature table (~193 billion rows):

| Stage | Table produced | Description |
|---|---|---|
| 1 | `hash_life` | First and last observed year for every distinct hash |
| 2 | `sample_stats` | Per-sample total hashes, new hashes, and fraction new |
| 3 | `cumulative_hashes_by_year`, `yearly_sample_stats`, `hash_life_counts` | Lightweight aggregates used by all plotting scripts |

The script auto-names the output database from the flags used
(e.g. `intermediate_full.duckdb`, `intermediate_wgs_metagenomics_illumina.duckdb`).

Stages 1 and 2 each require a full table scan and take several hours on the
complete dataset.  Stage 3 is fast.

**Key options:**

```
--test [--num-test-samples N]   Smoke test on N samples (default 500, ~minutes)
--metadata-filter SQL           Extra AND clause on metadata columns
--output-db PATH                Override the auto-named output path
--stage {1,2,3}                 Run only a single stage
--resume                        Skip stages whose tables already exist
--threads N                     DuckDB thread count (default 768)
--memory-limit SIZE             DuckDB memory limit (default 2500GB)
--no-csv                        Skip CSV export of aggregate tables
```

### Steps 02–04 — plotting

Each plotting script reads the relevant aggregate tables from the intermediate
database and writes PNGs to `data/plots/`.  They accept `--db PATH` to target
a specific database, or auto-detect the most recent `intermediate_full.duckdb`.

### compute_application0_stats.py

Queries both the full and WGS-filtered intermediate databases and prints every
statistic cited in the paper's Application 0 section, with labelled output for
traceability.  Requires both `intermediate_full.duckdb` and
`intermediate_wgs_metagenomics_illumina.duckdb` to already exist.

## Running the pipeline

### Smoke test (~minutes)

```bash
conda activate logan
cd scripts/
bash run_pipeline.sh --test --num-test-samples 500
```

### Full run (expect several hours; use nohup or tmux)

```bash
conda activate logan
cd scripts/
nohup bash run_pipeline.sh > run_pipeline.log 2>&1 &
```

### WGS-filtered stratum

```bash
conda activate logan
cd scripts/
nohup bash run_pipeline.sh \
  --metadata-filter "librarysource='METAGENOMIC' AND libraryselection IN ('RANDOM', 'RANDOM PCR') AND mbases>500 AND assay_type='WGS' AND platform='ILLUMINA'" \
  --output-db ../data/intermediate_wgs_metagenomics_illumina.duckdb \
  > run_WGS.log 2>&1 &
```

### Resume after interruption

```bash
bash run_pipeline.sh --resume
```

### Reproduce paper statistics only (after databases are built)

```bash
conda activate logan
cd scripts/
python compute_application0_stats.py
```

## Dependencies

| Package | Role |
|---|---|
| `duckdb` | Analytical query engine for large-scale hash computations |
| `matplotlib` | Figure generation |
| `numpy` | Array operations and Kaplan-Meier computation |
| `pandas` | DataFrame manipulation |

Install via the `logan` conda environment.

## Data sources

The analysis reads from two pre-existing databases:

| Variable | Path | Description |
|---|---|---|
| `SRC_DB` | `/scratch/shared_data/Logan_yacht_data/processed_data/database_all.db` | Logan FracMinHash signature table |
| `META_DB` | `/scratch/shared_data/Logan_yacht_data/metadata/aws_sra_metadata/metadata_geo_joined_5M.duckdb` | SRA metadata with release dates |

These are read-only inputs; all outputs are written under `data/`.

## Output figures

See `data/plots/figure_captions.md` for full figure descriptions.

| File | Content |
|---|---|
| `*__cumulative_hashes.png` | Cumulative and annual new distinct hashes (log scale) |
| `*__new_hashes_per_sample_median.png` | Per-sample novelty over time, median summary |
| `*__new_hashes_per_sample_mean.png` | Per-sample novelty over time, mean summary |
| `*__survival_km_overall.png` | Kaplan-Meier survival curve for all hashes |
| `*__survival_km_by_birth.png` | KM curves stratified by birth year |
| `*__birth_death_heatmap.png` | 2-D heatmap of hash birth year vs. death year |
