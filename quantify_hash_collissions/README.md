# Quantifying FracMinHash Hash Collision Rate

## Overview

This project measures how often MurmurHash3 — the hash function used in
FracMinHash (FMH) sketching — produces **hash collisions**: two distinct
canonical k-mers that map to the same 64-bit hash value.  We ran FMH
sketching on a large collection of metagenomic assemblies from the Logan
dataset and counted (a) how many collisions occur within individual sketches
and (b) how many distinct hash values appear across all sketches combined.

---

## Background: FracMinHash sketching

FracMinHash is a k-mer sketching technique that represents a sequence by
the subset of its k-mers whose hash value falls below a threshold:

```
threshold = scale × 2^64
```

With `scale = 0.001`, roughly 0.1 % of all k-mers are retained.  Because
the threshold is fixed (not data-dependent), sketches from different samples
are directly comparable using Jaccard similarity.

**Canonical k-mers.**  For each k-mer, only the lexicographically smaller of
the k-mer and its reverse complement is hashed.  This makes the sketch
strand-independent (a sequence and its reverse complement yield the same
sketch).

**Hash function.**  MurmurHash3_x64_128 with seed 42, retaining the low 64
bits (`MurmurHash3_x64_128_low64`).

**Parameters used throughout this project:**

| Parameter | Value |
|-----------|-------|
| k-mer size | 31 |
| Scale | 0.001 |
| Hash seed | 42 |
| Canonical | yes |
| Skip ambiguous bases | yes |

---

## Data source: Logan

Logan is a large-scale public metagenomics resource.  Assembled unitig
sequences (compressed `.fa.zst` files) are available from S3:

```
https://s3.amazonaws.com/logan-pub/u/<accession>/<accession>.unitigs.fa.zst
```

The accession list used is `data/sorted_WGS_accessions.csv`, which contains
**171,032 whole-genome-sequencing (WGS) accessions** sorted from largest to
smallest by (approximate) assembly size. These were the largest metagenomes (using the metadata filter of a WGS assay type, metagenome sample type, random selectio method, and mbases>1000).  When running a subset for testing,
take rows from the **bottom** of the file (the smallest assemblies).

Pre-computed FMH sketches for the Logan dataset are stored in a DuckDB
database at:

```
/scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db
```

Relevant table: `sigs_dna.signature_mins`, columns `sample_id`, `min_hash`,
`ksize`.  Open this database in **read-only mode only**.

---

## Repository layout

```
quantify_hash_collissions/
├── data/
│   ├── sorted_WGS_accessions.csv        # 171,032 accession IDs, large→small
│   └── collision_results/               # output from run_fmh_collisions.sh
│       ├── summary.tsv                  # one row per accession (see below)
│       └── <acc>.collisions.tsv         # per-accession collision details
├── kmer-sketch/
│   ├── include/                         # FMH implementation (headers-only C++)
│   │   ├── Hash.hpp                     # MurmurHash3 + SplitMix64
│   │   ├── KmerScanner.hpp              # k-mer extraction + canonicalisation
│   │   ├── FastxReader.hpp              # FASTA/FASTQ reader
│   │   └── Sketches.hpp                 # FracMinHash class (and others)
│   ├── src/
│   │   └── fmh_collisions.cpp           # NEW: collision-counting executable
│   └── bin/
│       └── fmh_collisions               # compiled binary
└── scripts/
    ├── run_fmh_collisions.sh            # parallel download + process pipeline
    ├── count_distinct_hashes_parallel.py  # count distinct hashes in DuckDB
    └── count_distinct_hashes.py           # earlier single-threaded attempt (slow)
```

---

## Step 1 — The C++ collision detector (`fmh_collisions`)

**Source:** `kmer-sketch/src/fmh_collisions.cpp`  
**Binary:** `kmer-sketch/bin/fmh_collisions`

### What it does

For a given FASTA file the tool:

1. Scans every k-mer position in every sequence.
2. Computes the canonical form of each k-mer (lex-min of k-mer and reverse
   complement), skipping any k-mer containing a non-ACGT base.
3. Hashes the canonical k-mer with MurmurHash3.
4. If `hash ≤ threshold` (i.e. the k-mer would be selected into the FMH
   sketch), records the mapping `hash → set of distinct canonical k-mers`.
5. After scanning the whole file, writes two output files:
   - **Summary file** — key/value pairs including `sketch_size` (number of
     distinct hash values selected), `n_collision_hashes` (hashes with ≥ 2
     distinct k-mers), and `n_kmers_in_collisions`.
   - **Collisions file** — one TSV line per colliding hash:
     `hash_value \t kmer1 \t kmer2 \t ...`

### Build

```bash
cd kmer-sketch
make bin/fmh_collisions
```

### Usage

```bash
kmer-sketch/bin/fmh_collisions \
    --input    sequences.fa \
    --kmer     31 \
    --scale    0.001 \
    --seed     42 \
    --summary  out_summary.tsv \
    --collisions out_collisions.tsv
```

Always uses canonical k-mers and skips ambiguous bases (no flags needed).

---

## Step 2 — Parallel pipeline (`run_fmh_collisions.sh`)

**Script:** `scripts/run_fmh_collisions.sh`

### What it does

For each accession in the input list:

1. Downloads `<accession>.unitigs.fa.zst` from the Logan S3 bucket via
   `wget`.
2. Decompresses on the fly with `zstd -d`.
3. Runs `fmh_collisions` on the decompressed FASTA.
4. Parses the per-accession summary and emits one TSV row to stdout.
5. Deletes all temporary files.

All accessions are processed in parallel using **GNU parallel**.  The master
process writes a combined `summary.tsv` with one row per accession.

### Usage

```bash
# From the repo root:
bash scripts/run_fmh_collisions.sh \
    -N 10          \   # number of accessions (omit for all 171,032)
    -t 16          \   # parallel jobs
    -o data/collision_results \   # output directory
    -d /dev/shm        # temp dir for decompressed FASTA files
```

**Flags:**

| Flag | Meaning | Default |
|------|---------|---------|
| `-N` | Number of accessions to process (taken from the **bottom** of the list — smallest first, good for testing) | all |
| `-t` | Parallel jobs | 4 |
| `-o` | Output directory | `<repo>/results` |
| `-d` | Temporary file directory (use `/dev/shm` or a large scratch path to avoid filling system `/tmp`) | `/tmp` |

### Output files

`<output_dir>/summary.tsv` — tab-separated, one row per accession:

| Column | Description |
|--------|-------------|
| `acc` | Accession ID |
| `sketch_size` | Number of distinct hash values ≤ threshold |
| `n_collision_hashes` | Hashes with ≥ 2 distinct canonical k-mers |
| `n_kmers_in_collisions` | Total k-mers involved in collisions |
| `total_kmer_positions` | All k-mer positions scanned (counting duplicates) |
| `total_distinct_kmers_selected` | Distinct canonical k-mers selected into sketch |

`<output_dir>/<acc>.collisions.tsv` — collision detail for each accession
(empty header only if no collisions).

### Disk space / temp files

At any point there are at most `-t` decompressed FASTA files alive
simultaneously.  Each is deleted immediately after `fmh_collisions`
finishes.  If the script is killed mid-run, stale files can be cleaned up
with:

```bash
rm -f /dev/shm/SRR*.fa /dev/shm/ERR*.fa \
      /dev/shm/DRR*.fa /dev/shm/SRR*_sum.*.tsv
```

### Summing sketch sizes across all accessions

```bash
awk -F'\t' 'NR>1 { s += $2 } END { print s }' data/collision_results/summary.tsv
```

---

## Step 3 — Count distinct hashes across all sketches

**Script:** `scripts/count_distinct_hashes_parallel.py`

### Motivation

`sum(sketch_size)` over all accessions **over-counts** the distinct hashes
because the same hash value can (and does) appear in multiple sketches.
This script counts the number of **globally distinct** `min_hash` values
across all the accessions listed in `summary.tsv`, using the pre-computed
sketches in the Logan DuckDB.

### Why not a simple `COUNT(DISTINCT min_hash)`?

The `signature_mins` table contains ~193 billion rows for `ksize = 31`.
A single `COUNT(DISTINCT)` query attempted to hold all distinct hash values
in memory at once and was killed (OOM / excessive runtime).

### Algorithm

Two phases:

**Phase 1 (parallel).** T worker processes are launched, each responsible
for a disjoint batch of accessions.  Each worker:
- Opens its own read-only DuckDB connection (1 internal thread each).
- Queries all `min_hash` values for its accessions via a single JOIN.
- Fetches the result as a numpy `int64` array (zero-copy, no Python loop).
- Partitions the hashes by `uint64(hash) % B` into B binary temp files.

Because all workers read the **same DB file simultaneously**, the OS page
cache reads each physical page from SSD once and serves subsequent requests
from RAM — effective I/O equals one logical scan regardless of T.

**Phase 2 (sequential).** For each of the B buckets:
- Load all T workers' slice for that bucket (a few hundred MB at most).
- `numpy.unique()` to count distinct values.
- Delete the files.

Because equal hash values always land in the same bucket, the sum of
per-bucket distinct counts equals the global distinct count.

### Usage

```bash
python3 scripts/count_distinct_hashes_parallel.py \
    --db      /scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db \
    --summary data/collision_results/summary.tsv \
    -t 32          \   # worker processes
    -b 256         \   # hash-space buckets
    -d /dev/shm        # temp dir (needs ~total_hashes × 8 bytes free)
```

**Flags:**

| Flag | Meaning | Default |
|------|---------|---------|
| `-t` | Worker processes | `min(32, cpu_count)` |
| `-b` | Partitions for dedup phase | 256 |
| `-d` | Temp dir | `/dev/shm` |
| `-k` | k-mer size filter | 31 |

Prints a single integer (the distinct count) to stdout.

If interrupted, clean up with: `rm -f /dev/shm/b????_w????.bin`

---

## Results (as of April 2026)

The pipeline was run on **1,541 accessions** from `sorted_WGS_accessions.csv`.

### Collision statistics (from `data/collision_results/summary.tsv`)

| Metric | Value |
|--------|-------|
| Accessions processed | 1,541 |
| Total k-mer positions scanned | 4,663,627,680,883 (~4.7 trillion) |
| Total hashes selected (sum of sketch sizes, with cross-sketch duplicates) | 4,663,778,741 |
| Total distinct canonical k-mers selected | 4,663,778,742 |
| **Hash collisions detected** | **1** |
| K-mers involved in that collision | 2 |

The single observed collision was in accession `ERR1726685`:

```
Hash value:  13086129273328152
k-mer 1:     AGTAAAAATAAAGAATTTTCTGACACTTCAG
k-mer 2:     AGTGGCTTCTCTTCAGAAGCCATCTCTTTCA
```

These two distinct canonical 31-mers produced identical MurmurHash3 values.

### Distinct hash count (from `count_distinct_hashes_parallel.py`)

Across all 1,541 accessions, the number of **globally distinct** hash values
(i.e. the size of the union of all sketches) is:

> **~2.6 billion**

### Implied collision rate

Over ~4.66 trillion k-mer positions scanned, exactly **1 hash collision** was
observed among ~2.6 billion distinct hash values that passed the scale
threshold.

The correct space to reason about is not the full 2^64 hash space but the
region *below the FMH threshold*, since only hashes in that region are
retained.  With scale = 0.001 the effective space has size:

```
M = 2^64 x 0.001 ~= 1.84 x 10^16
```

The birthday-paradox formula for the expected number of collisions among N
items drawn uniformly from a space of size M is:

```
E[collisions] ~= N^2 / (2M)
```

With N = 2.6 x 10^9 distinct hash values and M = 1.84 x 10^16:

```
E[collisions] ~= (2.6 x 10^9)^2 / (2 x 1.84 x 10^16)
               = 6.76 x 10^18 / (3.68 x 10^16)
              ~= 184
```

Observing only 1 collision against an expectation of ~184 means the 2.6
billion distinct hash values are **not** 2.6 billion independent draws from
that range.  The most likely explanation is that many k-mers are shared
across samples: the distinct hash count reflects the size of the *union* of
all sketches, not the total number of independent hashing events.  The true
number of independent draws is much smaller, and the observed collision count
of 1 is consistent with a near-ideal hash function operating on that smaller
effective input.

---

## Reproducing the results

```bash
# 0. Build the C++ tool (only needed once)
cd kmer-sketch && make bin/fmh_collisions && cd ..

# 1. Run the collision pipeline (adjust -N and -t as needed)
bash scripts/run_fmh_collisions.sh \
    -t 32 \
    -o data/collision_results \
    -d /dev/shm

# 2. Aggregate results
awk -F'\t' 'NR>1 { sketch+=$2; col+=$3 } END {
    print "Total sketch size:", sketch
    print "Total collisions:", col
}' data/collision_results/summary.tsv

# 3. Count globally distinct hashes
python3 scripts/count_distinct_hashes_parallel.py \
    --db      /scratch/shared_data_new/Logan_yacht_data/processed_data/database_all.db \
    --summary data/collision_results/summary.tsv \
    -t 32 -b 256 -d /dev/shm
```
