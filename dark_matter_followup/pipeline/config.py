"""
Shared configuration for the blood dark-matter follow-up.

Standalone / portable: every path here is absolute so the enclosing folder can be
copied elsewhere and still run (the input data paths stay where they are).

Run everything in the `logan` conda env.
"""

import os
from pathlib import Path

# --------------------------------------------------------------------------
# Inputs (absolute, read-only)
# --------------------------------------------------------------------------
INPUT_PARQUET = (
    "/scratch/dmk333/Logan_Analyses/verify_diversity_correlation/data/"
    "hash_diversity_results_full_cov_0.015625/analysis/"
    "filtered_analysis_min_mbases_1620_min_diversity_10/"
    "filtered_data_with_dmi_host_removed.parquet"
)

# 606 GB DuckDB: sample_hashes(sample_id, min_hash) + reference_hashes(hash)
DMI_DB = (
    "/scratch/dmk333/Logan_Analyses/dark_matter_index/data/"
    "samples_with_reference_hashes_cov_0.015625_min_mbases_1620_min_diversity_10.db"
)
# Human host hashes (HPRC Y1/Y2 + T2T CHM13 + GRCh38), k=31
HOST_DB = "/scratch/dmk333/Logan_Analyses/dark_matter_index/data/human_hashes_k31.db"

# --------------------------------------------------------------------------
# Sketch parameters -- verified empirically against the data, not assumed.
#   - max_hash was confirmed as sourmash's own value for scaled=1000
#   - 158,186,253 unitig k-mers / 1000 = 158,186 expected vs 157,772 observed
# NOTE: the "cov_0.015625" in the paths above is an unrelated coverage
#       parameter, NOT the FracMinHash scaled factor.
# --------------------------------------------------------------------------
KSIZE = 31
SCALED = 1000
SEED = 42
MAX_HASH = 18446744073709552  # int(round(2**64 / 1000)), == MinHash(...)._max_hash

# --------------------------------------------------------------------------
# Cohort selection -- must match analyze_body_sites.py exactly
#
# Which body site is processed is set by the DMI_SITE environment variable
# (default "blood"), so every script below is site-agnostic:
#
#     python 01_select_samples.py                 # blood
#     DMI_SITE=skin python 01_select_samples.py   # skin
#
# Blood outputs stay at the top level for backwards compatibility; every other
# site gets its own subdirectory.
# --------------------------------------------------------------------------
SITES = {
    "blood": "human blood metagenome",
    "skin": "human skin metagenome",
    "oral": "human oral metagenome",
    "vaginal": "human vaginal metagenome",
    "nasopharyngeal": "human nasopharyngeal metagenome",
    "lung": "human lung metagenome",
    "gut": "human gut metagenome",
    "feces": "human feces metagenome",
    "saliva": "human saliva metagenome",
    "sputum": "human sputum metagenome",
    "urinary": "human urinary tract metagenome",
}

SITE = os.environ.get("DMI_SITE", "blood").strip().lower()
if SITE not in SITES:
    raise SystemExit(
        f"DMI_SITE={SITE!r} not recognised; choose one of {sorted(SITES)}"
    )
ORGANISM = SITES[SITE]

MIN_MBASES = 100
MAX_HOST_FRACTION = 0.95
MIN_TOTAL_HASHES = 1000
PLATFORMS = ["ILLUMINA"]

# --------------------------------------------------------------------------
# Outputs
# --------------------------------------------------------------------------
# Layout:
#   <root>/pipeline/   this file and the numbered scripts
#   <root>/results/<site>/   small, committable outputs (reports, CSVs, BLAST)
#   <root>/data/<site>/      bulk artifacts (signatures, hashes, sequences) -- gitignored
PIPELINE_DIR = Path(__file__).resolve().parent
ROOT = PIPELINE_DIR.parent
BASE = ROOT / "results" / SITE       # small outputs; "BASE" kept for brevity below
DATA = ROOT / "data" / SITE          # bulk outputs, excluded from git
BLAST_DIR = BASE / "blast"
for _d in (BASE, DATA, BLAST_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# --- small, version-controlled ---
ACCESSIONS_TXT = BASE / f"{SITE}_accessions.txt"
SAMPLES_CSV = BASE / f"{SITE}_samples.csv"
COHORT_NOTES = BASE / f"{SITE}_cohort_notes.txt"
REPORT = BASE / f"{SITE}_dark_matter_report.txt"
SIG_MANIFEST = BASE / "sig_manifest.csv"
LOG_DIR = BASE / "extract_log"

# --- bulk, gitignored (see README for sizes and how to regenerate) ---
SIG_DIR = DATA / "sigs"
HASH_DIR = DATA / "hashes"
SEQ_DIR = DATA / "seqs"
KMER_DIR = DATA / "kmers"
BLASTABLE_FA = DATA / "dark_contigs_blastable.fa"
DARK_DENSITY_CSV = DATA / "contig_dark_density.csv"

# Logan S3 (public, no credentials)
S3_CONTIGS = "s3://logan-pub/c/{acc}/{acc}.contigs.fa.zst"
S3_UNITIGS = "s3://logan-pub/u/{acc}/{acc}.unitigs.fa.zst"


def apply_filters(df):
    """The exact cohort filter used by analyze_body_sites.py."""
    return df[
        (df["mbases_x"] >= MIN_MBASES)
        & (df["host_fraction"] <= MAX_HOST_FRACTION)
        & (df["total_hashes_host_filtered"] >= MIN_TOTAL_HASHES)
        & (df["platform"].isin(PLATFORMS))
        & (df["dmi_host_filtered"].notna())
    ]
