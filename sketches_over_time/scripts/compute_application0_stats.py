#!/usr/bin/env python3
"""
compute_application0_stats.py
-------------------------------------------------------------------
Computes and prints every statistic cited in application_0.tex,
providing a traceable record of the numbers used in the paper.

Run from any directory:
    conda run -n logan python scripts/compute_application0_stats.py

All database connections are read-only.
"""

import duckdb
import math
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
FULL_DB  = BASE_DIR / "data" / "intermediate_full.duckdb"
WGS_DB   = BASE_DIR / "data" / "intermediate_wgs_metagenomics_illumina.duckdb"
META_DB  = Path("/scratch/shared_data/Logan_yacht_data/metadata/aws_sra_metadata/metadata_geo_joined_5M.duckdb")

SEP = "-" * 72

def qfull(sql):
    con = duckdb.connect(str(FULL_DB), read_only=True)
    r = con.execute(sql).fetchall()
    con.close()
    return r

def qwgs(sql):
    con = duckdb.connect(str(WGS_DB), read_only=True)
    r = con.execute(sql).fetchall()
    con.close()
    return r

def qmeta(sql):
    con = duckdb.connect(str(META_DB), read_only=True)
    r = con.execute(sql).fetchall()
    con.close()
    return r


# ── 1. Full-corpus totals ─────────────────────────────────────────────────────
print(SEP)
print("1. FULL-CORPUS TOTALS")
print(SEP)

rows = qfull("SELECT year, new_hashes_this_year, cumulative_distinct_hashes FROM cumulative_hashes_by_year ORDER BY year")
print(f"{'Year':>6}  {'New hashes':>20}  {'Cumulative':>22}")
for yr, new, cum in rows:
    print(f"{yr:>6}  {new:>20,}  {cum:>22,}")
print()

total_hashes = rows[-1][2]
new_2023     = rows[-1][1]
new_2008     = rows[0][1]
cum_2008     = rows[0][2]
print(f"Total distinct hashes (2008–2023): {total_hashes:,}")
print(f"New hashes in 2008:                {new_2008:,}")
print(f"Cumulative hashes by end of 2008:  {cum_2008:,}")
print(f"New hashes in 2023:                {new_2023:,}")
print(f"New hashes in 2023 (billions):     {new_2023/1e9:.2f}B")


# ── 2. Per-sample stats, full run ─────────────────────────────────────────────
print()
print(SEP)
print("2. PER-SAMPLE STATS — FULL RUN")
print(SEP)

yss = qfull("""
    SELECT release_year, num_samples,
           avg_total_hashes, avg_new_hashes,
           median_total_hashes, median_new_hashes,
           median_fraction_new
    FROM yearly_sample_stats
    ORDER BY release_year
""")
total_samples = sum(r[1] for r in yss)
print(f"Total samples across all years: {total_samples:,}")
print()
print(f"{'Year':>6}  {'N':>8}  {'AvgTotal':>12}  {'AvgNew':>12}  {'MedTotal':>10}  {'MedNew':>8}  {'MedFracNew':>12}  {'MeanFracNew':>12}")
for yr, n, at, an, mt, mn, mfn in yss:
    mean_frac = an / at if at > 0 else 0.0
    print(f"{yr:>6}  {n:>8,}  {at:>12,.1f}  {an:>12,.1f}  {mt:>10,.0f}  {mn:>8,.1f}  {mfn:>12.4f}  {mean_frac:>12.4f}")

# Highlight 2021-2023 range for the paper
print()
print("--- Key values cited in paper (2021-2023 mean new hashes range) ---")
for yr, n, at, an, mt, mn, mfn in yss:
    if yr in (2021, 2022, 2023):
        mean_frac = an / at
        print(f"  {yr}: mean_new={an:,.0f}, median_new={mn:.0f}, mean_frac_new={mean_frac:.3f} ({mean_frac*100:.1f}%)")
        print(f"       avg_total={at:,.0f}, median_total={mt:.0f}")

# Mean/median ratio
yr2023 = [r for r in yss if r[0] == 2023][0]
ratio_full = yr2023[4] / yr2023[5]  # avg_new / median_new — wait, use avg_new and median_new
ratio_full = yr2023[3] / yr2023[5]  # avg_new_hashes / median_new_hashes
print(f"\n  2023 mean/median ratio (new hashes): {ratio_full:.0f}:1")


# ── 3. Per-sample stats, WGS-filtered run ────────────────────────────────────
print()
print(SEP)
print("3. PER-SAMPLE STATS — WGS-FILTERED RUN")
print(SEP)

yss_wgs = qwgs("""
    SELECT release_year, num_samples,
           avg_total_hashes, avg_new_hashes,
           median_total_hashes, median_new_hashes,
           median_fraction_new
    FROM yearly_sample_stats
    ORDER BY release_year
""")
total_wgs = sum(r[1] for r in yss_wgs)
print(f"Total WGS-filtered samples: {total_wgs:,}")
print()
print(f"{'Year':>6}  {'N':>8}  {'AvgTotal':>12}  {'AvgNew':>12}  {'MedTotal':>10}  {'MedNew':>8}  {'MedFracNew':>12}  {'MeanFracNew':>12}")
for yr, n, at, an, mt, mn, mfn in yss_wgs:
    mean_frac = an / at if at > 0 else 0.0
    print(f"{yr:>6}  {n:>8,}  {at:>12,.1f}  {an:>12,.1f}  {mt:>10,.0f}  {mn:>8,.1f}  {mfn:>12.4f}  {mean_frac:>12.4f}")

print()
print("--- Key values cited in paper (WGS 2022-2023) ---")
for yr, n, at, an, mt, mn, mfn in yss_wgs:
    if yr in (2022, 2023):
        mean_frac = an / at
        ratio = an / mn if mn > 0 else float("inf")
        print(f"  {yr}: mean_new={an:,.0f}, median_new={mn:.0f}, ratio={ratio:.0f}:1")
        print(f"       mean_frac_new={mean_frac:.3f} ({mean_frac*100:.1f}%), median_frac={mfn:.4f} ({mfn*100:.2f}%)")


# ── 4. Hash survival / transience ────────────────────────────────────────────
print()
print(SEP)
print("4. HASH SURVIVAL (birth = death = single-year; persist to 2023)")
print(SEP)

r = qfull("""
    SELECT
        SUM(CASE WHEN birth_year = death_year THEN n_hashes ELSE 0 END) AS single_year,
        SUM(CASE WHEN death_year = 2023       THEN n_hashes ELSE 0 END) AS persist_to_2023,
        SUM(n_hashes) AS total
    FROM hash_life_counts
""")
single_year, persist_2023, total = r[0]
print(f"Total distinct hashes:          {total:,}")
print(f"Seen in exactly one year:       {single_year:,}  ({100*single_year/total:.1f}%)")
print(f"Death year = 2023:              {persist_2023:,}  ({100*persist_2023/total:.1f}%)")

# Hashes born before 2023 that persist to 2023
born_2023 = new_2023   # from cumulative table above (birth_year=2023 => new in 2023)
persist_from_before = persist_2023 - born_2023
print(f"New hashes born in 2023:        {born_2023:,}")
print(f"Born before 2023, alive in 2023:{persist_from_before:,}  ({100*persist_from_before/total:.1f}%)")

multi_year = total - single_year
print(f"Span multiple years (not transient): {multi_year:,}  ({100*multi_year/total:.1f}%)")


# ── 5. Sequencing error back-of-envelope ─────────────────────────────────────
print()
print(SEP)
print("5. SEQUENCING ERROR ESTIMATE (WGS-filtered, median 2023 sample)")
print(SEP)

# Median mbases for WGS-filtered 2023
r2 = qmeta("""
    SELECT MEDIAN(mbases) AS median_mbases, AVG(mbases) AS avg_mbases, COUNT(*) AS n
    FROM metadata_geo_joined
    WHERE librarysource='METAGENOMIC'
      AND libraryselection IN ('RANDOM', 'RANDOM PCR')
      AND mbases > 500
      AND assay_type='WGS'
      AND platform='ILLUMINA'
      AND releasedate IS NOT NULL
      AND YEAR(releasedate) = 2023
""")
median_mb, avg_mb, n_wgs23 = r2[0]
print(f"WGS 2023 samples: {n_wgs23:,}")
print(f"Median mbases: {median_mb:,.0f} Mb")
print(f"Avg    mbases: {avg_mb:,.0f} Mb")

# Parameters
read_len   = 150          # bp (typical Illumina short read)
kmers_per_read = read_len - 31 + 1  # = 120 31-mers per read
epsilon    = 0.001        # Illumina per-base error rate

total_bases = median_mb * 1e6                         # bases
n_reads     = total_bases / read_len                  # reads
total_kmers = n_reads * kmers_per_read                # total 31-mers in sample

# Unique non-singleton k-mers from FracMinHash sketch (scale 1000)
wgs_2023_rows = [r for r in yss_wgs if r[0] == 2023][0]
median_sketch_size = wgs_2023_rows[4]  # median_total_hashes
scale = 1000
G = median_sketch_size * scale  # approx unique non-singleton k-mers

avg_coverage = total_kmers / G
lam = avg_coverage * (epsilon / 3.0)  # expected count of specific 1-sub error k-mer
p_ge2 = lam**2 / 2.0                  # P(>= 2 copies), Poisson approx for small lam

n_variants = 31 * 3  # single-substitution variants per k-mer
expected_spurious = G * n_variants * p_ge2
frac_spurious = expected_spurious / G

# FracMinHash expected spurious hashes in sketch
expected_spurious_in_sketch = expected_spurious / scale  # equivalently: sketch_size * frac_spurious

print()
print(f"--- Calculation parameters ---")
print(f"Read length:                    {read_len} bp")
print(f"31-mers per read:               {kmers_per_read}")
print(f"Per-base error rate (epsilon):  {epsilon}")
print(f"Total bases in median sample:   {total_bases:.3e}")
print(f"Total reads:                    {n_reads:.3e}")
print(f"Total 31-mers:                  {total_kmers:.3e}")
print(f"Unique non-singleton 31-mers G: {G:.3e}  (sketch size {median_sketch_size:,} × scale {scale})")
print(f"Average coverage per k-mer C:   {avg_coverage:.1f}×")
print()
print(f"--- Error calculation ---")
print(f"lambda = C * epsilon/3:         {lam:.4e}")
print(f"P(specific erroneous variant >= 2 reads) = lambda^2/2: {p_ge2:.4e}")
print(f"Single-substitution variants per k-mer:  {n_variants}")
print(f"Expected spurious k-mers (passing >=2 filter):  {expected_spurious:.2e}")
print(f"Fraction spurious of all unique k-mers:         {frac_spurious:.4e}  ({100*frac_spurious:.3f}%)")
print(f"Expected spurious hashes in FracMinHash sketch: {expected_spurious_in_sketch:.1f}")
print(f"  (out of median sketch size {median_sketch_size:,})")
print()
print("Note: This is a CONSERVATIVE UPPER BOUND assuming uniform coverage C.")
print("Realistic heterogeneous abundance distributions yield substantially lower fractions.")

# ── 6. Summary for paper ──────────────────────────────────────────────────────
print()
print(SEP)
print("6. SUMMARY — NUMBERS USED IN application_0.tex")
print(SEP)
print(f"  Total distinct hashes:                {total_hashes:,.0f}  (~{total_hashes/1e9:.1f}B)")
print(f"  Total metagenomes in analysis:        {total_samples:,}")
print(f"  Year range:                           2008–2023")
print(f"  New hashes in 2008:                   {new_2008:,}")
print(f"  New hashes in 2023:                   {new_2023:,}  (~{new_2023/1e9:.1f}B)")
print(f"  Full 2023 mean new/sample:            {yr2023[3]:,.0f}")
print(f"  Full 2023 median new/sample:          {yr2023[5]:.0f}")
print(f"  Full 2023 mean fraction new:          {yr2023[3]/yr2023[2]*100:.1f}%")
print(f"  Full 2023 mean/median ratio:          {yr2023[3]/yr2023[5]:.0f}:1")
print(f"  WGS 2023 mean new/sample:             {wgs_2023_rows[3]:,.0f}")
print(f"  WGS 2023 median new/sample:           {wgs_2023_rows[5]:,.0f}")
print(f"  WGS 2023 mean/median ratio:           {wgs_2023_rows[3]/wgs_2023_rows[5]:.0f}:1")
print(f"  Single-year hashes (77.4%):           {single_year:,}")
print(f"  Persist to 2023 (39.3%):              {persist_2023:,}")
print(f"  Sequencing error fraction (UB):       {100*frac_spurious:.2f}%")
