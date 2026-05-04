#!/usr/bin/env python3
"""
SRA Metadata Statistics Analysis
Generates publication-ready LaTeX tables summarizing ~5M SRA accessions.

Usage:
    python analyze_sra_metadata.py

Output:
    sra_metadata_tables.tex  -- LaTeX file with all tables
    sra_stats.json           -- Raw stats in JSON for reproducibility checks
"""

import duckdb
import json
from pathlib import Path

DB_PATH = "/scratch/shared_data_new/Logan_yacht_data/metadata/aws_sra_metadata/metadata_geo_joined_5M.duckdb"
TEX_OUT = "sra_metadata_tables.tex"
JSON_OUT = "sra_stats.json"

con = duckdb.connect(DB_PATH, read_only=True)

# ---------------------------------------------------------------------------
# Query helpers
# ---------------------------------------------------------------------------

def q(sql):
    return con.execute(sql).fetchall()

def qdf(sql):
    return con.execute(sql).df()


# ---------------------------------------------------------------------------
# 1. Grand totals
# ---------------------------------------------------------------------------

total_runs = q("SELECT COUNT(*) FROM metadata_geo_joined")[0][0]
# mbases/mbytes columns are in units of megabases/megabytes;
# divide by 1e6 to get terabases/terabytes, or by 1e9 for peta-.
total_tbases = q("SELECT ROUND(SUM(mbases)/1e6, 1) FROM metadata_geo_joined")[0][0]
total_tbytes = q("SELECT ROUND(SUM(mbytes)/1e6, 1) FROM metadata_geo_joined")[0][0]
n_bioprojects = q("SELECT COUNT(DISTINCT bioproject) FROM metadata_geo_joined WHERE bioproject IS NOT NULL")[0][0]

# ---------------------------------------------------------------------------
# 2. Assay type (sequencing strategy)
# ---------------------------------------------------------------------------

assay_rows = q("""
    SELECT assay_type, COUNT(*) as n
    FROM metadata_geo_joined
    GROUP BY assay_type
    ORDER BY n DESC
""")
# Collapse rare assay types into 'Other'
assay_threshold = 1000
assay_data = []
other_n = 0
for row in assay_rows:
    if row[1] >= assay_threshold:
        assay_data.append(row)
    else:
        other_n += row[1]
if other_n > 0:
    assay_data.append(("Other", other_n))

# ---------------------------------------------------------------------------
# 3. Amplicon target gene breakdown (for AMPLICON assay type)
# ---------------------------------------------------------------------------

amplicon_total = q("SELECT COUNT(*) FROM metadata_geo_joined WHERE assay_type='AMPLICON'")[0][0]

amplicon_gene_rows = q("""
    SELECT
      CASE
        WHEN (jattr ILIKE '%16S%' OR jattr ILIKE '%16s rRNA%' OR jattr ILIKE '%16S rDNA%' OR jattr ILIKE '%16SrRNA%')
          AND jattr NOT ILIKE '%18S%' AND (jattr NOT ILIKE '%ITS%' OR jattr ILIKE '%target_gene%16S%') THEN '16S rRNA'
        WHEN jattr ILIKE '%ITS%' AND jattr NOT ILIKE '%16S%' AND jattr NOT ILIKE '%18S%' THEN 'ITS (fungal/eukaryotic)'
        WHEN (jattr ILIKE '%18S%' OR jattr ILIKE '%ribosomal RNA 18S%')
          AND jattr NOT ILIKE '%16S%' THEN '18S rRNA'
        WHEN jattr ILIKE '%target_gene%' THEN 'Other specified target'
        ELSE 'Target unspecified'
      END as amplicon_type,
      COUNT(*) as n
    FROM metadata_geo_joined
    WHERE assay_type = 'AMPLICON'
    GROUP BY amplicon_type
    ORDER BY n DESC
""")

# ---------------------------------------------------------------------------
# 4. Sequencing platform
# ---------------------------------------------------------------------------

platform_rows = q("""
    SELECT platform, COUNT(*) as n
    FROM metadata_geo_joined
    GROUP BY platform
    ORDER BY n DESC
""")

# ---------------------------------------------------------------------------
# 5. Top 10 organisms (verbatim from the organism field)
# ---------------------------------------------------------------------------

organism_rows = q("""
    SELECT organism, COUNT(*) as n
    FROM metadata_geo_joined
    WHERE organism IS NOT NULL
    GROUP BY organism
    ORDER BY n DESC
    LIMIT 10
""")
n_distinct_organisms = q("SELECT COUNT(DISTINCT organism) FROM metadata_geo_joined WHERE organism IS NOT NULL")[0][0]

# ---------------------------------------------------------------------------
# 6. Geographic distribution (continent)
# ---------------------------------------------------------------------------

continent_rows = q("""
    SELECT geo_loc_name_country_continent_calc as continent, COUNT(*) as n
    FROM metadata_geo_joined
    WHERE geo_loc_name_country_continent_calc IS NOT NULL
      AND geo_loc_name_country_continent_calc != 'uncalculated'
    GROUP BY continent
    ORDER BY n DESC
""")
n_with_geo = q("""
    SELECT COUNT(*) FROM metadata_geo_joined
    WHERE geo_loc_name_country_calc IS NOT NULL
      AND geo_loc_name_country_calc != 'uncalculated'
""")[0][0]

# Top countries
country_rows = q("""
    SELECT geo_loc_name_country_calc as country, COUNT(*) as n
    FROM metadata_geo_joined
    WHERE geo_loc_name_country_calc IS NOT NULL
      AND geo_loc_name_country_calc != 'uncalculated'
    GROUP BY geo_loc_name_country_calc
    ORDER BY n DESC
    LIMIT 15
""")

# ---------------------------------------------------------------------------
# 7. Temporal distribution (by year)
# ---------------------------------------------------------------------------

year_rows = q("""
    SELECT YEAR(releasedate) as year, COUNT(*) as n
    FROM metadata_geo_joined
    WHERE releasedate IS NOT NULL
    GROUP BY year
    ORDER BY year
""")

# ---------------------------------------------------------------------------
# Serialize raw stats to JSON
# ---------------------------------------------------------------------------

stats = {
    "total_runs": total_runs,
    "total_tbases": total_tbases,
    "total_tbytes": total_tbytes,
    "n_bioprojects": n_bioprojects,
    "assay_types": [{"assay_type": r[0], "n": r[1]} for r in assay_data],
    "amplicon_target_genes": [{"target": r[0], "n": r[1]} for r in amplicon_gene_rows],
    "platforms": [{"platform": r[0], "n": r[1]} for r in platform_rows],
    "top_organisms": [{"organism": r[0], "n": r[1]} for r in organism_rows],
    "n_distinct_organisms": n_distinct_organisms,
    "continents": [{"continent": r[0], "n": r[1]} for r in continent_rows],
    "top_countries": [{"country": r[0], "n": r[1]} for r in country_rows],
    "by_year": [{"year": r[0], "n": r[1]} for r in year_rows],
}
Path(JSON_OUT).write_text(json.dumps(stats, indent=2))
print(f"Wrote {JSON_OUT}")

# ---------------------------------------------------------------------------
# LaTeX helpers
# ---------------------------------------------------------------------------

def pct(n, total):
    return f"{100.0 * n / total:.1f}\\%"

def fmt(n):
    return f"{n:,}"


# ---------------------------------------------------------------------------
# Pretty-print helpers
# ---------------------------------------------------------------------------

PLATFORM_LABELS = {
    "ILLUMINA": "Illumina",
    "LS454": "Roche 454",
    "ION_TORRENT": "Ion Torrent",
    "PACBIO_SMRT": "PacBio SMRT",
    "OXFORD_NANOPORE": "Oxford Nanopore",
    "BGISEQ": "BGI-Seq",
    "DNBSEQ": "DNB-Seq",
    "ABI_SOLID": "ABI SOLiD",
    "COMPLETE_GENOMICS": "Complete Genomics",
    "CAPILLARY": "Sanger (capillary)",
    "HELICOS": "Helicos",
}

ASSAY_LABELS = {
    "AMPLICON": "Amplicon sequencing",
    "WGS": "Whole-genome shotgun (WGS)",
    "OTHER": "Other (unspecified)",
    "RNA-Seq": "RNA-Seq",
    "WGA": "Whole-genome amplification (WGA)",
    "Targeted-Capture": "Targeted capture",
    "WCS": "Whole-chromosome shotgun (WCS)",
    "RAD-Seq": "RAD-Seq",
    "POOLCLONE": "Pool clone",
    "CLONE": "Clone",
    "miRNA-Seq": "miRNA-Seq",
    "WXS": "Whole-exome sequencing (WXS)",
    "Tn-Seq": "Tn-Seq",
    "CLONEEND": "Clone end",
    "Other": "Other (rare strategies)",
}

# ---------------------------------------------------------------------------
# Build LaTeX
# ---------------------------------------------------------------------------

lines = []

def emit(s=""):
    lines.append(s)

emit(r"\documentclass[11pt]{article}")
emit(r"\usepackage{booktabs}")
emit(r"\usepackage{longtable}")
emit(r"\usepackage{array}")
emit(r"\usepackage{geometry}")
emit(r"\usepackage{caption}")
emit(r"\usepackage{pdflscape}")
emit(r"\geometry{margin=1in}")
emit(r"\begin{document}")
emit()

# ---- Table 1: Dataset overview ----
emit(r"\begin{table}[h!]")
emit(r"\centering")
emit(r"\caption{Summary of the SRA metadata corpus. Terabases (Tbases) and terabytes (Tbytes) are computed from the \texttt{mbases} and \texttt{mbytes} fields (both stored in megabase/megabyte units) by dividing by $10^6$.}")
emit(r"\label{tab:overview}")
emit(r"\begin{tabular}{lr}")
emit(r"\toprule")
emit(r"\textbf{Metric} & \textbf{Value} \\")
emit(r"\midrule")
emit(f"Total SRA runs & {fmt(total_runs)} \\\\")
emit(f"Total sequencing data (Tbases) & {fmt(int(total_tbases))} \\\\")
emit(f"Total data size (Tbytes) & {fmt(int(total_tbytes))} \\\\")
emit(f"Unique BioProjects & {fmt(n_bioprojects)} \\\\")
emit(f"Runs with geographic annotation & {fmt(n_with_geo)} ({pct(n_with_geo, total_runs)}) \\\\")
emit(r"\bottomrule")
emit(r"\end{tabular}")
emit(r"\end{table}")
emit()

# ---- Table 2: Assay type (sequencing strategy) ----
emit(r"\begin{table}[h!]")
emit(r"\centering")
emit(r"\caption{Distribution of SRA runs by sequencing strategy (\texttt{assay\_type}). Strategies with fewer than 1{,}000 runs are aggregated into `Other'.}")
emit(r"\label{tab:assay_type}")
emit(r"\begin{tabular}{lrr}")
emit(r"\toprule")
emit(r"\textbf{Sequencing strategy} & \textbf{Runs} & \textbf{\%} \\")
emit(r"\midrule")
for row in assay_data:
    name, n = row[0], row[1]
    label = ASSAY_LABELS.get(name, name)
    emit(f"{label} & {fmt(n)} & {pct(n, total_runs)} \\\\")
emit(r"\midrule")
emit(f"\\textbf{{Total}} & \\textbf{{{fmt(total_runs)}}} & 100.0\\% \\\\")
emit(r"\bottomrule")
emit(r"\end{tabular}")
emit(r"\end{table}")
emit()

# ---- Table 3: Amplicon target gene ----
emit(r"\begin{table}[h!]")
emit(r"\centering")
emit(r"\caption{Classification of the {fmt(fmt(amplicon_total))} AMPLICON runs by inferred target gene. Target gene was derived by keyword matching against the \texttt{jattr} metadata field.}".replace("{fmt(fmt(amplicon_total))}", fmt(amplicon_total)))
emit(r"\label{tab:amplicon}")
emit(r"\begin{tabular}{lrr}")
emit(r"\toprule")
emit(r"\textbf{Target gene} & \textbf{Runs} & \textbf{\% of AMPLICON} \\")
emit(r"\midrule")
for row in amplicon_gene_rows:
    name, n = row[0], row[1]
    emit(f"{name} & {fmt(n)} & {pct(n, amplicon_total)} \\\\")
emit(r"\midrule")
emit(f"\\textbf{{Total}} & \\textbf{{{fmt(amplicon_total)}}} & 100.0\\% \\\\")
emit(r"\bottomrule")
emit(r"\end{tabular}")
emit(r"\end{table}")
emit()

# ---- Table 4: Sequencing platform ----
emit(r"\begin{table}[h!]")
emit(r"\centering")
emit(r"\caption{Distribution of SRA runs by sequencing platform.}")
emit(r"\label{tab:platform}")
emit(r"\begin{tabular}{lrr}")
emit(r"\toprule")
emit(r"\textbf{Platform} & \textbf{Runs} & \textbf{\%} \\")
emit(r"\midrule")
for row in platform_rows:
    name, n = row[0], row[1]
    label = PLATFORM_LABELS.get(name, name)
    emit(f"{label} & {fmt(n)} & {pct(n, total_runs)} \\\\")
emit(r"\midrule")
emit(f"\\textbf{{Total}} & \\textbf{{{fmt(total_runs)}}} & 100.0\\% \\\\")
emit(r"\bottomrule")
emit(r"\end{tabular}")
emit(r"\end{table}")
emit()

# ---- Table 5: Top 10 organisms (verbatim) ----
emit(r"\begin{table}[h!]")
emit(r"\centering")
emit(f"\\caption{{Top 10 most frequent values of the \\texttt{{organism}} metadata field across all {fmt(total_runs)} SRA runs. Values are reported verbatim as deposited; {fmt(n_distinct_organisms)} distinct values exist in total.}}")
emit(r"\label{tab:organisms}")
emit(r"\begin{tabular}{lrr}")
emit(r"\toprule")
emit(r"\textbf{Organism} & \textbf{Runs} & \textbf{\%} \\")
emit(r"\midrule")
for row in organism_rows:
    name, n = row[0], row[1]
    safe_name = name.replace("&", r"\&").replace("_", r"\_")
    emit(f"\\textit{{{safe_name}}} & {fmt(n)} & {pct(n, total_runs)} \\\\")
emit(r"\bottomrule")
emit(r"\end{tabular}")
emit(r"\end{table}")
emit()

# ---- Table 6: Geographic distribution by continent ----
emit(r"\begin{table}[h!]")
emit(r"\centering")
emit(r"\caption{Geographic distribution of SRA runs by continent. Percentages are of runs with a calculable geographic annotation ($n = $" + fmt(n_with_geo) + r").}")
emit(r"\label{tab:continent}")
emit(r"\begin{tabular}{lrr}")
emit(r"\toprule")
emit(r"\textbf{Continent} & \textbf{Runs} & \textbf{\% of geo-annotated} \\")
emit(r"\midrule")
for row in continent_rows:
    name, n = row[0], row[1]
    emit(f"{name} & {fmt(n)} & {pct(n, n_with_geo)} \\\\")
emit(r"\midrule")
emit(f"\\textbf{{Total (geo-annotated)}} & \\textbf{{{fmt(n_with_geo)}}} & 100.0\\% \\\\")
emit(r"\bottomrule")
emit(r"\end{tabular}")
emit(r"\end{table}")
emit()

# ---- Table 7: Top countries ----
emit(r"\begin{table}[h!]")
emit(r"\centering")
emit(r"\caption{Top 15 countries of origin for SRA runs. Percentages are of runs with a calculable country annotation ($n = $" + fmt(n_with_geo) + r").}")
emit(r"\label{tab:countries}")
emit(r"\begin{tabular}{lrr}")
emit(r"\toprule")
emit(r"\textbf{Country} & \textbf{Runs} & \textbf{\% of geo-annotated} \\")
emit(r"\midrule")
for row in country_rows:
    name, n = row[0], row[1]
    emit(f"{name} & {fmt(n)} & {pct(n, n_with_geo)} \\\\")
emit(r"\bottomrule")
emit(r"\end{tabular}")
emit(r"\end{table}")
emit()

# ---- Table 8: Temporal distribution ----
emit(r"\begin{table}[h!]")
emit(r"\centering")
emit(r"\caption{Annual deposition of SRA runs in this corpus by public release year.}")
emit(r"\label{tab:temporal}")
emit(r"\begin{tabular}{rrr}")
emit(r"\toprule")
emit(r"\textbf{Year} & \textbf{Runs deposited} & \textbf{\%} \\")
emit(r"\midrule")
for row in year_rows:
    year, n = row[0], row[1]
    emit(f"{year} & {fmt(n)} & {pct(n, total_runs)} \\\\")
emit(r"\midrule")
emit(f"\\textbf{{Total}} & \\textbf{{{fmt(total_runs)}}} & 100.0\\% \\\\")
emit(r"\bottomrule")
emit(r"\end{tabular}")
emit(r"\end{table}")
emit()

emit(r"\end{document}")

tex_content = "\n".join(lines)
Path(TEX_OUT).write_text(tex_content)
print(f"Wrote {TEX_OUT}")

# ---------------------------------------------------------------------------
# Print a human-readable summary to stdout
# ---------------------------------------------------------------------------

print()
print("=" * 60)
print("DATASET SUMMARY")
print("=" * 60)
print(f"  Total runs        : {fmt(total_runs)}")
print(f"  Total Tbases      : {fmt(int(total_tbases))}")
print(f"  Total Tbytes      : {fmt(int(total_tbytes))}")
print(f"  Unique BioProjects: {fmt(n_bioprojects)}")
print()
print("TOP ASSAY TYPES")
for r in assay_data[:8]:
    print(f"  {r[0]:<25s} {fmt(r[1]):>10s}  ({pct(r[1], total_runs)})")
print()
print("AMPLICON TARGET GENES")
for r in amplicon_gene_rows:
    print(f"  {r[0]:<35s} {fmt(r[1]):>10s}  ({pct(r[1], amplicon_total)} of AMPLICON)")
print()
print("TOP 10 ORGANISMS")
for r in organism_rows:
    print(f"  {r[0]:<35s} {fmt(r[1]):>10s}  ({pct(r[1], total_runs)})")
print()
print("Done.")
