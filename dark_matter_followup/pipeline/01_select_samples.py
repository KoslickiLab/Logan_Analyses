#!/usr/bin/env python3
"""
Step 1: select the samples of one body site that appear in body_site_analysis.png.

Re-applies the exact analyze_body_sites.py filters, then writes:
  {site}_accessions.txt  -- one accession per line
  {site}_samples.csv     -- metrics + flattened SRA attributes
  {site}_cohort_notes.txt -- bioproject breakdown and provenance warnings

Usage:  conda activate logan && python 01_select_samples.py  (set DMI_SITE=skin for skin)
"""

import sys

import pandas as pd

import config as C

# Attribute keys worth promoting to their own CSV columns: they carry the
# library-prep / host-depletion provenance that decides how to read the DMI.
PROVENANCE_KEYS = [
    "isolation_source_sam",
    "env_medium_sam",
    "env_broad_scale_sam",
    "sample_material_processing_sam",
    "data_preprocessing_sam",
    "host_sam",
    "collection_date_sam",
    "lat_lon_sam_s_dpl34",
]


def flatten_attributes(attrs):
    """SRA `attributes` is a list of {'k':..., 'v':...}; drop the noisy
    primary_search entries and collapse to a dict."""
    out = {}
    if attrs is None:
        return out
    for entry in attrs:
        k, v = entry.get("k"), entry.get("v")
        if k and k != "primary_search" and k not in out:
            out[k] = v
    return out


def main():
    df = pd.read_parquet(C.INPUT_PARQUET)
    n_loaded = len(df)

    cohort = C.apply_filters(df)
    cohort = cohort[cohort["organism"] == C.ORGANISM].copy()
    n_site_raw = int((df["organism"] == C.ORGANISM).sum())

    if cohort.empty:
        sys.exit("ERROR: no samples survived the filters")

    cohort = cohort.sort_values("dmi_host_filtered", ascending=False)

    # ---- accessions -------------------------------------------------------
    C.ACCESSIONS_TXT.write_text("\n".join(cohort["accession"]) + "\n")

    # ---- per-sample table -------------------------------------------------
    flat = cohort["attributes"].apply(flatten_attributes)
    cols = [
        "accession", "bioproject", "sra_study", "center_name", "instrument",
        "librarylayout", "libraryselection", "librarysource", "releasedate",
        "mbases_x", "total_hashes_original", "host_hashes", "host_fraction",
        "total_hashes_host_filtered", "unmapped_hashes_host_filtered",
        "mapped_hashes_host_filtered", "dmi_host_filtered", "alpha_diversity",
        "hashes_per_mb_host_filtered", "diversity_per_mb",
    ]
    out = cohort[cols].copy()
    for key in PROVENANCE_KEYS:
        out[key] = flat.apply(lambda d, k=key: d.get(k, ""))
    out.to_csv(C.SAMPLES_CSV, index=False)

    # ---- cohort notes -----------------------------------------------------
    grp = (
        cohort.groupby(["bioproject", "center_name"])
        .agg(
            n=("accession", "size"),
            dmi_median=("dmi_host_filtered", "median"),
            alpha_div_median=("alpha_diversity", "median"),
            hashes_median=("total_hashes_original", "median"),
            host_frac_median=("host_fraction", "median"),
        )
        .sort_values("n", ascending=False)
        .reset_index()
    )

    lines = [
        f"{C.SITE.capitalize()} cohort from body_site_analysis.png",
        "=" * 62,
        f"Source parquet: {C.INPUT_PARQUET}",
        "",
        "FILTERS (identical to analyze_body_sites.py)",
        "-" * 62,
        f"  organism == '{C.ORGANISM}'",
        f"  mbases_x >= {C.MIN_MBASES}; host_fraction <= {C.MAX_HOST_FRACTION}; "
        f"total_hashes_host_filtered >= {C.MIN_TOTAL_HASHES}",
        f"  platform in {C.PLATFORMS}; dmi_host_filtered not null",
        f"  {n_loaded:,} rows loaded -> {n_site_raw} cohort -> {len(cohort)} analysed",
        "",
        "BIOPROJECT BREAKDOWN",
        "-" * 62,
        grp.to_string(index=False),
        "",
        "PROVENANCE WARNINGS",
        "-" * 62,
    ]

    # Surface library-prep provenance that changes how the DMI should be read.
    for bp, sub in out.groupby("bioproject", sort=False):
        proc = sorted({v for v in sub["sample_material_processing_sam"] if v})
        prep = sorted({v for v in sub["data_preprocessing_sam"] if v})
        src = sorted({v for v in sub["isolation_source_sam"] if v}) or sorted(
            {v for v in sub["env_medium_sam"] if v}
        )
        lines.append(f"  {bp}  (n={len(sub)})")
        if src:
            lines.append(f"     isolation source        : {'; '.join(src)}")
        if proc:
            lines.append(f"     material processing     : {'; '.join(proc)}")
        if prep:
            lines.append(f"     submitter preprocessing : {'; '.join(prep)}")
        lines.append("")

    # Concentration check: a site carried by one study cannot be read as a
    # property of the body site.
    top_share = grp["n"].iloc[0] / len(cohort)
    lines += [
        f"READ THIS BEFORE INTERPRETING THE {C.SITE.upper()} DMI",
        "-" * 62,
        f"  bioprojects contributing : {cohort['bioproject'].nunique()}",
        f"  largest single project   : {grp['bioproject'].iloc[0]} "
        f"({grp['n'].iloc[0]} runs, {100 * top_share:.1f}% of the cohort)",
        "",
    ]
    if top_share >= 0.5:
        lines += [
            "  WARNING: over half of this cohort is one bioproject, so its DMI is a",
            "  property of that study's protocol as much as of the body site. Check",
            "  the material-processing and preprocessing fields above: VLP/virome",
            "  enrichment, submitter-side host removal, or whole-genome amplification",
            "  all inflate the apparent dark fraction of a low-biomass sample.",
            "",
        ]
    else:
        lines += [
            "  No single bioproject dominates this cohort, so per-study protocol",
            "  artefacts are less likely to drive the site-level DMI. Per-project",
            "  breakdown is still worth checking before drawing conclusions.",
            "",
        ]
    C.COHORT_NOTES.write_text("\n".join(lines))

    print(f"cohort samples analysed : {len(cohort)} (of {n_site_raw} in parquet)")
    print(grp.to_string(index=False))
    print(f"\nwrote {C.ACCESSIONS_TXT}")
    print(f"wrote {C.SAMPLES_CSV}")
    print(f"wrote {C.COHORT_NOTES}")


if __name__ == "__main__":
    main()
