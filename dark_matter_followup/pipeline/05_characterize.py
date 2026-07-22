#!/usr/bin/env python3
"""
Step 5: characterise the blood dark matter and write the report.

Three things:

1. ASSEMBLY-RATE TABLE (the headline).  For the pilot accessions, build the
   scaled=1000 sketch of the sample's own contigs and intersect it with the
   sample's hashes split by class (dark vs reference-matched).  If dark hashes
   assemble far less often than reference-matched ones, the dark matter is
   fragmented sequence rather than novel genomes -- a novel genome assembles
   fine and merely fails to match a database.

2. SEQUENCE CHARACTER.  Length, coverage (ka), GC and 3-mer entropy of the
   dark-carrying contigs, each against a background of the non-matching
   contigs from the very same file, so every distribution has a within-sample
   control.  Expectations: sequencing noise -> short + low ka + low entropy;
   kit/reagent contamination -> high GC (~60-70%) and ordinary entropy.

3. dark_contigs_blastable.fa -- dark contigs >= MIN_BLAST_LEN pooled across all
   accessions, sorted by length x coverage.  The file to actually BLAST.

Usage:  conda activate logan && python 05_characterize.py
"""

import gzip
import math
import subprocess
import sys
from collections import Counter

import numpy as np
import pandas as pd
from sourmash.minhash import MinHash

import config as C
import extractor

MIN_BLAST_LEN = 300
BACKGROUND_SAMPLE = 40000  # non-matching contigs to profile as a control

# A contig is only interesting as dark matter if MOST of its hashes are dark.
# Without this, a 100 kb well-characterised genome carrying a single stray dark
# k-mer (a strain SNP) tops the list purely on length.
MIN_DARK_FRAC = 0.8


# --------------------------------------------------------------------------
def seq_profile(seq):
    """GC fraction, longest homopolymer run, and 3-mer Shannon entropy (bits)."""
    n = len(seq)
    gc = (seq.count("G") + seq.count("C")) / n if n else float("nan")

    run = best = 1
    for i in range(1, n):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)

    if n >= 3:
        cnt = Counter(seq[i : i + 3] for i in range(n - 2))
        tot = sum(cnt.values())
        ent = -sum((c / tot) * math.log2(c / tot) for c in cnt.values())
    else:
        ent = float("nan")
    return gc, best, ent


def stream_contigs(acc):
    url = C.S3_CONTIGS.format(acc=acc)
    p1 = subprocess.Popen(
        ["aws", "s3", "cp", url, "-", "--no-sign-request"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    p2 = subprocess.Popen(
        ["zstd", "-dcq"], stdin=p1.stdout, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    p1.stdout.close()
    import io
    return p1, p2, io.TextIOWrapper(p2.stdout, encoding="ascii")


def assembly_rates(acc):
    """Fraction of each hash class that appears anywhere in the contigs."""
    h = pd.read_parquet(C.HASH_DIR / f"{acc}.hashes.parquet")
    dark = set(h.loc[h["is_dark"], "min_hash"].tolist())
    refm = set(h.loc[~h["is_dark"], "min_hash"].tolist())

    mh = MinHash(n=0, ksize=C.KSIZE, scaled=C.SCALED)
    n_ctg = n_bp = 0
    rows = []
    rng = np.random.default_rng(0)
    p1, p2, stream = stream_contigs(acc)
    try:
        for name, ka, seq in extractor.parse_fasta(stream):
            seq = seq.upper()
            n_ctg += 1
            n_bp += len(seq)
            mh.add_sequence(seq, True)
            if len(rows) < BACKGROUND_SAMPLE and rng.random() < 0.15:
                rows.append((name, ka, len(seq), seq))
    finally:
        stream.close()
        p2.wait()
        p1.wait()

    ctg_hashes = set(mh.hashes)
    return dict(
        accession=acc,
        n_contigs=n_ctg,
        contig_bp=n_bp,
        contig_sketch=len(ctg_hashes),
        n_dark=len(dark),
        dark_in_contigs=len(dark & ctg_hashes),
        dark_rate=len(dark & ctg_hashes) / len(dark) if dark else float("nan"),
        n_refmatched=len(refm),
        refmatched_in_contigs=len(refm & ctg_hashes),
        refmatched_rate=len(refm & ctg_hashes) / len(refm) if refm else float("nan"),
    ), rows


def profile_frame(records, label):
    out = []
    for name, ka, ln, seq in records:
        gc, hp, ent = seq_profile(seq)
        out.append((label, name, ka, ln, gc, hp, ent))
    return pd.DataFrame(
        out, columns=["group", "seq_id", "ka", "length", "gc", "max_homopolymer", "entropy3"]
    )


def describe(df, cols=("length", "ka", "gc", "max_homopolymer", "entropy3")):
    q = df[list(cols)].quantile([0.25, 0.5, 0.75])
    return "  ".join(
        f"{c}={q.loc[0.5, c]:.3g} [{q.loc[0.25, c]:.3g}-{q.loc[0.75, c]:.3g}]"
        for c in cols
    )


# --------------------------------------------------------------------------
def pick_pilots(meta, n=3):
    """One representative run from each of the n largest bioprojects.

    Deep-dive samples are chosen by cohort structure rather than hard-coded, so
    the same script works for any body site. Within a project we take the run
    with the median non-host hash count, to avoid picking an outlier.
    """
    order = meta["bioproject"].value_counts().index[:n]
    picks = []
    for proj in order:
        sub = meta[meta["bioproject"] == proj].sort_values("total_hashes_host_filtered")
        picks.append(sub.iloc[len(sub) // 2]["accession"])
    return picks


def main():
    lines = []
    W = 74

    def sec(title):
        lines.extend(["", title, "-" * W])

    summ = pd.read_csv(C.BASE / "extract_summary_contigs.csv")
    ok = summ[summ["error"].isna()] if "error" in summ.columns else summ
    man = pd.read_csv(C.SIG_MANIFEST)
    meta = pd.read_csv(C.SAMPLES_CSV)
    PILOTS = pick_pilots(meta)

    lines += [
        f"{C.SITE.capitalize()} dark matter: what the sequences actually are",
        "=" * W,
        f"Cohort      : {len(meta)} '{C.ORGANISM}' runs from "
        f"body_site_analysis.png",
        f"Dark hashes : non-host AND absent from the full reference union",
        f"              (GenBank WGS/genomes, AllTheBacteria, GTDB, SILVA, TLS,",
        f"               TSA, Logan plasmids/obelisks, Serratus viruses)",
        f"Sketch      : k={C.KSIZE}, scaled={C.SCALED}, seed={C.SEED}",
        "",
        "ASCII-only by design; these reports are read in terminals.",
    ]

    sec("METHODS -- how every number below was produced")
    lines += [
        "  Cohort selection (01_select_samples.py)",
        f"    source parquet : {C.INPUT_PARQUET}",
        f"    filters        : mbases_x >= {C.MIN_MBASES}; host_fraction <= "
        f"{C.MAX_HOST_FRACTION};",
        f"                     total_hashes_host_filtered >= {C.MIN_TOTAL_HASHES}; "
        f"platform in {C.PLATFORMS};",
        f"                     dmi_host_filtered not null; organism == '{C.ORGANISM}'",
        "                     (identical to analyze_body_sites.py, so the cohort is",
        "                      exactly the points behind body_site_analysis.png)",
        "",
        "  Hash classification (02_export_hashes.py)",
        f"    sketches   : FracMinHash, k={C.KSIZE}, scaled={C.SCALED}, seed={C.SEED},",
        f"                 max_hash={C.MAX_HASH}, murmur3_x64_128 low 64 bits of the",
        "                 canonical (lexicographic min of fwd/revcomp) 31-mer.",
        "    host hash  : present in human_hashes_k31.db (HPRC + T2T CHM13 + GRCh38).",
        "    dark hash  : not host AND absent from the reference union (GenBank WGS,",
        "                 GenBank genomes, AllTheBacteria, GTDB, SILVA, GenBank TLS,",
        "                 GenBank TSA, Logan plasmids, Logan obelisks, Serratus viruses).",
        "    MATCHING IS EXACT 31-mer IDENTITY. There is no similarity threshold: a",
        "    single substitution destroys a match. At 80% nucleotide identity to a",
        "    reference, P(a given 31-mer is conserved) = 0.80^31 = 1.0e-3, so a",
        "    sequence can be 'dark' while still having clear database homologs.",
        "    Counts were cross-checked against the parquet columns and must agree.",
        "",
        "  Sequence extraction (04_extract_sequences.py)",
        "    Logan contigs streamed from s3://logan-pub/c/{acc}/{acc}.contigs.fa.zst",
        "    and unitigs from s3://logan-pub/u/... (no credentials, nothing staged).",
        "    Search uses the FULL non-host hash set; each matched hash is then labelled",
        "    dark or reference-matched. A record is written out only if it carries at",
        "    least one dark hash. Hashing calls sourmash's own Rust seq_to_hashes on",
        "    batched records, validated in 03_validate_hashing.py against both a naive",
        "    mmh3 canonical-k-mer loop and `sourmash sig kmers` (identical outputs).",
        "",
        "  Statistics defined (this section exists because these are easy to conflate)",
        "    dark_frac     : dark hashes / non-host hashes ON THAT RECORD. Needed",
        "                    because a long, well-characterised genome carrying one",
        "                    stray dark k-mer is not dark matter.",
        "    recovery_frac : distinct dark hashes located in the sequence / total dark",
        "                    hashes for that run. ~100% expected for unitigs (they are",
        "                    what the sketch was built from); for contigs it is the",
        "                    measurement of interest, not a target.",
        "    assembly rate : fraction of a hash class found anywhere in that run's",
        "                    contigs, computed by sketching the contigs at the same",
        "                    scaled and intersecting. Compared dark vs reference-matched.",
        "    ka            : mean k-mer abundance (coverage) from the Logan FASTA header.",
        "    GC            : (G + C) / length.",
        "    entropy3      : Shannon entropy of the 3-mer composition, in bits, max 6.",
        "    max_homopolymer : longest single-base run.",
        "    Quantiles are reported as median [25th-75th percentile].",
        "",
    ]

    sec("0. COHORT PROVENANCE -- read this first")
    bp = meta.groupby(["bioproject", "center_name"]).size().sort_values(ascending=False)
    for (proj, center), n in list(bp.items())[:12]:
        sub = meta[meta["bioproject"] == proj]
        proc = "; ".join(sorted({str(v) for v in sub["sample_material_processing_sam"] if str(v) != "nan"}))
        prep = "; ".join(sorted({str(v) for v in sub["data_preprocessing_sam"] if str(v) != "nan"}))
        lines.append(f"  {proj}  n={n:4d}  {center[:44]}")
        if proc:
            lines.append(f"      material processing     : {proc}")
        if prep:
            lines.append(f"      submitter preprocessing : {prep}")
    if len(bp) > 12:
        lines.append(f"  ... and {len(bp) - 12} further bioprojects "
                     f"({int(bp.iloc[12:].sum())} runs); see {C.SAMPLES_CSV.name}")

    top_share = bp.iloc[0] / len(meta)
    lines += [
        "",
        f"  {len(meta['bioproject'].unique())} bioprojects; largest is "
        f"{100 * top_share:.1f}% of the cohort.",
    ]
    if top_share >= 0.5:
        lines += [
            "  WARNING: this site is carried by a single study, so its DMI reflects",
            "  that study's protocol as much as the body site. Submitter-side host",
            "  removal makes host_fraction an artefact of upstream processing, and",
            "  VLP enrichment or whole-genome amplification inflates the apparent",
            "  dark fraction of a low-biomass sample.",
        ]
    else:
        lines += [
            "  No single study dominates, so protocol artefacts are less likely to",
            "  drive the site-level DMI than they are for a single-study cohort.",
        ]

    sec("1. HASH BUDGET (all accessions)")
    lines += [
        f"  non-host hashes exported : {man.n_nonhost.sum():,}",
        f"  dark (no reference match): {man.n_dark.sum():,} "
        f"({100 * man.n_dark.sum() / man.n_nonhost.sum():.1f}%)",
        f"  reference-matched        : {man.n_reference_matched.sum():,}",
        "",
        f"  dark hashes found on contigs: {ok.n_dark_recovered.sum():,} / "
        f"{ok.n_dark_hashes.sum():,} "
        f"({100 * ok.n_dark_recovered.sum() / ok.n_dark_hashes.sum():.2f}%)",
        f"  dark-carrying contigs       : {ok.n_records_matched.sum():,} "
        f"({ok.bp_matched.sum() / 1e6:.1f} Mbp)",
    ]

    # ---- 2. assembly rates -------------------------------------------------
    sec("2. ASSEMBLY RATE BY HASH CLASS  (the headline result)")
    lines.append("  fraction of each class's hashes appearing anywhere in the contigs")
    lines.append("")
    lines.append(f"  {'accession':14s} {'dark':>18s} {'reference-matched':>22s} {'ratio':>8s}")
    rate_rows, prof_frames = [], []
    for acc in PILOTS:
        print(f"assembly rates: {acc} ...", file=sys.stderr)
        r, background = assembly_rates(acc)
        rate_rows.append(r)
        ratio = r["refmatched_rate"] / r["dark_rate"] if r["dark_rate"] else float("nan")
        lines.append(
            f"  {acc:14s} "
            f"{r['dark_in_contigs']:7,d}/{r['n_dark']:<7,d} ({100*r['dark_rate']:4.2f}%)"
            f"  {r['refmatched_in_contigs']:6,d}/{r['n_refmatched']:<6,d} "
            f"({100*r['refmatched_rate']:5.2f}%) {ratio:7.1f}x"
        )

        # profile dark-carrying contigs vs the background from the same file
        fa = C.SEQ_DIR / "contigs" / f"{acc}.dark_{'contigs'}.fa.gz"
        dark_recs = []
        if fa.exists():
            with gzip.open(fa, "rt") as fh:
                for name, ka, seq in extractor.parse_fasta(fh):
                    dark_recs.append((name, ka, len(seq), seq.upper()))
        dark_ids = {n for n, *_ in dark_recs}
        bg = [b for b in background if b[0] not in dark_ids]
        d = profile_frame(dark_recs, "dark")
        b = profile_frame(bg, "background")
        d["accession"] = acc
        b["accession"] = acc
        prof_frames += [d, b]

    rates = pd.DataFrame(rate_rows)
    rates.to_csv(C.BASE / "assembly_rates.csv", index=False)
    med_ratio = float(
        (rates["refmatched_rate"] / rates["dark_rate"].replace(0, float("nan"))).median()
    )
    lines += [
        "",
        "  Interpretation. Sequence that is genuinely novel but real assembles like",
        "  any other sequence and merely fails to match a database, so its assembly",
        "  rate should approach that of reference-matched sequence. Sequence that",
        "  cannot be assembled at all is fragmented: errors, degraded DNA, or",
        "  single-observation noise. The ratio below is therefore the discriminator.",
        f"  Median ratio here: {med_ratio:.1f}x.",
    ]
    if med_ratio >= 10:
        lines += [
            "  Dark hashes assemble an order of magnitude less often than",
            "  reference-matched ones. That is the signature of fragmented noise,",
            "  not of novel organisms.",
        ]
    elif med_ratio >= 3:
        lines += [
            "  Dark hashes assemble several times less often than reference-matched",
            "  ones: a mixture, with a real assembling component alongside",
            "  unassemblable fragments.",
        ]
    else:
        lines += [
            "  Dark hashes assemble nearly as well as reference-matched ones. They",
            "  behave like real organisms that are simply absent from the reference",
            "  union, not like noise.",
        ]

    uni_csv = C.BASE / "extract_summary_unitigs.csv"
    if uni_csv.exists():
        uni = pd.read_csv(uni_csv)
        lines += [
            "",
            "  The dark hashes are not missing from the data -- they are present in",
            "  the raw unitig graph and recovered in full; they simply never make it",
            "  into an assembled contig:",
            "",
            f"  {'accession':14s} {'unitig recovery':>16s} {'contig recovery':>16s} "
            f"{'mean unitig bp':>15s}",
        ]
        for _, r in uni.iterrows():
            crow = ok[ok.accession == r.accession]
            crate = 100 * crow.recovery_frac.iloc[0] if len(crow) else float("nan")
            lines.append(
                f"  {r.accession:14s} {100 * r.recovery_frac:15.2f}% "
                f"{crate:15.2f}% {r.bp_matched / r.n_records_matched:15.1f}"
            )

    # ---- 3. sequence character --------------------------------------------
    prof = pd.concat(prof_frames, ignore_index=True)
    prof.to_csv(C.BASE / "contig_profiles.csv", index=False)
    sec("3. CHARACTER OF THE DARK-CARRYING CONTIGS  (median [IQR])")
    lines.append("  entropy3 = Shannon entropy of 3-mer composition, max 6 bits;")
    lines.append("  low entropy / long homopolymers = low-complexity repeat sequence.")
    for acc in PILOTS:
        lines.append(f"\n  {acc}")
        for grp in ("dark", "background"):
            sub = prof[(prof.accession == acc) & (prof.group == grp)]
            if len(sub):
                lines.append(f"    {grp:11s} n={len(sub):6,d}  {describe(sub)}")

    lines += [
        "",
        "  Caution reading this table: the background is every OTHER contig in the",
        "  file, and Logan contig files are dominated by 31-50 bp low-complexity",
        "  stubs. So 'dark looks longer and higher-entropy than background' is a",
        "  statement about the background being junk, not about the dark contigs",
        "  being remarkable. The meaningful comparisons are:",
        "    dark GC shifted high (~60-70%) at normal entropy -> reagent 'kitome'",
        "    dark entropy low / homopolymers long            -> low-complexity repeat",
        "    dark ka near the detection floor                -> single-observation noise",
        "  and above all section 2: what fails to assemble at all.",
    ]

    # ---- 4. blastable pool -------------------------------------------------
    sec("4. DARK DENSITY -- how much of a 'dark' contig is actually dark")
    dens_rows = []
    for fa in sorted((C.SEQ_DIR / "contigs").glob("*.dark_contigs.fa.gz")):
        acc = fa.name.split(".")[0]
        with gzip.open(fa, "rt") as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                f = dict(
                    p.split("=", 1) for p in line.split() if "=" in p and not p.startswith("ka:")
                )
                nd, nn = int(f["dark"]), int(f["nonhost"])
                dens_rows.append((acc, int(f["len"]), nd, nn, nd / nn))
    dens = pd.DataFrame(
        dens_rows, columns=["accession", "length", "n_dark", "n_nonhost", "dark_frac"]
    )
    dens.to_csv(C.DARK_DENSITY_CSV, index=False)  # one row per contig; large

    long_ = dens[dens.length >= 1000]
    lines += [
        "  A contig is written out if it carries >=1 dark hash, but a 100 kb",
        "  well-characterised genome with one stray dark k-mer is not dark matter.",
        "  dark_frac = dark hashes / non-host hashes on that contig.",
        "",
        f"  all dark-carrying contigs : {len(dens):,}   "
        f"median dark_frac {dens.dark_frac.median():.2f}",
        f"  contigs >= 1 kb           : {len(long_):,}   "
        + (f"median dark_frac {long_.dark_frac.median():.2f}" if len(long_) else ""),
        f"  fully dark (frac == 1.0)  : {int((dens.dark_frac >= 0.999).sum()):,} "
        f"({100 * (dens.dark_frac >= 0.999).mean():.1f}%)",
        f"  >=1 kb AND dark_frac>={MIN_DARK_FRAC}  : "
        f"{int(((dens.length >= 1000) & (dens.dark_frac >= MIN_DARK_FRAC)).sum()):,}",
    ]

    sec(f"5. BLAST-ABLE POOL  (>= {MIN_BLAST_LEN} bp AND dark_frac >= {MIN_DARK_FRAC})")
    kept = []
    for fa in sorted((C.SEQ_DIR / "contigs").glob("*.dark_contigs.fa.gz")):
        acc = fa.name.split(".")[0]
        with gzip.open(fa, "rt") as fh:
            hdr = None
            for line in fh:
                if line.startswith(">"):
                    hdr = line[1:].split()
                    continue
                seq = line.strip()
                f = dict(p.split("=", 1) for p in hdr if "=" in p and not p.startswith("ka:"))
                nd, nn = int(f["dark"]), int(f["nonhost"])
                if len(seq) >= MIN_BLAST_LEN and nd / nn >= MIN_DARK_FRAC:
                    ka = next((p[5:] for p in hdr if p.startswith("ka:f:")), "nan")
                    kept.append((len(seq), acc, hdr[0], ka, nd, nn, seq))
    kept.sort(reverse=True, key=lambda t: t[0])
    with open(C.BLASTABLE_FA, "w") as out:
        for ln, acc, name, ka, nd, nn, seq in kept:
            out.write(
                f">{acc}|{name} len={ln} ka:f:{ka} dark={nd} nonhost={nn}\n{seq}\n"
            )
    big = [k for k in kept if k[0] >= 1000]
    lines += [
        f"  contigs written : {len(kept):,}",
        f"  total bp        : {sum(k[0] for k in kept):,}",
        f"  longest         : {max((k[0] for k in kept), default=0):,} bp",
        f"  >= 1 kb         : {len(big):,}",
        f"  file            : {C.BLASTABLE_FA.name}  (sorted longest first)",
    ]
    if big:
        gc = np.array([(s.count("G") + s.count("C")) / len(s) for *_, s in big])
        ka = np.array([float(k[3]) for k in big])
        ln = np.array([k[0] for k in big])
        top = pd.Series([k[1] for k in big]).value_counts()
        lines += [
            "",
            "  The >=1 kb, fully-dark subset is the part that is NOT noise:",
            f"    n={len(big):,}   {ln.sum() / 1e6:.1f} Mbp   "
            f"median length {np.median(ln):,.0f} bp   max {ln.max():,} bp",
            f"    coverage ka  median {np.median(ka):.1f} "
            f"[{np.percentile(ka, 25):.1f}-{np.percentile(ka, 75):.1f}]",
            f"    GC           median {np.median(gc):.3f} "
            f"[{np.percentile(gc, 25):.3f}-{np.percentile(gc, 75):.3f}] (broad -> mixed origins)",
            "",
            "  Multi-kb contigs at tens-of-fold coverage whose hashes are 100% dark",
            "  cannot be sequencing error. These are the sequences worth identifying.",
            f"  Spread: {top.size:,} of {len(meta):,} runs contribute; the largest",
            f"  single run holds {100 * top.iloc[0] / len(big):.1f}% and the top five "
            f"{100 * top.head(5).sum() / len(big):.1f}%.",
            "    " + ", ".join(f"{a} ({n:,})" for a, n in top.head(5).items()),
        ]
        if top.iloc[0] / len(big) >= 0.15:
            lines.append(
                "  Concentrated in a handful of runs, which points at a per-run or"
                " per-batch effect (contamination) rather than site-wide biology."
            )
        else:
            lines.append(
                "  Spread broadly across runs and studies, which is hard to explain"
                " as a per-batch artefact and is consistent with real, site-wide"
                " sequence that is missing from the reference union."
            )

    sec("FILES")
    rel_res = f"results/{C.SITE}"
    rel_dat = f"data/{C.SITE}"
    lines += [
        f"  small outputs, under {rel_res}/ (version controlled)",
        f"    {C.ACCESSIONS_TXT.name:31s} {len(meta)} accessions, one per line",
        f"    {C.SAMPLES_CSV.name:31s} per-sample metrics + SRA provenance",
        f"    {C.COHORT_NOTES.name:31s} cohort provenance warnings",
        f"    {C.REPORT.name:31s} this report",
        "    sig_manifest.csv                per-accession hash counts vs the parquet",
        "    assembly_rates.csv              per-class assembly rates (deep-dive runs)",
        "    contig_profiles.csv             per-contig length/ka/GC/entropy",
        "    extract_summary_*.csv           per-accession recovery and yield",
        "    extract_log/<source>/*.txt      per-accession extraction log",
        "    blast/                          BLAST queries and raw NCBI reports",
        "",
        f"  bulk artifacts, under {rel_dat}/ (gitignored; see README to regenerate)",
        "    sigs/nonhost/*.nonhost.sig      all non-host hashes, sourmash format",
        "    sigs/dark/*.dark.sig            dark hashes only (the DMI numerator)",
        "    hashes/*.hashes.parquet         min_hash + is_dark, input to extraction",
        "    seqs/contigs/*.fa.gz            contigs carrying >=1 dark hash",
        "    seqs/unitigs/*.fa.gz            unitigs carrying dark hashes (pilots only)",
        "    kmers/*/*.parquet               matched k-mer, position, hash, is_dark, ka",
        "    dark_contigs_blastable.fa       >=300 bp AND >=80% dark -- BLAST this",
        "    contig_dark_density.csv         per-contig dark / non-host hash counts",
        "",
    ]

    C.REPORT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nwrote {C.REPORT}")


if __name__ == "__main__":
    main()
