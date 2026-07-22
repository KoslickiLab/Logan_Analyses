#!/usr/bin/env python3
"""
Step 4: pull the Logan sequences that carry each sample's DARK hashes
(non-host AND matching nothing in the reference union -- the DMI numerator).

Two sources, same extractor:

  --source contigs   s3://logan-pub/c/{acc}/{acc}.contigs.fa.zst   (~11 MB/sample)
      The BLAST-able deliverable. Assembled, up to ~8.6 kb, but reaches only
      a small fraction of the dark hashes.

  --source unitigs   s3://logan-pub/u/{acc}/{acc}.unitigs.fa.zst   (~400 MB/sample)
      Complete (recovery must be ~100%) but mean record length is ~41 bp.
      Used to quantify how much dark matter fails to assemble at all.

Nothing is staged to disk: each file is streamed straight from S3 through zstd.

Outputs per accession
  seqs/{source}/{acc}.dark_{source}.fa.gz
  kmers/{source}/{acc}.matched_kmers.parquet
  extract_log/{source}/{acc}.txt

Usage:
  python 04_extract_sequences.py --source contigs --jobs 32
  python 04_extract_sequences.py --source unitigs --accessions SRR14209595 ...
"""

import argparse
import gzip
import io
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

import config as C
import extractor

SOURCES = {
    "contigs": (C.S3_CONTIGS, "contigs"),
    "unitigs": (C.S3_UNITIGS, "unitigs"),
}


def open_logan_stream(url):
    """aws s3 cp <url> - | zstd -dc, as a text stream."""
    p1 = subprocess.Popen(
        ["aws", "s3", "cp", url, "-", "--no-sign-request"],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    p2 = subprocess.Popen(
        ["zstd", "-dcq"],
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    p1.stdout.close()
    return p1, p2, io.TextIOWrapper(p2.stdout, encoding="ascii")


def process_one(acc, source):
    url, tag = SOURCES[source]
    url = url.format(acc=acc)

    # Search with the FULL non-host set, then classify each matched hash as
    # dark or reference-matched. A contig is only written out if it carries at
    # least one dark hash, but knowing its reference-matched hashes too is what
    # lets step 5 compute dark DENSITY -- without it, a 100 kb well-known genome
    # carrying one stray dark k-mer looks like a headline discovery.
    hashes = pd.read_parquet(C.HASH_DIR / f"{acc}.hashes.parquet")
    nonhost = hashes["min_hash"].to_numpy(dtype=np.uint64)
    dark_set = set(hashes.loc[hashes["is_dark"], "min_hash"].tolist())
    n_dark = len(dark_set)
    n_nonhost = len(nonhost)

    fa_path = C.SEQ_DIR / tag / f"{acc}.dark_{tag}.fa.gz"
    km_path = C.KMER_DIR / tag / f"{acc}.matched_kmers.parquet"
    log_path = C.LOG_DIR / tag / f"{acc}.txt"
    for p in (fa_path, km_path, log_path):
        p.parent.mkdir(parents=True, exist_ok=True)

    seen_dark = set()
    seen_nonhost = set()
    rows = []
    n_seq_out = 0
    bp_out = 0
    scanned = [0, 0]

    def progress(n_rec, n_bp):
        scanned[0], scanned[1] = n_rec, n_bp

    t0 = time.time()
    p1, p2, stream = open_logan_stream(url)
    try:
        with gzip.open(fa_path, "wt") as out:
            for name, ka, seq, matches in extractor.extract(
                stream, nonhost, C.MAX_HASH, progress=progress
            ):
                flagged = [(pos, kmer, hv, hv in dark_set) for pos, kmer, hv in matches]
                seen_nonhost.update(hv for _, _, hv, _ in flagged)
                n_dk = sum(1 for *_, d in flagged if d)
                if n_dk == 0:
                    continue  # reference-matched only; not dark matter
                seen_dark.update(hv for _, _, hv, d in flagged if d)
                n_seq_out += 1
                bp_out += len(seq)
                ka_s = "nan" if ka != ka else f"{ka:g}"
                out.write(
                    f">{name} ka:f:{ka_s} len={len(seq)} dark={n_dk} "
                    f"nonhost={len(flagged)}\n{seq}\n"
                )
                for pos, kmer, hv, is_dark in flagged:
                    rows.append((name, pos, kmer, hv, is_dark, len(seq), ka))
    finally:
        stream.close()
        p2.wait()
        p1.wait()

    if p1.returncode != 0 or p2.returncode != 0:
        raise RuntimeError(
            f"{acc}: stream failed (aws={p1.returncode} zstd={p2.returncode})"
        )
    if scanned[0] == 0:
        raise RuntimeError(f"{acc}: no records read from {url}")

    pd.DataFrame(
        rows,
        columns=["seq_id", "pos", "kmer", "hashval", "is_dark", "seq_len", "seq_ka"],
    ).to_parquet(km_path, index=False)

    recovery = len(seen_dark) / n_dark if n_dark else float("nan")
    elapsed = time.time() - t0
    stats = dict(
        accession=acc,
        source=tag,
        n_dark_hashes=n_dark,
        n_dark_recovered=len(seen_dark),
        recovery_frac=recovery,
        n_nonhost_hashes=n_nonhost,
        n_nonhost_recovered=len(seen_nonhost),
        nonhost_recovery_frac=len(seen_nonhost) / n_nonhost if n_nonhost else float("nan"),
        n_records_scanned=scanned[0],
        bp_scanned=scanned[1],
        n_records_matched=n_seq_out,
        bp_matched=bp_out,
        seconds=round(elapsed, 1),
    )
    log_path.write_text(
        "\n".join(f"{k:20s} {v}" for k, v in stats.items()) + "\n"
    )
    return stats


def _worker(args):
    acc, source = args
    try:
        return process_one(acc, source)
    except Exception as exc:  # keep one bad accession from killing the run
        return dict(accession=acc, source=source, error=f"{type(exc).__name__}: {exc}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=sorted(SOURCES), required=True)
    ap.add_argument("--accessions", nargs="*", help="default: all in blood_accessions.txt")
    ap.add_argument("--jobs", type=int, default=16)
    args = ap.parse_args()

    accs = args.accessions or C.ACCESSIONS_TXT.read_text().split()
    missing = [a for a in accs if not (C.HASH_DIR / f"{a}.hashes.parquet").exists()]
    if missing:
        sys.exit(f"ERROR: run step 2 first; missing hashes for {missing[:5]}")

    print(f"extracting dark {args.source} for {len(accs)} accessions, {args.jobs} jobs")
    results = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futs = {pool.submit(_worker, (a, args.source)): a for a in accs}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            results.append(r)
            if "error" in r:
                print(f"  [{i}/{len(accs)}] {r['accession']} ERROR {r['error']}",
                      file=sys.stderr)
            else:
                print(f"  [{i}/{len(accs)}] {r['accession']} "
                      f"recovered {r['n_dark_recovered']:,}/{r['n_dark_hashes']:,} "
                      f"({100 * r['recovery_frac']:.2f}%) "
                      f"-> {r['n_records_matched']:,} records, "
                      f"{r['bp_matched'] / 1e6:.2f} Mbp, {r['seconds']}s")

    df = pd.DataFrame(results)
    summary = C.BASE / f"extract_summary_{args.source}.csv"
    df.to_csv(summary, index=False)

    if "error" in df.columns:
        errs, ok = df[df["error"].notna()], df[df["error"].isna()]
    else:
        errs, ok = df.iloc[0:0], df
    print(f"\ndone in {time.time() - t0:.0f}s; wrote {summary}")
    if len(ok):
        print(f"  accessions ok        : {len(ok)}")
        print(f"  dark hashes recovered: {ok.n_dark_recovered.sum():,} / "
              f"{ok.n_dark_hashes.sum():,} "
              f"({100 * ok.n_dark_recovered.sum() / ok.n_dark_hashes.sum():.2f}%)")
        print(f"  sequence written     : {ok.n_records_matched.sum():,} records, "
              f"{ok.bp_matched.sum() / 1e6:.1f} Mbp")
    if len(errs):
        print(f"  FAILED               : {len(errs)}", file=sys.stderr)

    # Unitigs are the complete source: anything much below 100% means the
    # sketches and the unitig files disagree, which would invalidate the run.
    if args.source == "unitigs" and len(ok):
        worst = ok["recovery_frac"].min()
        if worst < 0.999:
            sys.exit(
                f"ERROR: unitig recovery dipped to {100 * worst:.2f}% "
                "-- sketch and unitigs disagree"
            )
        print("  unitig recovery ~100% as required")


if __name__ == "__main__":
    main()
