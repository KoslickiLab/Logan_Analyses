#!/usr/bin/env python3
"""
Count distinct min_hash values for a set of accessions in the Logan DuckDB.

Two-phase parallel approach:

  Phase 1 (parallel) — T worker processes each query a disjoint batch of
    accessions.  Each worker fetches its hashes as a numpy array, partitions
    them by (unsigned_hash % B) into B binary temp files, and exits.
    Because all workers read the same DB file simultaneously, the OS page
    cache serves each physical page only once from SSD and then satisfies
    all parallel reads from RAM -- effective I/O is one logical scan.

  Phase 2 (sequential) — For each of the B buckets, load all T workers'
    slice for that bucket, run numpy.unique, accumulate the count, and
    delete the files.  Peak memory per bucket is (total_hashes / B * 8)
    bytes -- tiny even for billions of hashes.

Cleanup: all temp files are deleted in Phase 2.  On interrupt, leftover
  files in --tmp-dir can be removed with: rm -f <tmp-dir>/b????_w????.bin
"""

import argparse
import csv
import os
import sys
import numpy as np
import duckdb
from multiprocessing import Pool


def worker_fn(args):
    wid, accs_batch, tmp_dir, db_path, ksize, B = args

    # 1 DuckDB thread per worker; parallelism comes from T workers, not
    # from each worker spawning more threads.
    con = duckdb.connect(config={"threads": 1})
    con.execute(f"ATTACH '{db_path}' AS logan (READ_ONLY)")
    con.execute("CREATE TEMP TABLE accs AS SELECT unnest(?::VARCHAR[]) AS acc",
                [accs_batch])

    result = con.execute(f"""
        SELECT s.min_hash AS min_hash
        FROM   logan.sigs_dna.signature_mins s
        JOIN   accs ON s.sample_id = accs.acc
        WHERE  s.ksize = {ksize}
    """).fetchnumpy()

    hashes = result["min_hash"]   # numpy int64 array -- zero-copy from DuckDB
    n = len(hashes)
    print(f"  worker {wid:>4d}: {n:>13,} hashes, partitioning into {B} buckets ...",
          flush=True, file=sys.stderr)

    # Partition hashes into B buckets using a counting-sort-style approach:
    #   1. compute bucket index for every hash (O(N), vectorised)
    #   2. argsort by bucket (O(N log N) at numpy C speed)
    #   3. write each contiguous bucket slice to its binary file
    ubuckets = hashes.view(np.uint64) % np.uint64(B)   # uint64 for correct modulo
    order    = np.argsort(ubuckets, kind="stable")
    shashes  = hashes[order]
    sbuckets = ubuckets[order]

    counts  = np.bincount(sbuckets.astype(np.intp), minlength=B)
    offsets = np.empty(B + 1, dtype=np.intp)
    offsets[0] = 0
    np.cumsum(counts, out=offsets[1:])

    for b in range(B):
        lo, hi = int(offsets[b]), int(offsets[b + 1])
        if hi > lo:
            shashes[lo:hi].tofile(
                os.path.join(tmp_dir, f"b{b:04d}_w{wid:04d}.bin")
            )

    print(f"  worker {wid:>4d}: done.", flush=True, file=sys.stderr)
    return n


def main():
    ap = argparse.ArgumentParser(
        description="Count distinct min_hash values (parallel extraction + bucketed dedup)."
    )
    ap.add_argument("--db",      required=True, help="Path to DuckDB database (read-only)")
    ap.add_argument("--summary", required=True, help="summary.tsv with 'acc' column")
    ap.add_argument("-t", "--threads", type=int, default=min(32, os.cpu_count()),
                    help="Parallel worker processes (default: min(32, cpu_count))")
    ap.add_argument("-b", "--buckets", type=int, default=256,
                    help="Hash-space partitions for dedup phase (default: 256)")
    ap.add_argument("-d", "--tmp-dir", default="/dev/shm",
                    help="Temp dir for intermediate files (default: /dev/shm)")
    ap.add_argument("-k", "--ksize", type=int, default=31)
    args = ap.parse_args()

    T, B = args.threads, args.buckets

    with open(args.summary) as f:
        accs = [row["acc"] for row in csv.DictReader(f, delimiter="\t")]
    print(f"Accessions: {len(accs):,}   Workers: {T}   Buckets: {B}   Tmp: {args.tmp_dir}",
          file=sys.stderr)
    os.makedirs(args.tmp_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Phase 1: parallel extraction and partitioning
    # ------------------------------------------------------------------
    print("Phase 1: extracting hashes (parallel) ...", file=sys.stderr)

    batches   = [accs[i::T] for i in range(T)]
    work_args = [(i, batches[i], args.tmp_dir, args.db, args.ksize, B)
                 for i in range(T)]

    with Pool(T) as pool:
        counts = pool.map(worker_fn, work_args)

    total_with_dups = sum(counts)
    print(f"Phase 1 done.  Total hashes (with duplicates): {total_with_dups:,}",
          file=sys.stderr)

    # ------------------------------------------------------------------
    # Phase 2: per-bucket deduplication (sequential; each bucket is small)
    # ------------------------------------------------------------------
    print("Phase 2: deduplicating per bucket ...", file=sys.stderr)
    total_distinct = 0

    for b in range(B):
        parts = []
        for wid in range(T):
            p = os.path.join(args.tmp_dir, f"b{b:04d}_w{wid:04d}.bin")
            if os.path.exists(p):
                parts.append(np.fromfile(p, dtype=np.int64))
                os.remove(p)

        if parts:
            merged = np.concatenate(parts) if len(parts) > 1 else parts[0]
            total_distinct += len(np.unique(merged))

        if (b + 1) % 16 == 0 or b == B - 1:
            pct = 100 * (b + 1) / B
            print(f"\r  {b+1:>{len(str(B))}}/{B} buckets ({pct:5.1f}%)  "
                  f"distinct so far: {total_distinct:,}   ",
                  end="", flush=True, file=sys.stderr)

    print(f"\r{'':80}", file=sys.stderr)
    print(total_distinct)


if __name__ == "__main__":
    main()
