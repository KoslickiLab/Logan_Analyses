#!/usr/bin/env python3
"""
Count distinct min_hash values for a set of accessions in the Logan DuckDB.

Streams all matching min_hash values in a single table scan and deduplicates
them in a Python set in memory.  At ~50 bytes/element, a billion distinct
hashes needs ~50 GB RAM -- well within the available headroom.
"""

import argparse
import csv
import sys
import duckdb

CHUNK = 10_000_000   # rows fetched per iteration (~80 MB at 8 bytes/hash)


def main():
    p = argparse.ArgumentParser(
        description="Count distinct min_hash values across Logan sketches (single-pass)."
    )
    p.add_argument("--db",      required=True, help="Path to DuckDB database (read-only)")
    p.add_argument("--summary", required=True, help="summary.tsv produced by run_fmh_collisions.sh")
    p.add_argument("-k", "--ksize", type=int, default=31, help="k-mer size filter (default: 31)")
    args = p.parse_args()

    with open(args.summary) as f:
        accs = [row["acc"] for row in csv.DictReader(f, delimiter="\t")]
    print(f"Accessions: {len(accs):,}", file=sys.stderr)

    con = duckdb.connect()
    con.execute(f"ATTACH '{args.db}' AS logan (READ_ONLY)")
    con.execute("CREATE TEMP TABLE accs AS SELECT unnest(?::VARCHAR[]) AS acc", [accs])

    cursor = con.execute(f"""
        SELECT s.min_hash
        FROM   logan.sigs_dna.signature_mins s
        JOIN   accs ON s.sample_id = accs.acc
        WHERE  s.ksize = {args.ksize}
    """)

    seen   = set()
    n_rows = 0
    while chunk := cursor.fetchmany(CHUNK):
        n_rows += len(chunk)
        seen.update(r[0] for r in chunk)
        print(f"\r  rows processed: {n_rows:>15,}   distinct hashes: {len(seen):>15,}",
              end="", flush=True, file=sys.stderr)

    print(f"\r{'':80}", file=sys.stderr)
    print(len(seen))


if __name__ == "__main__":
    main()
