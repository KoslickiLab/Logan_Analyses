#!/usr/bin/env python3
"""
Step 2: export the non-host hashes of the blood cohort as sourmash signatures.

Does ONE pass over the 606 GB sample-hash DuckDB (a scan costs the same whether
you ask for 2 accessions or 125), left-joining the human host DB and the
reference-union DB, and emits per accession:

  sigs/nonhost/{acc}.nonhost.sig  -- every hash surviving the human host filter
  sigs/dark/{acc}.dark.sig        -- non-host AND absent from the reference union
                                     (the DMI numerator; the interesting subset)
  hashes/{acc}.hashes.parquet     -- min_hash + is_dark, for the extractor
  sig_manifest.csv                -- counts, cross-checked against the parquet

Signatures are built with the sourmash Python API so that `max_hash` and
`md5sum` are computed by sourmash itself -- hand-writing the JSON is exactly
where the max_hash / md5sum gotchas bite.

Both DuckDB files are opened READ-ONLY; they are shared assets.

Usage:  conda activate logan && python 02_export_hashes.py
"""

import sys
import time

import duckdb
import pandas as pd
import sourmash
from sourmash.minhash import MinHash

import config as C


def write_sig(hashes, acc, kind, out_dir):
    """Build a MinHash via the sourmash API and save it as a .sig."""
    mh = MinHash(n=0, ksize=C.KSIZE, scaled=C.SCALED)
    if mh._max_hash != C.MAX_HASH:  # guard against a sourmash behaviour change
        sys.exit(f"ERROR: max_hash {mh._max_hash} != expected {C.MAX_HASH}")
    mh.add_many(hashes)
    if len(mh) != len(hashes):
        sys.exit(
            f"ERROR: {acc}/{kind}: sourmash kept {len(mh)} of {len(hashes)} hashes "
            "(some exceeded max_hash -- wrong scaled?)"
        )
    sig = sourmash.SourmashSignature(mh, name=f"{acc}.{kind}", filename=acc)
    path = out_dir / f"{acc}.{kind}.sig"
    with open(path, "w") as fp:
        sourmash.save_signatures([sig], fp)
    return path


def main():
    accessions = C.ACCESSIONS_TXT.read_text().split()
    if not accessions:
        sys.exit(f"ERROR: no accessions in {C.ACCESSIONS_TXT}; run step 1 first")
    print(f"{len(accessions)} accessions")

    for d in (C.SIG_DIR / "nonhost", C.SIG_DIR / "dark", C.HASH_DIR):
        d.mkdir(parents=True, exist_ok=True)

    # ---- one scan of the big DB ------------------------------------------
    con = duckdb.connect(C.DMI_DB, read_only=True)
    con.execute("SET threads=64")
    con.execute(f"ATTACH '{C.HOST_DB}' AS host (READ_ONLY)")
    con.execute("CREATE TEMP TABLE wanted(acc VARCHAR)")
    con.executemany("INSERT INTO wanted VALUES (?)", [(a,) for a in accessions])

    print("scanning sample_hashes (606 GB, expect a few minutes) ...")
    t0 = time.time()
    df = con.execute(
        """
        SELECT s.sample_id, s.min_hash, (r.hash IS NULL) AS is_dark
        FROM sample_hashes s
        JOIN wanted w              ON s.sample_id = w.acc
        LEFT JOIN host.human_hashes h ON s.min_hash = h.hash
        LEFT JOIN reference_hashes  r ON s.min_hash = r.hash
        WHERE h.hash IS NULL          -- non-host only
        """
    ).fetchdf()
    con.close()
    print(f"  {len(df):,} non-host hashes in {time.time() - t0:.0f}s")

    # ---- expected counts from the parquet, for cross-checking -------------
    par = pd.read_parquet(
        C.INPUT_PARQUET,
        columns=[
            "accession",
            "total_hashes_host_filtered",
            "unmapped_hashes_host_filtered",
        ],
    ).set_index("accession")

    # ---- per-accession outputs -------------------------------------------
    rows, problems = [], []
    for acc, sub in df.groupby("sample_id", sort=False):
        nonhost = sub["min_hash"].to_numpy()
        dark = sub.loc[sub["is_dark"], "min_hash"].to_numpy()

        write_sig(nonhost, acc, "nonhost", C.SIG_DIR / "nonhost")
        write_sig(dark, acc, "dark", C.SIG_DIR / "dark")
        sub[["min_hash", "is_dark"]].to_parquet(
            C.HASH_DIR / f"{acc}.hashes.parquet", index=False
        )

        exp_nh = int(par.at[acc, "total_hashes_host_filtered"])
        exp_dk = int(par.at[acc, "unmapped_hashes_host_filtered"])
        ok = len(nonhost) == exp_nh and len(dark) == exp_dk
        if not ok:
            problems.append(
                f"{acc}: nonhost {len(nonhost)} vs {exp_nh}, dark {len(dark)} vs {exp_dk}"
            )
        rows.append(
            dict(
                accession=acc,
                n_nonhost=len(nonhost),
                n_dark=len(dark),
                n_reference_matched=len(nonhost) - len(dark),
                expected_nonhost=exp_nh,
                expected_dark=exp_dk,
                matches_parquet=ok,
            )
        )

    man = pd.DataFrame(rows).sort_values("accession")
    man.to_csv(C.SIG_MANIFEST, index=False)

    missing = sorted(set(accessions) - set(man["accession"]))
    if missing:
        problems.append(f"{len(missing)} accessions absent from sample_hashes: {missing[:5]}")

    print(f"\nwrote {len(man)} signature pairs to {C.SIG_DIR}")
    print(f"  total non-host hashes : {man.n_nonhost.sum():,}")
    print(f"  total dark hashes     : {man.n_dark.sum():,} "
          f"({100 * man.n_dark.sum() / man.n_nonhost.sum():.1f}%)")
    print(f"wrote {C.SIG_MANIFEST}")

    if problems:
        print("\nCOUNT MISMATCHES:", file=sys.stderr)
        for p in problems:
            print("  " + p, file=sys.stderr)
        sys.exit("ERROR: exported counts disagree with the parquet")
    print("\nall counts reconcile with the parquet")


if __name__ == "__main__":
    main()
