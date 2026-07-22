#!/usr/bin/env python3
"""
Step 3: correctness gate. Step 4 must not run unless every check here passes.

Checks
  1. sourmash's max_hash for scaled=1000 is the value we hard-coded.
  2. A written .sig round-trips, and its md5sum matches an independent
     re-implementation (md5 of ksize followed by the sorted hashes as strings).
  3. On synthetic sequence, our batched seq_to_hashes extractor finds exactly
     the same (position, hash) pairs as a naive per-k-mer
     mmh3.hash64(min(kmer, revcomp(kmer)), 42) loop.
  4. On one real sample's contigs (23 Mbp, small enough that `sourmash sig
     kmers` is tractable), our extractor and `sourmash sig kmers` select an
     identical set of sequence IDs.

Usage:  conda activate logan && python 03_validate_hashing.py
"""

import hashlib
import io
import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path

import mmh3
import numpy as np
import sourmash
from sourmash.minhash import MinHash

import config as C
import extractor

FAILURES = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def naive_hashes(seq):
    """Reference implementation: canonical k-mer -> murmur3_x64_128 lo 64 bits."""
    out = []
    for i in range(len(seq) - C.KSIZE + 1):
        kmer = seq[i : i + C.KSIZE]
        if set(kmer) - set("ACGT"):
            continue
        canon = min(kmer, revcomp(kmer))
        out.append((i, mmh3.hash64(canon, C.SEED, signed=False)[0]))
    return out


# --------------------------------------------------------------------------
def test_max_hash():
    mh = MinHash(n=0, ksize=C.KSIZE, scaled=C.SCALED)
    check(
        "max_hash for scaled=1000",
        mh._max_hash == C.MAX_HASH and mh.seed == C.SEED and mh.ksize == C.KSIZE,
        f"max_hash={mh._max_hash} seed={mh.seed} k={mh.ksize}",
    )


def test_sig_roundtrip():
    sigs = sorted((C.SIG_DIR / "dark").glob("*.dark.sig"))
    if not sigs:
        check("signature round-trip", False, "no .sig files; run step 2 first")
        return
    path = sigs[0]
    raw = json.loads(path.read_text())[0]["signatures"][0]

    # independent md5: md5(str(ksize) + each hash as decimal, in sorted order)
    m = hashlib.md5()
    m.update(str(raw["ksize"]).encode())
    for h in raw["mins"]:
        m.update(str(h).encode())
    check(
        "md5sum recomputed independently",
        m.hexdigest() == raw["md5sum"],
        f"{path.name} {raw['md5sum'][:12]}",
    )
    check(
        "sig header fields",
        raw["ksize"] == C.KSIZE
        and raw["seed"] == C.SEED
        and raw["max_hash"] == C.MAX_HASH
        and raw["mins"] == sorted(raw["mins"]),
        f"k={raw['ksize']} seed={raw['seed']} max_hash={raw['max_hash']} n={len(raw['mins'])}",
    )

    loaded = sourmash.load_file_as_signatures(str(path))
    sig = next(iter(loaded))
    check(
        "sourmash re-loads the sig",
        sig.minhash.scaled == C.SCALED and len(sig.minhash) == len(raw["mins"]),
        f"scaled={sig.minhash.scaled} n={len(sig.minhash)}",
    )


def test_extractor_vs_naive():
    rnd = random.Random(0)
    recs = []
    for i in range(50):
        n = rnd.randint(31, 400)
        recs.append((f"rec{i}", 3.5, "".join(rnd.choice("ACGT") for _ in range(n))))
    # records that must contribute nothing
    recs.append(("short", 1.0, "ACGT" * 5))  # < k
    recs.append(("ambig", 1.0, "ACGTN" * 40))  # every k-mer touches an N

    fasta = "".join(f">{n} ka:f:{k}\n{s}\n" for n, k, s in recs)

    truth = {}  # hashval -> set of (record, position)
    for name, _, seq in recs:
        for pos, h in naive_hashes(seq):
            if h < C.MAX_HASH:
                truth.setdefault(h, set()).add((name, pos))
    targets = np.array(sorted(truth), dtype=np.uint64)

    got = {}
    for name, _ka, _seq, matches in extractor.extract(
        io.StringIO(fasta), targets, C.MAX_HASH, batch_bp=2000
    ):
        for pos, kmer, h in matches:
            got.setdefault(h, set()).add((name, pos))

    check(
        "batched extractor == naive mmh3 loop",
        got == truth,
        f"{len(truth)} hashes, {sum(len(v) for v in truth.values())} placements",
    )

    # k-mer strings reported must actually hash to the reported value
    bad = 0
    for name, _ka, seq, matches in extractor.extract(
        io.StringIO(fasta), targets, C.MAX_HASH, batch_bp=2000
    ):
        for pos, kmer, h in matches:
            if kmer != seq[pos : pos + C.KSIZE]:
                bad += 1
            elif mmh3.hash64(min(kmer, revcomp(kmer)), C.SEED, signed=False)[0] != h:
                bad += 1
    check("reported k-mer strings hash to reported values", bad == 0, f"{bad} bad")


def test_against_sourmash_sig_kmers(acc):
    sig_path = C.SIG_DIR / "dark" / f"{acc}.dark.sig"
    if not sig_path.exists():
        check("vs `sourmash sig kmers`", False, f"missing {sig_path}")
        return

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        fa = td / f"{acc}.contigs.fa"
        print(f"    downloading {acc} contigs ...")
        url = C.S3_CONTIGS.format(acc=acc)
        with open(fa, "wb") as out:
            p1 = subprocess.Popen(
                ["aws", "s3", "cp", url, "-", "--no-sign-request"],
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            )
            p2 = subprocess.Popen(["zstd", "-dc"], stdin=p1.stdout, stdout=out)
            p1.stdout.close()
            p2.wait()
            p1.wait()
        if p2.returncode != 0 or not fa.stat().st_size:
            check("vs `sourmash sig kmers`", False, "download/decompress failed")
            return

        sig = next(iter(sourmash.load_file_as_signatures(str(sig_path))))
        targets = np.array(sorted(sig.minhash.hashes), dtype=np.uint64)

        print("    running our extractor ...")
        with open(fa) as fh:
            ours = {
                name
                for name, _ka, _seq, _m in extractor.extract(fh, targets, C.MAX_HASH)
            }

        print("    running `sourmash sig kmers` (the slow oracle) ...")
        out_fa = td / "oracle.fa"
        r = subprocess.run(
            ["sourmash", "sig", "kmers", "--signatures", str(sig_path),
             "--sequences", str(fa), "--save-sequences", str(out_fa), "-q"],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            check("vs `sourmash sig kmers`", False, r.stderr.strip()[:200])
            return
        theirs = {
            ln[1:].split()[0]
            for ln in out_fa.read_text().splitlines()
            if ln.startswith(">")
        }

    check(
        "our extractor == `sourmash sig kmers` sequence set",
        ours == theirs,
        f"ours={len(ours)} sourmash={len(theirs)} "
        f"only_ours={len(ours - theirs)} only_sourmash={len(theirs - ours)}",
    )


def main():
    acc = sys.argv[1] if len(sys.argv) > 1 else "SRR14209595"
    print("1. sketch parameters")
    test_max_hash()
    print("2. signature format")
    test_sig_roundtrip()
    print("3. extractor vs naive mmh3")
    test_extractor_vs_naive()
    print(f"4. extractor vs `sourmash sig kmers` on real contigs ({acc})")
    test_against_sourmash_sig_kmers(acc)

    print()
    if FAILURES:
        sys.exit(f"GATE FAILED: {len(FAILURES)} check(s) failed: {FAILURES}")
    print("GATE PASSED -- step 4 may run")


if __name__ == "__main__":
    main()
