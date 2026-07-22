#!/usr/bin/env python3
"""
Predict protein structures with ESMFold, then look for remote homologs by
STRUCTURE with Foldseek.

Why: for proteins this divergent, sequence search (blastp, CD-Search) often finds
nothing, because sequence identity decays far faster than fold. Structure is the
last resort that still carries homology signal at <20% identity.

Services used (no local install, no GPU):
  ESMFold   POST https://api.esmatlas.com/foldSequence/v1/pdb/   (raw sequence body)
  Foldseek  https://search.foldseek.com/api/ticket -> /ticket/<id> -> /result/<id>/0

ESMFold is single-sequence (no MSA), so confidence is lower than AlphaFold2 for
shallow families. Mean pLDDT is reported per model; treat anything below ~70 as
low confidence and anything below ~50 as unreliable.

Usage:
  python fold_and_search.py --query fold_targets.faa --outdir structures \
      [--databases afdb50 pdb100] [--skip-foldseek]
"""

import argparse
import json
import os
import time
import urllib.request
import urllib.error

ESMFOLD = "https://api.esmatlas.com/foldSequence/v1/pdb/"
FOLDSEEK = "https://search.foldseek.com/api"
UA = "logan-dmi-followup"


def read_fasta(path):
    out, name, buf = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def esmfold(seq, retries=3):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                ESMFOLD, data=seq.encode(),
                headers={"User-Agent": UA, "Content-Type": "text/plain"},
            )
            with urllib.request.urlopen(req, timeout=600) as r:
                return r.read().decode()
        except Exception as exc:
            if attempt == retries - 1:
                return f"ERROR: {exc}"
            time.sleep(20)


def mean_plddt(pdb_text):
    """ESMFold writes per-residue pLDDT into the B-factor column."""
    vals, seen = [], set()
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            resi = line[22:27]
            if resi in seen:
                continue
            seen.add(resi)
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                pass
    return sum(vals) / len(vals) if vals else float("nan")


def foldseek_submit(pdb_path, databases):
    """Multipart upload; returns a ticket id."""
    boundary = "----logandmi"
    parts = []
    for db in databases:
        parts.append(
            f"--{boundary}\r\nContent-Disposition: form-data; name=\"database[]\"\r\n\r\n{db}\r\n"
        )
    parts.append(
        f"--{boundary}\r\nContent-Disposition: form-data; name=\"mode\"\r\n\r\n3diaa\r\n"
    )
    with open(pdb_path) as fh:
        content = fh.read()
    parts.append(
        f"--{boundary}\r\nContent-Disposition: form-data; name=\"q\"; filename=\"q.pdb\"\r\n"
        f"Content-Type: chemical/x-pdb\r\n\r\n{content}\r\n"
    )
    parts.append(f"--{boundary}--\r\n")
    body = "".join(parts).encode()
    req = urllib.request.Request(
        f"{FOLDSEEK}/ticket", data=body,
        headers={"User-Agent": UA,
                 "Content-Type": f"multipart/form-data; boundary={boundary}"},
    )
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.loads(r.read().decode())


def foldseek_wait(ticket, interval=15, limit=1800):
    t0 = time.time()
    while time.time() - t0 < limit:
        req = urllib.request.Request(f"{FOLDSEEK}/ticket/{ticket}",
                                     headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=120) as r:
            st = json.loads(r.read().decode()).get("status")
        if st == "COMPLETE":
            return True
        if st == "ERROR":
            return False
        time.sleep(interval)
    return False


def foldseek_results(ticket):
    req = urllib.request.Request(f"{FOLDSEEK}/result/{ticket}/0",
                                 headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=300) as r:
        return json.loads(r.read().decode())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--databases", nargs="*",
                    default=["afdb50", "pdb100", "afdb-proteome"])
    ap.add_argument("--skip-foldseek", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    summary = []
    for name, seq in read_fasta(args.query):
        pdb_path = os.path.join(args.outdir, f"{name}.pdb")
        if os.path.exists(pdb_path) and os.path.getsize(pdb_path) > 1000:
            pdb = open(pdb_path).read()
            print(f"{name}: cached")
        else:
            print(f"{name}: folding {len(seq)} aa ...", flush=True)
            pdb = esmfold(seq)
            if pdb.startswith("ERROR"):
                print(f"   {pdb}")
                summary.append((name, len(seq), float("nan"), "fold failed"))
                continue
            open(pdb_path, "w").write(pdb)
        plddt = mean_plddt(pdb)
        print(f"   mean pLDDT {plddt:.1f}")

        if args.skip_foldseek:
            summary.append((name, len(seq), plddt, "not searched"))
            continue
        try:
            tk = foldseek_submit(pdb_path, args.databases)["id"]
            ok = foldseek_wait(tk)
            if not ok:
                summary.append((name, len(seq), plddt, "foldseek failed"))
                continue
            res = foldseek_results(tk)
            top = []
            for db in res.get("results", []):
                for aln in db.get("alignments", [])[:3]:
                    if isinstance(aln, list):
                        aln = aln[0] if aln else {}
                    top.append((db.get("db"), aln.get("target"),
                                aln.get("eval"), aln.get("prob"),
                                aln.get("taxName")))
            json.dump(res, open(os.path.join(args.outdir, f"{name}.foldseek.json"), "w"))
            best = "; ".join(
                f"{d}:{t} E={e} prob={p}" for d, t, e, p, _ in top[:3]) or "no hits"
            print(f"   foldseek: {best}")
            summary.append((name, len(seq), plddt, best))
        except Exception as exc:
            print(f"   foldseek error: {exc}")
            summary.append((name, len(seq), plddt, f"error {exc}"))

    with open(os.path.join(args.outdir, "summary.tsv"), "w") as o:
        o.write("protein\tlength_aa\tmean_plddt\tfoldseek_top\n")
        for row in summary:
            o.write("\t".join(str(x) for x in row) + "\n")
    print(f"\nwrote {args.outdir}/summary.tsv")


if __name__ == "__main__":
    main()
