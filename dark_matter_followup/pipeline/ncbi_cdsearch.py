#!/usr/bin/env python3
"""
Minimal client for NCBI Batch CD-Search (Conserved Domain Database).

Why this and not BLAST: for sequence this divergent, position-specific scoring
matrices built from curated domain alignments detect homology that pairwise
blastp misses. CD-Search is the cheapest route to "what protein family is this"
when the nearest relative in nr is 30% identical.

API: https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi
  POST queries=<fasta> db=cdd ...   -> returns a "#cdsid  QM-..." handle
  POST cdsid=<handle> ...           -> "#status 3" while running, "#status 0" when done

Usage:
  python ncbi_cdsearch.py --query proteins.faa --out cdd_hits.tsv [--evalue 0.01]
"""

import argparse
import sys
import time
import urllib.parse
import urllib.request

URL = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
TOOL = "logan-dmi-followup"
EMAIL = "dmkoslicki@gmail.com"


def post(params):
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(URL, data=data, headers={"User-Agent": TOOL})
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.read().decode("utf-8", "replace")


def submit(fasta, evalue, db, maxhit):
    resp = post({
        "queries": fasta,
        "db": db,                # cdd = CDD + Pfam + SMART + COG + TIGRFAM + NCBIfam
        "smode": "auto",
        "useid1": "true",
        "compbasedadj": "1",
        "filter": "true",
        "evalue": str(evalue),
        "maxhit": str(maxhit),
        "dmode": "full",         # all hits, not just the representative one
        "tdata": "hits",
        "tool": TOOL,
        "email": EMAIL,
    })
    for line in resp.splitlines():
        if line.startswith("#cdsid"):
            return line.split()[-1].strip()
    sys.exit(f"no cdsid returned:\n{resp[:1200]}")


def poll(cdsid, interval=20, limit=3600):
    t0 = time.time()
    while time.time() - t0 < limit:
        resp = post({
            "cdsid": cdsid, "tdata": "hits", "dmode": "full",
            "cddefl": "true", "qdefl": "true", "tool": TOOL, "email": EMAIL,
        })
        status = None
        for line in resp.splitlines():
            if line.startswith("#status"):
                parts = line.split("\t")
                status = parts[1].strip() if len(parts) > 1 else line.split()[-1]
                break
        print(f"    [{time.time() - t0:5.0f}s] {cdsid} status={status}", flush=True)
        if status == "0":
            return resp
        if status not in ("3", None):   # 3 = still running
            sys.exit(f"CD-Search failed, status={status}\n{resp[:1200]}")
        time.sleep(interval)
    sys.exit("CD-Search timed out")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True, help="protein FASTA")
    ap.add_argument("--out", required=True)
    ap.add_argument("--evalue", type=float, default=0.01)
    ap.add_argument("--db", default="cdd")
    ap.add_argument("--maxhit", type=int, default=10)
    args = ap.parse_args()

    fasta = open(args.query).read()
    n = fasta.count(">")
    print(f"submitting {n} proteins ({len(fasta):,} bytes) to CD-Search db={args.db}")
    cdsid = submit(fasta, args.evalue, args.db, args.maxhit)
    print(f"  cdsid={cdsid}")
    resp = poll(cdsid)
    open(args.out, "w").write(resp)
    hits = [l for l in resp.splitlines() if l and not l.startswith("#") and "\t" in l]
    print(f"  wrote {args.out} ({len(hits)} hit rows)")


if __name__ == "__main__":
    main()
