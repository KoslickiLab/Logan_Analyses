#!/usr/bin/env python3
"""
Minimal NCBI BLAST URL-API client (no local BLAST install / no nt download).

Follows NCBI etiquette: identifies tool+email, submits one job at a time, and
polls no more often than every 60 s.

Usage:
  python ncbi_blast.py --query top_dark_contigs.fa --program blastn \
      --megablast --db nt --out blast_megablast.tsv --max-seqs 15
"""

import argparse
import sys
import time
import urllib.parse
import urllib.request

URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
TOOL = "logan-dmi-blood-followup"
EMAIL = "dmkoslicki@gmail.com"

TAB_FIELDS = (
    "qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
    "evalue bitscore staxids sscinames scomnames stitle"
)


def post(params):
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(URL, data=data, headers={"User-Agent": TOOL})
    with urllib.request.urlopen(req, timeout=180) as r:
        return r.read().decode("utf-8", "replace")


def get(params):
    q = urllib.parse.urlencode(params)
    req = urllib.request.Request(f"{URL}?{q}", headers={"User-Agent": TOOL})
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.read().decode("utf-8", "replace")


def read_fasta(path, max_seqs, max_len):
    out, name, chunks = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                out.append((name, "".join(chunks)))
            name, chunks = line[1:].strip(), []
        else:
            chunks.append(line.strip())
    if name:
        out.append((name, "".join(chunks)))
    out = out[:max_seqs]
    return [(n, s[:max_len]) for n, s in out]


def submit(fasta, program, db, megablast, hitlist, evalue, task):
    params = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": db,
        "QUERY": fasta,
        "HITLIST_SIZE": str(hitlist),
        "EXPECT": str(evalue),
        "FORMAT_TYPE": "Text",
        "tool": TOOL,
        "email": EMAIL,
    }
    if megablast:
        params["MEGABLAST"] = "on"
    if task:
        params["BLAST_PROGRAMS"] = task
    resp = post(params)
    rid = rtoe = None
    for line in resp.splitlines():
        line = line.strip()
        if line.startswith("RID = "):
            rid = line[6:].strip()
        elif line.startswith("RTOE = "):
            rtoe = int(line[7:].strip())
    if not rid:
        sys.exit(f"no RID returned; response head:\n{resp[:1500]}")
    return rid, rtoe or 30


def wait(rid, poll=60, limit=5400):
    t0 = time.time()
    while time.time() - t0 < limit:
        time.sleep(poll)
        info = get({"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid,
                    "tool": TOOL, "email": EMAIL})
        status = ""
        for line in info.splitlines():
            if "Status=" in line:
                status = line.split("Status=")[1].strip()
        print(f"    [{time.time() - t0:5.0f}s] {rid} {status or '?'}", flush=True)
        if status == "READY":
            return True
        if status in ("FAILED", "UNKNOWN"):
            print(f"    job {status}", file=sys.stderr)
            return False
    print("    timed out", file=sys.stderr)
    return False


def fetch(rid, fmt="Tabular"):
    params = {"CMD": "Get", "RID": rid, "FORMAT_TYPE": fmt,
              "tool": TOOL, "email": EMAIL}
    if fmt == "Tabular":
        params["ALIGNMENT_VIEW"] = "Tabular"
    return get(params)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True)
    ap.add_argument("--program", default="blastn")
    ap.add_argument("--db", default="nt")
    ap.add_argument("--megablast", action="store_true")
    ap.add_argument("--task", default=None, help="e.g. discoMegablast, blastn")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-seqs", type=int, default=15)
    ap.add_argument("--max-len", type=int, default=30000)
    ap.add_argument("--hitlist", type=int, default=10)
    ap.add_argument("--evalue", type=float, default=1e-3)
    args = ap.parse_args()

    seqs = read_fasta(args.query, args.max_seqs, args.max_len)
    fasta = "".join(f">{n}\n{s}\n" for n, s in seqs)
    print(f"submitting {len(seqs)} sequences ({len(fasta):,} bytes) "
          f"to {args.program}/{args.db}"
          f"{' megablast' if args.megablast else ''}{' ' + args.task if args.task else ''}")

    rid, rtoe = submit(fasta, args.program, args.db, args.megablast,
                       args.hitlist, args.evalue, args.task)
    print(f"  RID={rid}  estimated {rtoe}s")
    if not wait(rid):
        sys.exit(1)
    text = fetch(rid, "Text")
    open(args.out, "w").write(text)
    print(f"  wrote {args.out} ({len(text):,} bytes)")


if __name__ == "__main__":
    main()
