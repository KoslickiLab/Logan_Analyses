#!/usr/bin/env python3
"""
Parse an NCBI BLAST plain-text report and quantify HOW homology is distributed
along each query.

Motivation: a best-HSP alignment covering 0.6% of a 22 kb contig does not mean
"0.6% similar". It means homology is confined to a short window -- typically one
conserved gene -- while the rest of the contig has no detectable relative. That
pattern is what you would see for a conserved/housekeeping locus, or for a
horizontally transferred element, sitting in otherwise unplaceable sequence.
Distinguishing those requires knowing whether the homologous windows are few and
localised or many and scattered, which is what this script measures.

Definitions used here (stated explicitly because BLAST reports several):
  HSP            a single local alignment (one "Score = ..." block).
  pct_identity   identities / alignment_length * 100 for that HSP, exactly as
                 printed on the "Identities = a/b (c%)" line. The denominator
                 includes gap positions, so it is gapped identity.
  best_hsp_cov   alignment_length of the single highest-scoring HSP / query
                 length. This is NOT NCBI's "Query Cover" column.
  union_cov      fraction of query positions covered by at least one HSP to the
                 single best subject, merging overlaps. Comparable to NCBI's
                 per-subject "Query Cover".
  any_cov        fraction of query positions covered by at least one HSP to ANY
                 subject in the report, merging overlaps across all subjects.

Usage:  python analyze_blast_hsps.py blast_megablast_nt.txt
"""

import re
import sys
from collections import defaultdict


def merge(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    out = [list(intervals[0])]
    for a, b in intervals[1:]:
        if a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def covered(intervals):
    return sum(b - a + 1 for a, b in merge(intervals))


def parse(path):
    text = open(path).read()
    queries = []
    for block in text.split("\nQuery= ")[1:]:
        qname = block.split("\n")[0].strip()
        m = re.search(r"\nLength=(\d+)", block)
        qlen = int(m.group(1)) if m else 0

        # Split into per-subject sections (alignment part only)
        subj_sections = re.split(r"\n> ?", block)[1:]
        subjects = []
        for sec in subj_sections:
            title = sec.split("\n")[0].strip()
            hsps = []
            # each HSP starts at " Score = "
            for hsp in re.split(r"\n\s*Score = ", sec)[1:]:
                idm = re.search(r"Identities = (\d+)/(\d+) \((\d+)%\)", hsp)
                coords = re.findall(r"^Query\s+(\d+)\s+\S+\s+(\d+)", hsp, re.M)
                if not idm or not coords:
                    continue
                starts = [int(a) for a, _ in coords]
                ends = [int(b) for _, b in coords]
                lo, hi = min(starts + ends), max(starts + ends)
                hsps.append(
                    dict(ident=int(idm.group(3)), alen=int(idm.group(2)),
                         nident=int(idm.group(1)), qstart=lo, qend=hi)
                )
            if hsps:
                subjects.append(dict(title=title, hsps=hsps))
        queries.append(dict(name=qname, qlen=qlen, subjects=subjects))
    return queries


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "blast_megablast_nt.txt"
    queries = parse(path)

    print(f"# {path}")
    print(f"{'query':22s} {'qlen':>7s} {'subj':>5s} {'HSPs':>5s} {'bestID':>7s} "
          f"{'bestHSPcov':>10s} {'unionCov':>9s} {'anyCov':>7s} {'nWindows':>8s}")
    for q in queries:
        if not q["subjects"]:
            print(f"{q['name'].split('|')[-1]:22s} {q['qlen']:7d} "
                  f"{0:5d} {0:5d} {'-':>7s} {'-':>10s} {'-':>9s} {'-':>7s} {'-':>8s}")
            continue
        all_hsps = [h for s in q["subjects"] for h in s["hsps"]]
        best = max(all_hsps, key=lambda h: h["alen"])
        best_subj = max(
            q["subjects"], key=lambda s: sum(h["alen"] for h in s["hsps"])
        )
        union = covered([(h["qstart"], h["qend"]) for h in best_subj["hsps"]])
        anyc = merge([(h["qstart"], h["qend"]) for h in all_hsps])
        print(
            f"{q['name'].split('|')[-1]:22s} {q['qlen']:7d} "
            f"{len(q['subjects']):5d} {len(all_hsps):5d} {best['ident']:6d}% "
            f"{100 * best['alen'] / q['qlen']:9.2f}% "
            f"{100 * union / q['qlen']:8.2f}% "
            f"{100 * sum(b - a + 1 for a, b in anyc) / q['qlen']:6.2f}% "
            f"{len(anyc):8d}"
        )

    print("\n# homology windows on the least-covered queries "
          "(query coordinates, merged across all subjects)")
    for q in queries:
        all_hsps = [h for s in q["subjects"] for h in s["hsps"]]
        if not all_hsps:
            continue
        anyc = merge([(h["qstart"], h["qend"]) for h in all_hsps])
        frac = sum(b - a + 1 for a, b in anyc) / q["qlen"]
        if frac < 0.10:
            wins = ", ".join(f"{a}-{b} ({b - a + 1} nt)" for a, b in anyc[:8])
            print(f"  {q['name'].split('|')[-1]:22s} qlen={q['qlen']:6d} "
                  f"covered={100 * frac:.2f}% in {len(anyc)} window(s): {wins}")


if __name__ == "__main__":
    main()
