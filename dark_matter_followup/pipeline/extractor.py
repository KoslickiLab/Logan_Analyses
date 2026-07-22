"""
Fast extraction of the FASTA records carrying a given set of FracMinHash hashes.

Why not `sourmash sig kmers`
----------------------------
`sig kmers` does copy_and_clear() + add_sequence() + two intersection() calls
*per FASTA record* (sourmash/sig/__main__.py), and with --save-kmers it tests
`hashval in query_mh.hashes` per k-mer, rebuilding a dict wrapper each time.
A Logan unitig file has ~14.1M records, so that is ~1.8e9 Python round-trips.
It also cannot read .zst and cannot annotate the matches.

How this works instead
----------------------
Concatenate ~5 Mbp of records into one string joined by runs of 'N', and make a
single Rust call:

    MinHash(ksize=31, scaled=1).seq_to_hashes(batch, force=True,
                                              bad_kmers_as_zeroes=True)

That returns exactly one hash per k-mer start position, with invalid k-mers
(anything touching an 'N', i.e. every junction k-mer) reported as 0. Positions
therefore map back to records by pure arithmetic, and everything after the Rust
call is vectorised numpy. ~120 FFI calls per unitig file instead of 14.1M.

The hashing itself is sourmash's own, so it cannot drift from the sketches in
the database. (Verified equivalent to
`mmh3.hash64(min(kmer, revcomp(kmer)), 42, signed=False)[0]` in 03_validate_hashing.py.)
"""

import numpy as np
from sourmash.minhash import MinHash

KSIZE = 31
SEPARATOR = "N" * KSIZE  # >=1 would do; any spanning k-mer must contain an N
DEFAULT_BATCH_BP = 5_000_000


def parse_fasta(fh):
    """Stream (name, ka, sequence) from a Logan unitig/contig FASTA.

    Headers look like:  >SRR14209595_0 ka:f:32.0   L:+:1:+ L:-:9:-
    `ka` is the mean k-mer abundance (coverage) of the record.
    """
    name, ka, chunks = None, float("nan"), []
    for line in fh:
        if line.startswith(">"):
            if name is not None:
                yield name, ka, "".join(chunks)
            chunks = []
            head = line[1:].rstrip()
            parts = head.split()
            name = parts[0]
            ka = float("nan")
            for p in parts[1:]:
                if p.startswith("ka:f:"):
                    try:
                        ka = float(p[5:])
                    except ValueError:
                        pass
                    break
        else:
            chunks.append(line.strip())
    if name is not None:
        yield name, ka, "".join(chunks)


def _isin_sorted(values, sorted_targets):
    """Vectorised membership test against a sorted uint64 array."""
    if len(sorted_targets) == 0:
        return np.zeros(len(values), dtype=bool)
    idx = np.searchsorted(sorted_targets, values)
    np.clip(idx, 0, len(sorted_targets) - 1, out=idx)
    return sorted_targets[idx] == values


def _flush_batch(names, kas, seqs, offsets, hasher, sorted_targets, max_hash):
    """Hash one concatenated batch and return the records that matched.

    Yields (name, ka, seq, [(pos_in_seq, kmer, hashval), ...]).
    """
    batch = SEPARATOR.join(seqs)
    raw = hasher.seq_to_hashes(
        batch, force=True, bad_kmers_as_zeroes=True, is_protein=False
    )
    arr = np.array(raw, dtype=np.uint64)

    # Invalid k-mers are 0 and would otherwise pass the max_hash test.
    cand = np.nonzero((arr != 0) & (arr < max_hash))[0]
    if cand.size == 0:
        return
    vals = arr[cand]
    keep = _isin_sorted(vals, sorted_targets)
    if not keep.any():
        return
    pos = cand[keep]
    vals = vals[keep]

    starts = np.asarray(offsets, dtype=np.int64)
    rec = np.searchsorted(starts, pos, side="right") - 1

    out = {}
    for p, v, r in zip(pos.tolist(), vals.tolist(), rec.tolist()):
        out.setdefault(r, []).append((p - starts[r], batch[p : p + KSIZE], v))
    for r in sorted(out):
        yield names[r], kas[r], seqs[r], out[r]


def extract(fh, target_hashes, max_hash, batch_bp=DEFAULT_BATCH_BP, progress=None):
    """Stream a FASTA and yield records containing at least one target hash.

    Yields (name, ka, seq, matches) where matches is a list of
    (position_in_record, kmer, hashval).
    """
    sorted_targets = np.sort(np.asarray(target_hashes, dtype=np.uint64))
    hasher = MinHash(n=0, ksize=KSIZE, scaled=1)

    names, kas, seqs, offsets = [], [], [], []
    cursor, n_rec, n_bp = 0, 0, 0

    for name, ka, seq in parse_fasta(fh):
        seq = seq.upper()
        n_rec += 1
        n_bp += len(seq)
        names.append(name)
        kas.append(ka)
        seqs.append(seq)
        offsets.append(cursor)
        cursor += len(seq) + len(SEPARATOR)

        if cursor >= batch_bp:
            yield from _flush_batch(
                names, kas, seqs, offsets, hasher, sorted_targets, max_hash
            )
            names, kas, seqs, offsets = [], [], [], []
            cursor = 0
            if progress:
                progress(n_rec, n_bp)

    if names:
        yield from _flush_batch(
            names, kas, seqs, offsets, hasher, sorted_targets, max_hash
        )
    if progress:
        progress(n_rec, n_bp)


def scan_stats(fh):
    """Count records/bp without hashing (used for logging denominators)."""
    n, bp = 0, 0
    for _, _, seq in parse_fasta(fh):
        n += 1
        bp += len(seq)
    return n, bp
