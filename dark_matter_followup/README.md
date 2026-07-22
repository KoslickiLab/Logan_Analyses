# Dark Matter Index follow-up: what is the "uncharacterized" sequence, really?

ASCII-only by convention -- these files are read in terminals and plain text viewers.

## 1. What this is

`../body_site_analysis/body_site_analysis.png` ranks human body sites by **Dark Matter
Index (DMI)**: the fraction of a sample's non-host FracMinHash hashes that match nothing
in a large reference union. Blood came out top (median DMI 0.81), skin fourth (0.46).

This directory answers the obvious next question: **is that dark sequence real novel
biology, or an artefact?** It exports the dark hashes as sourmash signatures, pulls the
actual Logan contigs/unitigs carrying them, characterises them, and BLASTs the longest.

**Headline result: the two sites give opposite answers, for the same DMI.**

| | blood | skin |
|---|---|---|
| runs analysed | 125 | 1,633 |
| bioprojects | 3 | 26 |
| largest project's share of the cohort | 92.8% | 30.9% |
| dark fraction of non-host hashes | 81.6% | 41.4% |
| **dark vs reference-matched assembly rate** | **18.5x worse** | **2.1x worse** |
| fully-dark contigs >= 1 kb | 11,492 (19.7 Mbp) | 152,025 (355.1 Mbp) |
| longest fully-dark contig | 22,182 bp | 124,998 bp |
| fully-dark contigs >= 50 kb | 0 | 49 |
| dark hashes shared across >= 2 bioprojects | 12 (0.002%) | 115,792 (5.35%) |
| what the longest ones BLAST to | soil/freshwater taxa (cyanobacteria, soil actinomycetes) | human commensals + their phages |
| verdict | protocol artefact of one virome study | genuine under-characterisation |

The interpretation, with full methods and caveats, is in
`results/blood/blast_findings.md` and `results/skin/blast_findings.md`.

## 2. Layout

```
dark_matter_followup/
  README.md                  this file
  .gitignore                 excludes data/ and __pycache__
  pipeline/                  all code; site-agnostic
    config.py                paths, sketch parameters, cohort filters, SITES map
    01_select_samples.py     parquet -> cohort accession list + provenance notes
    02_export_hashes.py      DuckDB -> per-run .sig files (nonhost + dark)
    03_validate_hashing.py   correctness gate; must pass before 04
    04_extract_sequences.py  S3 Logan contigs/unitigs -> dark-carrying sequences
    05_characterize.py       assembly rates, sequence character, BLAST-able pool, report
    extractor.py             fast hash->sequence extractor (see section 7)
    ncbi_blast.py            minimal NCBI BLAST URL-API client
    analyze_blast_hsps.py    parses BLAST text reports; HSP coverage statistics
  results/<site>/            small outputs, version controlled (~27 MB total)
    <site>_accessions.txt    the cohort, one accession per line
    <site>_samples.csv       per-sample metrics + flattened SRA attributes
    <site>_cohort_notes.txt  bioproject breakdown + provenance warnings
    <site>_dark_matter_report.txt   generated report, incl. its own METHODS section
    blast_findings.md        hand-written interpretation of the BLAST results
    sig_manifest.csv         per-accession hash counts, cross-checked vs the parquet
    assembly_rates.csv       per-class assembly rates for the deep-dive runs
    contig_profiles.csv      per-contig length/ka/GC/entropy (dark vs background)
    extract_summary_*.csv    per-accession recovery and yield
    extract_log/<source>/    per-accession extraction logs
    blast/                   BLAST query FASTAs and raw NCBI reports
  analyses/                  one-off deep dives that are not part of the per-site pipeline
    tools/                   pip --target installs (pyrodigal); GITIGNORED
    unplaced_skin_contigs/   the two skin contigs that matched nothing in nt
      FINDINGS.md            what they turned out to be, with methods and caveats
      gene_map.csv           per-gene coords, strand, length, best CDD domain
      unplaced_proteins.faa  103 predicted proteins
      cdd_hits.tsv           raw CD-Search output
      blastp_markers.txt     blastp vs nr for the marker domains
      structures/            ESMFold models + gzipped Foldseek results + summary.tsv
  data/<site>/               bulk artefacts, GITIGNORED (~11 GB) -- see section 8
    sigs/{nonhost,dark}/     sourmash signatures
    hashes/                  min_hash + is_dark parquet, input to extraction
    seqs/{contigs,unitigs}/  the dark-carrying sequences themselves
    kmers/                   which k-mer matched where, with is_dark and coverage
    dark_contigs_blastable.fa   the filtered pool actually worth BLASTing
    contig_dark_density.csv     one row per dark-carrying contig
```

`blood` is not special: it is just the first site processed. Every script picks its
cohort from the `DMI_SITE` environment variable (default `blood`).

## 3. Prerequisites

- **conda env `logan`** (`conda activate logan`). Used here as
  `/scratch/dmk333/conda/envs/logan/bin/python`. Needs: `duckdb` (1.5.3),
  `pandas`, `pyarrow`, `numpy`, `sourmash` (4.8.14), `mmh3`, `screed`.
- **aws-cli v2** on PATH, for anonymous S3 reads (`--no-sign-request`). No credentials.
- **zstd** on PATH (the conda env provides it).
- **Outbound HTTPS** to `blast.ncbi.nlm.nih.gov` and `logan-pub.s3.amazonaws.com`.
- A big machine helps but is not required: the reference box has 768 cores / 3 TB RAM.
  Peak RAM observed was ~36 GB (skin step 2). Free disk needed: ~15 GB.

## 4. External inputs (NOT in this repo)

| path | size | what |
|---|---|---|
| `/scratch/dmk333/Logan_Analyses/verify_diversity_correlation/data/hash_diversity_results_full_cov_0.015625/analysis/filtered_analysis_min_mbases_1620_min_diversity_10/filtered_data_with_dmi_host_removed.parquet` | 19 MB | per-run metrics + DMI. The same file behind `body_site_analysis.png`. |
| `/scratch/dmk333/Logan_Analyses/dark_matter_index/data/samples_with_reference_hashes_cov_0.015625_min_mbases_1620_min_diversity_10.db` | 606 GB | DuckDB: `sample_hashes(sample_id, min_hash)` + `reference_hashes(hash)` |
| `/scratch/dmk333/Logan_Analyses/dark_matter_index/data/human_hashes_k31.db` | 346 MB | DuckDB: `human_hashes(hash)`, 3,701,086 host hashes (HPRC + T2T CHM13 + GRCh38) |
| `s3://logan-pub/c/{acc}/{acc}.contigs.fa.zst` | ~11 MB/run | Logan contigs |
| `s3://logan-pub/u/{acc}/{acc}.unitigs.fa.zst` | ~400 MB/run | Logan unitigs |

**Both DuckDB files are shared, read-only assets. Always open with
`duckdb.connect(path, read_only=True)`.** The pipeline does; do not change that.

`sample_hashes` has **no index**, so any query costs one full 606 GB scan (~4-13 min).
Batch all accessions into a single query -- never loop per accession.

## 5. Parameters, and the one gotcha that will bite you

```
k = 31, scaled = 1000, seed = 42, max_hash = 18446744073709552
hash = murmur3_x64_128 low 64 bits of the canonical 31-mer
canonical = lexicographic min(kmer, reverse_complement(kmer))
```

**The `cov_0.015625` in the input paths is NOT the scaled factor.** It is an unrelated
coverage parameter. Reading it as scaled=64 gives a wrong `max_hash` and silently
corrupts every signature you build. Verified: one sample's 158,186,253 unitig 31-mers
divided by 1000 gives 158,186 expected hashes versus 157,772 actually stored.

Cohort filters (identical to `../body_site_analysis/analyze_body_sites.py`, so the
cohort is exactly the points behind the figure):
`mbases_x >= 100`, `host_fraction <= 0.95`, `total_hashes_host_filtered >= 1000`,
`platform == ILLUMINA`, `dmi_host_filtered` not null, `organism == <site organism>`.

Hash classes:
- **host** -- present in `human_hashes_k31.db`.
- **dark** -- not host AND absent from the reference union: GenBank WGS, GenBank
  genomes, AllTheBacteria, GTDB, SILVA, GenBank TLS, GenBank TSA, Logan plasmids,
  Logan obelisks, Serratus viruses.

Matching is **exact 31-mer identity**; there is no similarity threshold. At 80%
nucleotide identity to a reference, the chance a given 31-mer is conserved is
0.80^31 = 1.0e-3. **So "dark" means "no close relative (>~95% identity) in the
databases", not "no homolog exists."** Every long dark contig BLASTed here had
identifiable bacterial relatives at 66-86% identity. This is the single most important
thing to remember about the DMI.

## 6. Reproducing a site, end to end

From `pipeline/`, with `logan` activated. Substitute any key from `SITES` in
`config.py` (blood, skin, oral, vaginal, nasopharyngeal, lung, gut, feces, saliva,
sputum, urinary).

```bash
cd pipeline
export DMI_SITE=skin          # omit for blood

python 01_select_samples.py                              # seconds
python 02_export_hashes.py                               # see runtimes below
python 03_validate_hashing.py                            # ~5 min; MUST pass
python 04_extract_sequences.py --source contigs --jobs 48
python 04_extract_sequences.py --source unitigs --jobs 3 --accessions ACC1 ACC2 ACC3
python 05_characterize.py
```

Observed runtimes:

| step | blood (125 runs) | skin (1,633 runs) | dominated by |
|---|---|---|---|
| 01 | seconds | seconds | reading the parquet |
| 02 | ~5 min | ~90 min | 606 GB scan (248 s / 754 s) then writing signatures |
| 03 | ~5 min | ~5 min | downloading one contig file + the `sourmash sig kmers` oracle |
| 04 contigs | 58 s @ 32 jobs | 647 s @ 48 jobs | S3 download |
| 04 unitigs | 172 s for 3 runs | not run | 400 MB/run download |
| 05 | ~2 min | ~12 min | re-reading every extracted FASTA twice |

Step 02 is dominated by JSON signature serialisation, not the scan; for a large cohort
most of the wall clock is writing `data/<site>/sigs/`.

**Step 03 is a hard gate and is worth understanding**, because the fast extractor is the
one place a subtle bug would silently corrupt everything downstream. It asserts:
1. sourmash's own `max_hash` for scaled=1000 equals the hard-coded constant;
2. a written `.sig` round-trips, and its `md5sum` matches an independent
   reimplementation (md5 of ksize followed by the sorted hashes as decimal strings);
3. on synthetic sequence, the batched extractor returns **identical** (position, hash)
   pairs to a naive per-k-mer `mmh3.hash64(min(kmer, revcomp), 42)` loop;
4. on one real contig file, the extractor and `sourmash sig kmers` select an
   **identical** set of sequence IDs.
Independently, unitig extraction recovered **100.00%** of dark hashes on all 3 blood
pilots, which is the strongest end-to-end check available.

### Adding a new body site

Add `"<key>": "<organism string>"` to `SITES` in `config.py`, then run the steps above
with `DMI_SITE=<key>`. Nothing else is site-specific: deep-dive runs are chosen
automatically (`pick_pilots` takes the median-hash run from each of the three largest
bioprojects), and the report's conclusions are generated from thresholds on the measured
values rather than hard-coded.

## 7. Why the extractor is written the way it is

`sourmash sig kmers` is the obvious tool and is the wrong one here. It does
`copy_and_clear()` + `add_sequence()` + two `intersection()` calls **per FASTA record**;
a Logan unitig file has ~14.1 million records, so that is ~1.8e9 Python round-trips per
sample. With `--save-kmers` it also tests `hashval in query_mh.hashes` per k-mer,
rebuilding a dict wrapper each time. It cannot read `.zst` and cannot annotate matches.

`extractor.py` instead concatenates ~5 Mbp of records joined by runs of `N` and makes a
single Rust call: `MinHash(ksize=31, scaled=1).seq_to_hashes(batch, force=True,
bad_kmers_as_zeroes=True)`. That returns exactly one hash per k-mer start position, with
every junction k-mer reported as 0, so positions map back to records by arithmetic and
everything after the call is vectorised numpy. About 120 FFI calls per unitig file
instead of 14.1 million. The hashing is still sourmash's own, so it cannot drift from the
sketches in the database -- and check 3 above proves it matches a naive `mmh3` loop.

## 8. Bulk artefacts (gitignored) and how to regenerate them

Run from `pipeline/` with `DMI_SITE` set. Sizes as built.

| artefact | blood | skin | regenerate with |
|---|---|---|---|
| `data/<site>/sigs/` | 452 MB | 5.6 GB | `python 02_export_hashes.py` |
| `data/<site>/hashes/` | 141 MB | 2.1 GB | same |
| `data/<site>/seqs/` | 40 MB | 845 MB | `python 04_extract_sequences.py --source contigs` |
| `data/<site>/kmers/` | 27 MB | 416 MB | same |
| `data/<site>/dark_contigs_blastable.fa` | 53 MB | 771 MB | `python 05_characterize.py` |
| `data/<site>/contig_dark_density.csv` | 7.7 MB | 201 MB | same |

Everything under `results/` is small (largest single file 8.3 MB) and is committed.

## 9. BLAST: how it was run, and two traps

All searches went to NCBI over the URL API (`ncbi_blast.py`), no local BLAST install and
no local `nt`. Full parameter tables are in each `blast_findings.md`. Summary: megablast
vs `nt` for first-pass identification; discontiguous megablast vs `nt` for anything
megablast could not place; blastx vs SwissProt for protein-level identity (blood only).

**Trap 1: NCBI's plain-text report silently omits queries that return zero hits.** They
appear neither with results nor with a "No hits found" marker. The skin megablast run
returned 9 of 15 queries. Do not infer "no hits" from absence -- resubmit the missing
queries alone and confirm they appear with zero subjects. Both skin "no hit" results
were established that way.

**Trap 2: NCBI enforces per-job CPU and size limits.** A blastx-vs-`nr` job on 5 x 6 kb
died with `CPU usage limit was exceeded ... SIGXCPU (24)`; a megablast job on 15 x 30 kb
died with `Process size limit exceeded ... SIGXFSZ (25)`. Keep total query bytes near
150 kb and truncate long contigs. Note that **coverage percentages are then relative to
the truncated query**, not the full contig.

**Coverage means three different things.** `analyze_blast_hsps.py` reports them
separately because conflating them is easy and changes conclusions: `bestHSPcov` (single
longest HSP / query length -- NOT NCBI's "Query Cover"), `unionCov` (merged HSPs to the
best subject), `anyCov` (merged across all subjects). Quoting `bestHSPcov` once made a
contig with 96.9% union coverage look like an isolated gene hit.

## 10. Deep dive: the two contigs that matched nothing

`analyses/unplaced_skin_contigs/FINDINGS.md` follows up the two skin contigs that
returned zero hits from both megablast and discontiguous megablast against all of `nt`.
Both turned out to be **mobile genetic elements**, which is the cleanest explanation of
why they are dark: MGEs are the fastest-evolving and least-sampled part of microbial
genomes, so they sit well outside the ~95%-identity radius exact 31-mer matching needs.

- `SRR15899555_383` (74,940 bp, 93.5% coding) is an **integrative and conjugative
  element**: MobL relaxase, VirD4 coupling protein, site-specific recombinase,
  transposase, phage lysozyme, LysM. Relaxase is 47% aa identical to *Arthrobacter*
  (Actinobacteria). It is detected at >=10% of its hashes in **14 runs across 6
  bioprojects on two continents** -- a widespread skin element absent from every
  reference database.
- `ERR7738256_16833` (62,324 bp, 97.2% coding, all 50 genes on one strand) carries its
  own replication and transcription machinery including a **fused rpoB-rpoC** RNA
  polymerase (2,606 aa; beta domain at residues 462-1142, beta-prime at 1248-2326), at
  65-72% aa identity to uncultured *Clostridium* / Lachnospiraceae. No capsid or
  terminase was found, so it is not called a phage. It is a singleton.

Extra tooling used only here, all documented in that file: pyrodigal (gene calling, via
`pip install --target analyses/tools`, so the shared env is untouched), NCBI Batch
CD-Search (`pipeline/ncbi_cdsearch.py`), and ESMFold + Foldseek
(`pipeline/fold_and_search.py`). Structure prediction added little: mean pLDDT was
0.33-0.71 and only two of eight searches found anything significant, so those are leads,
not assignments.

## 11. Caveats that apply to everything here

- The BLASTed contigs are the **longest** ones, which is a length-biased convenience
  sample (15 of 152,025 for skin). Nothing here supports claims about the composition of
  the dark fraction as a whole.
- Assembly rates come from **3 runs per site**, chosen as project medians. Consistent in
  direction, not a survey.
- There are **no negative controls** in either dataset, so reagent contamination and
  genuine environmental DNA cannot be formally separated.
- Blood's cross-bioproject comparison is **underpowered**: 3 projects, two of which
  contain 8 and 1 runs.
- Skin's 26 bioprojects include 7 from NISC totalling 52% of runs, so "many
  bioprojects" overstates independence.
- Taxonomic labels are best-BLAST-hit at 66-86% nucleotide identity, i.e. below species
  and often below genus. They indicate approximate relatedness, not identification. No
  synteny or phylogenetic analysis was done, so horizontal transfer versus ordinary
  orthology is unresolved.

## 12. Provenance

The parent figure and its cohort definition come from
`../body_site_analysis/analyze_body_sites.py`. The DMI itself and the host-filtered
parquet are produced upstream by
`../../scripts/compute_dmi_host_filtered.py`; the reference union membership is
defined in `../../scripts/DB_info.json`.
