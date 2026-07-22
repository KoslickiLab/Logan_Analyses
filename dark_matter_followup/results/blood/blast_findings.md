# What the blood "dark matter" actually is -- BLAST results

ASCII-only by design (these reports get read in terminals and plain text viewers).

## 1. METHODS

### 1.1 What was BLASTed

Input file: `dark_contigs_blastable.fa`, produced by `05_characterize.py`. Membership
rule: a Logan contig is included if it is >= 300 bp AND at least 80% of the sample's
non-host sketch hashes falling on that contig are "dark" (see 1.2). That pool holds
71,067 contigs (49.8 Mbp), of which 11,492 are >= 1 kb.

From that pool I took the **30 longest contigs** (`top_dark_contigs.fa`) and submitted
the **top 15** to nucleotide BLAST. This is a length-biased convenience sample of
0.13% of the >= 1 kb pool -- see the caveats in section 5.

### 1.2 How "dark" is defined upstream

- Sketches: FracMinHash, k = 31, scaled = 1000, seed = 42,
  max_hash = 18446744073709552, hash function murmur3_x64_128 lower 64 bits of the
  canonical (lexicographically smaller of forward / reverse-complement) 31-mer.
- host hash: present in `human_hashes_k31.db` (HPRC + T2T CHM13 + GRCh38).
- dark hash: NOT a host hash AND absent from the reference union
  (GenBank WGS, GenBank genomes, AllTheBacteria, GTDB, SILVA, GenBank TLS,
  GenBank TSA, Logan plasmids, Logan obelisks, Serratus viruses).
- Matching is **exact 31-mer identity**. There is no similarity threshold; a single
  substitution destroys the match. This is the crux of section 3.

### 1.3 BLAST runs (all via the NCBI URL API, script `ncbi_blast.py`)

No local BLAST install and no local nt copy; every search ran on NCBI servers through
`https://blast.ncbi.nlm.nih.gov/Blast.cgi` (CMD=Put then CMD=Get), identifying
tool and email per NCBI etiquette, polling at 60 s intervals.

| # | program | database | task / flags | E-value | hitlist | queries | output |
|---|---------|----------|--------------|---------|---------|---------|--------|
| 1 | blastn  | nt        | `MEGABLAST=on` (megablast, default word size 28) | 1e-3 | 10 | top 15 contigs, full length (139,280 nt total) | `blast_megablast_nt.txt` |
| 2 | blastx  | nr        | default | 1e-3 | 8 | 5 contigs truncated to first 6,000 nt | FAILED |
| 3 | blastx  | swissprot | default | 1e-3 | 6 | same 5 contigs, first 6,000 nt | `blast_blastx_swissprot.txt` |
| 4 | blastn  | nt        | `BLAST_PROGRAMS=discoMegablast` (discontiguous megablast) | 1e-3 | 20 | `SRR14209604_503` only, full 22,182 nt | `blast_disco_503.txt` |

Run 2 returned `[blastsrv4.REAL] Error: CPU usage limit was exceeded, resulting in
SIGXCPU (24)` and produced no results. Run 3 (SwissProt, a much smaller protein
database) was substituted. **No blastx-vs-nr result exists**, so protein-level
identification rests on SwissProt only, which is curated and small and will miss
uncultured-organism proteins.

All other BLAST parameters were NCBI web defaults: for megablast, word size 28,
match/mismatch 1/-2, gap costs linear; for discontiguous megablast, word size 11 with
template length 16, match/mismatch 2/-3, gap open 5 extend 2; for blastx, word size 3,
BLOSUM62, gap open 11 extend 1, standard genetic code, all 6 reading frames. Low
complexity filtering (DUST for nucleotide, SEG for protein) was left at its default
(on). No taxonomic restriction was applied. Databases were whatever NCBI served on
2026-07-21.

### 1.4 How the numbers in section 2 are defined

BLAST reports several different "identity" and "coverage" quantities and they are
easy to conflate, so each is spelled out. Computed by `analyze_blast_hsps.py`, which
parses the plain-text report.

- **HSP**: one local alignment, i.e. one `Score = ...` block.
- **pct identity**: taken verbatim from the `Identities = a/b (c%)` line of an HSP,
  so it is `a / b` where `b` is the **gapped alignment length** (gap positions count
  in the denominator). For run 1 and 4 this is **nucleotide** identity. For run 3
  (blastx) it is **amino acid identity** -- yes, "aa" in this document means amino
  acid identity of the translated query against the protein subject. Nucleotide and
  amino acid identity are not comparable numbers.
- **bestHSPcov**: alignment length of the single longest HSP divided by query length.
  This is *not* NCBI's "Query Cover" column. It was the number in my first summary and
  it understates homology whenever it is split across several HSPs.
- **unionCov**: fraction of query positions covered by at least one HSP **to the single
  best subject**, merging overlapping HSPs. This is the closest analogue to NCBI's
  per-subject "Query Cover".
- **anyCov**: fraction of query positions covered by at least one HSP to **any**
  subject in the report, merging across all subjects.
- **nWindows**: number of disjoint intervals making up anyCov. A small nWindows with
  a low anyCov means homology is confined to a few discrete loci.

Coverage denominators are always the full query contig length.

### 1.5 Non-BLAST checks (local, `05_characterize.py` and ad hoc)

- **Concatemer test**: count of distinct 31-mers divided by (contig length - 30). A
  rolling-circle or otherwise tiled amplification product would repeat sequence and
  score well below 1.0.
- **Coding test**: longest stop-free run, in nucleotides, over all 6 reading frames,
  using the standard genetic code stop codons TAA/TAG/TGA. Compared against the
  expectation for random sequence at the observed GC content.
- **GC content**: (G + C) / length, on the contig as written.

## 2. RESULTS

### 2.1 Cohort context (from the literature)

**PRJNA634526** contributes 11,491 of the 11,492 >= 1 kb fully-dark contigs. It is
Zhang et al., mSphere 2021, "HIV-1 Infection Alters the Viral Composition of Plasma in
Men Who Have Sex with Men" (https://journals.asm.org/doi/10.1128/msphere.00081-21),
Institut Pasteur of Shanghai + Shenzhen CDC. The SRA `lat_lon` of 22.54 N 114.06 E is
Shenzhen. Protocol per the paper:

- plasma, **0.45 um filtration plus nuclease digestion** to enrich virus-like particles
- **MALBAC whole-genome amplification** from very low input
- human and bacterial reads removed *bioinformatically* after sequencing
- NovaSeq, 2 x 150 bp; 101 MSM across four HIV/ART/CD4 strata plus 20 controls
- reported virome dominated by anelloviruses, then pegivirus and HBV

**PRJNA787952** (8 runs, 1 contig) -- Stanford, Severyn et al. JCI Insight 2022, gut
decontamination trial in paediatric allogeneic HCT (NCT02641236); these are its blood
samples. **PRJNA544865** (1 run) -- UCSF / CZ Biohub, `isolation_source = platelet bag`.

### 2.2 Every long dark contig has a bacterial relative, at 71-86% nucleotide identity

megablast vs nt (run 1). Identity is of the best HSP; coverages as defined in 1.4.

| contig | len | GC | best subject | ident | bestHSPcov | unionCov | anyCov | windows |
|---|---|---|---|---|---|---|---|---|
| SRR14209604_503 | 22,182 | 0.380 | Ca. Woesebacteria MAG | 81% | 0.65% | 0.63% | 1.32% | 2 |
| SRR14209591_1830 | 10,211 | 0.694 | uncultured Actinomycetota MAG | 78% | 5.46% | 7.29% | 7.29% | 2 |
| SRR14209591_3125 | 9,533 | 0.544 | Flavisolibacter sp. MAG | 75% | 16.03% | 15.91% | 45.03% | 3 |
| SRR14209591_1447 | 8,537 | 0.703 | Microlunatus phosphovorus | 77% | 41.70% | 40.92% | 41.20% | 1 |
| SRR14209591_579 | 8,184 | 0.488 | Leptodesmis sichuanensis | 79% | 15.26% | 15.49% | 80.03% | 7 |
| SRR14209591_7030 | 8,090 | 0.626 | Microbacterium sp. | 83% | 5.22% | 5.13% | 5.57% | 1 |
| SRR14209591_3634 | 7,788 | 0.565 | Ilyomonas sp. MAG | 71% | 11.31% | 11.12% | 17.72% | 3 |
| SRR14209591_2852 | 7,608 | 0.497 | Leptolyngbya sp. NIES-3755 | 82% | 55.28% | **96.90%** | 96.90% | 3 |
| SRR14209591_166 | 7,562 | 0.492 | Leptolyngbya sp. NIES-3755 | 79% | 30.63% | 63.12% | 68.95% | 3 |
| SRR14209591_870 | 7,085 | 0.720 | Nocardioides aquaticus | 78% | 15.46% | 15.26% | 15.30% | 1 |
| SRR14209633_152403 | 7,029 | 0.450 | uncultured bacterium MAG | 72% | 17.47% | 17.24% | 31.03% | 3 |
| SRR14209591_1379 | 7,005 | 0.680 | Qipengyuania sediminis | 80% | 73.02% | 72.46% | 72.46% | 1 |
| SRR14209591_1306 | 7,002 | 0.541 | Flavisolibacter ginsenosidimutans | 79% | 54.54% | **86.76%** | 88.47% | 3 |

Genera across all hits in run 1: Leptolyngbya (21 HSPs), Streptomyces (8),
Qipengyuania (8), Microbacterium (8), Hymenobacter (6), Flavisolibacter (5),
Cellulomonas (5), Georgenia (4), Mucilaginibacter (3), Bradyrhizobium (3),
Chitinophaga, Kribbella, Nonomuraea, Nocardioides, Pseudocnuella,
Candidatus Woesebacteria.

Robustness note on the assembly-rate result quoted elsewhere: `blood_dark_matter_report.txt`
now auto-selects its deep-dive runs (median non-host hash count within each of the three
bioprojects) rather than using a hard-coded list. That swapped SRR14209595 for
SRR14209668 in PRJNA634526 and moved the dark-vs-reference-matched assembly-rate ratio
from 21.2x to 18.5x. The conclusion does not depend on which run is picked.

### 2.3 Why these are dark: the arithmetic

At 80% nucleotide identity, the probability that any given 31-mer is exactly conserved
is 0.80^31 = 1.0e-3. Each of these contigs contributes only about 7-23 hashes to a
scaled=1000 sketch (one hash per ~1000 k-mers), so the expected number of *matching
sketch hashes* is about 0.02. Observing zero is exactly what the model predicts.

> DMI "dark" means **no close relative (roughly >= 95% nucleotide identity) in the
> reference union**. It does not mean no homolog exists. These sequences are
> recognisably bacterial; they are 15-30% divergent from anything sequenced.

### 2.4 Coverage is often high once HSPs are merged -- and the low ones were partly an artefact

My first summary quoted bestHSPcov, which is the least generous of the three coverage
measures. Merging HSPs changes the reading substantially: SRR14209591_2852 covers
**96.9%** of its length against a single Leptolyngbya genome, and _1306 covers 86.8%.
Those are not "one conserved gene" cases; they are divergent genome fragments aligned
along nearly their whole length.

The genuinely low-coverage cases were then tested directly. Re-running
`SRR14209604_503` with **discontiguous megablast** (run 4), which is designed for
cross-species sensitivity, raised coverage from **1.32% to 34.16%** (7 windows, 84 HSPs
across 20 subjects) with best identity dropping to **66%**. So most of the apparent
"no homology" was megablast's word-size-28 seeding failing on diverged sequence, not
a real absence of relatives.

This matters for your HGT question. With only megablast, two 150 nt islands in 22 kb
looked like an isolated conserved gene in otherwise alien sequence -- the classic
signature you would want to invoke HGT for. With a sensitive search, a third of the
contig aligns at ~66% identity. **These data cannot distinguish HGT from ordinary
orthology in a deeply divergent lineage**, and the megablast-only view was actively
misleading on the point. Deciding it would need synteny (are the homologous windows
co-linear with the subject, or shuffled?) and phylogenetics per gene, neither of which
was done here.

### 2.5 What SRR14209604_503 is: probably CPR / Patescibacteria

The 20 discontiguous-megablast subjects are almost entirely Candidate Phyla Radiation
MAGs, at 65-71% identity:

```
OZ246394.1  uncultured Gammaproteobacteria MAG          721 bits  0.0     66%
AP024681.1  Candidatus Parcubacteria bacterium KatS3mg097 685     0.0     69%
CP055307.1  Candidatus Woesebacteria bacterium           581     5e-159  67%
OY970230.1  Candidatus Gottesmanbacteria bacterium       498     4e-134  66%
OZ490593.1  Candidatus Levyibacteriota bacterium         352     4e-90   70%
OZ490654.1  Minisyncoccia bacterium                      331     4e-84   70%
CP066688.1  Candidatus Daviesbacteria bacterium          319     2e-80   71%
OZ246395.1  uncultured Candidatus Saccharibacteria       267     1e-64   68%
AP024679.1  Candidatus Dojkabacteria bacterium           241     5e-57   67%
```

CPR / Patescibacteria are ultra-small (roughly 0.2-0.4 um), reduced-genome bacteria
known almost exclusively from MAGs rather than isolates. Two things line up:

1. they are **under-represented in reference databases** by construction, which is why
   nothing gets closer than ~70% identity, which is why the contig is dark; and
2. they are **small enough to pass a 0.45 um filter**, so a VLP-enrichment protocol
   would concentrate rather than remove them.

The SwissProt blastx hit for this contig was competence protein ComM at 37% **amino
acid** identity, consistent with a divergent bacterium and not informative about phylum.

This is the strongest genuine-novelty candidate in the set. It is a hypothesis, not an
identification: all supporting subjects are themselves uncultured MAGs, identity is
only 65-71%, and the single best hit is labelled Gammaproteobacteria rather than CPR.

### 2.6 They are real coding DNA, not noise or amplification artefact

- **Not concatemers**: distinct-31-mer fraction = 1.0000 for all 15 top contigs, so
  MALBAC did not tile a small circle into a fake long contig.
- **Not random sequence**: longest stop-free ORFs run 1,302-5,370 nt. At GC = 0.38 a
  stop codon is expected every ~15 codons, so a 718-codon stop-free run has
  p ~ 1e-17 even after correcting for 6 frames and all start positions.
- **Protein hits are bacterial housekeeping / metabolism** (SwissProt, amino acid
  identity): uronate isomerase 76%, ABC transporters YdiF / YheS 26-29%, type-4
  uracil-DNA glycosylase 46-48%, competence protein ComM 37%, and a traX-finO
  intergenic ORF at 25% (a conjugative-plasmid region).
- **No viral hallmark genes, no anellovirus, nothing eukaryotic** among the top hits --
  notable given this is a virome preparation.

### 2.7 The taxa are environmental, not human-associated

Leptolyngbya and Woronichinia are obligately photoautotrophic cyanobacteria and cannot
replicate in plasma. The remainder are soil / freshwater Actinomycetota (Streptomyces,
Microbacterium, Cellulomonas, Nocardioides, Kribbella, Microlunatus),
Chitinophagaceae / Bacteroidota (Flavisolibacter, Chitinophaga, Mucilaginibacter,
Ilyomonas), and Alphaproteobacteria (Qipengyuania, Bradyrhizobium -- the latter a
canonical reagent contaminant). None is a bloodstream pathogen or a human-microbiome
taxon.

With 0.45 um filtration (intact ordinary bacteria removed, leaving free DNA and
ultramicrobacteria), very low input, and MALBAC amplification (which amplifies trace
template enormously), the parsimonious reading is **environmental and/or reagent DNA
amplified during library preparation**, plus a CPR component that the filtration step
would actively enrich.

### 2.8 It is systematic, not one bad tube

The 11,492 >= 1 kb fully-dark contigs come from **117 of 125 runs**. The single worst
run holds 25.4%, the top five 40.8%, the top ten 47.3%. Being spread across the cohort
with a heavy tail implicates a shared reagent or protocol batch rather than one
contaminated sample.

## 3. WHAT I AM NOT CLAIMING

1. **15 contigs is not a sample of 11,492**, and they were selected by *length*, which
   biases toward whatever assembled best. The other 99.87% may look different.
2. Those 15 come from **3 runs** (SRR14209591 alone supplies 10 of them), so run-level
   and contig-level effects are not separable here.
3. **There are no negative controls in this dataset.** I cannot distinguish reagent
   contamination from environmental DNA genuinely circulating in donor plasma. The
   protocol makes the former more likely; it does not prove it.
4. The CPR assignment for SRR14209604_503 is a **hypothesis from MAG hits at 65-71%
   identity**, not an identification.
5. **No blastx-vs-nr was completed** (NCBI CPU limit). Protein evidence is SwissProt
   only, which is small and curated and will miss uncultured-organism proteins.
6. Low coverage in run 1 was substantially a **megablast sensitivity artefact**
   (section 2.4). Any statement of the form "no better hit exists" applies to the
   specific program and database used, not to nr or to sensitive search generally.
7. Nothing here was checked for **synteny or phylogeny**, so HGT versus orthology is
   open.
8. This says nothing about other body sites. See `../skin/` for the skin comparison.

## 4. BEARING ON THE TALK

Blood's number-one DMI rank is best explained as a **protocol artefact of one virome
study** rather than evidence that blood harbours the most uncharacterized biology.

The finding that survives is methodological, and it is a better talk: there exists a
large, real, protein-coding, 15-30%-divergent bacterial fraction -- including probable
CPR -- that exact 31-mer matching cannot see at all. That is a quantitative statement
about how far reference databases sit from environmental diversity, and the DMI is
measuring exactly that distance.
