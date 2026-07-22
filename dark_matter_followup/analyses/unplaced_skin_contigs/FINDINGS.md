# The two unplaceable skin contigs: what they are

ASCII-only by convention.

These are the two contigs from `results/skin/blast_findings.md` section 2.7 that
returned **zero hits** from both megablast and discontiguous megablast against the whole
of `nt`. This directory pulls them apart by gene prediction, conserved-domain search,
protein BLAST, cross-cohort prevalence, and structure prediction.

## Bottom line

Both are **mobile genetic elements**, which is the single best explanation for why they
are dark: MGEs are the fastest-evolving, most mosaic, and least well-sampled part of
microbial genomes, so they sit far outside the 95%-identity radius that exact 31-mer
matching requires.

| | SRR15899555_383 | ERR7738256_16833 |
|---|---|---|
| length | 74,940 bp | 62,324 bp |
| coverage (ka) | 76.5x | 11.3x |
| GC | 45.65% | 29.86% |
| genes (Prodigal meta) | 53 | 50 |
| coding density | 93.5% | 97.2% |
| strand organisation | 30 fwd / 23 rev, 5 switches | **all 50 on one strand** |
| proteins with a CDD domain | 25 / 53 | 14 / 50 |
| **identity** | **integrative and conjugative element (ICE)** | **large mobile element, Clostridia-associated** |
| nearest relatives | MobL relaxase 47% aa to *Arthrobacter* (Actinobacteria) | fused rpoB-rpoC 65-72% aa to uncultured *Clostridium* / Lachnospiraceae |
| prevalence in 1,633 skin runs | **14 runs, 6 bioprojects at >=10%** | 1 run (singleton) |
| source | PRJNA763232, U. Wisconsin | PRJEB49206, Stanford |

## Methods

Run from `analyses/unplaced_skin_contigs/`. `pyrodigal` was installed into
`analyses/tools/` with `pip install --target` so the shared `logan` env was not modified;
add it to `sys.path` to reuse.

| step | tool | how |
|---|---|---|
| full-length extraction | `pipeline/extractor.py` | pulled from `data/skin/seqs/contigs/*.fa.gz` (the BLAST queries had been truncated to 8-10 kb; these are full length) |
| gene prediction | pyrodigal 3.7.1, `GeneFinder(meta=True)` | metagenomic mode, no training; standard code |
| domain search | NCBI Batch CD-Search, `pipeline/ncbi_cdsearch.py` | db=cdd (CDD+Pfam+SMART+COG+TIGRFAM+NCBIfam), E<=0.01, dmode=full, all 103 proteins in one job |
| protein BLAST | NCBI blastp vs `nr`, `pipeline/ncbi_blast.py` | E<=1e-5, hitlist 10. Whole-protein jobs failed on NCBI's CPU cap, so single domains of <=700 aa were submitted |
| prevalence | sourmash | contig sketched at k=31/scaled=1000, intersected against every run's `data/skin/hashes/*.parquet` |
| structure | ESMFold API + Foldseek server, `pipeline/fold_and_search.py` | 8 domain-less proteins of 150-400 aa; Foldseek vs afdb50, pdb100, afdb-proteome |

**Percent identity** below is amino acid identity over the gapped alignment length of the
best HSP, taken verbatim from BLAST's `Identities = a/b (c%)` line.

## SRR15899555_383 -- an integrative and conjugative element

Annotated gene order (full table in `gene_map.csv`):

```
gene   coords          aa    domain                E-value
   2   1763-3571      602    Phage_lysozyme2       3.3e-30   endolysin
   7   8198-8710      170    MobL                  1.2e-11   relaxase
   8   8701-11379     892    MobL                  8.3e-47   relaxase
  13   15587-16186    199    COG2452               1.6e-49   site-specific recombinase
  14   16176-17654    492    InsQ                  1.0e-31   transposase
  17   19788-23717   1309    ClpA                  4.8e-37   AAA+ ATPase
  20   27849-29960    703    VirD4                 1.9e-12   T4SS coupling protein
  22   32678-34978    766    VirD4                 4.2e-27   T4SS coupling protein
  32   49322-52858   1178    RecB                  4.6e-27   nuclease
  34   53790-55808    672    clpC                  1.9e-84   AAA+ ATPase
  43   62504-63625    373    RseP                  9.6e-52   site-2 protease
  46   66028-66663    211    LysM                  9.2e-09   peptidoglycan binding
  49   67270-69189    639    DEXQc_bact_SNF2       1.7e-63   helicase
  51   69747-70154    135    VRR_NUC               2.0e-06   nuclease
```

Relaxase (MobL) + coupling protein (VirD4) + site-specific recombinase + transposase +
peptidoglycan-degrading lysozyme + LysM is the **canonical ICE / conjugative transposon
gene set**. The relaxase is 47% aa identical to MobL relaxases of *Arthrobacter* and
*Saccharopolyspora* (Actinobacteria) -- consistent with an actinobacterial host, which
fits skin, where Actinomycetota (*Cutibacterium*, *Corynebacterium*, *Micrococcus*)
dominate. 47% identity over 601 residues is a genuine but distant homolog.

**It is not rare.** Sketching the contig and intersecting against all 1,633 skin runs:

| detection threshold | runs | bioprojects |
|---|---|---|
| >=1 shared hash | 211 | 11 |
| >=3 hashes | 73 | 8 |
| >=10% of its 85 hashes | **14** | **6** |
| >=25% | 6 | 2 |
| 100% (source) | 1 | 1 |

Strong detections: SRR15899555 (100%), SRR15899489 (48%), SRR3188290 (46%),
SRR15899460 (45%), SRR15899511 (41%), SRR3189141 (39%), SRR17534316 (24%),
SRR18292695 (21%) -- spanning PRJNA763232 (Wisconsin), PRJNA46333 and PRJNA604820
(NISC), PRJNA814343 (A*STAR Singapore), PRJNA471898 (NISC).

So this is a **conjugative element circulating in skin communities on at least two
continents, at meaningful abundance in ~1% of runs and detectably in ~13%, that matches
nothing in the reference union.** Lead with the >=10% tier (14 runs, 6 projects); the
211-run figure is mostly single-hash detections that could be one shared conserved
fragment.

## ERR7738256_16833 -- a large Clostridia-associated element

```
gene   coords          aa    domain                E-value
   4   1694-3844      716    FtsK_SpoIIIE          1.1e-20   DNA translocase
   9   9058-9717      219    DnaQ                  3.2e-21   3'-5' proofreading exonuclease
  29   25341-29063   1240    PolA                  1.4e-28   DNA polymerase I
  33   32854-33363    169    McrA                  6.5e-15   HNH endonuclease
  34   33368-41188   2606    PRK09603/rpoB/rpoC    1.7e-116  RNA polymerase, see below
  40   43357-44574    405    DNA_ligase_aden       4.7e-06   DNA ligase
  43   45590-46333    247    COG4509               2.0e-56
  44   46348-47937    529    SpaA                  4.2e-04   pilin/adhesin
  47   50219-57238   2339    ClfA                  1.6e-48   MSCRAMM surface adhesin
  48   57429-61058   1209    ClfA                  1.0e-42   MSCRAMM surface adhesin
```

**Gene 34 is a fused rpoB-rpoC.** The domain hits partition cleanly along the 2,606
residue protein with no overlap:

```
   462-1142   rpoB / RNA_pol_B_RPB2 / RNAP beta      E=4.6e-99
  1248-2326   rpoC / RNAP_beta'_N / RPOLA_N          E=1.2e-93
```

i.e. an RNA polymerase beta subunit followed by a beta-prime subunit in a single open
reading frame. blastp of the two halves separately gives **72%** (rpoB half, 681 aa) and
**65%** (rpoC half, 683 aa) amino acid identity, in both cases to the *same* subject
proteins -- "DNA-directed RNA polymerase subunit beta family protein" from **uncultured
*Clostridium* sp.**, *Paraclostridium bifermentans*, *Anaerobutyricum hallii/soehngenii*,
*Eubacterium* sp. and a Lachnospiraceae bacterium. Those subjects are therefore also
fused, so the fusion is a feature of this divergent family rather than a one-off.

The element carries its own replication and transcription machinery (RNAP, DNA
polymerase I, ligase, proofreading exonuclease, FtsK translocase, HNH endonuclease) plus
two very large LPXTG-type MSCRAMM surface adhesins (2,339 and 1,209 aa) and a pilin.
GC of 29.9% is a good fit for Clostridia.

**I stop short of calling it a phage.** It has the replication and transcription modules
of a large phage, unidirectional transcription across all 62 kb, and monotonic GC skew.
But CD-Search found **no capsid, terminase, portal or tail genes**, so a virion-forming
phage is not supported by the evidence. A chromosomal island, a phage-plasmid, or a
fragment of a larger element are all consistent. It is a singleton -- present in exactly
1 of 1,633 skin runs -- so there is no cross-sample evidence to help.

An earlier guess of mine that the rpoB-rpoC fusion pointed to Campylobacterota
(*Helicobacter*-like) was **wrong**; blastp puts it firmly in Clostridia/Lachnospiraceae.

## Neither is cellular, and they are unrelated to each other

- **No ribosomal proteins, tRNA synthetases or elongation factors** in either contig.
  The apparent "RecA-like" CD-Search hits are the RecA-like ATPase fold of ClpA/ClpB,
  not RecA. A 62-75 kb fragment of a real chromosome would not be expected to contain
  many of these by chance, so this is suggestive rather than decisive -- but combined
  with the density and organisation it points away from cellular genome fragments.
- **Zero shared 31-mers** between the two contigs: unrelated elements.
- **No terminal direct repeats** >= 31 bp in either, so neither is a circularly permuted
  or DTR-complete phage genome as assembled.

## Structure prediction added little

Eight domain-less proteins of 150-400 aa were folded with ESMFold and searched with
Foldseek against afdb50, pdb100 and afdb-proteome (`structures/summary.tsv`).

- Mean pLDDT was **0.33-0.71** (i.e. 33-71 on the usual 0-100 scale). ESMFold is
  single-sequence, and for families with no relatives it is expected to be poor. Most
  models are low confidence and one failed outright.
- Only two searches produced anything significant: `SRR15899555_383_18` matched
  AF-Q58446 *DNA-directed RNA polymerase subunit Rpo1C* (E=6.6e-4, prob=1), and
  `ERR7738256_16833_45` matched an uncharacterised AFDB entry (E=3.1e-7, prob=1).
  Everything else had best E-values of 0.4-8, i.e. no detectable structural homolog.
- **Given the low pLDDT, even the two significant hits should be treated as leads, not
  assignments.** A structure search on a 40-pLDDT model is not reliable evidence.

The honest summary is that roughly half of `SRR15899555_383` (28/53 proteins) and nearly
three quarters of `ERR7738256_16833` (36/50) have no assignable domain, and structure
prediction did not rescue them.

## What I am not claiming

1. Host assignment is by **best BLAST hit at 47-72% amino acid identity**, which is
   genus-level at best and often above it. "Actinobacterial" and "Clostridia-associated"
   are the strength of claim the data supports; naming a species would not be.
2. `ERR7738256_16833` is a **singleton**. Nothing here excludes it being a rare organism,
   an assembly chimera, or a fragment of something much larger. It was not independently
   verified.
3. No **tRNA or ncRNA scan** was run (no tRNAscan-SE/Aragorn/Infernal available), so
   "no tRNAs" is not established -- it was simply not tested.
4. Prodigal in metagenomic mode over-calls short ORFs in dense sequence; gene counts and
   coding density should be read as approximate.
5. Contigs are **as assembled by Logan**; neither is a closed or validated element, and
   assembly artefacts (chimeras) cannot be excluded without read-level checking.
6. Structure conclusions are weak for the reason given above.
7. The prevalence scan uses scaled=1000 sketches, so a contig contributes only 59-85
   hashes; single-hash detections are weak evidence and are reported separately for that
   reason.

## Files

```
unplaced_full.fa            the two contigs, full length
unplaced_proteins.faa       103 predicted proteins
unplaced_genes.fna          the corresponding nucleotide CDS
gene_map.csv                gene, coords, strand, length, best domain, E-value
cdd_hits.tsv               raw CD-Search output (526 hit rows)
blastp_markers.txt          blastp vs nr for rpoB, rpoC and MobL domains
blastp_nr_longest.txt       FAILED job (NCBI CPU limit) -- kept for the record
longest_proteins.faa        query for the failed job
marker_proteins.faa         the three marker domains that were searched
nodomain_top.faa            longest domain-less proteins
fold_targets.faa            the 8 sent to ESMFold
structures/*.pdb            ESMFold models (B-factor column = pLDDT)
structures/*.foldseek.json  raw Foldseek results
structures/summary.tsv      length, mean pLDDT, top structural hit
```
