# Cross-species expansion: references and datasets (2026-08-24)

Goal: test whether the per-nt Ribo-seq signal model is broadly applicable beyond human/mouse.
Two tiers: (A) the three non-human primates already sitting in the Ruiz-Orera study, and (B) four
distant model organisms (C. elegans, zebrafish, Drosophila, yeast).

## A. References -- ALL SEVEN LOCATED, none missing

| species | assembly | annotation | date | source | status |
|---|---|---|---|---|---|
| Gorilla | GCF_029281585.2 NHGRI_mGorGor1-v2.1 (T2T) | RefSeq | 2025-04 | NCBI | fetching |
| Chimp | GCF_028858775.2 NHGRI_mPanTro3-v2.1 (T2T) | RefSeq | 2026-05 | NCBI | fetching |
| Macaque | GCF_049350105.2 T2T-MMU8v2.0 (complete) | RefSeq | 2025-08 | NCBI | fetching |
| C. elegans | GCF_000002985.6 WBcel235 | WormBase WS298 | 2025-12 | NCBI | fetching |
| Zebrafish | GCF_000002035.6 GRCz11 | RefSeq | 2024-08 | NCBI | fetching |
| Yeast | GCF_000146045.2 R64 | SGD R64-5-1 | 2026-07 | NCBI | fetching |
| Drosophila | dmel r6.67 | FlyBase r6.67 | -- | **ALREADY LOCAL** | `annotations/fly/` |

**Use NCBI, not Ensembl, for the primates.** Ensembl 116 still ships gorGor4 (2011-era) and
Pan_tro_3.0 for gorilla and chimp; the T2T assemblies above are far better and, importantly, carry
CURRENT RefSeq annotation, which is what RiboCode needs. Ensembl is fine for the other species but
NCBI was used uniformly for consistency.

**Drosophila needs no download.** The local FlyBase r6.67 is NEWER than both Ensembl BDGP6.54 and
NCBI's FlyBase 6.54 mirror.

Downloads are chained SERIALLY behind the human RNA fetch (head-node rule) via
`scripts/fetch_species_refs.sh`. Each writes to the SHARED tree
`genomes/xspecies_refs/<label>/{genome.fna,genomic.gtf,PROVENANCE.json}` so other projects reuse it.

**Still to build per species:** STAR index (`STAR_indexes/star_index_<species>_<annver>/`) and,
where quantification is wanted, a decoy-aware Salmon index. Rough index cost: ~30 GB per primate,
~25 GB zebrafish, ~4 GB fly, ~3 GB worm, ~1 GB yeast. Quota has 2.24 TB free, so this is not a
constraint.

## B. Datasets -- the metadata trap that shapes this whole search

**`library_strategy` is USELESS for finding Ribo-seq in these organisms.** The `Ribo-Seq` value was
added to the SRA/ENA vocabulary relatively recently, and essentially all worm/fish/fly/yeast
ribosome profiling predates it. Measured distributions (5,000-run samples per taxon): ZERO runs
labelled `Ribo-Seq` for any of the four. They are deposited as `RNA-Seq` or `OTHER`.

Consequences, all verified rather than assumed:
1. Search **study titles** for "ribosome profiling"/"Ribo-seq"/"translatome", not `library_strategy`.
2. Within a study, separate Ribo from RNA by **sample_title**, not strategy.
3. A matched pair is often **split across BioProjects** (GEO SuperSeries), so a project whose title
   says "matched" may still contain only one arm.

### Candidates, with verification status

| species | study | composition | size | status |
|---|---|---|---|---|
| C. elegans | **PRJNA1183817** (low-input Ribo-seq, early embryogenesis) | 28 Ribo (`OTHER`, "ribosomes profiling") + 21 RNA (`RNA-Seq`, "Total RNA-seq") | 56.4 + 4.1 GB | **matched IN ONE PROJECT, verified** -- but see poly(A) caveat |
| Zebrafish | PRJNA288987 ("RNA-Seq matched to ribosome profiling over a developmental timecourse") | 24 runs, 8 stages x 3 reps, all `RNA-Seq`, titles are stage names only | ~40 GB | **UNRESOLVED**: looks like the RNA arm only; the Ribo half is probably a companion BioProject. Must locate before use |
| Drosophila | PRJNA307183 (enterocyte translatome, dietary restriction) | 8 runs, all `RNA-Seq`, titles "Ad Libitum rep2" etc | small | **UNRESOLVED**: cannot tell Ribo from RNA by title; needs GEO lookup |
| Yeast | PRJNA231536 and others | PRJNA231536 has only 1 run | -- | **INSUFFICIENT**: need a study with replicates + matched RNA |

### The poly(A) caveat, which may disqualify the C. elegans candidate

Its RNA arm is labelled **"Total RNA-seq"**. The standing rule
([[feedback_rnaseq_must_be_polya]]) exists because Ribo-Zero total RNA retains tRNA and 7SL, which
are typed as lncRNA and survive the biotype filter -- this crushed the GSE243134 liver universe to
3,321 tx. A total-RNA arm is a red flag for universe construction even when the library is
otherwise fine. Verify the protocol text before committing, and consider whether the universe can
be built from a different source for that organism.

## C. Recommended order

1. **Primates first.** The Ribo AND RNA already exist in PRJEB65856, accessions in hand
   (gg 6 RNA runs/38 GB, pt 11/58.9 GB, rm 10/49.3 GB, plus their Ribo arms). No dataset hunting
   required, and it is the cleanest generalisation test: same study, same protocol, same tissue,
   only the species differs. That controls for batch effects in a way cross-study comparisons cannot.
2. **Yeast next** -- smallest genome, fastest index, densest ORF annotation, and the organism where
   Ribo-seq periodicity is cleanest. A good smoke test of the whole pipeline on a non-vertebrate.
3. **Worm / fly / fish** once their matched pairs are actually confirmed.

## D. Open questions before committing compute

- Zebrafish and Drosophila matched Ribo arms are **not yet located**. Do not queue downloads for
  these until the companion BioProject is confirmed.
- No yeast study with replicates + matched RNA has been vetted yet.
- The model's ORF-candidate track assumes ATG-initiated ORFs with a 5'UTR. Yeast 5'UTRs are short
  and its genome is intron-poor; worth checking the track is not degenerate there before drawing
  conclusions from a yeast result.
