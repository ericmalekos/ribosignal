# Are the model's extra ORF calls hidden by multimap filtering?

Mouse liver, 2026-08-21. **Question:** the model calls ORFs that real Ribo-seq misses. Are those
false positives, or is the real data blind to them because `--outFilterMultimapNmax 1` discarded the
reads?

**Status: three orthogonal correlational lines all say the blindness is real. A direct test
(mm25 re-alignment) is in flight and is the only one that can settle it.**

---

## 0. A 4x counting error, corrected

The first pass reported **9,181** attn-specific calls. The real figure is **2,284**. The error was
loading `*_collapsed.txt` by hand instead of through `scripts/compare_dropin_calls.build_loader`,
which applies four filters:

- restriction to the model's test transcripts
- `pval <= 0.05`
- ORF length `>= 90 nt`
- **predicted density over the ORF `>= 0.5x` uniform** -- predicted side only

That last one exists precisely to drop weak predicted calls, which is the population being counted.
`build_loader`'s own docstring warned about this: *"Two tables that disagree because one of them
quietly skipped the enrichment filter is a manuscript-level hazard."*

**Lesson: use the shared loader. Do not reimplement the load.** Two further keying subtleties
matter and are handled by it: the poster keys on `(gene_id, ORF_gstop)` (genomic), not `ORF_ID`,
because one genomic ORF appears once per transcript it sits on -- the same isoform-multiplicity
inflation that shows up everywhere in this project.

A conclusion flipped on correction: unfiltered, model-model Jaccard on model-specific calls was
0.415 vs 0.682 overall, which read as "the extras are architecture-specific noise". Filtered it is
**0.568**, and attn-only shrank 5.6x. The architectures agree substantially more than first reported.

### Corrected counts (gse243134 RNA arm, genomic keying, filtered)

| | calls |
|---|--:|
| REAL union over all 3 liver datasets | 13,170 |
| attn predicted / model-specific | 14,035 / **2,284** |
| mamba4 predicted / model-specific | 13,198 / **1,805** |
| **both models, no real support** | **1,481** |
| attn only | 803 |
| mamba4 only | 324 |

---

## 1. Three orthogonal measures, all pointing the same way

Enrichment vs REAL-supported genes. **No measure derives from another**: one is a mapping outcome,
one is sequence content, one is annotation coordinates.

| measure | BOTH | mamba4-only | attn-only | REAL |
|---|--:|--:|--:|--:|
| **pseudogene locus overlap >10%** | **17.4x** | 16.9x | 8.5x | 1.0x (0.5% of genes) |
| **multimap depletion** (expr-controlled) | **6.8x** | 7.3x | 4.0x | 1.0x |
| **exonic repeat >25%** | 3.7x | **5.4x** | 3.1x | 1.0x |
| **LTR bp per gene** | 3.1x | **5.2x** | 3.0x | 1.0x |

All p < 1e-16. Cross-architecture and mamba4-only calls rank high on all three; attn-only is
consistently about half.

**Composition:** model-specific calls are **53-63% lncRNA** against **4%** for real-supported
(13x), and 32% `Gm` predicted gene models against 2.5%.

### 1a. Multimap depletion, and how it is measured

`log2((pack coverage per nt + 0.1) / salmon TPM)`, median over the 3 liver datasets, residualised on
`log10(TPM)`. Salmon's EM **uses** multimappers; the pack's coverage is STAR `mm1`, which
**discards** them, so the ratio is the signal lost to filtering.

Validity: length confound Spearman **+0.014** (clean). Expression confound was **-0.368**, removed
by residualising (to +0.112) with the enrichments essentially unchanged -- and it ran *against* the
conclusion anyway, since the depleted set is the *more* expressed one (median TPM 12.5 vs 4.7) with
**53x less coverage** (0.75 vs 39.3 per nt). The mouse quants are decoy-aware, which biases
*against* the finding, so true depletion is likely larger.

Weakness: 21.9% of the extreme tail sits at the pseudocount floor (zero coverage), so its exact
value is set by the `+0.1`, not measured. The absolute score is meaningless; only between-set
comparison within the same datasets is valid.

### 1b. Repeat content is a THRESHOLD effect, not a linear one

| exonic repeat | genes | median depletion |
|---|--:|--:|
| zero | 2,562 | -0.029 |
| 0-5% | 3,139 | +0.092 |
| 5-15% | 2,238 | +0.003 |
| **15-30%** | 945 | **-0.287** |
| **>30%** | 688 | **-1.158** |

Below ~15% exonic repeat there is no depletion at all; above it, depletion appears and accelerates.
That is how multimapping should behave. **It also explains why the global Spearman is only -0.125**
(p=1.1e-34): 60% of genes have under 5% repeat and contribute noise. Quoting -0.125 alone would
understate a real effect concentrated in the top ~17% of genes.

### 1c. The pseudogene test I got wrong the first time

First attempt tested the **transcript biotype label** and found zero pseudogenes, concluding the
hypothesis was untestable because the universe is protein-coding + lncRNA only. **That was wrong.**
A pseudogene locus commonly persists in the annotation as an *overlapping lncRNA* with the same
sequence, so the label is absent while the locus is present. Testing by **genomic overlap** against
all 13,809 pseudogene loci finds it immediately, and it is the strongest of the three signals
(17.4x). Credit for the correction to the user.

---

## 2. The direct test (in flight)

`scripts/align_gse243134_mm25.sbatch` + `ribocode_gse243134_mm25.sbatch`, chained
`--dependency=afterok`. Re-aligns the 21 gse243134 RPF runs -- the **deepest** liver reference,
150.5M P-sites at mm1 -- with `--outFilterMultimapNmax 25`, everything else byte-identical, then
runs RiboCode with the identical invocation.

**Already measured: 23-27% of reads map to multiple loci**, against 47-49% unique. mm1 discarded
all of them.

**This DELIBERATELY violates the project's Ribo-seq rule.** RiboCode does not drop multimappers, so
mm25 P-site counts are inflated by construction. Output is confined to `data/mm25_diagnostic/` and
must never feed a pack, a training run, or a headline table. The question is not "are these counts
correct" but "do the missing ORFs APPEAR", which is robust to the inflation.

**The number to read: of the 1,481 both-models model-specific calls, how many appear in mm25?**

---

## 3. What this does and does not establish

Even a positive mm25 result shows the real data was **blind**, not that translation occurs. It
removes the strongest evidence *against* these calls; it does not supply evidence *for* them.
Immunopeptidome MS on the ~350 both-models severely-depleted genes is the arbiter, and that is a
far better-targeted question than testing all of them.

Scope: one RNA arm (gse243134), mouse liver only. The human held-outs are untested.

---

## 4. Infrastructure failure worth recording (2026-08-21)

Job 36834132 (gse243134-only mm25 RiboCode) died after **92 minutes** with
`OSError: truncated file` on `SRR26055549.sorted.bam`.

**Cause: one-writer-per-path violated.** The combined MM=25 job was launched while the
gse243134-only job was still running, and BOTH sort to the same `<run>.sorted.bam` paths. The
second truncated a file the first was mid-fetch on.

**The deeper defect, which was project-wide.** The sort-skip guard tested only the BAM header:

```bash
[[ "$(samtools view -H "$b" | sed -n 1p)" == *SO:coordinate* ]]
```

**A BAM truncated mid-write keeps a perfectly valid header**, so this guard accepts it and the
failure surfaces much later inside RiboCode as an opaque pysam error. Audit found the pattern in
5 scripts -- and inverted: the scripts that WRITE sorted BAMs (`sort_*_bams.sbatch`) already used
`samtools quickcheck`, while every CONSUMER checked the header only, including the live
`ribocode_gse243134.sbatch`.

Hardened with `samtools quickcheck` (verifies the EOF block): `ribocode_gse243134.sbatch`,
`human_01_pack.sbatch`, `riboseq_align.sbatch`, `ribocode_gse243134_mm25.sbatch`,
`ribocode_liver_combined_mm.sbatch`. `align_wang_ribo.sbatch` left alone -- already banner-marked
provenance-only.

**Generalisable rule: a header check is not a validity check.** Validate a BAM you did not just
write with `quickcheck`, not by reading its header.

The failed job was NOT resubmitted: the combined run supersedes it, covering all three liver
datasets rather than one, which is what the model-specific definition requires anyway.
