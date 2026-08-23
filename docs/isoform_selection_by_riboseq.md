# Can a gene's representative isoform be selected by Ribo-seq?

**Answer: no.** Under this project's multimapper convention, Ribo-seq is close to constant across a
gene's isoforms, so "pick the isoform with the most Ribo-seq" selects on noise and, in the worked
example, actively picks the wrong transcript. RNA-seq carries the information; Ribo-seq does not.

Written 2026-08-20 from analysis run that day. Every number below traces to an artifact named in
the "Where the numbers come from" section. Figures in `figures/X_gtf2i_isoforms/`.

---

## 1. Why this is even a question

The universe is **isoform-level, not gene-level**: 84,472 transcripts across 18,828 genes, a mean of
**4.49 isoforms per gene**, up to 69 for a single gene. Only 5,679 genes contribute one transcript.

Under **posture A** (methods.md 2.3) a footprint that is genomically unique still maps to every
compatible isoform of its gene, and is counted at weight 1 on all of them. Cross-gene multimappers
and non-coding transcripts are dropped first, so the ambiguity is always within one gene. This is
what produces the ~18x project-wide inflation of summed per-transcript P-site totals.

An obvious-looking fix is to keep one isoform per gene and pick it by the Ribo-seq itself. This
document is why that does not work.

---

## 2. The worked example: GTF2I in held-out hepatocytes

GTF2I has **35 isoforms** in the union universe, 10 to 35 exons each (median 33), and is
deliberately a hard case: its inflation is **34.2x**, against the ~18x project average.

| quantity | value |
|---|---|
| summed P-sites over 35 isoforms | 452,992 |
| per-isoform P-site range | 11,161 to 13,258 (**16% of max**) |
| RNA TPM range | 1.11 to 50.31 (**45x**) |
| Spearman(RNA TPM, P-sites) | **-0.248** (p = 0.151) |
| Pearson(log TPM, P-sites) | -0.135 |

**RNA separates these isoforms 45-fold. Ribo-seq separates them by 16%, with no usable rank order.**

### It picks the wrong isoform

| | transcript | P-sites | RNA TPM |
|---|---|--:|--:|
| argmax by **Ribo-seq** | `ENST00000901263.1` | 13,258 | **1.25** |
| argmax by **RNA** | `ENST00000620879.4` | 13,138 | **50.31** |

The genuinely dominant isoform ranks **14th of 35** by Ribo-seq count. Selecting on Ribo-seq
discards it in favour of one expressed **40x lower**.

### Periodicity does not rescue it

The natural next idea is that the real isoform should show cleaner 3-nt periodicity even if counts
are tied. It does not. Frame-0 fraction inside the CDS, computed per isoform relative to **that
isoform's own annotated CDS start**:

| min | median | max |
|--:|--:|--:|
| 85.3% | 85.6% | 85.8% |

A **0.5 percentage point** spread across isoforms differing 45x in abundance. The CDS starts differ
by hundreds of nucleotides (317 to 512 among the top ten) but those offsets are near-multiples of 3,
so the isoforms are in register with one another and the same reads land in frame 0 regardless of
which start you measure from.

### The information is not there to be extracted

How many observed P-sites sit at a position whose k-mer is **unique to one isoform**, i.e. could
ever discriminate:

| k | isoform-unique positions | P-sites there | isoforms with any |
|--:|--:|--:|--:|
| 19 | 0.73% | 123 / 452,992 (**0.027%**) | 2 / 35 |
| 25 | 0.79% | 498 (**0.110%**) | 3 / 35 |
| 29 | 0.83% | 527 (**0.116%**) | 3 / 35 |

At k=29, the longest window a 29-nt read supports, **32 of 35 isoforms have zero discriminating
reads**. No estimator recovers what is not in the data.

**Caveat on this measurement.** A read is counted as discriminating if the k-mer *starting at its
P-site* is unique. A real footprint spans roughly 12 nt upstream to 17 downstream of the P-site, so
the window is shifted rather than exact. This will not move 0.1% to the ~10% that would be needed,
but it is an approximation, not a precise figure.

---

## 3. The sub-question: run salmon on Ribo-seq?

Salmon's EM is the principled version of fractional assignment ("posture C"). Three obstacles:

1. **k-mer size, correctly identified.** Ribo-seq mode read length here is **26 to 29 nt**; both
   existing indexes are **k=31**, so they would match nothing. A smaller-k index is mandatory.
   But **lowering k makes discrimination worse, not better**: unique positions fall from 0.83% at
   k=29 to 0.73% at k=19. The k reduction buys mappability at the cost of the little specificity
   there is.
2. **The effective-length model assumes near-uniform coverage.** Ribo-seq is concentrated in the
   CDS, near-absent from UTRs, and carries 3-nt periodicity plus codon-specific pausing. The length
   correction would be systematically wrong in a way sequence-bias correction does not address.
3. **UTR-distinguishing isoforms are unresolvable in principle.** `ENST00000443166.5` differs
   largely by a long terminal UTR exon. Ribo-seq has essentially no UTR coverage, so no EM recovers
   that distinction at any k.

The shape that could work is the inverse: **let RNA-seq resolve the isoforms and apportion the
Ribo-seq accordingly.** Worth checking the Kingsford group's Ribomap, which I believe does exactly
this (verify before citing). The project's `make_representative_tx.py` is already a hard 100/0
version of the same idea.

---

## 4. The distinction that actually matters

**Filtering by signal strength is sound, and the pipeline already does it.** `min_train_signal = 50`
drops transcripts carrying too few P-sites to learn a shape from: 42,000 training transcripts
against 70,883 at test, where the filter is removed. That asks *"is there enough evidence here?"*,
which Ribo-seq answers well.

**Selecting among isoforms asks a different question** -- *"which of these near-identical sequences
is the real one?"* -- and Ribo-seq under posture A is constitutionally unable to answer it, because
the convention copies the same reads onto every compatible isoform.

This also disposes of an objection raised and then withdrawn during the discussion: that selecting
on Ribo-seq is circular and undeployable. It is **not**, provided the selection is applied to the
training set only and inference runs on everything passing an RNA threshold -- exactly the pattern
`min_train_signal` already uses. Using labels to curate training data is what labels are for. The
problem is not circularity. The problem is that the selection is uninformative.

---

## 5. What it costs when done properly

Task #100 measured train-restricted / predict-broad directly, scoring the posture-B checkpoints on
the **full** test set so training-composition is isolated from eval-composition:

| architecture | pc profile r (B - A) | pc count r (B - A) |
|---|--:|--:|
| baseline | -0.0157 | -0.0429 |
| orf_v2_attn (deployed) | -0.0155 | **-0.0682** |

**The count head takes the damage**, and more of it in the deployed architecture. Training only on
dominant isoforms never shows it the low-abundance end of the range it must predict at inference, so
it is asked to extrapolate below everything it saw.

That task also found that **Task 21's reported gains were eval-composition artifacts**: the claimed
count-r rise of +0.0158 is **-0.0682** on a common test set, and attn-B's uORF 0.6546 ("best of all
four") is **0.6123**, below A's 0.6258. Task 21's conclusion that posture A stands survives, and for
a stronger reason -- A now beats B on a common denominator rather than merely matching it.

---

## 6. Open: is GTF2I representative?

**This is one gene, and a worst case.** 35 isoforms against a project median of 4.49; 34.2x
inflation against ~18x. Genes with two or three isoforms may well have real discriminating coverage.

The genome-wide version is cheap, because uniqueness only has to be computed within each gene's own
isoform set -- 18,828 small independent problems. It would give the fraction of genes for which
isoform-level Ribo-seq assignment is even possible, and should precede any decision to build a
low-k index. **Not yet run.**

---

## Where the numbers come from

| claim | source |
|---|---|
| universe size, isoforms per gene | `data/union_universe.tsv` |
| per-isoform P-sites, GTF2I | `data/packed_union_Hepatocytes/{target_counts,offsets,lengths,tx_order}` |
| RNA TPM | `data/union_universe.tsv`, `max_tpm` |
| CDS starts, frame-0 fractions | `data/human_ribocode_annot_primary/transcripts_cds.txt` |
| exon structure, junctions | `genomes/gencode.v49.primary_assembly.annotation.gtf`; exon lengths verified to sum EXACTLY to packed transcript length for all 35 |
| k-mer uniqueness | `data/union_universe.fa` |
| read lengths | `results/dataset_census/dataset_census.tsv`, `mode_read_length` |
| index k | `Salmon_indexes/*/info.json` |
| posture-B cost | `results/postureB/*/extra_metrics.json` (task #100) |
| `min_train_signal`, split sizes | `figures/arch_attn/arch_spec_attn.json` |

Figures: `figures/X_gtf2i_isoforms/` -- `X_gtf2i_isoforms` (counts vs RNA),
`X_gtf2i_frames{,_top10}` (frame-coloured per isoform), `X_gtf2i_model_top10` (with exon structure).
All carry `model: null` in their values files -- observed Ribo-seq only, no prediction involved.

---

## 7. REVISION 2026-08-20: per-tissue TPM was recovered, and it changes two claims above

Everything in sections 1-6 stands. Two things need correcting, and one gets much stronger.

### 7.1 "RNA-seq carries the isoform information" is only half true

Section 2 contrasts Ribo-seq (uninformative) against RNA (45x spread). That 45x is **salmon TPM**.
The **RNA the model actually receives** -- the pack `coverage.npy` channel -- is built from the same
STAR transcriptome BAMs with within-gene isoform multimappers retained, so it inherits posture A and
is **also isoform-ambiguous.** Measured on GTF2I across all 8 tissues:

| | |
|---|---|
| \|Spearman(pack coverage, salmon TPM)\| | **<= 0.17 in every tissue** |
| Spearman(length, coverage per nt) | **~ -0.8** |
| total coverage spread across 35 isoforms | only **1.7-2.7x** |

**So BOTH the model's RNA input and its Ribo-seq label are isoform-ambiguous.** Salmon TPM is the
only isoform-discriminating quantity in the pipeline and is never a model input -- it is used for
universe definition and representative selection only. The model has no channel that could
distinguish a gene's isoforms even in principle. That is a stronger statement than section 2 makes,
and it strengthens the conclusion rather than weakening it.

### 7.2 Per-tissue TPM was recoverable after all -- no download needed

An earlier draft treated the per-tissue values as gone. They were recovered on 2026-08-20 from 55
surviving STAR transcriptome BAMs via alignment-mode salmon, with no re-download. See `results.md`
and `methods.md` (2026-08-20 sections). Decoy-free, so validated against the surviving decoy-aware
Fibroblast table: **Spearman 0.755 on expressed transcripts, median 0.84x, no length bias.** Good
for ranking, not for absolute TPM.

### 7.3 The GTF2I answer, per tissue

`ENST00000443166.5` ranks **1st in all 8 tissues**, carrying **66-90%** of the gene, with
`ENST00000901300.1` second everywhere at 10-33%. Only 2-4 of 35 isoforms clear 0.1 TPM. **GTF2I is a
two-isoform gene in every Chothani tissue.** Full table in `results.md`; machine-readable at
`data/tpm/chothani_alnmode/GTF2I_per_tissue.json`.

**`max_tpm` is a cross-tissue maximum and inverts within-tissue rankings.** It put 443166.5 4th in
hepatocytes behind three isoforms; per-tissue values put it 1st, and the isoform whose `max_tpm` is
50.31 does not appear in the hepatocyte top five at all. Never use `max_tpm` for a per-tissue
question.

---

## 8. Why picking an isoform by Ribo-seq alone is so hard (synthesis, 2026-08-21)

To tell two isoforms apart you need reads that land **where the isoforms actually differ**, and
enough of them. Ribo-seq fails that test on five independent counts, and they compound.

**1. The counting convention erases the difference before analysis starts.** Under posture A a
footprint compatible with 35 isoforms is counted on all 35. The per-transcript numbers are
near-identical *by construction*: GTF2I spans 16% across its isoforms while RNA spans 45x.

**2. The reads are too short.** Footprints are ~29 nt. Only **0.116%** of P-sites sit at a position
unique to one isoform, and **32 of 35 isoforms have zero** such reads. A 75 nt paired RNA read
covers more sequence and its mate brackets a whole fragment, so it constrains far more.

**3. It only covers coding sequence.** Isoforms frequently differ in UTRs -- alternative first
exons, alternative 3' ends. Ribo-seq has essentially no UTR signal, so those differences are
invisible *in principle*, not merely hard.

**4. Periodicity does not break the tie.** The natural hope is that the genuinely translated isoform
shows cleaner 3-nt rhythm. It does not: frame-0 fraction is **85.3-85.8% across all 35**, a
half-point spread. The isoforms' CDS starts differ by near-multiples of 3, so they sit in register
with one another and the same reads look in-frame whichever isoform you measure against.

**5. The discriminating reads are exactly the fragile ones.** At paralog-rich loci the reads that
separate isoforms are also the reads that multimap across genomic copies -- so a multimapper filter
destroys isoform resolution while leaving gene totals intact. This is not hypothetical: it is what
broke `chothani_alnmode` (see results.md 2026-08-21). And Ribo-seq pipelines apply that filter as
standard -- this project's own rule is `--outFilterMultimapNmax 1` for every Ribo-seq alignment.

**And the statistical machinery does not transfer.** Salmon's EM leans on a fragment-length
distribution and roughly uniform coverage. Ribo-seq has neither: density is concentrated in the CDS,
carries 3-nt structure, and shows codon-specific pausing. The tool that solves this for RNA cannot
simply be pointed at Ribo-seq.

**The constructive framing.** Isoform identity is a question about *where a transcript differs*;
translation level is a question about *how much ribosome occupies it*. RNA-seq can answer the first
and Ribo-seq answers the second. Asking Ribo-seq to do both is asking it to resolve differences it
is structurally blind to. That is why `make_representative_tx.py` ranks on RNA TPM, and why the
right architecture for a combined analysis is **RNA resolves the isoform, Ribo-seq quantifies the
translation on it.**
