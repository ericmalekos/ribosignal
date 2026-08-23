# Ribo-seq signal model: figure package for poster assembly

Everything needed to build a poster: every figure, the exact numbers behind it, the provenance note
for each, and the generator that drew it. 39 figure folders, 318 files, ~31 MB.

Assembled by `scripts/build_poster_package.py --model attn`. Re-run it to refresh.

---

## 0. How to use this package

```
figures/<name>/
  <name>.pdf / .png          the figure
  <name>_values.json         EVERY number plotted, written by the generator AFTER drawing
  <name>_values.tsv          the same as a flat table
  FIGURE_DATA_INPUTS.md      where the data came from + caveats a caption must carry
  make_<name>.py             the generator
DATASET_LABELS.md            internal key -> display label -> accession, for every dataset
dataset_registry.tsv         the source of truth that table is generated from
PROJECT_*.md                 project-level docs (layout, status, policy, results, methods)
MANIFEST.json                machine-readable index
```

**Use `DATASET_LABELS.md` when captioning.** Filenames and `*_values.json` keys carry INTERNAL keys
(`janich`, `cart`, `gse208041`, `Hepatocytes`); the poster should print display labels
(`mouse_liver_GSE67305`, `human_cart_GSE304796`, `human_thp1_GSE208041`,
`human_hepatocytes_GSE182371`). The internal keys are deliberately not renamed -- they are
load-bearing across paths, scripts and every values file.

**What is portable and what is not.** The `*_values.json` / `.tsv` files are fully portable: every
plotted number is in them, so any figure can be restyled, recombined or replotted from this package
alone. 31 of 39 folders ship a `*_values*` file; the other 8 (`benchmarks`, `depth_crossover`,
`kozak`, `o2_venn`, `prediction_examples`, `profile_exemplars`, `S_encoding`) write richer NAMED
TSVs instead, because their plotted data is per-nucleotide or per-region rather than a summary.
`main/` is the only folder with no data file -- it is a composed figure built from other panels'
values. The **generators are not** -- they read absolute paths into the project tree (packs, search
results, drop-in dumps) and need the full ~20 TB project to re-run end-to-end. Edit styling in the
generator, but replot from the values files if you are working off-cluster.

**Read `PROJECT_POSTER_SESSION_HANDOFF.md` first if you drafted panels before 2026-08-15.** A549 was
removed from every figure and the held-out human numbers moved, so some drafts are built on data that
no longer exists. Section 0 of that file lists exactly what changed.

**Then `PROJECT_POSTER_LAYOUT.md`.** It is a generated 4-column layout with live numbers, and every
panel carries a **MUST STATE** line. Those are not caption garnish: most exist because the broader,
more eye-catching version of the claim was measured and found false.

---

## 1. What the project is

**The problem.** Ribosome profiling (Ribo-seq) measures which parts of the transcriptome are being
translated, at nucleotide resolution. It is the only assay that directly observes translation, and it
is how small ORFs, upstream ORFs and translated lncRNAs get discovered. But it is expensive, needs
fresh material, and needs deep sequencing to see anything outside abundant coding genes.

**What this model does.** It predicts the per-nucleotide Ribo-seq P-site density of a transcript
**from sequence and RNA-seq alone** -- no ribosome profiling for the query. Feed the predicted
profile into a standard ORF caller (RiboCode) and you get ORF calls for a sample you never
Ribo-seq'd.

**Why that is useful.** RNA-seq is cheap and ubiquitous. If predicted profiles support ORF calling,
then any RNA-seq dataset becomes a candidate for translation analysis, and a proteogenomics search
database can be built for a sample where ribosome profiling was never done.

**Inputs and outputs.**

| | |
|---|---|
| input | one-hot sequence (4 ch) + ORF track (5 ch: frame 0/1/2 occupancy, start context, stop) + RNA-seq coverage (1 ch) |
| output | per-nucleotide P-site density, plus a total-count head |
| size | ~5M parameters, CPU-runnable at inference |
| the deployed arm | **standalone**: no Ribo-seq at inference, at all |

Two prediction arms are reported everywhere, and both must always be shown:

- `pred_obsdepth` -- predicted SHAPE scaled by the dataset's observed sequencing depth
- `pred_preddepth` -- **fully standalone**, depth predicted too. This is the arm the application uses.

---

## 2. Which model, and why the transformer

Two architectures were trained on identical data: a **transformer** (`attn`, 2 attention layers) and
a **Mamba state-space model** (`mamba4`, 4 attention layers). They are statistically
indistinguishable, replicated three times independently:

| benchmark | attn | mamba4 |
|---|--:|--:|
| B721.221 vs 327M measured footprints (F1) | **0.6661** | 0.6654 |
| macrophage cross-subtype proteomics | tie | tie |
| 4-dataset immunopeptidome panel | tie | tie |

**This package emphasises the transformer (`attn`).** It costs nothing scientifically, it is the
architecture most readers know, and it is the one that runs on CPU -- `mamba4` requires CUDA
(`mamba_ssm` dispatches to `causal_conv1d_cuda`), which is a real deployment constraint.

### Transformer data availability -- the gap is CLOSED as of 2026-08-15

| panel group | attn available? |
|---|---|
| held-out human validation (P5: 4 arms) | **yes** -- `P5_heldout_human_attn.*` |
| mouse cross-study 3x3 (E1, E2, E3, P9, B9) | **yes**, full 9-cell grid |
| ground-truth benchmark (A4) | **yes**, both calibration arms |
| interpretability (A1, A2, B4, B5, C11) | **yes** |
| **human** immunopeptidome proteomics (F2b, D12, D13) | **yes**, all 5 datasets searched with both |
| **mouse macrophage** proteomics (P11, P12, P13) | **yes** -- all 12 populations re-called and re-searched |

An earlier build of this package said the macrophage sweep was mamba-only and recommended dropping
it. **That is no longer true.** All 12 populations were re-called from attn predictions and
re-searched at the correct 30-aa tryptic floor, 7 arms per population against identical spectra. Both
architectures are in the same reports, so the comparison is like-for-like.

**The transformer is also the better database here**, on merit rather than familiarity:

| median across 12 populations | mamba4 | **attn** |
|---|--:|--:|
| database size | 1,252 | **1,062** |
| novel peptides | 27 | 26 |
| **discovery density (per 1k seqs)** | 17.9 | **21.5** |

attn finds the same peptides from a **15% smaller database**, so density rises 20%. Its 15-20% fewer
ORF calls (consistent across all 12 populations) were disproportionately the ones yielding no
confident peptide -- consistent with B9's finding that ~92% of experiment-unsupported calls are false
positives. Do NOT rank the two on delta-total peptides; that gap is inside the +-50 search-parameter
swing. Density is the robust axis.

**So the poster can be all-transformer with no gaps and no omissions.**

**P5 renders both.** `P5_heldout_human.*` is mamba4 (locked decision D1b, the main-figure default);
`P5_heldout_human_attn.*` is the transformer. **Use the `_attn` files for this poster.** Every panel
now carries its model in its own `*_values.json`, because an earlier version of this README printed
P5's mamba4 numbers under an "attn" heading.

---

## 3. Methods

### Training data

Human Ribo-seq across **7 cell types, 69 libraries** (Chothani et al.), pooled per cell type.
**Hepatocytes (5 libraries) held out entirely** as a leave-one-tissue-out (LOTO) test. Brain was
dropped upstream for low periodicity (period_obs 0.044) before any data was built.

- universe: **84,472 transcripts** -- protein-coding + lncRNA, salmon TPM >= 1, <= 10,000 nt,
  **chrM excluded** (the 13 mitochondrial mRNAs use a different genetic code and a different ribosome)
- **11x depth spread** across training cell types (Fibroblast 32 libraries, HUVEC 3). Training caps
  transcripts per tissue rather than pooling, and LOTO holds out a TISSUE -- a random transcript
  split would let the deepest tissue dominate both sides.

### Ribo-seq processing (one final recipe, enforced)

Every dataset in this package goes through one pipeline. Deviations are declared, never silent
(`PROJECT_PIPELINE_POLICY.md`).

1. adapter trim (`cutadapt`), footprints 20-40 nt
2. STAR with **`--alignEndsType EndToEnd`** and **`--outFilterMultimapNmax 1`** (unique mappers only)
3. drop rRNA / tRNA / miRNA / Mt_rRNA / Mt_tRNA loci, and cross-gene multimappers
4. UMI dedup where the protocol has UMIs
5. RiboCode `metaplots` for per-library P-site offsets, then per-nucleotide P-site counts

**Why EndToEnd matters:** soft-clipping displaces P-sites systematically, not randomly. Off-recipe
alignments inflated non-canonical ORF calls by ~22%, so off-recipe and final-recipe numbers are never
compared.

### ORF calling and scoring

Calls are made permissively (`--min_aa 5`) and filtered at **scoring**, so the model's calls and the
reference's calls are cut at the same place by one code path:

- ORF length >= **90 nt = exactly 30 amino acids** (`ORF_length` excludes the stop codon)
- p <= 0.05; predicted calls additionally need mean predicted density >= 0.5x uniform
- keyed **genomically** on `(gene_id, ORF_stop)`, so choosing a different representative isoform of
  the same ORF is not scored as a disagreement
- every table states **n_ref**, the size of the reference set it was scored against

### Proteogenomics

Predicted ORFs become a mass-spec search database, searched against real spectra alongside a
GENCODE-only baseline and naive-enumeration nulls.

- novel peptides at **1% class-specific FDR** (novel targets vs their own decoys); canonical columns
  at 1% global FDR. Mixing the two is how a bigger database looks better than it is.
- **ORF-length floor is set by ASSAY**: tryptic whole-cell lysate **30 aa**, MHC/HLA
  immunopeptidomics **7 aa** (HLA-I peptides are 8-11 aa). Enforced in code.

---

## 4. The figures

Grouped by the poster narrative. Headline numbers are **attn** where attn exists.

### Column 1 -- what it is and what it was trained on

**`arch_attn`** -- model architecture. One-hot sequence + 5-channel ORF track + RNA-seq coverage ->
per-nucleotide P-site density. ~5M parameters.
*Must state:* the standalone arm uses NO Ribo-seq at inference. "Predicts ribosome density" is
ambiguous about that, and the ambiguity is the whole application.

**`P3_training_data`** -- 7 cell types, 69 libraries, Hepatocytes held out, 84,472-tx universe, 11x
depth spread. Panel B lists the inference inputs explicitly.
*Must state:* Hepatocytes is coloured differently because it is the LOTO holdout, not because it was
trained on.

**`S_riboseq_qc`** -- library quality for all 15 dataset arms (8 training + 5 held-out + 2 cell
lines): frame-0 fraction 74.6-87.7% against a 33% no-periodicity floor.
*Must state:* establishes that a downstream failure is the model's, not the library's.

**`S_encoding`**, **`prediction_examples`**, **`profile_exemplars`** -- input encoding and observed
vs predicted per-nucleotide profile exemplars (uORF / CDS / lncRNA).

### Column 2 -- held-out human validation (the central claim)

**`P5_heldout_human`** -- four independent held-outs that fail in *different* ways: same-study
withheld, cross-study, cross-lab cell line, primary engineered cells. **attn, standalone arm:**

| arm | F1 | annotated CDS recall | non-canonical recall | n_ref |
|---|--:|--:|--:|--:|
| hepatocyte (LOTO) | **0.911** | 0.997 | **0.594** | **16,236** |
| iPSC-CM (Ruiz-Orera, cross-study) | **0.875** | 0.997 | **0.550** | **13,087** |
| THP-1 (GSE208041) | **0.831** | **0.997** | **0.525** | **13,225** |
| CAR-T (GSE304796) | **0.886** | 0.994 | 0.580 | **8,253** |

CAR-T is now on its **canonical** pack (it was 0.904 / n_ref 8,856 on a pack with ~29% PCR
duplicates). The reference lost duplicate-driven calls -- n_ref 8,856 -> 8,253 -- and the model
barely moved, which is the reassuring direction.

*Must state:* panel B is the honest half and must not be cropped. Overall F1 is dominated by
annotated CDS, so a single number reads ~0.9 everywhere and hides that non-canonical recall is
~0.55. **THP-1 is the one arm still on its original pack** -- see section 6.

**`A4_b721_ground_truth`** -- the only NON-null benchmark. B721.221: predictions from RNA-seq alone
vs **327 million measured** Ribo-seq footprints, 16,331 measured calls in scope. Scored against the
assay's own split-half ceiling (two disjoint halves of the same Ribo-seq agree at F1 0.889,
annotated 0.973, **non-canonical only 0.498**):

| | attn theta=1 | % of ceiling |
|---|--:|--:|
| F1 | 0.666 | **75%** |
| annotated CDS recall | 0.759 | 78% |
| non-canonical recall | 0.267 | **54%** |

*Must state:* never quote the 0.267 without the 0.498 denominator. Two independent measurements of
the same cells agree on only half of each other's non-canonical calls, so the model recovers about
half of what a replicate EXPERIMENT would -- not a quarter of perfect.

**`A1_localization_ceiling`** -- predicted translation localization beats the observed periodicity
ceiling on all-ORFs. *Must state:* on non-canonical ORFs it reaches 94-96% of ceiling but does NOT
beat it. Say both.

**`A2_ribocode_dropin`** -- the predicted profile dropped into RiboCode reproduces that dataset's own
ORF calls, with no Ribo-seq for the query.

**`A3_loto_spread`** -- cross-tissue generalization across 9 folds; the spread tracks held-out target
quality, not a gradient in model skill.

### Column 3 -- mouse cross-study: does it transfer, and where does it fail?

Three independent mouse-liver Ribo-seq experiments (Janich, GSE243134, Wang) on one pipeline and one
universe, so a difference is between EXPERIMENTS, not code paths.

**`E1_rna_provenance`** -- matched vs mismatched RNA-seq input across the 3x3 factorial.

**`E2_class_ceiling`** -- model vs the experiment-vs-experiment ceiling, per ORF class. Annotated CDS
sits near its ceiling; non-canonical classes do not. *Must state:* the ceiling is not 1.0; scoring
against 1.0 charges the model for irreproducibility in the reference.

**`E3_shallow_experiment_win`** -- the predicted profile calls ORFs *better* than measured Ribo-seq
below ~29M P-sites. This is the "why not just do Ribo-seq" answer, and it is depth-conditional.

**`P9_upset_mouse3x3`** (8 outputs: 2 model families x 4 ORF-class views) -- where the disagreement
lives. **`B9_overcall_validation`** -- and how much of it is real.

The model calls ~2,000 ORFs no single experiment supports. A **28-library merged** RiboCode reference
(all three datasets pooled, 26,388 calls) arbitrates:

| stratum | n | corroborated |
|---|--:|--:|
| model calls an experiment DID see (anchor) | 11,441 | **99.3%** |
| observed calls unique to ONE experiment (**control**) | 1,054 | **74.5%** |
| **model-only calls** (attn) | **2,056** | **8.6%** |
| called by BOTH architectures, no experiment | 1,504 | 8.7% |

*Must state:* **~92% of the model-only calls are false positives**, and the control must travel with
the number -- a genuinely real, singly-observed ORF corroborates only 74.5% of the time either.
**Cross-architecture agreement is NOT evidence** (8.7%, indistinguishable from one model); it was
tested as a confidence filter and failed.

### Column 4 -- proteogenomics: what a better ORF database is worth

**Use these three for an all-transformer poster.**

> **A549 WAS REMOVED FROM EVERY FIGURE ON 2026-08-15.** It was the project's only tryptic human
> dataset, so **this entire column is now HLA-I immunopeptidome only**. Nothing here may be described
> as a whole-proteome or tryptic result, and none of it generalises past HLA-I. A549 also carried the
> results least favourable to the model, so several claims got narrower, not stronger -- each is
> spelled out below. The archived data and the numbers it produced are in
> `proteogenomics/data/_archive_a549_2026_08_15/`.

**`F2b_discovery_forest`** -- 4 HLA-I datasets x 2 models x 2 arms = **16 points**, all above 1.0;
the model's database is a median **93.6x denser** in discoveries per sequence than naive enumeration
(range 11-2744x). Panel (b) is the absolute bottom line: change in TOTAL unique peptides vs searching
GENCODE alone -- the Poisson arm is positive in **8 of 8**.
*Must state:* panel (a) is a ratio; panel (b) is the absolute total and must be shown beside it.
**8/8 is not a strengthening.** The old 9/10 included the only two points where the arms disagreed
(A549: attn -35, mamba4 +11); removing that dataset removed the disagreement, not the doubt.

**`D12_cpat_cpc2`** -- against CPAT and CPC2, the standard sequence-only coding-potential selectors
(the honest competitor; "better than enumerating every ORF" is a weak claim on its own), on **3 HLA-I
immunopeptidomes**.
*Must state:* **the assay split is gone.** This panel used to show CPAT/CPC2 winning the tryptic
whole proteome and the model winning HLA-I; the tryptic half was A549 and is the half the model
LOST. Caption this as an HLA-I result and never as "the model beats CPAT/CPC2". The HLA-I advantage
also depends on sub-30-aa ORFs -- 27-33% of model discoveries there come only from short ORFs,
versus **0%** for CPAT/CPC2, whose databases are 0-2% short by construction.

**`D13_discovery_errorbars`** -- the database-design tradeoff in four panels: GAIN (density), COST
(canonical PSMs displaced), DRIVER (database size), TOTAL (net unique peptides). Both models drawn.
**3 datasets, and the `null_nc` (near-cognate null) ARM IS GONE** -- it existed only for A549,
because nonspecific digestion of a ~3M-sequence near-cognate null exceeds MSFragger's ~2e9 peptide
cap. No near-cognate null remains anywhere in the human proteogenomics.

**`S_fdr_rigor`** -- why class-specific FDR is required.
*Must state:* **substantially weakened by the A549 removal.** The headline was A549/null_atg, 375
novel peptides at global 1% FDR versus 10 at class-specific -- a **37.5x** inflation. The largest
discrepancy among the remaining HLA-I sets is **2.3x** (HBL-1/CPAT, 7 vs 3 peptides). "Global FDR
inflates novel discovery ~10-13x" was an A549 result and is **no longer supported**. What survives is
that both cuts differ in the same direction on every arm, on much smaller numbers.

**`S_rescore_substrate` is RETIRED** (`figures/_retired_S_rescore_substrate/`). Its claim was a
CONTRAST between substrates -- rescoring does nothing on the tryptic proteome and a lot on HLA-I --
and with A549 gone there is only one substrate, so there is no contrast to draw. The pipeline
decision it justified (report raw hyperscore class-specific FDR) still stands and is documented in
`PROJECT_PIPELINE_POLICY.md`; it simply has no figure any more.

**Mouse macrophage panels (now attn, 30-aa tryptic floor -- see section 2):**

**`P11_bmdm_proteomics`** -- 12 populations; the model's database is a median **78x smaller** and
**43x denser** than a naive AUG null. **Cost-neutral** on total unique peptides (above the
GENCODE-only baseline in **6/12**, median **+3** of ~55,000) while the null is **0/12**, median
~-1,100.
*Must state:* the ratios were **156x / 111x** in the 7-aa build and that was inflated -- ~65% of the
null's 235,686 sequences were sub-30-aa ORFs the ORF-call track excludes by rule. 78x/43x is the
honest figure. The advantage is about half what the 7-aa comparison implied, still large, still the
right direction.

**`P12_macrophage_xsubtype`** -- predicted vs a real Ribo-seq database. On 11 transfer populations a
fixed **506-sequence** real database matches the model on novel yield (26 vs 24), beats it on total
unique peptides (**+18 vs -3**; 8/11 vs 5/11) and is **2.4x denser**.
*Must state:* the model's case is **not needing** Ribo-seq, not beating it. Correcting the floor did
NOT rescue the model here -- it narrowed the density gap (2.9x -> 2.4x) but widened the
total-peptide gap. This panel is the honest counterweight to P11 and should be shown with it.

**`P13_macrophage_psm_venn`** -- the two databases share only **25%** of their novel PSMs (866 model
/ 1,065 Ribo-seq / 386 shared of a 1,545 union, overlap in 11/11 populations but Jaccard never above
0.34). They are **complementary, not redundant**.
*Must state:* this replicated unusually well -- 24.0% at 7 aa with mamba, 25.0% at 30 aa with attn.
Two floors and two architectures agreeing to within one point makes complementarity a property of
predicted-vs-measured databases, not of one build.

### Interpretability and rigor (supporting)

**`C10_saliency`** -- the start codon is among the most salient nucleotides in the whole transcript.
**`C11_channel_ablation`** -- ORF-track channel ablation, 6 arms x 2 models. The frame channels
*suppress* spurious internal ORFs; removing any one roughly triples internal F1 in the wrong
direction. Judge this on **precision, not F1**.
**`C12_codon_occupancy`** -- the model reproduces codon-level occupancy it was never given.
**`C13_memorisation`** -- how much is memorisation, at matched depth: the naive gap is ~4x inflated
by depth; matched, it is ~0.10.
**`B4_replicate_ceiling`**, **`B5_expression_independence`**, **`B6_input_ablation`** (cell-type
specificity enters through the COUNT head, not the shape), **`B7_multimap_posture`**,
**`S_reproducibility`** (can the published pipeline regenerate the training data? 3.04-11.66% of
signal moves, median 3.86%), **`kozak`**, **`depth_crossover`**, **`D11_vs_seq2ribo`**,
**`o2_venn`**, **`benchmarks`**.

---

## 5. Claims that were measured and found FALSE

Do not put these on a poster. Each was tested and failed.

- **"Two architectures agreeing is evidence an ORF is real."** 8.7% corroboration, indistinguishable
  from a single model.
- **"The predicted database beats Ribo-seq."** It does not, on 11 macrophage transfer populations.
  The claim is that it does not *require* Ribo-seq.
- **"The model finds more peptides."** On total unique peptides it is cost-neutral, not positive
  (6 of 12 macrophage populations; median +3 of ~55,000 -- indistinguishable from zero, and well
  inside the +-50 search-parameter swing).
- **"P11 shows the model is the only arm that increases total peptides."** That was BMDM alone; BMDM
  is the best of 12.
- **"The model's database is 156x smaller / 111x denser than the null."** Inflated at a 7-aa floor by
  sub-30-aa null sequences the ORF-call track excludes by rule. The honest figures are **78x / 43x**.
- **Any non-canonical F1 predating the final-recipe alignment redo** -- inflated by ~22%.
- **Any ORF-call precision quoted without its transcript space and n_ref.**

---

## 6. Known-provisional, as of this build

- **P5 is 4/4 ON-RECIPE as of 2026-08-16.** An earlier draft of this bullet said THP-1 was "the one
  remaining arm on its original pack", which was INVERTED -- THP-1 was already on-recipe and the
  off-recipe pair was hepatocyte + iPSC-CM. Both have since been fixed: hepatocyte re-dumped on
  packed_canon_Hepatocytes (#86), iPSC-CM rebuilt from re-downloaded ENA FASTQs (#92). No arm
  carries a substrate caveat now. Superseded text follows for the record: its final-recipe pack
  was completed 2026-08-15 and the drop-in dump is running for both architectures; until it lands,
  P5 prints `ORIGINAL PACK` for THP-1 on every run, so the mixed state is visible rather than
  implied. CAR-T is now canonical (F1 0.904 -> 0.886, n_ref 8,856 -> 8,253).
- **A549 is REMOVED, not provisional.** It was archived on 2026-08-15 and taken out of every figure
  rather than rebuilt at 30 aa. Consequence: the human proteogenomics column is **HLA-I only**, and
  the project has **no tryptic human dataset**. The remaining sets (HBL-1, SU-DHL-4, DoHH2, THP-1,
  B721.221) are `nonspecific` / `termini = 0`, where 7 aa is correct and permanent because HLA-I
  peptides are 8-11 aa. The 12 mouse macrophage populations remain the only tryptic proteogenomics
  in the project, and they are at 30 aa.
- **GSE39561** is a deliberate negative control: zero periodic read lengths under the canonical
  recipe, hence a coverage-only pack with no `real` arm. Brain, despite similar ~30% unique mapping,
  behaves differently -- it selects offsets for 5/5 libraries. Low mapping does NOT predict
  periodicity failure, and the two are not interchangeable as "bad data" examples.

`PROJECT_STATUS_CURRENT_VS_ARCHIVED.md` is the authoritative list of what is live versus superseded.
`DATASET_LABELS.md` maps every internal key to its display label and accession.
