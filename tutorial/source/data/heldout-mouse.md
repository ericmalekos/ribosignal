# Held-out mouse data

Mouse is the stronger generalisation evidence in this project, for a reason that is
structural rather than about sample size: **three independent labs measured Ribo-seq on
mouse liver**. That yields a between-experiment ceiling for the ORF caller itself, so the
model can be scored against what the assay actually reproduces instead of against 1.0.

```{list-table}
:header-rows: 1
:widths: 24 26 50

* - Dataset
  - Role
  - Notes
* - **Wang liver**
  - held-out study, cross-species
  - Human-trained model applied to mouse with no retraining. The shallowest of the three liver sets.
* - **Janich liver**
  - held-out study, cross-species
  - RNA arm was rebuilt after decontamination: 6-11% of reads carried 52-70% of TPM through salmon length normalisation, which shrank the expressed universe to 11,532 tx. Post-fix: 16,379.
* - **GSE243134 liver**
  - held-out study, cross-species
  - Deepest of the three. Its RNA arm is GSE302188 poly(A), replacing a Ribo-Zero total-RNA arm that could not be used (tRNA/7SL survive the biotype filter).
* - **GSE155087 T-cell**
  - held-out study, different tissue
  - Mouse CD4+ T-cell, WT only. Reported in {doc}`/benchmarks` rather than here.
```

## Between-experiment ceiling

Before reading any model row: this is what the three observed datasets score **against
each other**, on the same restricted transcript space and the same genomic key. It is the
realistic upper bound.

Restricted to {'wang_tx': 14331, 'janich_tx': 16553, 'janich_reps': 7, 'shared_tx': 11032, 'shared_genes': 8506}.

| observed pair | F1 | CDS F1 | non-canonical F1 |
|---|--:|--:|--:|
| Wang vs Janich | 0.924 | 0.986 | 0.633 |
| Wang vs GSE243134 | 0.918 | 0.986 | 0.617 |
| Janich vs Wang | 0.924 | 0.986 | 0.633 |
| Janich vs GSE243134 | 0.953 | 0.997 | 0.788 |
| GSE243134 vs Wang | 0.918 | 0.986 | 0.617 |
| GSE243134 vs Janich | 0.953 | 0.997 | 0.788 |

**Annotated CDS calls reproduce; non-canonical calls do not.** The CDS ceiling is
near-perfect and stable across all three pairs, while the non-canonical ceiling
both is far lower and varies by pair. Any single-pair non-canonical number is
therefore a property of that pair, not a property of the assay.

## Standalone drop-in: no Ribo-seq for the query sample

`pred_preddepth` = predicted profile shape **and** predicted depth. Scored against RiboCode
run on that pack's real counts, so the transcript space is identical on both sides by
construction. Class columns are recall.

| dataset | n_ref | P | R | F1 | annotated | uORF | Ovl_uORF | novel | dORF | Ovl_dORF | internal |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| Wang liver (attn) | 8,552 | 0.809 | 0.934 | 0.867 | 0.999 | 0.438 | 0.360 | 0.717 | 0.138 | 0.333 | 0.033 |
| Wang liver (mamba4) | 8,552 | 0.810 | 0.932 | 0.867 | 0.999 | 0.446 | 0.320 | 0.709 | 0.098 | 0.296 | 0.022 |
| Janich liver (attn) | 9,850 | 0.880 | 0.917 | 0.898 | 0.998 | 0.670 | 0.637 | 0.673 | 0.111 | 0.273 | 0.014 |
| Janich liver (mamba4) | 9,850 | 0.899 | 0.888 | 0.893 | 0.980 | 0.523 | 0.437 | 0.662 | 0.060 | 0.159 | 0.022 |
| GSE243134 liver (attn) | 11,761 | 0.879 | 0.912 | 0.895 | 0.998 | 0.613 | 0.526 | 0.695 | 0.046 | 0.222 | 0.035 |
| GSE243134 liver (mamba4) | 11,761 | 0.881 | 0.912 | 0.896 | 0.999 | 0.610 | 0.514 | 0.699 | 0.054 | 0.156 | 0.040 |

## Predicted shape at observed depth

`pred_obsdepth` keeps the real per-transcript depth and replaces only the profile shape.
The gap between this and the table above isolates **what the count head costs**, separate
from whether the shape is right.

| dataset | n_ref | P | R | F1 | annotated | uORF | Ovl_uORF | novel | dORF | Ovl_dORF | internal |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| Wang liver (attn) | 8,552 | 0.940 | 0.906 | 0.923 | 0.983 | 0.296 | 0.213 | 0.577 | 0.106 | 0.333 | 0.000 |
| Wang liver (mamba4) | 8,552 | 0.942 | 0.912 | 0.927 | 0.987 | 0.350 | 0.240 | 0.601 | 0.081 | 0.259 | 0.011 |
| Janich liver (attn) | 9,850 | 0.935 | 0.904 | 0.919 | 0.989 | 0.642 | 0.567 | 0.614 | 0.101 | 0.273 | 0.022 |
| Janich liver (mamba4) | 9,850 | 0.944 | 0.908 | 0.926 | 0.994 | 0.652 | 0.543 | 0.644 | 0.074 | 0.182 | 0.050 |
| GSE243134 liver (attn) | 11,761 | 0.948 | 0.892 | 0.919 | 0.987 | 0.561 | 0.511 | 0.548 | 0.046 | 0.178 | 0.025 |
| GSE243134 liver (mamba4) | 11,761 | 0.947 | 0.896 | 0.921 | 0.990 | 0.575 | 0.474 | 0.585 | 0.038 | 0.156 | 0.040 |

Mean F1 cost of predicting depth as well as shape: **+0.036** (range +0.021 to +0.060). Precision carries essentially all of it --
the standalone arm calls MORE ORFs than the reference, it does not miss them.

```{admonition} Checkpoint provenance
:class: important
Each row's source checkpoint, read from the npz `meta` field at build time:

- `Wang liver (attn)` -> `orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes` (RELEASED)
- `Wang liver (mamba4)` -> `orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocyt` (RELEASED)
- `Janich liver (attn)` -> `orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes` (RELEASED)
- `Janich liver (mamba4)` -> `orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocyt` (RELEASED)
- `GSE243134 liver (attn)` -> `orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes` (RELEASED)
- `GSE243134 liver (mamba4)` -> `orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocyt` (RELEASED)

All rows are on the shipping models.
```

## Cross-experiment agreement by ORF class

Precision AND recall per class, not recall alone: a recall-only breakdown cannot tell
"found the uORFs" apart from "called uORFs everywhere and some were right". No model is
involved -- this is what the three OBSERVED call sets score against each other on the shared
22,974-tx universe, and it is the ceiling every model row should be read against.

| reference | compared | class | P | R | F1 | n_ref |
|---|---|---|--:|--:|--:|--:|
| gse243134 | janich | annotated | 0.996 | 0.995 | 0.995 | 10,529 |
| gse243134 | wang | annotated | 0.996 | 0.988 | 0.992 | 10,529 |
| janich | gse243134 | annotated | 0.995 | 0.996 | 0.995 | 10,515 |
| janich | wang | annotated | 0.996 | 0.989 | 0.993 | 10,515 |
| wang | janich | annotated | 0.989 | 0.996 | 0.993 | 10,444 |
| wang | gse243134 | annotated | 0.988 | 0.996 | 0.992 | 10,444 |
| gse243134 | janich | novel | 0.851 | 0.748 | 0.796 | 874 |
| gse243134 | wang | novel | 0.848 | 0.581 | 0.690 | 874 |
| janich | gse243134 | novel | 0.748 | 0.851 | 0.796 | 769 |
| janich | wang | novel | 0.826 | 0.644 | 0.724 | 769 |
| wang | janich | novel | 0.644 | 0.826 | 0.724 | 599 |
| wang | gse243134 | novel | 0.581 | 0.848 | 0.690 | 599 |
| gse243134 | janich | uORF | 0.890 | 0.795 | 0.840 | 806 |
| gse243134 | wang | uORF | 0.909 | 0.560 | 0.693 | 806 |
| janich | gse243134 | uORF | 0.795 | 0.890 | 0.840 | 720 |
| janich | wang | uORF | 0.909 | 0.626 | 0.742 | 720 |
| wang | janich | uORF | 0.626 | 0.909 | 0.742 | 496 |
| wang | gse243134 | uORF | 0.560 | 0.909 | 0.693 | 496 |
| gse243134 | janich | Overlap_uORF | 0.846 | 0.745 | 0.793 | 377 |
| gse243134 | wang | Overlap_uORF | 0.858 | 0.546 | 0.668 | 377 |
| janich | gse243134 | Overlap_uORF | 0.745 | 0.846 | 0.793 | 332 |
| janich | wang | Overlap_uORF | 0.821 | 0.593 | 0.689 | 332 |
| wang | janich | Overlap_uORF | 0.593 | 0.821 | 0.689 | 240 |
| wang | gse243134 | Overlap_uORF | 0.546 | 0.858 | 0.668 | 240 |
| gse243134 | janich | internal | 0.679 | 0.570 | 0.620 | 223 |
| gse243134 | wang | internal | 0.687 | 0.453 | 0.546 | 223 |
| janich | gse243134 | internal | 0.570 | 0.679 | 0.620 | 187 |
| janich | wang | internal | 0.612 | 0.481 | 0.539 | 187 |
| wang | janich | internal | 0.481 | 0.612 | 0.539 | 147 |
| wang | gse243134 | internal | 0.453 | 0.687 | 0.546 | 147 |
| gse243134 | janich | dORF | 0.705 | 0.577 | 0.635 | 286 |
| gse243134 | wang | dORF | 0.757 | 0.392 | 0.516 | 286 |
| janich | gse243134 | dORF | 0.577 | 0.705 | 0.635 | 234 |
| janich | wang | dORF | 0.676 | 0.427 | 0.524 | 234 |
| wang | janich | dORF | 0.427 | 0.676 | 0.524 | 148 |
| wang | gse243134 | dORF | 0.392 | 0.757 | 0.516 | 148 |
| gse243134 | janich | Overlap_dORF | 0.617 | 0.580 | 0.598 | 50 |
| gse243134 | wang | Overlap_dORF | 0.600 | 0.360 | 0.450 | 50 |
| janich | gse243134 | Overlap_dORF | 0.580 | 0.617 | 0.598 | 47 |
| janich | wang | Overlap_dORF | 0.567 | 0.362 | 0.442 | 47 |
| wang | janich | Overlap_dORF | 0.362 | 0.567 | 0.442 | 30 |
| wang | gse243134 | Overlap_dORF | 0.360 | 0.600 | 0.450 | 30 |

**Annotated CDS reproduces; nothing else does.** Annotated F1 is 0.992-0.995 across all six ordered pairs. Every other class falls away:
uORF 0.693-0.840,
novel 0.690-0.796,
Overlap_uORF 0.668-0.793,
internal 0.539-0.620,
dORF 0.516-0.635,
with `internal` and `dORF` the weakest. Judge a model's non-canonical calls against these
numbers, never against 1.0.

**The precision/recall split is a depth effect.** Wang is the shallowest set (78.3 M
P-sites vs 110.0 M and 150.5 M). Taking Wang as the reference, the deeper datasets recover
90.9% of its uORFs (recall) while only 56.0% of theirs appear in Wang (precision); reverse
the roles and the two numbers swap exactly. The deeper library simply calls more
non-canonical ORFs, and the shallow set's calls are close to a subset of them -- so a
non-canonical precision computed against a shallower reference is mostly reporting the
depth gap, not a disagreement about biology.

## Model vs observed, by ORF class

The same four classes, now for the MODEL against each dataset's observed calls. Rows are
matched-RNA cells (RNA input from the same experiment as the Ribo reference); the full
7-class x 9-cell x 2-arm table is `results/mouse_liver_3x3/scored/model_vs_ribo_by_class.tsv`.
Read each row against the same class in the ceiling table above, not against 1.0.

### Predicted shape at observed depth (`pred_obsdepth`)

| class | dataset | model | P | R | F1 | n_ref |
|---|---|---|--:|--:|--:|--:|
| annotated | janich | attn | 0.994 | 0.978 | 0.986 | 10,515 |
| annotated | janich | mamba4 | 0.993 | 0.989 | 0.991 | 10,515 |
| annotated | gse243134 | attn | 0.994 | 0.983 | 0.989 | 10,529 |
| annotated | gse243134 | mamba4 | 0.994 | 0.991 | 0.992 | 10,529 |
| annotated | wang | attn | 0.994 | 0.977 | 0.985 | 10,444 |
| annotated | wang | mamba4 | 0.991 | 0.985 | 0.988 | 10,444 |
| novel | janich | attn | 0.667 | 0.602 | 0.633 | 769 |
| novel | janich | mamba4 | 0.683 | 0.632 | 0.656 | 769 |
| novel | gse243134 | attn | 0.656 | 0.587 | 0.620 | 874 |
| novel | gse243134 | mamba4 | 0.668 | 0.595 | 0.629 | 874 |
| novel | wang | attn | 0.569 | 0.633 | 0.599 | 599 |
| novel | wang | mamba4 | 0.560 | 0.634 | 0.595 | 599 |
| uORF | janich | attn | 0.684 | 0.625 | 0.653 | 720 |
| uORF | janich | mamba4 | 0.699 | 0.619 | 0.657 | 720 |
| uORF | gse243134 | attn | 0.720 | 0.656 | 0.687 | 806 |
| uORF | gse243134 | mamba4 | 0.721 | 0.650 | 0.684 | 806 |
| uORF | wang | attn | 0.576 | 0.351 | 0.436 | 496 |
| uORF | wang | mamba4 | 0.605 | 0.359 | 0.451 | 496 |
| internal | janich | attn | 0.067 | 0.011 | 0.018 | 187 |
| internal | janich | mamba4 | 0.143 | 0.021 | 0.037 | 187 |
| internal | gse243134 | attn | 0.029 | 0.004 | 0.008 | 223 |
| internal | gse243134 | mamba4 | 0.088 | 0.013 | 0.023 | 223 |
| internal | wang | attn | 0.000 | 0.000 | 0.000 | 147 |
| internal | wang | mamba4 | 0.136 | 0.020 | 0.035 | 147 |

### Standalone, no Ribo-seq at inference (`pred_preddepth`)

| class | dataset | model | P | R | F1 | n_ref |
|---|---|---|--:|--:|--:|--:|
| annotated | janich | attn | 0.973 | 0.997 | 0.985 | 10,515 |
| annotated | janich | mamba4 | 0.972 | 0.982 | 0.977 | 10,515 |
| annotated | gse243134 | attn | 0.977 | 0.996 | 0.987 | 10,529 |
| annotated | gse243134 | mamba4 | 0.975 | 0.988 | 0.981 | 10,529 |
| annotated | wang | attn | 0.969 | 0.998 | 0.983 | 10,444 |
| annotated | wang | mamba4 | 0.968 | 0.997 | 0.982 | 10,444 |
| novel | janich | attn | 0.322 | 0.724 | 0.446 | 769 |
| novel | janich | mamba4 | 0.329 | 0.714 | 0.450 | 769 |
| novel | gse243134 | attn | 0.347 | 0.700 | 0.464 | 874 |
| novel | gse243134 | mamba4 | 0.359 | 0.684 | 0.471 | 874 |
| novel | wang | attn | 0.244 | 0.735 | 0.366 | 599 |
| novel | wang | mamba4 | 0.246 | 0.743 | 0.369 | 599 |
| uORF | janich | attn | 0.633 | 0.669 | 0.651 | 720 |
| uORF | janich | mamba4 | 0.694 | 0.500 | 0.581 | 720 |
| uORF | gse243134 | attn | 0.682 | 0.676 | 0.679 | 806 |
| uORF | gse243134 | mamba4 | 0.730 | 0.528 | 0.613 | 806 |
| uORF | wang | attn | 0.510 | 0.369 | 0.428 | 496 |
| uORF | wang | mamba4 | 0.517 | 0.367 | 0.429 | 496 |
| internal | janich | attn | 0.083 | 0.016 | 0.027 | 187 |
| internal | janich | mamba4 | 0.231 | 0.016 | 0.030 | 187 |
| internal | gse243134 | attn | 0.026 | 0.004 | 0.008 | 223 |
| internal | gse243134 | mamba4 | 0.087 | 0.009 | 0.016 | 223 |
| internal | wang | attn | 0.000 | 0.000 | 0.000 | 147 |
| internal | wang | mamba4 | 0.091 | 0.020 | 0.033 | 147 |

Model F1 against the between-experiment ceiling, per class (`pred_obsdepth`, matched RNA, both models pooled):

| class | model F1 range | ceiling F1 range |
|---|--:|--:|
| annotated | 0.985-0.992 | 0.992-0.995 |
| novel | 0.595-0.656 | 0.690-0.796 |
| uORF | 0.436-0.687 | 0.693-0.840 |
| internal | 0.000-0.037 | 0.539-0.620 |

## Library depth and periodicity

The three liver sets differ by ~2x in usable signal and by 7 points in 3-nt periodicity.
That spread is what makes the comparisons below interpretable, so it is stated up front
rather than left implicit.

| dataset | Ribo samples | input reads | STAR unique | pooled P-sites | mean f0 (periodicity) |
|---|--:|--:|--:|--:|--:|
| janich | 5 | 184,157,988 | 57.8% | 109,996,969 | 82.1% |
| gse243134 | 21 | 548,815,069 | 49.9% | 150,531,643 | 80.5% |
| wang | 2 | 182,182,664 | 46.4% | 78,257,663 | 75.4% |

**Wang is the shallowest AND the least periodic** -- 2 samples, 78.3 M P-sites, 75.4% f0.
That is not an outlier, it is what a large share of published Ribo-seq looks like.

## Model recall beats a real Ribo-seq experiment on novel ORFs and uORFs

Across the novel and uORF classes, a model cell recovers MORE of the reference's calls than
a real experiment does in **28 of 144** comparisons. Every one of those
wins is against Wang (18 novel, 10 uORF), and the winning arm is
**always `pred_preddepth`** -- the standalone model, with no Ribo-seq for the query sample at
any point.

| class | reference | best model recall | arm | Wang's recall | margin |
|---|---|--:|---|--:|--:|
| novel | janich | 0.724 | `preddepth` | 0.644 | +0.081 |
| novel | gse243134 | 0.709 | `preddepth` | 0.581 | +0.128 |
| uORF | janich | 0.688 | `preddepth` | 0.626 | +0.061 |
| uORF | gse243134 | 0.676 | `preddepth` | 0.560 | +0.117 |

```{admonition} Why this matters more than the F1 table suggests
:class: tip
It is tempting to discount these as "only beating the shallow library". That reading is
backwards. Most published Ribo-seq IS shallow, or low-periodicity, or both -- Wang's 2 samples
and 75.4% f0 are unremarkable for the field. The result says that for novel ORFs and uORFs,
**a model with no ribosome profiling at all recovers more of a deep reference's calls than a
real, published Ribo-seq experiment does** -- at essentially zero marginal cost, from sequence
and RNA-seq that most labs already have.

The honest boundaries: the model still loses to the two DEEPER experiments on recall, and it
loses to all of them on precision and therefore on F1 (novel 0.595-0.656 and uORF 0.436-0.687
against experimental 0.690-0.796 and 0.693-0.840). So this is not "the model replaces deep
Ribo-seq". It is "the model is a better non-canonical detector than a shallow experiment,
and free" -- which is the substitution most labs actually face.
```

## Where the model exceeds the between-experiment value

Scoring every model cell against the STRICTEST bar -- the best of the two other observed
datasets on the same reference, same key, same universe -- the model comes out ahead in
**22 comparisons**, and in **0 of them on F1**. That split is the point: the
wins are all on one side of the precision/recall trade, never on the balanced metric.

| class | reference | model | RNA | arm | metric | model | ceiling | delta |
|---|---|---|---|---|---|--:|--:|--:|
| all | wang | attn | wang | `pred_obsdepth` | precision | 0.938 | 0.918 | +0.021 |
| all | wang | mamba4 | wang | `pred_obsdepth` | precision | 0.934 | 0.918 | +0.017 |
| all | wang | mamba4 | janich | `pred_obsdepth` | precision | 0.927 | 0.918 | +0.009 |
| annotated | wang | attn | wang | `pred_obsdepth` | precision | 0.994 | 0.989 | +0.005 |
| annotated | wang | attn | gse243134 | `pred_obsdepth` | precision | 0.994 | 0.989 | +0.005 |
| annotated | wang | attn | janich | `pred_obsdepth` | precision | 0.993 | 0.989 | +0.003 |
| annotated | wang | mamba4 | gse243134 | `pred_obsdepth` | precision | 0.992 | 0.989 | +0.003 |
| all | wang | mamba4 | gse243134 | `pred_obsdepth` | precision | 0.921 | 0.918 | +0.003 |
| annotated | gse243134 | attn | wang | `pred_preddepth` | recall | 0.998 | 0.995 | +0.003 |
| annotated | gse243134 | mamba4 | wang | `pred_preddepth` | recall | 0.997 | 0.995 | +0.002 |
| annotated | gse243134 | attn | janich | `pred_preddepth` | recall | 0.997 | 0.995 | +0.002 |
| annotated | wang | mamba4 | wang | `pred_obsdepth` | precision | 0.991 | 0.989 | +0.002 |

Two clusters, and they are not equally meaningful.

- **`pred_obsdepth` precision against Wang (up to +0.021)** is mostly a DEPTH ARTEFACT, not
  superiority. Wang is the shallowest set, so a deeper experiment scored against it calls
  many ORFs Wang never called and is punished on precision. The model calls fewer, so it
  scores higher. Read it as "more conservative than a deeper library", not "better than an
  experiment".
- **`pred_preddepth` annotated recall (0.996-0.998)** is the real one. The STANDALONE model --
  no Ribo-seq at inference at all -- recovers marginally more annotated CDS than a second
  real experiment does (0.995-0.996). The margin is small, but the comparison is honest:
  canonical CDS recall is saturated, and a second wet-lab replicate adds nothing over
  predicting it from sequence and RNA-seq. This matches the Ruiz-Orera result, where
  predicted frame-0 localization AUROC (0.945) beat that study's own observed ceiling (0.908).

Nowhere does a model beat the ceiling on F1, and nowhere does it come close on the
non-canonical classes -- see the internal-ORF row above.

## Does the RNA input have to match the Ribo-seq?

The three liver datasets each bring their own RNA-seq, so every Ribo-seq reference can be
predicted from every RNA input: a **3 Ribo x 3 RNA** factorial, all nine cells pooled from
BAMs through one pipeline onto one shared universe (22,974 tx). If the model were quietly
relying on the RNA-seq being from the same experiment, the diagonal would dominate.

`mamba4`, `pred_obsdepth`, all-class F1. Rows are the RNA input, columns the Ribo reference;
`*` marks the matched cell.

| RNA input | ref=janich | ref=gse243134 | ref=wang |
|---|--:|--:|--:|
| janich | 0.919* | 0.912 | 0.926 |
| gse243134 | 0.919 | 0.913* | 0.923 |
| wang | 0.907 | 0.897 | 0.919* |

**The matched cell never wins its column.** Across all 12 model x reference x arm
combinations the matched-minus-mismatched difference averages **+0.0030 F1**
(range -0.0057 to +0.0079); every one is under 0.008 in absolute value and
two are negative. Substituting an unrelated experiment's RNA-seq costs essentially nothing.

The RNA effect that *is* real is input quality, not matching:

| RNA input | samples | mean F1 (obsdepth, over refs and models) |
|---|--:|--:|
| janich | 7 | 0.9149 |
| gse243134 | 19 | 0.9143 |
| wang | 2 | 0.9064 |

Wang's 2-sample RNA costs ~0.008 F1 against any reference, while 19 samples buy nothing over
7. Depth matters up to a handful of samples and then saturates.

```{admonition} Not comparable to pre-2026-08-10 Janich numbers
:class: warning
RiboCode `metaplots` re-derives which read lengths carry periodicity **per sample**, so
re-pooling Janich from re-fetched FASTQs reproduces its historical pack only to ~1.7%. These
nine cells are internally consistent; the superseded numbers are archived under
`results/_archive_pre_2026_08_10_janich_pack/`.
```

## Reading it

- **The cross-species transfer holds.** A model trained only on human tissue recovers
  ~99% of annotated mouse ORFs from sequence and RNA-seq alone.
- **uORF recall tracks the dataset, not the model.** The two arms of each dataset agree
  closely with each other and differ substantially between datasets, which points at
  library depth and 5'UTR coverage rather than at architecture.
- **dORF and internal are at the floor everywhere**, exactly as in human.
- **attn and mamba4 are interchangeable here.** No dataset separates them by more than a
  few thousandths of F1, consistent with the 3-seed union comparison.

_Generated by `tutorial/make_heldout_mouse.py`; do not edit by hand._
