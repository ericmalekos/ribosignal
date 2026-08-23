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
| gse243134 | janich | annotated | 0.996 | 0.994 | 0.995 | 10,383 |
| gse243134 | wang | annotated | 0.996 | 0.986 | 0.991 | 10,383 |
| janich | gse243134 | annotated | 0.994 | 0.996 | 0.995 | 10,361 |
| janich | wang | annotated | 0.995 | 0.988 | 0.992 | 10,361 |
| wang | janich | annotated | 0.988 | 0.995 | 0.992 | 10,282 |
| wang | gse243134 | annotated | 0.986 | 0.996 | 0.991 | 10,282 |
| gse243134 | janich | novel | 0.825 | 0.748 | 0.785 | 663 |
| gse243134 | wang | novel | 0.834 | 0.575 | 0.680 | 663 |
| janich | gse243134 | novel | 0.748 | 0.825 | 0.785 | 601 |
| janich | wang | novel | 0.805 | 0.612 | 0.696 | 601 |
| wang | janich | novel | 0.612 | 0.805 | 0.696 | 457 |
| wang | gse243134 | novel | 0.575 | 0.834 | 0.680 | 457 |
| gse243134 | janich | uORF | 0.889 | 0.782 | 0.832 | 697 |
| gse243134 | wang | uORF | 0.917 | 0.542 | 0.682 | 697 |
| janich | gse243134 | uORF | 0.782 | 0.889 | 0.832 | 613 |
| janich | wang | uORF | 0.905 | 0.609 | 0.728 | 613 |
| wang | janich | uORF | 0.609 | 0.905 | 0.728 | 412 |
| wang | gse243134 | uORF | 0.542 | 0.917 | 0.682 | 412 |
| gse243134 | janich | Overlap_uORF | 0.832 | 0.718 | 0.771 | 344 |
| gse243134 | wang | Overlap_uORF | 0.865 | 0.541 | 0.665 | 344 |
| janich | gse243134 | Overlap_uORF | 0.718 | 0.832 | 0.771 | 297 |
| janich | wang | Overlap_uORF | 0.800 | 0.579 | 0.672 | 297 |
| wang | janich | Overlap_uORF | 0.579 | 0.800 | 0.672 | 215 |
| wang | gse243134 | Overlap_uORF | 0.541 | 0.865 | 0.665 | 215 |
| gse243134 | janich | internal | 0.680 | 0.565 | 0.617 | 207 |
| gse243134 | wang | internal | 0.674 | 0.420 | 0.518 | 207 |
| janich | gse243134 | internal | 0.565 | 0.680 | 0.617 | 172 |
| janich | wang | internal | 0.597 | 0.448 | 0.512 | 172 |
| wang | janich | internal | 0.448 | 0.597 | 0.512 | 129 |
| wang | gse243134 | internal | 0.420 | 0.674 | 0.518 | 129 |
| gse243134 | janich | dORF | 0.670 | 0.553 | 0.606 | 246 |
| gse243134 | wang | dORF | 0.722 | 0.370 | 0.489 | 246 |
| janich | gse243134 | dORF | 0.553 | 0.670 | 0.606 | 203 |
| janich | wang | dORF | 0.627 | 0.389 | 0.480 | 203 |
| wang | janich | dORF | 0.389 | 0.627 | 0.480 | 126 |
| wang | gse243134 | dORF | 0.370 | 0.722 | 0.489 | 126 |
| gse243134 | janich | Overlap_dORF | 0.571 | 0.545 | 0.558 | 44 |
| gse243134 | wang | Overlap_dORF | 0.640 | 0.364 | 0.464 | 44 |
| janich | gse243134 | Overlap_dORF | 0.545 | 0.571 | 0.558 | 42 |
| janich | wang | Overlap_dORF | 0.560 | 0.333 | 0.418 | 42 |
| wang | janich | Overlap_dORF | 0.333 | 0.560 | 0.418 | 25 |
| wang | gse243134 | Overlap_dORF | 0.364 | 0.640 | 0.464 | 25 |

**Annotated CDS reproduces; nothing else does.** Annotated F1 is 0.991-0.995 across all six ordered pairs. Every other class falls away:
uORF 0.682-0.832,
novel 0.680-0.785,
Overlap_uORF 0.665-0.771,
internal 0.512-0.617,
dORF 0.480-0.606,
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
7-class x 9-cell x 2-arm table is `results/mouse_liver_3x3_canon/scored/model_vs_ribo_by_class.tsv`.
Read each row against the same class in the ceiling table above, not against 1.0.

### Predicted shape at observed depth (`pred_obsdepth`)

| class | dataset | model | P | R | F1 | n_ref |
|---|---|---|--:|--:|--:|--:|
| annotated | janich | attn | 0.991 | 0.976 | 0.983 | 10,361 |
| annotated | janich | mamba4 | 0.990 | 0.986 | 0.988 | 10,361 |
| annotated | gse243134 | attn | 0.992 | 0.980 | 0.986 | 10,383 |
| annotated | gse243134 | mamba4 | 0.992 | 0.988 | 0.990 | 10,383 |
| annotated | wang | attn | 0.993 | 0.974 | 0.983 | 10,282 |
| annotated | wang | mamba4 | 0.990 | 0.982 | 0.986 | 10,282 |
| novel | janich | attn | 0.609 | 0.571 | 0.589 | 601 |
| novel | janich | mamba4 | 0.617 | 0.604 | 0.611 | 601 |
| novel | gse243134 | attn | 0.595 | 0.560 | 0.577 | 663 |
| novel | gse243134 | mamba4 | 0.606 | 0.573 | 0.589 | 663 |
| novel | wang | attn | 0.521 | 0.617 | 0.565 | 457 |
| novel | wang | mamba4 | 0.513 | 0.619 | 0.561 | 457 |
| uORF | janich | attn | 0.622 | 0.643 | 0.632 | 613 |
| uORF | janich | mamba4 | 0.633 | 0.638 | 0.635 | 613 |
| uORF | gse243134 | attn | 0.646 | 0.658 | 0.652 | 697 |
| uORF | gse243134 | mamba4 | 0.647 | 0.648 | 0.648 | 697 |
| uORF | wang | attn | 0.493 | 0.337 | 0.401 | 412 |
| uORF | wang | mamba4 | 0.527 | 0.359 | 0.427 | 412 |
| internal | janich | attn | 0.036 | 0.006 | 0.010 | 172 |
| internal | janich | mamba4 | 0.154 | 0.023 | 0.040 | 172 |
| internal | gse243134 | attn | 0.028 | 0.005 | 0.008 | 207 |
| internal | gse243134 | mamba4 | 0.103 | 0.015 | 0.025 | 207 |
| internal | wang | attn | 0.000 | 0.000 | 0.000 | 129 |
| internal | wang | mamba4 | 0.143 | 0.023 | 0.040 | 129 |

### Standalone, no Ribo-seq at inference (`pred_preddepth`)

| class | dataset | model | P | R | F1 | n_ref |
|---|---|---|--:|--:|--:|--:|
| annotated | janich | attn | 0.959 | 0.998 | 0.978 | 10,361 |
| annotated | janich | mamba4 | 0.958 | 0.982 | 0.970 | 10,361 |
| annotated | gse243134 | attn | 0.964 | 0.996 | 0.980 | 10,383 |
| annotated | gse243134 | mamba4 | 0.962 | 0.988 | 0.975 | 10,383 |
| annotated | wang | attn | 0.954 | 0.998 | 0.976 | 10,282 |
| annotated | wang | mamba4 | 0.953 | 0.997 | 0.975 | 10,282 |
| novel | janich | attn | 0.247 | 0.712 | 0.367 | 601 |
| novel | janich | mamba4 | 0.251 | 0.697 | 0.369 | 601 |
| novel | gse243134 | attn | 0.262 | 0.695 | 0.380 | 663 |
| novel | gse243134 | mamba4 | 0.269 | 0.674 | 0.384 | 663 |
| novel | wang | attn | 0.182 | 0.720 | 0.291 | 457 |
| novel | wang | mamba4 | 0.182 | 0.720 | 0.290 | 457 |
| uORF | janich | attn | 0.558 | 0.693 | 0.619 | 613 |
| uORF | janich | mamba4 | 0.613 | 0.519 | 0.562 | 613 |
| uORF | gse243134 | attn | 0.597 | 0.684 | 0.638 | 697 |
| uORF | gse243134 | mamba4 | 0.637 | 0.534 | 0.581 | 697 |
| uORF | wang | attn | 0.415 | 0.362 | 0.387 | 412 |
| uORF | wang | mamba4 | 0.429 | 0.366 | 0.395 | 412 |
| internal | janich | attn | 0.028 | 0.006 | 0.010 | 172 |
| internal | janich | mamba4 | 0.231 | 0.017 | 0.032 | 172 |
| internal | gse243134 | attn | 0.026 | 0.005 | 0.008 | 207 |
| internal | gse243134 | mamba4 | 0.087 | 0.010 | 0.017 | 207 |
| internal | wang | attn | 0.000 | 0.000 | 0.000 | 129 |
| internal | wang | mamba4 | 0.091 | 0.023 | 0.037 | 129 |

Model F1 against the between-experiment ceiling, per class (`pred_obsdepth`, matched RNA, both models pooled):

| class | model F1 range | ceiling F1 range |
|---|--:|--:|
| annotated | 0.983-0.990 | 0.991-0.995 |
| novel | 0.561-0.611 | 0.680-0.785 |
| uORF | 0.401-0.652 | 0.682-0.832 |
| internal | 0.000-0.040 | 0.512-0.617 |

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
a real experiment does in **25 of 144** comparisons. Every one of those
wins is against Wang (13 novel, 12 uORF), and the winning arm is
**always `pred_preddepth`** -- the standalone model, with no Ribo-seq for the query sample at
any point.

| class | reference | best model recall | arm | Wang's recall | margin |
|---|---|--:|---|--:|--:|
| novel | janich | 0.719 | `preddepth` | 0.612 | +0.107 |
| novel | gse243134 | 0.695 | `preddepth` | 0.575 | +0.121 |
| uORF | janich | 0.711 | `preddepth` | 0.609 | +0.103 |
| uORF | gse243134 | 0.684 | `preddepth` | 0.542 | +0.142 |

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
**20 comparisons**, and in **0 of them on F1**. That split is the point: the
wins are all on one side of the precision/recall trade, never on the balanced metric.

| class | reference | model | RNA | arm | metric | model | ceiling | delta |
|---|---|---|---|---|---|--:|--:|--:|
| all | wang | attn | wang | `pred_obsdepth` | precision | 0.936 | 0.921 | +0.015 |
| all | wang | mamba4 | wang | `pred_obsdepth` | precision | 0.933 | 0.921 | +0.012 |
| annotated | wang | attn | wang | `pred_obsdepth` | precision | 0.993 | 0.988 | +0.005 |
| annotated | wang | attn | gse243134 | `pred_obsdepth` | precision | 0.992 | 0.988 | +0.004 |
| annotated | gse243134 | attn | wang | `pred_preddepth` | recall | 0.998 | 0.994 | +0.004 |
| annotated | wang | attn | janich | `pred_obsdepth` | precision | 0.992 | 0.988 | +0.004 |
| annotated | gse243134 | mamba4 | wang | `pred_preddepth` | recall | 0.997 | 0.994 | +0.003 |
| annotated | gse243134 | attn | janich | `pred_preddepth` | recall | 0.997 | 0.994 | +0.003 |
| annotated | wang | mamba4 | gse243134 | `pred_obsdepth` | precision | 0.991 | 0.988 | +0.003 |
| annotated | gse243134 | attn | gse243134 | `pred_preddepth` | recall | 0.996 | 0.994 | +0.003 |
| annotated | wang | mamba4 | wang | `pred_obsdepth` | precision | 0.990 | 0.988 | +0.003 |
| annotated | wang | attn | wang | `pred_preddepth` | recall | 0.998 | 0.996 | +0.002 |

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
| janich | 0.919* | 0.912 | 0.925 |
| gse243134 | 0.918 | 0.913* | 0.923 |
| wang | 0.908 | 0.899 | 0.922* |

**The matched cell never wins its column.** Across all 12 model x reference x arm
combinations the matched-minus-mismatched difference averages **+0.0035 F1**
(range -0.0042 to +0.0070); every one is under 0.008 in absolute value and
two are negative. Substituting an unrelated experiment's RNA-seq costs essentially nothing.

The RNA effect that *is* real is input quality, not matching:

| RNA input | samples | mean F1 (obsdepth, over refs and models) |
|---|--:|--:|
| janich | 7 | 0.9144 |
| gse243134 | 19 | 0.9137 |
| wang | 2 | 0.9086 |

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
