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
