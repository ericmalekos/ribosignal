# Does "held out" mean the model has not seen it?

Short answer: **no, and it matters more than expected.**

## The split holds out a tissue, not transcripts

`dataset.py::loto_split` holds out a TISSUE. Train and val are chromosome-split against each other
for honest early stopping, but `test_tx` is "the held-out tissue's scorable tx ... **all
chromosomes**". The nine Chothani tissues share one fixed 36,668-tx universe, and per the LOTO header
"per-tissue expression enters via the coverage input, not the transcript set".

Measured on the Hepatocytes holdout:

| | count | share |
|---|--:|--:|
| test transcripts | 33,918 | |
| also in TRAIN gradients (via some held-in tissue) | 27,007 | 79.6% |
| also in VAL (early stopping) | 6,813 | 20.1% |
| **never seen in any held-in tissue** | **98** | **0.3%** |

For a given transcript the SEQUENCE and the ORF TRACK are byte-identical between its train and test
appearances. Only coverage and target differ.

## A three-tier exposure test

The 98 never-seen transcripts are confounded: they are unseen precisely BECAUSE they fall below the
signal threshold in every held-in tissue, so "never seen" is entangled with "Hepatocyte-specific".

The `VAL only` tier breaks that entanglement. Those transcripts ARE expressed in held-in tissues --
they passed the same signal filter -- but were excluded from gradient updates by the chromosome fold.
All tiers are depth-matched to the unseen count window (50-241 observed CDS counts), because depth
alone swings profile correlation by more than the effect being measured.

| tier | in gradients | expressed elsewhere | n | mean r (mamba4) | 95% CI |
|---|---|---|--:|--:|---|
| in TRAIN gradients | yes | yes | 1,382 | 0.301 | [0.289, 0.314] |
| **in VAL only** | **no** | **yes** | 435 | **0.233** | [0.213, 0.256] |
| never seen anywhere | no | no | 93 | 0.140 | [0.106, 0.176] |

`attn` is identical in shape: 0.310 / 0.241 / 0.144. **All three intervals are disjoint in both
models.**

## What this establishes

**Gradient exposure is worth about +0.068 profile correlation**, controlling for expression. That is
the TRAIN vs VAL-only contrast, where both tiers are expressed in held-in tissues and differ only in
whether the optimizer took gradients on them. It is a memorisation effect, not a tissue-specificity
effect.

The further drop to never-seen (-0.094) mixes remaining exposure with tissue-specificity and is not
resolved by this experiment.

## Reconciling with the cross-species result

Cross-species mouse shares NO transcripts with training and loses only ~0.02-0.04 **ORF-call F1**
against human Hepatocytes. That looks contradictory next to a -0.068 hit in **profile correlation**,
and the resolution is that these are different readouts:

- Profile correlation is fine-grained and clearly memorisation-sensitive.
- ORF calling is thresholded -- it needs the right ORF, not the exact profile -- and is largely
  robust to it.

Working model: **memorisation inflates held-out PROFILE metrics substantially while leaving ORF-CALL
metrics largely intact.** Stated as a hypothesis, not a conclusion; it predicts that any future
profile-level claim on Hepatocytes needs this caveat and any ORF-level claim needs it much less.

```{admonition} Consequence for reading the rest of this site
:class: warning
Hepatocytes PROFILE numbers (Pearson r, saliency) carry a memorisation component. Hepatocytes
ORF-CALL numbers (precision / recall / F1) appear not to. Where a clean substrate is required --
notably sequence attribution -- mouse is used as primary, with Hepatocytes retained deliberately as
the contaminated arm so the two can be compared.
```
