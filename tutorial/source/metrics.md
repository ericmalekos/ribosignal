# Metrics

We track two families of numbers: **profile metrics** (does the predicted P-site curve match the observed
one?) and **ORF-call metrics** (do the ORFs called from the prediction match the ORFs called from the real
experiment?). The central lesson of the project is that **these two do not track each other** -- so the
ORF-call metrics, not profile Pearson, decide which model ships.

## Profile metrics

`val_pearson` (profile fidelity)
: Per-transcript Pearson correlation between the predicted position distribution $p$ and the observed P-site
  distribution $q = c/\sum c$, reported as the **median over validation transcripts**. Scale-invariant (it
  compares shape, not depth). Used for **early stopping**.
: **Caveat -- it saturates.** Pearson is dominated by coverage magnitude, which the RNA input already
  supplies, so it plateaus at ~0.52 and a *better* ORF-caller can have a *lower* val_pearson. Do not rank
  models on it.

`frame0_pred_median` (ORF-aware periodicity)
: The fraction of predicted P-site mass in **frame 0** (the first nt of each codon) inside annotated CDS:
  $\text{frame0} = \sum_{p \equiv \text{CDS start} \,(\mathrm{mod}\,3)} \hat{y}_p \big/ \sum_{\text{CDS}} \hat{y}$.
  Random = 0.33; clean Ribo-seq observed ~0.86 (the ceiling). Unlike val_pearson, frame-0 **does** rank the
  architectures in the same order as ORF-call quality, so it is logged per epoch as a candidate
  checkpoint-selection signal (`best_frame0.pt` alongside `best.pt`).

## ORF-call metrics

The deliverable. RiboCode calls ORFs on the predicted density; each call is matched to the observed call set
on a genomic coordinate key (`gstart_gstop_len`, annotation-version robust), **per ORF class**: `annotated`
(canonical CDS), `uORF`, `novel`, `dORF`, and the `non-canonical` aggregate.

precision / recall / F1
: Standard set overlap of predicted vs observed calls, per class. **Precision** = fraction of predicted calls
  that are real; **recall** = fraction of observed calls recovered.

over-call ratio
: $n_\text{predicted} / n_\text{observed}$ per class. $>1$ means the model calls more ORFs than the experiment
  saw -- the failure mode for the non-canonical tail (a standalone prediction with no Ribo-seq tends to
  over-call novel ORFs ~2x).

## The two-arm protocol

Every ORF-call evaluation reports **two operating points**, because a raw standalone prediction over-calls:

```{list-table}
:header-rows: 1
:widths: 26 74

* - Arm
  - What it is
* - **standalone $\theta$=1**
  - deterministic: predicted shape x predicted total, no calibration. Exposes the over-calling honestly.
* - **Poisson CDS-anchored**
  - the count total is scaled by $\theta$ and re-sampled (Poisson), sweeping $\theta$ until canonical (CDS)
    recall hits a target (~0.90 of the **achievable** CDS, not of the whole annotation). Reports
    non-canonical precision at that fixed canonical-recall point. See {doc}`calling` for how $\theta$* is
    chosen and why the anchor is relative.
```

The Poisson arm is the **single most useful lever**: at matched CDS recall it raises novel-ORF precision by
~0.1 over the deterministic arm -- more than any architecture change buys.

```{warning}
The Poisson arm is **stochastic**: same model, same $\theta$, different seed, different call set. Counts are
reproducible (CV 0.21% over 5 seeds) but membership is not -- pairwise Jaccard **0.822**, and only 66.9% of
the union is called in all five draws. Always record the seed, never compare arms across different seeds,
and use `pgx.consensus` when the call list is an end product. Full numbers in {doc}`calling`.
```

## Minimum-ORF-length sweep

Non-canonical calls are dominated by short ORFs, so every model reports P/R/F1 restricted to ORFs
$\ge \{5, 10, 15, 20, 30, 40\}$ amino acids. **Primary floor = 20 aa** (RiboCode's own default). Calls are
made at a permissive 5 aa and length-filtered post hoc; the sweep shows precision rising monotonically as the
tiny-ORF tail is excluded, and confirms model *rankings* are stable across the floor.

```{admonition} Reading a benchmark row
:class: tip
`CDS 0.94/0.95/0.95` = canonical precision/recall/F1 (trusted). `novel 0.40/0.82` = novel precision/recall.
`non-canon F1 0.54` = F1 over all non-CDS ORFs. `over-call 2.0x` = twice as many novel calls as the
experiment saw at $\theta$=1. Higher non-canonical precision **at a fixed CDS recall** (Poisson arm) is the
number that matters.
```
