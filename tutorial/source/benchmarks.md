# Benchmarks

All numbers are **standalone $\theta$=1** ORF-call quality on the Hepatocytes hold-out (20-aa primary floor),
unless noted. The figure below is regenerated from `results/orf_call_metrics/*.json` on every build, so it
tracks the latest evals.

```{figure} /img/benchmarks/benchmark_comparison.png
:width: 100%

Per-model ORF-call quality. **Left:** non-canonical F1 (standalone $\theta$=1), grouped by evaluation
universe. **Right:** novel-ORF precision -- bar is the $\theta$=1 point, diamond is the calibrated Poisson
arm. Deeper transformer (`attn4`) leads the in-distribution fibroblast group; foundation-model embeddings
(orthrus/hydrarna/rinalmo) sit below the one-hot models; the mouse T-cell held-out (green) generalises cleanly.
```

## Architecture: depth helps, foundation-model embeddings do not

```{list-table}
:header-rows: 1
:widths: 34 16 16 16 18

* - Model (fibroblast universe)
  - CDS F1
  - novel P
  - non-canon F1
  - Poisson novel P
* - **attn4** (4 transformer layers)
  - 0.945
  - 0.437
  - **0.593**
  - **0.686**
* - mamba (2 SSM layers)
  - 0.938
  - 0.401
  - 0.541
  - 0.581
* - onehot 2-attn (deployed)
  - 0.945
  - 0.397
  - 0.535
  - 0.605
* - orthrus (FM, 512-d)
  - 0.942
  - 0.372
  - 0.523
  - 0.505
* - hydrarna (FM, 1024-d)
  - 0.940
  - 0.380
  - 0.516
  - 0.579
* - rinalmo (FM, 1280-d)
  - 0.936
  - 0.379
  - 0.511
  - 0.566
```

- **Depth does not reliably help -- it is within seed noise.** A single seed suggested `attn4` beat the
  deployed 2-layer model by +0.058 non-canonical F1. A multi-seed repeat overturns that: attn2 spans
  **0.513-0.581** across seeds (spread 0.068), and attn2's best seed (0.581) *beats* attn4's (0.549). The mean
  gap (~+0.028) is smaller than the attn2 seed spread, so the deployed 2-layer model sits well inside the band.
- **RNA foundation-model embeddings give no lift.** All three FMs land below both one-hot architectures.
  One-hot sequence + RNA coverage already carries the signal.

```{warning}
The single-seed comparisons below are shown for completeness, but **the architecture differences are
seed-fragile** -- rank models with error bars, not one run. Profile Pearson also does *not* track ORF-call
quality, so it is never the selection metric ({doc}`metrics`).
```

## Depth sweep: transformers saturate, Mamba scales

Sweeping mixer depth on the union universe (single-seed, same 17,284-call reference):

```{list-table}
:header-rows: 1
:widths: 20 26 26 28

* - Layers
  - Transformer (non-canon F1)
  - Mamba (non-canon F1)
  - Notes
* - 2
  - 0.547
  - 0.548
  - deployed depth
* - 4
  - 0.566
  - **0.570**
  - mamba4 = best union model
* - 6
  - **0.538** (regresses)
  - n/a
  - attn6 < attn4 < 2-layer
```

- **Transformer depth peaks at 4 then regresses** -- attn6 falls *below* even the 2-layer model. Adding
  attention layers is not a path forward.
- **Mamba keeps improving with depth** (mamba4 is the best union model) and scales where attention does not --
  the one genuinely promising direction. A **mamba seed confirm is running** to test whether that lead is
  robust or, like attn4's, a favorable seed.

## Universe: breadth trades canonical precision for discovery

The **union** universe (70,883 vs 33,918 held-out transcripts) carries **827 novel ORFs in its reference vs
454** for fibroblast -- ~80% more discoverable novel ORFs, because it includes Hepatocytes-expressed
transcripts the fibroblast universe omits. The cost is a small canonical-precision dip (CDS F1 0.92 vs 0.94).
Union-vs-fibroblast is not like-for-like (the reference sets differ), so it measures *coverage*, not a
head-to-head win.

## Two more results worth the figure

```{figure} /img/depth_crossover/depth_crossover.png
:width: 90%

Predict-vs-measure crossover: below a per-sample Ribo-seq depth, the model's **predicted** P-sites call ORFs
*better* than a shallow real experiment -- so for low-depth samples, prediction is not just a substitute but an
improvement.
```

```{figure} /img/prediction_examples/prediction_examples.png
:width: 100%

Example transcripts: predicted per-nt P-site profile (model) vs observed (experiment). The model recovers the
3-nt periodic CDS envelope and the ORF boundaries from sequence + RNA coverage alone.
```

## Ground truth: predicted calls vs *measured* translation in the same cells

Everything above scores the model against RiboCode run on the **same pack** it was evaluated on, or against a
null. Neither asks the sharper question: are the ORFs the model calls the ones actually being translated?

B721.221 answers it, because one cell line has both halves from independent labs. Sarkizova RNA-seq drives
the prediction; Ouspenskaia Ribo-seq (327 M unique footprints, 7 runs) supplies the measured calls. **The
model never sees the Ribo-seq.**

Scoring is genomic-keyed `(gene_id, ORF_gstop)` and restricted to the 11,527 genes in the model's universe.
The two call sets are *not* made over the same transcript space -- measured comes from the full prepared
annotation (21,974 calls, 16,331 in the model's gene space) -- and RiboCode's collapse keeps the longest ORF
per genomic stop, so restricting the transcript set reshuffles which isoform carries a call. Genomic keying is
invariant to that; transcript keying would score isoform bookkeeping as disagreement.

| arm | precision | recall | F1 | recall canonical | recall non-canonical |
|---|--:|--:|--:|--:|--:|
| mamba4 $\theta$=1 | 0.780 | 0.580 | 0.665 | 0.754 | 0.241 |
| attn $\theta$=1 | 0.761 | 0.592 | **0.666** | 0.759 | 0.267 |
| mamba4 Poisson | 0.887 | 0.506 | 0.644 | 0.714 | 0.099 |
| attn Poisson | **0.918** | 0.491 | 0.640 | 0.705 | 0.074 |

### The number is meaningless without a ceiling

An F1 of 0.67 reads very differently against a 0.95 bar than against a 0.70 one. So RiboCode was run
independently on **two disjoint halves of the same Ribo-seq**, split by HLA allele so each half holds complete
biological libraries (never two technical re-sequencings of one library split across halves, which would make
the halves more similar than two real measurements and understate the noise).

The ceiling is taken as the *shallower* half predicting the *deeper* one, because that is the geometry the
model faces -- one measurement scored against a larger reference:

| | F1 | recall canonical | recall non-canonical |
|---|--:|--:|--:|
| **ceiling** (half A vs half B) | **0.889** | 0.973 | **0.498** |
| attn $\theta$=1 | 0.666 (**75%**) | 0.759 (**78%**) | 0.267 (**54%**) |

```{admonition} Two independent measurements of the same cells agree on only half their non-canonical calls
:class: important
The 0.24-0.27 non-canonical recall looks alarming against a denominator of 1.0. The real denominator is
**0.498**. Non-canonical ORF calling is intrinsically ~50% reproducible at this depth, so the model recovers
roughly *half what a replicate Ribo-seq experiment would* -- not a quarter of perfect.

The same ~0.50 shows up in the between-experiment mouse-liver ceiling (Wang vs Janich, non-canonical F1
0.505). Canonical, by contrast, drops from 0.973 within-experiment to 0.735 between-experiment -- that gap is
where the batch effect actually lives.
```

Robust across the one free parameter: sweeping minimum ORF length 0 / 30 / 90 / 150 nt leaves every
conclusion intact. Canonical recall is essentially invariant (0.754-0.759 at all four cuts, since annotated
CDSs are long), so the F1 gain with longer cuts comes entirely from removing short non-canonical ORFs from the
reference, not from the model improving. Non-canonical recall actually *falls* as the cut rises
(0.261 -> 0.228), so the model is relatively better on short non-canonical ORFs than long ones.

Two caveats that cap this comparison: the RNA-seq and Ribo-seq come from *differently HLA-transduced*
B721.221 sub-lines, so it is same-parental-line but not same-sample; and B721.221 is a lymphoblastoid line,
far outside the Chothani training tissues.
