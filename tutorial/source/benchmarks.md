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
