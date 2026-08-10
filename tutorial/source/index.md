# Ribo-seq signal model: predicting translation from sequence + RNA-seq

A **code-first** walkthrough of a model that predicts the **per-nucleotide ribosome P-site profile** of a
transcript from its **sequence and RNA-seq coverage alone** -- no ribosome-profiling experiment required.
Feed the predicted profile to a standard ORF caller (RiboCode) and you get **de novo ORF calls** -- canonical
CDS plus uORFs, novel ORFs, and dORFs -- for any transcript in any tissue you have RNA-seq for. Prose is
short on purpose; the data, the code in `scripts/`, and the numbers are the point.

## What the model does

```{mermaid}
flowchart LR
    S["Transcript sequence<br/><small>one-hot or FM embedding</small>"]
    R["RNA-seq coverage<br/><small>per-nt, unique-mappers</small>"]
    O["ORF track<br/><small>reading-frame prior from sequence</small>"]
    M["Dilated CNN + mixer<br/><small>transformer / Mamba, ~5M params</small>"]
    P["Predicted P-site profile<br/><small>per-nt distribution + total</small>"]
    C["RiboCode drop-in<br/><small>de novo ORF calls</small>"]
    S --> M
    R --> M
    O --> M
    M --> P --> C
    classDef s fill:#eaf2fb,stroke:#3a6ea5,color:#16314a;
    class S,R,O,M,P,C s;
```

The model is a dual-head dilated-CNN with a global-context mixer (`scripts/model.py`). It is trained to
reproduce the observed P-site profile of held-in tissues, then applied to a held-out tissue or an entirely
new dataset. Because it learns the **shape** of translation from sequence + expression, it can stand in for a
Ribo-seq experiment where none exists.

## The headline

- **It works cross-tissue and cross-species.** The deployed model recovers canonical ORFs at **F1 ~0.94-0.97**
  and the non-canonical (uORF/novel/dORF) tail at F1 ~0.5, on a held-out human tissue *and* on mouse cell
  types it never saw (liver, BMDM, CD4+ T cells).
- **It matches a real experiment where it counts, and the ceiling says how much that is worth.** Against
  *measured* translation in the same cell line (B721.221, 327 M footprints, model never sees the Ribo-seq),
  the standalone model reaches **75% of the assay's own split-half reproducibility ceiling** on F1 and 78% on
  canonical recall. Non-canonical recall is 54% of ceiling -- and the ceiling itself is only **0.498**, because
  two independent measurements of the same cells agree on barely half their non-canonical calls. Judge
  non-canonical performance against ~0.5, never against 1.0.
- **Sequence encoding: one-hot is enough.** Swapping the one-hot input for RNA foundation-model embeddings
  (RiNALMo, Orthrus, HydraRNA) gives **no lift** for ORF calling.
- **The RNA-seq does not have to come from the same experiment.** A 3 Ribo-seq x 3 RNA-seq factorial on
  mouse liver -- all nine cells pooled through one pipeline onto one shared universe -- finds that matching
  the RNA input to the Ribo-seq reference is worth **+0.003 F1** on average, and is *negative* in 2 of 12
  comparisons. What matters is RNA library quality, not provenance: a 2-sample RNA arm costs ~0.008 F1,
  while 19 samples buy nothing over 7 ({doc}`data/heldout-mouse`).
- **Two models ship: `mamba4` (primary) and `attn` (CPU-compatible).** Across 3 seeds on the broad universe,
  `mamba4` beats `attn` by **+0.021** held-out test Pearson (0.680 vs 0.659) with **non-overlapping** seed
  ranges, and wins **5/5** mouse datasets on cross-species CDS F1. But `mamba4` **cannot run on CPU**
  (`causal_conv1d` is CUDA-only), so `attn` is maintained and released for CPU-only inference. Within the
  transformer family, depth is seed-fragile: 4 layers looked best on one seed but the gap vanishes across
  seeds, and 6 layers *regresses*. Profile Pearson does **not** predict ORF-call quality, so we select on a
  frame-aware metric.
- **Calibration is the real lever.** A Poisson count-scaling arm trades non-canonical recall for precision at
  a fixed canonical-recall target -- more useful than any architecture change. It is also **stochastic**:
  counts reproduce to CV 0.21% across seeds, but call *membership* only to Jaccard 0.822, so a single draw is
  one sample from a distribution ({doc}`calling`).
- **Predicted ORFs carry independent peptide evidence.** Feeding the calls to a mass-spec search finds novel
  peptides at **123x** the discovery density of a naive enumeration and at **28x** lower cost to canonical
  detection -- and recovers peptides a 154x larger database contains but cannot report at 1% FDR
  ({doc}`proteogenomics`).

## How to read this site

Each page is short and focused. **Setup** and **The model** cover the machinery; **Training data** and
**Test sets** define exactly what the model learned from and is judged on; **Calling** covers turning a
predicted profile into ORF calls, including the calibration dial and its randomness; **Metrics** defines
every number we track; **Benchmarks** is the living results table + figure, regenerated on each build;
**Proteogenomics** is the independent mass-spec test of the predictions.

```{toctree}
:maxdepth: 2
:caption: Setup

setup/environments
```

```{toctree}
:maxdepth: 2
:caption: The model

model/architecture
```

```{toctree}
:maxdepth: 2
:caption: Data

data/training-data
data/test-sets
data/heldout-human
data/heldout-mouse
```

```{toctree}
:maxdepth: 2
:caption: Calling ORFs

calling
```

```{toctree}
:maxdepth: 2
:caption: Evaluation

metrics
benchmarks
proteogenomics
```
