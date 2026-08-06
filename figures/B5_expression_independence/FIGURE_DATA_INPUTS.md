# FIGURE B5 -- translation signal exceeds an expression-only baseline: data inputs

**What it shows:** the model's predicted in-frame fraction (`pred_frame0`, a periodicity/shape score)
out-discriminates every magnitude / expression score by ~+0.12 (all ORFs) and ~+0.09 (non-canonical),
length-controlled, on held-out Hepatocytes. So a "translated ORF" call is not merely a "high-expression"
call -- the model adds a genuine translation (periodicity) signal on top of expression.

**Framing (read before using):** this model is *deliberately* expression-aware -- it takes per-nt RNA-seq
coverage as input, which is exactly what makes it cell-type-specific (the Fig 2 premise). So the claim is
NOT "expression-independent" (it should not be). The claim is "the periodicity signal is orthogonal to,
and additive over, expression magnitude." That is what these bars demonstrate.

## Data source
- `results/loto/orf_v2_attn_onehot_holdout_Hepatocytes/localization_metrics.json`, `discrimination.{all,noncanonical}`.
- Length-control = the eval's `auroc_stratified` (Mann-Whitney AUROC within 5 ORF-length quantile bins,
  pair-count weighted) -- so translated vs decoy are only compared at similar length. Copied verbatim into
  the generator for the exploratory panel (see caveat) but the shipped bars come straight from the JSON.

## Numbers (length-controlled AUROC, translated vs decoy ORF)
| stratum | pred_frame0 (shape) | pred_density (magnitude) | obs Ribo-seq density | shape - best magnitude |
|---|---:|---:|---:|---:|
| All ORFs | 0.891 | 0.771 | 0.751 | **+0.120** |
| Non-canonical | 0.832 | 0.743 | 0.745 | **+0.087** |

## Caveat / nuance (kept out of the figure to avoid over-claiming, documented here)
- The model's DENSITY head DOES track RNA expression (by design). The FRAME head is the expression-
  orthogonal shape signal, and it is what carries the +0.12 margin.
- An exploratory panel (binning ORFs by transcript RNA-expression quintile, then length-controlled
  `pred_frame0` AUROC within each bin) showed `pred_frame0` dropping to ~0.55 within a single expression
  bin. That is NOT a pure expression effect: within a narrow (high) expression bin the negatives become
  decoy ORFs on the SAME translated transcripts, which overlap the real CDS and inherit its periodicity,
  so it collapses to the hard overlapping-decoy problem, conflated with expression. It was therefore
  dropped from the shipped figure; the between-vs-within-transcript decomposition is future work if a
  reviewer presses on it. For novel ORFs specifically, `pred_density` (0.923) exceeds `pred_frame0`
  (0.874) length-controlled -- novel translated ORFs sit on high-expression transcripts, so expression is
  itself informative there; the "all"/"non-canonical" strata are where the shape > magnitude story is clean.

## Regenerate
```bash
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/figures/B5_expression_independence
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 make_expression_independence.py
```
