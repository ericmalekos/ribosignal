# S-REPRO -- training-data regeneration across 8 tissues: data inputs

Can the published pipeline reproduce the arrays the model trained on? Measured on all 8 Chothani
tissues, not one.

## Panel -> data

| element | source | produced by |
|---|---|---|
| both panels, every point | `results/chothani_regeneration/divergence_by_tissue.tsv` | `scripts/summarise_chothani_divergence.py` |
| per-tissue divergence JSONs behind it | `results/chothani_regeneration/_canonical_<T>/divergence_canonical.json` | `scripts/chothani_canonical_diverge.sbatch` (array 0-7) |
| the canonical re-alignments | `data/ribo_bam_e2e_human/<T>/*.Aligned.toTranscriptome.out.bam` (74 runs) | `scripts/riboseq_align.sbatch` + `data/chothani_ribo_samplesheet.tsv` |
| the training packs compared against | `data/packed_union` (Fibroblast) and `data/packed_union_<T>` | the 2026-07 build |

## The metric

`sum|new - old| / sum(old)` over per-nt P-sites, on each pack's OWN universe and row order (84,472
tx). Reported alongside per-nt Pearson r because the two answer different questions -- how much
signal moved, versus whether it moved coherently. Net change is in the TSV as well.

## Validation before trusting the other 7

The generalised script reproduces the original single-tissue HUVEC result EXACTLY: 3.04% of signal
moved, r = 0.9567, against the previously reported 3.04% / 0.957. That positive control is why the
remaining seven tissues are trustworthy.

## Which pack belongs to Fibroblast

`data/packed_union`, NOT `data/packed`. The latter is the pre-union single-tissue Fibroblast pack
(36,668 tx vs 84,472) whose `coverage_norm.json` records `n_ribo_samples: 0` -- it is coverage-only,
so diffing P-sites against it would compare reads to a pack that has none. All 8 tissues share the
84,472-tx union universe, which is what makes the numbers comparable to each other.

## Caveats a caption must carry

- **HUVEC is the best of the eight, not a typical one.** 3.04% vs a median of 3.86% and a worst case
  of 11.66%. Any single-tissue quotation of this claim understates it.
- **HA_EC (r = 0.432) is an unexplained outlier.** Four hypotheses were tested and rejected: uniform
  P-site shift (best-lag correlation peaks at lag 0), defective old pack (all 8 old packs' start-codon
  metagenes peak at +0 nt with frame0 0.732-0.825; HA_EC is mid-range at 0.744), changed P-site
  calling (old vs new frame0 agree within tissue; HA_EC +0.028 vs HUVEC +0.027, and HUVEC has
  r = 0.957), and metaplots offset disagreement (recorded as a prediction, then falsified -- ES has
  offset disagreement and the HIGHEST r at 0.980). **No causal story about HA_EC belongs in a
  caption.**
- **No depth axis, deliberately.** Hepatocytes (1.16B P-sites) and HUVEC (96M) both diverge at ~3.1%
  across a 12x depth range, and the worst tissue sits mid-range on depth. A depth panel would imply
  a trend the data does not support.
- RNA-seq is not part of this test. Only Ribo-seq was re-aligned; coverage is unchanged.

## Regenerate

```bash
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
cd figures/S_reproducibility && $PY make_reproducibility.py
```

Upstream, if the TSV is missing:

```bash
sbatch --array=0-7 scripts/chothani_canonical_diverge.sbatch
$PY scripts/summarise_chothani_divergence.py --dir results/chothani_regeneration \
  --out-tsv results/chothani_regeneration/divergence_by_tissue.tsv \
  --out-md  results/chothani_regeneration/divergence_by_tissue.md
```
