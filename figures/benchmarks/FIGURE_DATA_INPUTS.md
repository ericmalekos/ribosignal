# FIGURE benchmarks -- standalone ORF-call quality per model, by evaluation universe

**What it shows:** per-model standalone ORF-call quality at the uncalibrated operating point
(theta = 1) -- non-canonical F1 and novel-ORF precision -- grouped by evaluation universe. It answers
"which trained model actually calls better ORFs", separately from profile Pearson.

Universes plotted (models sorted within each by non-canonical F1):

| key | display | role |
|---|---|---|
| `fibroblast` | Fibroblast universe | in-distribution |
| `union` | Union universe | broad |
| `mouse-tcell` | Mouse T-cell | cross-species held-out |

## This is a LIVING figure

Unlike the other figures here, it is not a one-off render. `plot_benchmarks.py` globs whatever
per-model JSONs exist and is re-run on every `make html`, so it always tracks the latest evaluations.
Consequences:

- The set of bars CHANGES as models are added or removed. A panel from an older tutorial build is not
  necessarily reproducible from the current JSON directory.
- It no-ops cleanly when no JSONs are present, so a missing figure means "no evals landed", not "the
  build failed".
- Do not hand-edit the PNG or treat a specific bar ordering as stable.

## Data inputs

- `results/orf_call_metrics/*.json` -- one per model, written by
  `scripts/heldout/orf_call_metrics.py`. Each carries the model label, evaluation universe key, and
  the standalone theta = 1 metrics.

## Regeneration

```
# any python with matplotlib, e.g. conda_envs/cas12a
python scripts/tutorial/plot_benchmarks.py
# or implicitly:
cd tutorial && make html
```

Output: `figures/benchmarks/benchmark_comparison.png`, consumed by `tutorial/source/benchmarks.md`
(which references it as `/img/benchmarks/benchmark_comparison.png` after the tutorial's image copy).

## Caveat

Bars are the **theta = 1 uncalibrated** arm only. Per the standing two-arm rule
(`feedback_two_arm_orf_calling`), an uncalibrated number alone overstates novel-ORF calling; the
calibrated (Poisson, CDS-anchored theta) arm lives in the ORF-call eval tables, not in this figure.
Read it as a relative comparison BETWEEN models at a fixed operating point, not as the achievable
precision of any one of them.
