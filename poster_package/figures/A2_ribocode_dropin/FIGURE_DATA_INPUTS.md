# FIGURE A2 -- RiboCode drop-in concordance: data inputs

> **SUPERSEDED NUMBERS BELOW (flagged 2026-08-16).** The F1 0.923 / annotated n=10,947 figures in
> this file describe an OLDER checkpoint than the one its own `A2_values*.json` now holds. The JSON
> is current; this prose is not. A2 now also reads the FINAL-RECIPE hepatocyte dropin
> (`dropin_canonpack`) where present, and records `model`, `on_final_recipe`, `pack` and `source`,
> so provenance no longer has to be inferred by float comparison. Re-derive any number here from
> the JSON before quoting it.


**What it shows:** the model's PREDICTED per-nt profile, fed into RiboCode in place of real Ribo-seq,
reproduces RiboCode's OWN ORF calls at F1 0.923 (no Ribo-seq for the query). Panel a = precision/recall/F1;
panel b = recall by ORF type.

**Generator:** `make_ribocode_dropin.py` (cas12a). Reusable -- reads the JSON below; regenerates if retrained.

## Data source
- `results/loto/orf_v2_attn_onehot_holdout_Hepatocytes/dropin/dropin_metrics.json` (Task 14, held-out Hepatocytes).
- Comparison used = `pred_preddepth_vs_real` (fully predicted: predicted profile AND predicted depth, vs
  RiboCode's calls on the REAL data). `pred_obsdepth_vs_real` (predicted profile, observed depth) shown for reference.

## Numbers (current model, multimap-20 RNA coverage)
- Fully predicted vs real RiboCode calls: **P 0.925 / R 0.920 / F1 0.923**.
- Recall per ORF type: annotated 0.996 (n=10,947), novel 0.785 (n=377), uORF 0.618, Overlap_uORF 0.608,
  Overlap_dORF 0.255, dORF 0.101, internal 0.070.

## Caveats
- Low dORF/internal recall is expected -- the model's predicted profile strips spurious 3'UTR dORFs and
  internal in-frame fragments the real caller over-calls (this is a feature, not a miss; precision stays 0.925).
- Current model = multimap-20 RNA coverage; re-run vs the mm1-ablation model (task #25) to update if it retrains.

## Regenerate
```bash
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/figures/A2_ribocode_dropin
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 make_ribocode_dropin.py
```
