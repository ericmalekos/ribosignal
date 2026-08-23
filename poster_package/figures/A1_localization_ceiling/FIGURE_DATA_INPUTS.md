# FIGURE A1 -- localization-beats-ceiling: data inputs

**What it shows:** the model's PREDICTED translation-localization (`pred_frame0` AUROC) matches or beats the
OBSERVED Ribo-seq periodicity ceiling on held-out data. Grouped bars, length-controlled AUROC.

**Generator:** `make_localization_ceiling.py` (cas12a env). Reusable -- reads the JSONs below, so it
regenerates automatically if the model is retrained (e.g. the mm1 RNA-coverage ablation).

## Panels / bars -> data

| bar group | metric | source file | JSON path |
|---|---|---|---|
| Held-out Hepatocytes, all ORFs | `auroc_pred_frame0_lengthctrl` (model), `auroc_obs_frame0_ceiling_lengthctrl` (observed) | `results/loto/orf_v2_attn_onehot_holdout_Hepatocytes/localization_metrics.json` | `discrimination.all` |
| Held-out Hepatocytes, non-canonical | same keys | same file | `discrimination.noncanonical` |
| Cross-study Ruiz-Orera, all ORFs | same keys | `results/heldout/human_ruizorera/onehot/localization_metrics_human_ruizorera.json` | `discrimination.all` |
| Cross-study Ruiz-Orera, non-canonical | same keys | same file | `discrimination.noncanonical` |

`A1_values.json` (written by the generator) caches the exact numbers used, incl. the raw (non-length-ctrl)
AUROCs for the caption.

## Numbers (length-controlled), current model (multimap-20 RNA coverage)

| dataset / stratum | model pred | observed ceiling | verdict |
|---|---:|---:|---|
| Hepatocytes, all | 0.891 | 0.828 | BEATS ceiling |
| Hepatocytes, non-canonical | 0.832 | 0.879 | 95% of ceiling |
| Ruiz-Orera, all | 0.894 | 0.825 | BEATS ceiling |
| Ruiz-Orera, non-canonical | 0.842 | 0.874 | 96% of ceiling |

## Provenance / caveats
- Model = deployed `orf_v2_attn` one-hot, LOTO Hepatocytes fold. Ruiz-Orera = cross-study human held-out.
- These metrics are from the CURRENT model (RNA coverage built multimap-20). If the mm1-coverage ablation
  (task #25) produces a retrained model, re-run the generator against the new run dir to update.
- "Observed ceiling" = how well the OBSERVED Ribo-seq periodicity itself separates translated vs not --
  the model beating it on the all-stratum means the prediction is cleaner than the measurement on those ORFs.

## Regenerate
```bash
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/figures/A1_localization_ceiling
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 make_localization_ceiling.py
```
