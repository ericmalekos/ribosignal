# FIGURE A3 -- even cross-tissue generalization (9-fold LOTO): data inputs

**What it shows:** across 9 leave-one-tissue-out folds, whole-tx profile Pearson looks spread out
(0.23 Brain to 0.64 Hepatocytes), but that spread is 0.81-correlated with the held-out tissue's OWN
observed 3-nt periodicity (period_obs) -- i.e. it is held-out data QUALITY, not a generalization
gradient. Abundance/count transfer is even across tissues (pc count Pearson 0.59-0.88, CV 0.11).
Panel a = profile Pearson vs held-out period_obs (r=0.81); panel b = count Pearson across 9 tissues.

**Generator:** `make_loto_spread.py` (cas12a). Reads the 9 per-fold JSONs, regenerates on re-eval.

## Data source
- `results/loto_9fold/onehot_orf_v2_attn_holdout_<tissue>/test_metrics.json` -> `protein_coding.pearson_median`,
  `protein_coding.period_obs_median`, `protein_coding.period_pred_median`.
- `.../extra_metrics.json` -> `protein_coding.count_pearson`.
- Tissues (9): Brain, ES, Fat, Fibroblast, HA_EC, HCAEC, Hepatocytes, HUVEC, VSMC. Deployed model =
  onehot orf_v2_attn, trained on the other 8 tissues per fold.

## Numbers
- pc profile Pearson: 0.233 (Brain) to 0.638 (Hepatocytes), CV 0.24.
- pc count Pearson: 0.587 (Brain) to 0.882 (Fibroblast), CV 0.11, median 0.827.
- corr(period_obs, pc profile Pearson) = **0.811** (r^2 = 0.66): two-thirds of the profile-Pearson spread
  is explained by target periodicity. Brain (period_obs 0.044, lowest-quality target) is the low outlier.
- The model's period_pred stays 0.30-0.45 across folds regardless of target quality -- the model always
  predicts periodicity; low-Pearson folds are where the TARGET lacks it.

## Caveats / upgrade path
- This uses per-fold WHOLE-TX profile + count Pearson (what loto_9fold stores). The plan's ideal Fig 1a
  is per-fold drop-in F1 / localization AUROC (target-quality-robust). Those were only computed for the
  Hepatocytes fold; per-fold F1/AUROC would require re-running `eval_localization.py` + the RiboCode
  drop-in on each of the 9 `loto_9fold/*/best.pt` checkpoints (checkpoints are present; per-tissue packed
  targets under `data/packed_<Tissue>/` needed). Flagged for the final Main Fig 1a.
- period_obs is the median over pc transcripts of the observed in-frame fraction; it is the cleanest
  available proxy for held-out Ribo-seq quality.

## Regenerate
```bash
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/figures/A3_loto_spread
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 make_loto_spread.py
```
