# FIGURE B4 -- replicate-concordance ceiling: data inputs
**What it shows:** the Ribo-seq profile's reproducible ceiling (replicate concordance, Spearman-Brown-
corrected) vs the model's achieved per-tx Pearson. Honest "how much signal exists / how much captured."
**Generator:** `make_replicate_ceiling.py` (cas12a). Ceiling read from `results/replicate_concordance.json`
(32 Fibroblast samples, two 16-sample half-pools, 5 seeds, Spearman-Brown to full depth). MODEL dict = current
orf_v2_attn per-tx Pearson (results.md Task 22) -- UPDATE on mm1 retrain.
| stratum | ceiling | model | % |
|---|---:|---:|---:|
| PC whole-tx | 0.952 | 0.638 | 67% |
| lncRNA whole-tx | 0.890 | 0.364 | 41% |
| PC CDS per-codon | 0.950 | 0.526 | 55% |

## Regenerate
```bash
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/figures/B4_replicate_ceiling
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 make_replicate_ceiling.py
```
