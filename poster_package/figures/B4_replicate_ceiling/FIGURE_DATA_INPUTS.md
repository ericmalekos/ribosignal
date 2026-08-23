# FIGURE B4 -- replicate-concordance ceiling: data inputs

**What it shows:** the Ribo-seq profile's reproducible ceiling (replicate concordance,
Spearman-Brown-corrected to full depth) against the released model's achieved per-transcript Pearson,
over four matched strata. The honest "how much signal exists / how much of it is captured" panel.

**Generator:** `make_replicate_ceiling.py` (cas12a env). Both sides are now READ FROM FILES; nothing is
hardcoded, so a stale figure is detectable rather than plausible.

| side | source file | what is read |
|---|---|---|
| ceiling | `results/replicate_concordance.json` | per-seed median half-pool Pearson per stratum, then `2r/(1+r)` and averaged over the 5 seeds |
| model | `results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/extra_metrics.json` | `<stratum>.profile_pearson_median` |

The run name is a module constant `RUN` in the generator and is written into `B4_values.json`, so the
figure's provenance can be checked without re-reading the code.

## Panel values (released attn model, union universe)

| stratum | ceiling key | model key | ceiling | model | % of ceiling |
|---|---|---|---:|---:|---:|
| PC whole-tx | `pc_whole` | `protein_coding` | 0.952 | 0.670 | 70% |
| lncRNA whole-tx | `lncRNA_whole` | `lncRNA` | 0.890 | 0.458 | 51% |
| PC 5'UTR (uORF) | `pc_uorf5` | `uorf_5utr` | 0.953 | 0.582 | 61% |
| PC 3'UTR (dORF) | `pc_dorf3` | `dorf_3utr` | 0.749 | 0.255 | 34% |

Note the dORF ceiling is itself low (0.749): 3'UTR Ribo-seq signal does not reproduce well between
half-pools, so the model's weak dORF number is partly a property of the assay.

## Two corrections applied 2026-08-08

1. **A bar compared two different quantities.** The old third bar was `pc_cds_codon`: ceiling 0.9501
   (per-codon CDS *elongation shape*, periodicity deliberately removed) against a model value of 0.526
   that was a *periodicity* score. `scripts/replicate_concordance.py` records the model side of that
   stratum as `None` with the comment `our per-codon TBD` -- it was never computed. The derived
   "55% of ceiling" was therefore meaningless. That bar is replaced by the two strata the concordance
   script was written to pair with the eval (`pc_uorf5` / `pc_dorf3` against `uorf_5utr` / `dorf_3utr`,
   same windows, same transcripts). The per-codon CDS comparison can return once the model side is
   actually computed.
2. **The model values were pre-union.** `MODEL` was a hardcoded dict (pc 0.638, lncRNA 0.364) from the
   pre-union eval, while the project ships the union models (pc 0.670, lncRNA 0.458). It is now read
   from `extra_metrics.json`.

## Regenerate

```bash
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/figures/B4_replicate_ceiling
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 make_replicate_ceiling.py
```

Re-run after any retrain. If a stratum is missing from either source JSON the generator exits with a
message naming the file, rather than silently dropping a bar.
