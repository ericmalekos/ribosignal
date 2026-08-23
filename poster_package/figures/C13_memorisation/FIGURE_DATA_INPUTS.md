# C13 -- memorisation control at matched depth: data inputs

Three tiers of held-out Hepatocytes transcripts, per model. The point of the figure is that the
naive comparison (all-seen vs never-seen) is inflated ~4x by depth.

## Panel -> data

| element | source | produced by |
|---|---|---|
| all three bars, CIs, and n per tier | `results/memorisation_control/hepatocytes_{attn,mamba4}.json` | `scripts/memorisation_control.py` |
| the predicted profiles behind them | `results/loto/orf_v2_{attn,mamba4}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin/pred_profiles.npz` | the deployed LOTO runs |
| which transcripts count as "never seen" | transcripts absent from EVERY training tissue's pack at `--min_train_signal` | `scripts/memorisation_control.py` |

## The three tiers

| tier | n | what it is |
|---|--:|---|
| all seen | 70,785 | every held-out transcript the model also trained on. Unmatched on depth -- this is the misleading number the figure exists to intercept |
| seen, depth-matched | 980 | 10 seen transcripts per unseen one, matched on log10 total P-sites and on length (median counts 67 vs 66) |
| never seen | 98 | absent from all 7 training tissues |

## Why depth matching is the whole experiment

Profile Pearson rises steeply with counts, and the never-seen transcripts are shallow. Without
matching, "memorisation" absorbs the entire depth confound: attn 0.468, mamba4 0.492. With matching
it is 0.100 / 0.108, with disjoint 95% CIs. Both numbers are in `C13_values.json` (`naive_delta` and
`matched_delta`) so a caption cannot quote the inflated one by accident.

## Substrate

The DEPLOYED models on their own LOTO packs. NOT rebuilt on final-recipe alignments, deliberately:
this measures the deployed model against its own training data, so the training packs are the
correct reference and re-pooling them would change the question.

## Caveats a caption must carry

- n=98 never-seen transcripts. The CI is wide (0.121-0.177 attn) and the figure shows it.
- LOTO holds out a TISSUE, not transcripts -- ~99.7% transcript overlap is by design, not an error.
- This measures PROFILE CORRELATION, which is the memorisation-sensitive metric. Thresholded ORF
  calling is much less so (cross-species drop-in F1 0.929-0.931), which is why the generalization
  claim rests on the ORF-call results and not on these numbers.

## Regenerate

```bash
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
cd figures/C13_memorisation && $PY make_memorisation.py
```
