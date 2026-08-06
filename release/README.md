# Released models

The canonical answer to "which of the 201 checkpoints under `results/` is the model?"

Two ship. They differ only in the sequence mixer; the body, heads, inputs, training data, and ORF
track are identical.

| | `orf_v2_mamba4` | `orf_v2_attn` |
|---|---|---|
| role | **primary / headline** | supplemental, maintained CPU release |
| mixer | dilated CNN + 4 Bi-Mamba blocks | dilated CNN + 2 transformer layers |
| params | 7,521,026 | 5,071,106 |
| device | **GPU only** (mamba-ssm CUDA kernels) | **CPU or GPU** |
| held-out Pearson, 3 seeds | **0.6799** (0.6753-0.6851) | 0.6595 (0.6585-0.6603) |
| seed-0 Pearson (the checkpoint here) | 0.6851 | 0.6585 |
| checkpoint | 29 MB | 20 MB |

**Quote the 3-seed mean and range, not the seed-0 number.** mamba4's seed spread is 0.0098 (~5x
attn's 0.0018), so a single-seed mamba4 figure overstates its precision. The checkpoints shipped here
are seed 0 of each.

Pick by device. `--device cpu` against the mamba4 run dir fails inside the mamba-ssm CUDA kernels;
that is a hard constraint of the dependency, not a configuration problem, and it is the entire reason
attn stays released.

## Where the weights actually are

Checkpoints are NOT in git (`*.pt` is ignored; 46 GB of run output lives outside the repo). Canonical
paths on the group filesystem:

```
results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/best.pt
results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/best.pt
```

`SHA256SUMS` in each subdirectory pins the exact file. Verify with
`sha256sum -c release/orf_v2_<mixer>/SHA256SUMS` from inside the run directory. For external release
the weights need a separate host (Zenodo / HuggingFace) -- they cannot go in the repo.

`args.json` and `test_metrics.json` are copied here verbatim so the config and headline metrics are
version-controlled even though the weights are not.

## Companion artifacts required at inference

A checkpoint alone is not enough. All three must match what the model trained on:

| env var | path | note |
|---|---|---|
| `RIBO_PACK_DIR` | `data/packed_union/` | union universe, 84,472 tx |
| `RIBO_ORF_TRACK` | `data/packed_union/orf_track_v2_nokozak.npy` | **see the warning below** |
| `RIBO_ONEHOT_FASTA` | `data/union_universe.fa` | |

### ORF-track naming hazard -- read this

`data/packed_union/orf_track_v2.npy` carries the BARE name but contains `--kozak none` content. The
pack predates the 2026-07-31 naming fix, under which bare means *heuristic Kozak*. Its
`orf_track_v2_meta.json` is authoritative and says `"kozak": "none"`.

Use `orf_track_v2_nokozak.npy` (a symlink to the same bytes, added for exactly this reason) so the
name states the content. Misreading the bare name as heuristic is the mistake that previously fed a
nokozak-trained model the Kozak track. Full detail in `data/packed_union/NAMING_WARNING.md`.

The released models are `--kozak none` models. That is the correct inference setting for both.

## Training configuration

Identical across the two except `mixer`:

- Held-out tissue: Hepatocytes. Trained on Fibroblast, VSMC, ES, Fat, HA_EC, HCAEC, HUVEC (**Brain
  dropped**). 42,000 train / 8,400 val / 70,883 test transcripts.
- One-hot encoding. FM embeddings give no lift (Task 15), so there is **no foundation-model
  dependency** in either release.
- Inputs: sequence + RNA-seq coverage (`input_mode: both`), mm1 coverage, `cov_norm: global_mean`.
- 256 channels, 10 blocks, dropout 0.1, budget 16,000 nt, `count_weight` 0.1, `peakiness_weight` 0.0
  (the anti-smoothing loss is off -- O3 track 2 was a negative result).
- 30 epochs, lr 3e-4, weight decay 0.01, warmup 2, patience 6.

## Metrics caveat

`test_metrics.json` for attn predates the addition of `frame0_pred_median`, so that field is present
for mamba4 (0.8196) and absent for attn. Not a failure -- an older eval schema. Re-run
`eval_localization.py` on the attn run dir if the field is needed for a figure.

Measured on the same 70,883 held-out transcripts:

| | Pearson | period_pred | period_obs | frame0_pred |
|---|--:|--:|--:|--:|
| mamba4 | 0.6851 | 0.2632 | 0.1457 | 0.8196 |
| attn | 0.6585 | 0.3145 | 0.1457 | (not in schema) |

## Provenance warning for downstream results

Not every result in this project was produced with these models. As of 2026-08-04 the macrophage
proteogenomics (12 populations) uses both; the A549 and immunopeptidome results still use the OLD
`orf_v2_attn_onehot_holdout_Hepatocytes` with the heuristic Kozak track. Check the `RUN=` line of the
relevant `dump_*.sbatch` before attributing any downstream number to a released model.
