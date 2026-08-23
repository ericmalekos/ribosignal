# FIGURE B6 -- input-modality ablation: data inputs

**What it shows:** what each input modality contributes, measured by retraining the deployed recipe
with one input zeroed. Two panels: profile SHAPE (left), per-transcript MAGNITUDE (right).

**Why it exists:** it is the direct test of the claim the proteogenomics argument rests on. If the
model reached the same accuracy with the RNA-seq coverage channel zeroed, it would be a sequence-only
ORF prior emitting the same output for every cell type, and every cell-type-specific discovery claim
would collapse. This measures that rather than assuming it -- and the answer turned out to be
narrower than the claim originally assumed. See "Reading it" below.

**Generator:** `make_input_ablation.py` (cas12a env). Every value is read from a result JSON at build
time; the run directory paths are written into `B6_values.json`.

| arm | run directory | source files |
|---|---|---|
| both (deployed) | `results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes` | `extra_metrics.json`, `test_metrics.json` |
| sequence only | `results/ablation/orf_v2_attn_onehot_union_nokozak_input_emb_holdout_Hepatocytes` | same |
| RNA-seq only | `results/ablation/orf_v2_attn_onehot_union_nokozak_input_cov_holdout_Hepatocytes` | same |

Fields used: `protein_coding.profile_pearson_median`, `lncRNA.profile_pearson_median`,
`protein_coding.count_pearson` (from `extra_metrics.json`), `protein_coding.period_pred_median` and
`n` (from `test_metrics.json`). The generator asserts all three arms scored the same `n` and exits if
they did not, since arms scored over different transcript sets are not comparable.

## Panel values (held-out Hepatocytes, n = 70,883)

| input_mode | pc profile r | lncRNA profile r | pc count r | pc periodicity (pred) |
|---|--:|--:|--:|--:|
| both (deployed) | 0.6699 | 0.4583 | 0.8990 | 0.3359 |
| sequence only (RNA-seq zeroed) | 0.6601 | 0.4462 | 0.7063 | 0.3704 |
| RNA-seq only (sequence zeroed) | 0.1226 | 0.0838 | 0.8078 | -0.0155 |

Derived and drawn on the figure: sequence-only keeps **98.5%** of the deployed profile Pearson; the
count head loses **-0.193** Pearson without RNA-seq.

## How the arms were produced

`scripts/train_union_inputablation.sbatch --array=0-1`. It clones `train_loto_union.sbatch` line for
line and changes only `--input_mode`, including keeping `RIBO_PACK_SUFFIX=union`, the union ORF track,
and the union one-hot FASTA. The `both` arm is the deployed run itself rather than a re-train, so the
contrast carries no extra seed noise.

The older `scripts/train_loto_ablation.sbatch` also accepts `INPUT_MODE`, but it is pinned to the
pre-union Fibroblast config (`data/packed/orf_track_v2.npy`, no `--kozak none`, no
`RIBO_PACK_SUFFIX`). Ablating there contrasts against a model that is not the shipped one, so it
cannot speak to the deployed claim. That is the whole reason a second script exists.

## Reading it

- **SHAPE is sequence-borne.** Zeroing RNA-seq costs 0.0098 of pc profile Pearson. Periodicity is not
  merely preserved, it sharpens (0.3704 vs 0.3359): the coverage channel adds a smooth envelope that
  slightly dilutes fine periodicity.
- **MAGNITUDE is RNA-seq-borne.** Count Pearson 0.8990 -> 0.7063 without it, and RNA-seq alone (0.8078)
  beats sequence alone at per-transcript depth.
- **RNA-seq alone cannot do shape at all** (0.1226, periodicity -0.0155).
- **Consequence for the manuscript.** Cell-type specificity enters through the COUNT HEAD -- a
  transcript not expressed in the query cell type gets low predicted depth, fails the caller's
  significance test, and is never called. It does NOT enter through the profile shape, which is nearly
  cell-type-invariant. Write "predicts cell-type-specific translation" and cite -0.193; do not write
  "predicts a cell-type-specific profile shape".

## Regenerate

```bash
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/figures/B6_input_ablation
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 make_input_ablation.py
```
