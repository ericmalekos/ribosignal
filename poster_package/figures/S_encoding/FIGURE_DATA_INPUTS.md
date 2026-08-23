# SUPPLEMENTAL: encoding -- one-hot vs RNA foundation-model embeddings

**What it shows:** the simplest, dependency-free sequence encoding matches 650M-parameter RNA language
models on this task. That is the justification for the released models carrying no foundation-model
dependency and running on CPU. Framed as accessibility, not as a horse race.

**Generator:** `make_encoding_panel.py` (cas12a). Reads `test_metrics.json` from four runs; writes
`S_encoding.tsv` with the plotted values.

## No retraining was needed

Task 41 already trained all four encodings on the **identical** recipe, holdout and transcript set,
differing only in `--emb_backend`. The panel is a read of existing results.

| encoding | run | pc profile r | lncRNA profile r | predicted 3-nt periodicity |
|---|---|--:|--:|--:|
| **one-hot** (shipped) | `orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes` | **0.6250** | 0.4275 | 0.4427 |
| RiNALMo (650M) | `..._rinalmo_...` | 0.6089 | **0.4499** | **0.4564** |
| Orthrus (Mamba) | `..._orthrus_...` | 0.6056 | 0.4234 | 0.4481 |
| HydraRNA (full-length) | `..._hydrarna_...` | 0.6217 | 0.4348 | 0.4348 |

n = 33,918 scored transcripts in all four arms; the generator asserts this and exits if the arms ever
diverge, since arms scored over different transcript sets are not comparable.

## Reading it honestly

- **On protein-coding profile shape, one-hot is best**, by 0.003 to 0.020. The FM embeddings do not
  help, and two of the three are measurably worse.
- **RiNALMo does win two of the three metrics**: lncRNA profile Pearson (+0.022) and predicted
  periodicity (+0.014). The panel should not claim one-hot wins everywhere, because it does not. The
  claim is that the gap never justifies a 650M-parameter dependency.
- Effect sizes are small in absolute terms. This is a "no meaningful difference" result, which is
  exactly what makes the dependency-free choice defensible.

## Recipe caveat -- state it on the panel

These are `noBrain_nokozak_mm1` runs on the **Fibroblast universe**, not the union universe the
released models use. Union FM embeddings were never extracted, so a union version of this panel would
mean three more 12-17 h trainings. The four arms share one recipe, so the ENCODING comparison is
internally clean; only the absolute Pearson values are recipe-specific and must not be quoted beside
released-model numbers (released one-hot pc is 0.6699, mamba4 0.6799).

## Regenerate

```bash
cd /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/figures/S_encoding
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 make_encoding_panel.py
```
