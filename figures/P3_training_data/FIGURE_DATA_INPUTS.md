# P3 -- training composition and model inputs: data inputs

Poster panel. A: per-cell-type Ribo-seq depth with the LOTO holdout marked. B: what the model
receives at inference.

## Panel -> data

| element | source | produced by |
|---|---|---|
| P-sites per cell type | `data/packed_canon{,_<T>}/target_counts.npy` | `scripts/build_canonical_chothani_packs.sbatch` |
| Ribo library count per cell type | `data/packed_canon{,_<T>}/coverage_norm.json` -> `n_ribo_samples` | same |
| universe size | `data/packed_canon/tx_order.txt` | same |
| input channel description (B) | fixed by the architecture, not read from data | `scripts/model.py`, `scripts/build_orf_track.py` |

## Counts

7 training cell types, **69** Ribo-seq libraries. Hepatocytes (**5** libraries) is held out entirely.
74 total. The suptitle states training and holdout separately on purpose -- "74 libraries" reads as
if all 74 were trained on.

## Caveats a caption must carry

- **Depth is very uneven**: Fibroblast 32 libraries / 1,070M P-sites vs HUVEC 3 / 94M, an 11x spread.
  This is why training caps at `--max_tx_per_tissue 6000` and why LOTO holds out a TISSUE rather than
  random transcripts.
- **Brain is not shown** because no pack exists for it. It was dropped from the source study's
  9-tissue panel for low periodicity (period_obs 0.044) before packs were built.
- Panel A reads the CANONICAL packs, so its P-site totals are 2.1-3.0% below the deployed models'
  training packs (`packed_union*`). The deployed models were trained on the latter. The composition,
  ordering and depth spread are unaffected; only the absolute totals differ.

## Regenerate

```bash
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
cd figures/P3_training_data && $PY make_training_data.py
```
