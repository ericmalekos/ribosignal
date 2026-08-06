# Leave-one-tissue-out (LOTO) plan

Objective: test whether the per-nt Ribo-seq signal model, trained on several Chothani tissues,
predicts the P-site profile of a held-out tissue it never saw. This is the project's real
generalization objective (cross-tissue, then cross-dataset). It differs from the existing
held-out-chromosome test: chromosome folds test gene transfer within one tissue; LOTO tests
tissue transfer (gene overlap across tissues is expected and realistic).

## Tissues (9 Chothani, GENCODE v49 / GRCh38)

| Tissue      | Ribo-seq samples (target) | RNA-seq samples (coverage input) |
|-------------|---------------------------|----------------------------------|
| Fibroblast  | 32 | 32 |
| VSMC        | 11 | 5  |
| ES          | 6  | 2  |
| Fat         | 6  | 2  |
| HA_EC       | 6  | 6  |
| Brain       | 5  | 2  |
| HCAEC       | 5  | 3  |
| Hepatocytes | 5  | 2  |
| HUVEC       | 3  | 3  |

## Design decisions

1. **Fixed universe (36,668 tx), embeddings + ORF track reused.** The FM per-token embeddings
   (RiNALMo, Orthrus) and the ORF-candidate track are sequence-based and tissue-independent.
   Per-tissue expression enters through the RNA-seq coverage input, NOT by changing the
   transcript set. So the same 36,668-tx Fibroblast-expressed universe is the fixed feature set
   for every tissue; zero embedding re-extraction. A tissue-silent gene simply has ~0 coverage
   and ~0 target in that tissue and is filtered by the scorable threshold.
2. **Per-tissue packs on the shared axis.** `data/packed_<Tissue>/` holds that tissue's
   `target_counts.npy` + `coverage.npy` + per-tissue `coverage_norm.json`; `tx_order`/`offsets`/
   `lengths`/`sum_L` are identical across tissues (same universe), so the shared
   `data/packed/orf_track*.npy` and embeddings align to every pack.
3. **Depth normalization per tissue** (`--cov_norm global_mean`): each tissue divides its coverage
   by its OWN global-mean per-nt depth before log1p, so a shallow-RNA-seq tissue (e.g. Hepatocytes,
   2 samples) is comparable to a deep one. Depth-norm is exactly the transferability fix for LOTO.
4. **Mitochondrial exclusion** (project rule, methods sec 7.1): chrM dropped from every split via
   `dataset.excluded_tx()`.
5. **Tissue split, not chromosome split.** Train = held-in tissues (ConcatDataset of their per-tissue
   RiboDatasets), eval = held-out tissue. Validation for early stopping = a held-out fraction of the
   training tissues (stratified). No gene-disjointness required (LOTO holds out the tissue axis); a
   stricter tissue-AND-chromosome hold-out is a possible follow-up.
6. **Scorable / training-tx filter (per tissue):** train on universe tx with >= threshold pooled
   P-sites in that tissue (learn shape where shape exists); eval with the same profile-Pearson
   scorable windows used for Fibroblast (pc/lncRNA whole-tx, uORF/dORF UTR windows, >= 20 window
   P-sites). Threshold TBD at pack time from the per-tissue signal distribution.
7. **Configs:** baseline + orf_v2_attn (the established best all-around config from the Fibroblast
   matrix). RiNALMo embeddings first; Orthrus as a follow-up if warranted.

## Pipeline + SLURM chain (input builds are hold-out-independent -- all 9 tissues built)

- `build_psite_target.py <Tissue>` -> `build_target_loto.sbatch` (array 0-7)      job 35179263
- `align_rnaseq_sample.sh` per SRR -> `align_rnaseq_loto.sbatch` (array 1-25%10)   job 35179264
- `pool_rnaseq_coverage_tissue.py <Tissue>` -> `pool_coverage_loto.sbatch` (afterok align)  job 35179286
- `pack_target_coverage_tissue.py <Tissue>` -> `pack_loto.sbatch` (afterok target+pool)      job 35179287
  => produces `data/packed_<Tissue>/` for all 8 non-Fibroblast tissues (Fibroblast pack exists).

## Open decisions (confirm before the expensive multi-tissue training)

- **Held-out tissue for the first LOTO run.** Recommend Hepatocytes (biologically distinct from
  the endothelial/fibroblast/VSMC/brain/ES training mix, 5 Ribo-seq samples for a usable target;
  shallow 2-sample RNA-seq matters less since profile shape comes from the tissue-independent FM
  embeddings and coverage is depth-normalized). Alternative: VSMC (11 Ribo-seq / 5 RNA-seq = the
  deepest, least-noisy test target, but less distinct from the vascular training tissues).
- **Scope: single held-out (pilot) vs full 9-fold LOTO.** Start with one held-out tissue to validate
  the whole multi-tissue pipeline + get the first cross-tissue number; extend to the full 9-fold
  (each tissue held out once) if the pilot is sound. Training cost scales ~8x the single-tissue run
  per fold.

## Not yet built (pending the decisions above)

- Multi-tissue LOTO trainer (`train_loto.py` or a `--loto_holdout` mode on train.py): builds one
  PackedStore per tissue, ConcatDataset over held-in tissues, evaluates on the held-out tissue with
  `eval_extra.py`'s pc/lncRNA/uORF/dORF metrics. Reuses the existing model, loss, and eval.
- Within-held-out-tissue reference baseline (train on the held-out tissue's own genes, gene-disjoint)
  to contextualize the cross-tissue number against the tissue's own ceiling.
