# Plan: GSE120762 (mouse +/- LPS) held-out eval + Chothani mm1 retrain (2026-07-22)

Two parallel streams requested by the user. Downloads are SERIAL on the head node (ENA per-IP rate limit);
compute fans out on SLURM.

## Stream A -- GSE120762 mouse macrophage +/- LPS (primary interest)
GEO GSE120762 = SRP163147 / PRJNA494404 (Jackson 2018, "Translation of Non-Canonical ORFs Controls Mucosal
Immunity"). BMDM +/- LPS. Downloading 10 samples (~25 GB): 6 RNA-seq (WT/NT rep1-3 + LPS rep1-3) + 4
ribosome-profiling (CHX; NT rep1-2 + LPS rep1-2). Skipping the 6 RiboTag (different modality).

Goal: "check how the two conditions work in the model" -- predict per-nt Ribo-seq P-site signal for NT and
for LPS SEPARATELY (each fed its own condition RNA-seq coverage + the shared sequence), compare to the
observed CHX Ribo-seq footprints, and ask whether the model captures the LPS-induced translational program.

Pipeline (mouse GRCm39 / GENCODE vM38, reusing the Wang cross-species mouse infra):
1. RNA-seq: STAR `--outFilterMultimapNmax 1` (mm1) -> transcriptome -> rnaseq_coverage.py, pooled per
   condition (NT = WTrep1-3, LPS = LPSrep1-3) + salmon TPM for the expressed universe.
2. Ribo-seq: adapter/UMI check -> STAR single-mapper + ncRNA-locus drop -> RiboCode P-sites, per condition
   (NT = CHX_NT rep1-2, LPS = CHX_LPS rep1-2).
3. Mouse universe (expressed, per condition or shared) + build_heldout_pack entries `mouse_gse120762_nt` /
   `_lps` (coverage_only or with the Ribo target for the observed-ceiling comparison).
4. Predict with the deployed one-hot orf_v2_attn (shared RIBO_ORF_TRACK mouse track from Wang) -> per-nt
   profiles; eval per-tx Pearson + localization vs observed, and NT-vs-LPS differential (does predicted
   translation change track observed LPS induction?).

Status: DOWNLOADING (bqsozw0jn, ~1h at 7 MB/s, md5-verified).

## Stream B -- Chothani mm1 retrain (consistency + drop Brain)
The deployed model trained on 8 tissues (all except Hepatocytes); user wants (a) mm1 coverage everywhere and
(b) Brain dropped (noise; period_obs 0.044, the A3 low outlier) -> retrain on **7 tissues**: Fibroblast,
VSMC, ES, Fat, HA_EC, HCAEC, HUVEC. The Ribo-seq TARGET is unchanged by mm1; only the RNA-seq COVERAGE input
changes, but the training FASTQ were deleted -> must re-download + re-align mm1.

Samplesheet: 55 RNA-seq samples (Fibroblast 32 + VSMC 5 + HA_EC 6 + HCAEC 3 + HUVEC 3 + ES 2 + Fat 2 =
53 train; + Hepatocytes 2 for the held-out eval coverage). Brain excluded entirely. ENA PRJNA756023, all 55
matched (fastq_ftp + md5).

Pipeline:
1. Re-download 55 RNA-seq (~150 GB, ~6h; STAGED, launches after Stream A to keep ENA serial).
2. Re-align mm1 (align_rnaseq_sample.sh, now `--outFilterMultimapNmax 1`) -> per-sample coverage hd5 (SLURM array).
3. Pool per-tissue (pool_rnaseq_coverage_tissue.py) -> rebuild packed_<Tissue>/coverage.npy with mm1
   (pack_target_coverage_tissue.py; target_counts + geometry unchanged).
4. Retrain: `train_loto.py --holdout Hepatocytes --train_tissues Fibroblast,VSMC,ES,Fat,HA_EC,HCAEC,HUVEC
   --emb_backend onehot --cov_norm global_mean --use_orf_track --n_attn_layers 2` with
   RIBO_ORF_TRACK=data/packed/orf_track_v2.npy -> results/loto/orf_v2_attn_onehot_mm1_noBrain_holdout_Hepatocytes.
5. Eval on held-out Hepatocytes (mm1 coverage) + compare to the mm20 deployed model (expect ~equal per the
   proxy; the real change is dropping Brain).

Status: STAGED (download list + dl_chothani_mm1.sh ready; launches on Stream A completion).

## Notes / decisions
- mm1 is now the standard coverage posture (already applied to all align scripts). This retrain extends it to
  the training data for full provenance consistency; the proxy (Task 25) predicts the model itself barely
  changes, so the SUBSTANTIVE change is dropping the noisy Brain tissue.
- The deployed model is NOT discarded; the mm1+noBrain model is a parallel, cleaner version for comparison.
- Serial downloads (A then B) per the ENA rate-limit rule; aligns/retrain fan out on SLURM.
