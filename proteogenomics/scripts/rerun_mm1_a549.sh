#!/usr/bin/env bash
# Strict-uniform mm1 re-run of A549 (tryptic, fresh universe), fully SLURM-native (afterok chain).
# A549's canonical RNA input is the ENCODE FASTQ (ENCSR000CON: ENCFF000EJJ/EJW/EJV/EKB); the "fresh"
# lineage packs that coverage onto A549's own salmon universe (A549_universe.fa). Chain:
#   align_a549_encode(mm1) -> pack(human_a549_fresh) -> dump_a549_fresh(GPU) ->
#   post_predict_a549 (enumerate -> db_fresh -> tryptic search -> MS2Rescore -> class-aware compare).
# Result: db_comparison_rescored_classaware_mm1.md (mm20 preserved as _mm20.md).
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
S=$NEW/proteogenomics/scripts
RNA=$NEW/proteogenomics/data/A549_pilot/encode_rna
for f in ENCFF000EJJ ENCFF000EJW ENCFF000EJV ENCFF000EKB; do
  [ -s "$RNA/$f.fastq.gz" ] || { echo "missing A549 ENCODE FASTQ $f (re-download first)" >&2; exit 1; }
done

AJ=$(sbatch --parsable $S/align_a549_encode.sbatch)
PKJ=$(sbatch --parsable --dependency=afterok:$AJ --export=ALL,DS=human_a549_fresh $NEW/scripts/heldout/build_heldout_pack.sbatch)
DJ=$(sbatch --parsable --dependency=afterok:$PKJ $S/dump_a549_fresh.sbatch)
XJ=$(sbatch --parsable --dependency=afterok:$DJ $S/post_predict_a549.sbatch)
echo "A549 (tryptic, fresh) mm1 chain:  align=$AJ -> pack=$PKJ -> predict=$DJ -> postpred=$XJ"
