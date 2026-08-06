#!/usr/bin/env bash
#SBATCH --job-name=dump_gse_nk
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G
#SBATCH --time=02:00:00
#SBATCH --array=0-1
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval/dump_gse_nk_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval/dump_gse_nk_%A_%a.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
G=$NEW/proteogenomics/data/gse120762
RUN=$NEW/results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes   # mm1 no-kozak model
CONDS=(nt lps); COND=${CONDS[$SLURM_ARRAY_TASK_ID]}
PACK=$NEW/data/packed_heldout_mouse_gse120762_${COND}
# predict with the mm1 no-kozak model + the NO-KOZAK GSE track; output to a *_nokozak dir
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$PACK/orf_track_v2_nokozak.npy" \
  --env RIBO_ONEHOT_FASTA="$G/gse120762_${COND}_universe.fa" \
  -B /private/groups/carpenterlab/emalekos/RNAZoo_meta "/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img" \
  python3 $NEW/scripts/dump_pred_profiles.py --run "$RUN" \
    --heldout mouse_gse120762_${COND} --tx_list "$G/gse120762_${COND}_universe_tx.txt" \
    --fasta "$G/gse120762_${COND}_universe.fa" --out "$G/pred_${COND}_nokozak" --device cuda
echo "dump $COND (no-kozak model) done -> $G/pred_${COND}_nokozak/pred_profiles.npz"
