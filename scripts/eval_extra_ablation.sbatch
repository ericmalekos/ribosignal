#!/usr/bin/env bash
#SBATCH --job-name=eval_extra_abl
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval_extra_abl_%A_%a.out
set -euo pipefail
# The input ablation's test_metrics.json reports PROFILE Pearson only, which is shape and is
# scale-invariant per transcript. That measures the wrong thing for the cell-type-specificity claim:
# RNA-seq coverage feeds the COUNT HEAD, i.e. per-transcript MAGNITUDE -- which transcripts are
# expressed and how much. eval_extra adds count_pearson, so this is the arm of the ablation that
# actually tests the claim.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
MODES=(emb cov)
M=${MODES[$SLURM_ARRAY_TASK_ID]}
RUN=$NEW/results/ablation/orf_v2_attn_onehot_union_nokozak_input_${M}_holdout_Hepatocytes
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_PACK_SUFFIX=union \
  --env RIBO_ORF_TRACK="$NEW/data/packed_union/orf_track_v2.npy" \
  --env RIBO_ONEHOT_FASTA="$NEW/data/union_universe.fa" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$RUN" --eval_cap 0 --device cuda
echo "done $M=$(date -Iseconds)"
