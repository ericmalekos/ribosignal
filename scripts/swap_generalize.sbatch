#!/usr/bin/env bash
#SBATCH --job-name=swap_gen
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=8:00:00
#SBATCH --array=0-15
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/swap/gen_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/swap/gen_%A_%a.err
set -euo pipefail
# Generalise the GTF2I input-swap result: same mm1-trained model, mm1 vs mm10 coverage, over
# 5,653 HIGH-multimap test tx (frac_discarded>0.5) + 5,653 LOW-multimap controls (<0.10).
#
# The control group is the point. mm10 changes coverage everywhere (mean RNA multimap rate 3.8%),
# so predictions shifting under mm10 proves nothing on its own; the claim is that the shift is
# CONCENTRATED in high-multimap transcripts. Without the low-multimap arm there is no way to tell
# a targeted correction from a global one.
#
# CPU, not GPU: the attn model is CPU-deployable and the GPU queue is full of the A/B retrains.
# 16 array tasks = 8 shards x 2 coverage postures.
R=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
RUN=$R/results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
G=$R/data/mm25_diagnostic/gtf2i_inference
NSH=8
SH=$(( SLURM_ARRAY_TASK_ID % NSH ))
IDX=$(( SLURM_ARRAY_TASK_ID / NSH ))
SUF=$([ "$IDX" = 0 ] && echo union || echo mm10cov)
O=$G/swap_${SUF}; mkdir -p "$O"
[ -s "$O/pred_profiles_shard${SH}of${NSH}.npz" ] && { echo "SKIP shard $SH $SUF (exists)"; exit 0; }
echo "shard $SH/$NSH coverage=$SUF host=$(hostname) start=$(date -Iseconds)"
singularity exec --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_PACK_SUFFIX=$SUF \
  --env RIBO_ORF_TRACK="$R/data/packed_union/orf_track_v2.npy" \
  --env RIBO_ONEHOT_FASTA="$R/data/union_universe.fa" \
  -B "$BASE" "$SIF" \
  python3 "$R/scripts/dump_pred_profiles.py" --run "$RUN" --out "$O" \
    --tx_list "$G/swap_tx_all.txt" --fasta "$R/data/union_universe.fa" \
    --min_signal 0 --nshards "$NSH" --shard "$SH" --device cpu
echo "done shard $SH $SUF $(date -Iseconds)"
