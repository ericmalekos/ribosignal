#!/usr/bin/env bash
#SBATCH --job-name=ruiz_loc_rel
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/ruiz_loc_rel_%A_%a.out
set -euo pipefail
# Ruiz-Orera localization on the RELEASED models -- the second half of figure A1.
#
# A1 plots predicted translation-localization against the observed Ribo-seq periodicity ceiling on TWO
# datasets: held-out Hepatocytes and cross-study Ruiz-Orera. The Hepatocytes half is produced by
# eval_localization_released.sbatch. This is the Ruiz-Orera half. Both were previously on
# `orf_v2_attn_*_holdout_Hepatocytes` -- pre-nokozak / pre-mm1 / pre-union -- and the existing
# Ruiz-Orera outputs under results/heldout/human_ruizorera/{onehot,orthrus,rinalmo}/ are all from that
# generation of backends.
#
#   task 0  attn     -> results/heldout/human_ruizorera/released_attn/
#   task 1  mamba4   -> results/heldout/human_ruizorera/released_mamba4/
#
# mamba4 is the one A1 needs (locked decision D1b makes it the MAIN Fig 1 model); attn is run too
# because D1b keeps attn shipped as the CPU inference path and supplemental figures use it.
#
#   sbatch --array=0-1 scripts/eval_ruizorera_localization_released.sbatch
#
# HELD-OUT PACK, NOT THE TRAINING PACK. Unlike eval_localization_released.sbatch, this evaluates a
# held-out dataset, so there is NO RIBO_PACK_SUFFIX: `heldout_pack_dir()` resolves --heldout directly
# to data/packed_heldout_<name> and ignores the suffix. The ORF track and one-hot FASTA are therefore
# the ones built for THAT pack, not the union training ones. Mixing the two is the mistake that
# eval_localization_released.sbatch documents: a same-kozak track from a different universe passes a
# kozak check and silently mis-aligns every ORF.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
DATASET=human_ruizorera
PACK=$NEW/data/packed_heldout_$DATASET
# Built by ruizorera_released_dumps.sbatch, row-aligned to this pack, --kozak none.
ORFTRACK=$PACK/orf_track_v2_nokozak.npy
FASTA=$NEW/data/union_universe.fa                        # chrM-free one-hot input
TXLIST=$NEW/data/heldout_refs/ruizorera_nochrM_tx.txt

MODELS=(attn mamba4)
M=${MODELS[$SLURM_ARRAY_TASK_ID]}
case "$M" in
  attn)   RUN=results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba4) RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img ;;
esac
OUT=$NEW/results/heldout/$DATASET/released_${M}
mkdir -p "$OUT"

for f in "$ORFTRACK" "$FASTA" "$TXLIST" "$NEW/$RUN/best.pt"; do
  [[ -s "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
META=${ORFTRACK%.npy}_meta.json
[[ -s "$META" ]] && { grep -q '"kozak": *"none"' "$META" || \
  { echo "ABORT: $META is not kozak=none" >&2; exit 1; }; }
echo "model=$M dataset=$DATASET run=$RUN host=$(hostname) start=$(date -Iseconds)"

singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_ORF_TRACK="$ORFTRACK" --env RIBO_ONEHOT_FASTA="$FASTA" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_localization.py" --run "$NEW/$RUN" \
    --heldout "$DATASET" --out "$OUT" --min_orf_nt 30 --eval_cap 0 \
    --fasta "$FASTA" --device cuda

RES=$OUT/localization_metrics_${DATASET}.json
[[ -s "$RES" ]] || { echo "ABORT: no $RES written" >&2; exit 1; }
echo "done=$(date -Iseconds) -> $RES"
