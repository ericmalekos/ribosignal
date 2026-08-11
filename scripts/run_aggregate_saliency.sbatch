#!/usr/bin/env bash
#SBATCH --job-name=agg_sal
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=11:45:00   # medium partition caps at 12 h; gpu allows more, override per-submit if needed
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/agg_sal_%A_%a.out
set -euo pipefail
# Phase 3.1: aggregate sequence saliency, 2 substrates x 2 models.
#
#   task 0  attn   x mouse GSE243134   (PRIMARY -- transcripts never seen in training)
#   task 1  attn   x human Hepatocytes (CONTAMINATED arm -- 99.7% of test tx seen, kept on purpose)
#   task 2  mamba4 x mouse GSE243134
#   task 3  mamba4 x human Hepatocytes
#
# ARCHITECTURE SPLIT. `attn` runs on CPU. `mamba4` CANNOT: mamba_ssm has no CPU build
# (ModuleNotFoundError in the cas12a env, and its kernels are CUDA-only anyway), so tasks 2-3 must go
# to the gpu partition inside the orthrus SIF. Same reason mamba4 cannot do CPU dumps.
#   CPU:  sbatch --partition=medium --array=0-1 scripts/run_aggregate_saliency.sbatch
#   GPU:  sbatch --partition=gpu --gres=gpu:A5500:1 --array=2-3 scripts/run_aggregate_saliency.sbatch
#
# Substrate pairing is explicit rather than inferred: each pack goes with the ORF track and one-hot
# FASTA that DEFINE its universe, and with a collapsed call table from the same universe. Crossing
# those would attribute saliency at coordinates that do not correspond to the ORFs being named.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
OUT=$NEW/results/aggregate_saliency
mkdir -p "$OUT"

MODELS=(attn attn mamba4 mamba4)
SUBS=(mouse human mouse human)
M=${MODELS[$SLURM_ARRAY_TASK_ID]}
S=${SUBS[$SLURM_ARRAY_TASK_ID]}

case "$M" in
  attn)   RUN=$NEW/results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba4) RUN=$NEW/results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
          SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img ;;
esac

case "$S" in
  mouse)  PACK=$NEW/data/packed_heldout_l3x3_ribo-gse243134_rna-gse243134
          TRACK=$NEW/data/liver3x3/orf_track_union3_nokozak.npy
          FA=$NEW/data/heldout_refs/mouse_liver_union3_universe.fa
          COLL=$NEW/results/mouse_liver_3x3/mamba4_ribo-gse243134_rna-gse243134/real_collapsed.txt ;;
  human)  PACK=$NEW/data/packed_union_Hepatocytes
          TRACK=$NEW/data/packed_union/orf_track_v2_nokozak.npy
          FA=$NEW/data/union_universe.fa
          COLL=$NEW/results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin/real_collapsed.txt ;;
esac

for f in "$RUN/best.pt" "$PACK/tx_order.txt" "$TRACK" "$FA" "$COLL"; do
  [[ -e "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
# Name must match content (feedback_kozak_never_set).
META=${TRACK%.npy}_meta.json
[[ -s "$META" ]] && { grep -q '"kozak": *"none"' "$META" || { echo "ABORT: $META not kozak=none" >&2; exit 1; }; }

LAB="${M}_${S}"
echo "task=$SLURM_ARRAY_TASK_ID model=$M substrate=$S host=$(hostname) start=$(date -Iseconds)"
ARGS=(--run "$RUN" --pack "$PACK" --orf-track "$TRACK" --fasta "$FA" --collapsed "$COLL"
      --label "$LAB" --out "$OUT/${LAB}.tsv" --per-class 300)
if [[ "$M" == "mamba4" ]]; then ARGS+=(--device cuda); fi

if [[ "$M" == "attn" ]]; then
  "$PY" "$NEW/scripts/aggregate_saliency.py" "${ARGS[@]}"
else
  singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 -B "$BASE" "$SIF" \
    python3 "$NEW/scripts/aggregate_saliency.py" "${ARGS[@]}"
fi
echo "done=$(date -Iseconds) -> $OUT/${LAB}.tsv"
