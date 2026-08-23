#!/usr/bin/env bash
#SBATCH --job-name=eval_arch
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=140G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/eval_arch_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/eval_arch_%j.err
set -euo pipefail
# Predict on the Hepatocytes holdout for an arch variant: dump per-nt profiles (GPU, correct SIF) then run
# RiboCode detectORF on real / pred_obsdepth / pred_preddepth (ribocode env). Same as the standard model's eval.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
ORFTRACK=$NEW/data/packed/orf_track_v2_nokozak.npy
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
ARCH=${ARCH:?set ARCH=attn2|attn4|mamba}
SEED=${SEED:-0}
case "$ARCH" in
  attn2) NAME=orf_v2_attn;  SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  attn4) NAME=orf_v2_attn4; SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba) NAME=orf_v2_mamba;  SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img ;;
  mamba4) NAME=orf_v2_mamba4; SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img ;;
  *) echo "bad ARCH=$ARCH" >&2; exit 1 ;;
esac
SUF=""; SEEDLBL=""; [ "$SEED" = 0 ] || { SUF="_seed${SEED}"; SEEDLBL=" seed${SEED}"; }
RUN=results/loto/${NAME}_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes${SUF}
RIBOPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
ANNOT=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/human_ribocode_annot_primary
echo "eval ARCH=$ARCH run=$RUN host=$(hostname) start=$(date -Iseconds)"

# 1. dump predicted profiles on the Hepatocytes holdout (GPU)
# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" --device cuda

# 2. RiboCode drop-in: real (truth) + pred_obsdepth (shape @ real depth) + pred_preddepth (standalone)
export PYTHONNOUSERSITE=1
for V in real pred_obsdepth pred_preddepth; do
  echo "=== dropin $V ==="
  "$RIBOPY" "$NEW/scripts/ribocode_dropin.py" \
    --profiles "$NEW/$RUN/dropin/pred_profiles.npz" --annot "$ANNOT" \
    --variant "$V" --out "$NEW/$RUN/dropin" --min_aa 5 --pval 0.05
done

# --- two-arm Poisson CDS-anchored sweep + full standing-rule ORF-call metrics -> manuscript doc ---
SW=$NEW/$RUN/dropin_sweep_poisson
for TH in 0.02 0.05 0.1 0.2; do
  O=$SW/theta_$TH; rm -rf "$O"; mkdir -p "$O"
  "$RIBOPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$NEW/$RUN/dropin/pred_profiles.npz" \
    --annot "$ANNOT" --variant pred_preddepth --pred_scale "$TH" --pred_poisson --seed 0 \
    --out "$O" --min_aa 5 --pval 0.05
done
CASPY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
mkdir -p "$NEW/results/orf_call_metrics"
"$CASPY" "$NEW/scripts/heldout/orf_call_metrics.py" --run "$NEW/$RUN" --label "onehot ${ARCH}${SEEDLBL}" \
  --universe fibroblast --poisson-sweep "$SW" --cds-recall 0.90 \
  --out-json "$NEW/results/orf_call_metrics/onehot_${ARCH}${SUF}_fibroblast.json" \
  --append-md "$NEW/docs/manuscript_orf_call_eval.md"
echo "eval ARCH=$ARCH done=$(date -Iseconds)"
