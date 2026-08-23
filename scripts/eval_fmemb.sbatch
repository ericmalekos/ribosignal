#!/usr/bin/env bash
#SBATCH --job-name=eval_fmemb
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=140G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/eval_fmemb_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/eval_fmemb_%j.err
set -euo pipefail
# Predict on Hepatocytes holdout for an FM-embedding variant: dump per-nt profiles (GPU; the trained model
# reads the pre-computed FM embeddings via tx_index, so the rinalmo SIF/torch works for every backend) then
# RiboCode detectORF on real / pred_obsdepth / pred_preddepth (ribocode env).
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
ORFTRACK=$NEW/data/packed/orf_track_v2_nokozak.npy
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
BACKEND=${BACKEND:?set BACKEND=rinalmo|orthrus|hydrarna}
RUN=results/loto/orf_v2_attn_${BACKEND}_noBrain_nokozak_mm1_holdout_Hepatocytes
RIBOPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
ANNOT=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/human_ribocode_annot_primary
mkdir -p "$NEW/logs/dropin"
echo "eval BACKEND=$BACKEND run=$RUN host=$(hostname) start=$(date -Iseconds)"

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True --env RIBO_ORF_TRACK="$ORFTRACK" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" --device cuda

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
"$CASPY" "$NEW/scripts/heldout/orf_call_metrics.py" --run "$NEW/$RUN" --label "$BACKEND" \
  --universe fibroblast --poisson-sweep "$SW" --cds-recall 0.90 \
  --out-json "$NEW/results/orf_call_metrics/${BACKEND}_fibroblast.json" \
  --append-md "$NEW/docs/manuscript_orf_call_eval.md"
echo "eval BACKEND=$BACKEND done=$(date -Iseconds)"
