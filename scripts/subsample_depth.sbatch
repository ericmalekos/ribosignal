#!/usr/bin/env bash
#SBATCH --job-name=ribo_subsamp
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --array=0-5
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/subsamp_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/subsamp_%A_%a.err
set -euo pipefail
# Depth-crossover experiment: binomially thin the REAL Hepatocytes P-site profile to a fraction of
# depth, run RiboCode detectORF on the thinned counts, so the "measure at depth f" ORF-call F1 can be
# compared against the depth-independent "predict from RNA-seq (no Ribo-seq)" line. One array task per
# (depth, seed); override DEPTHS/SEEDS + --array to change the grid.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
RUN=${RUN:-results/loto/orf_v2_attn_onehot_holdout_Hepatocytes}
PROFILES=$NEW/$RUN/dropin/pred_profiles.npz
ANNOT=$ECH/data/ribocode_annot
OUT=$NEW/$RUN/dropin/subsample
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
mkdir -p "$NEW/logs/dropin" "$OUT"

# depth grid x seeds; task = depth_idx * n_seeds + seed_idx
IFS=' ' read -r -a DEPTHS <<< "${DEPTHS:-0.02 0.01 0.005 0.002 0.001 0.0005}"
IFS=' ' read -r -a SEEDS <<< "${SEEDS:-0}"
NS=${#SEEDS[@]}
di=$(( SLURM_ARRAY_TASK_ID / NS )); si=$(( SLURM_ARRAY_TASK_ID % NS ))
F=${DEPTHS[$di]}; S=${SEEDS[$si]}
echo "task=$SLURM_ARRAY_TASK_ID depth=$F seed=$S host=$(hostname) start=$(date -Iseconds)"
export PYTHONNOUSERSITE=1
"$PY" "$NEW/scripts/ribocode_dropin.py" \
  --profiles "$PROFILES" --annot "$ANNOT" --variant real \
  --subsample "$F" --seed "$S" --out "$OUT" --min_aa 5 --pval 0.05
echo "done depth=$F seed=$S $(date -Iseconds)"
