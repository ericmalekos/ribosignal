#!/usr/bin/env bash
#SBATCH --job-name=wang_pois
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --array=0-7
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/o2_liver/wang_pois_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/o2_liver/wang_pois_%A_%a.err
set -euo pipefail
# O2 two-arm: Poisson-calibrated arm for the Wang de-novo calls. CPU-only on the gpu partition's idle cores
# (no --gres). Re-runs RiboCode detectORF on the SAVED predicted profiles, Poisson-sampled at depth*theta.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
THETAS=(0.02 0.05 0.1 0.2 0.35 0.5 0.7 1.0)
TH=${THETAS[$SLURM_ARRAY_TASK_ID]}
OUT=$NEW/results/o2_liver/wang_nokozak/dropin_sweep_poisson/theta_${TH}
rm -rf "$OUT"; mkdir -p "$OUT" "$NEW/logs/o2_liver"
echo "wang Poisson theta=$TH -> $OUT"
"$PY" "$NEW/scripts/ribocode_dropin.py" \
  --profiles "$NEW/results/o2_liver/wang_nokozak/pred_profiles.npz" \
  --annot "$ECH/data/ribocode_annot_mouse" --variant pred_preddepth \
  --pred_scale "$TH" --pred_poisson --seed 0 \
  --out "$OUT" --min_aa 5 --pval 0.05
echo "done theta=$TH: $(($(wc -l < "$OUT/pred_preddepth_collapsed.txt")-1)) calls"
