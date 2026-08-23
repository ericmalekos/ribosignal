#!/usr/bin/env bash
#SBATCH --job-name=ruiz_dropin
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=110G
#SBATCH --time=8:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/ruiz_dropin_%A_%a.out
set -euo pipefail
# RiboCode drop-in on the RELEASED-model Ruiz-Orera dumps (36,655 tx, chrM excluded).
# 2 models x 3 density variants. `real` is the reference (RiboCode on the pack's real counts);
# pred_obsdepth isolates profile SHAPE; pred_preddepth is the standalone no-Ribo-seq arm.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
RPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python
MODELS=(attn attn attn mamba4 mamba4 mamba4)
VARS=(real pred_obsdepth pred_preddepth real pred_obsdepth pred_preddepth)
M=${MODELS[$SLURM_ARRAY_TASK_ID]}; V=${VARS[$SLURM_ARRAY_TASK_ID]}
D=$NEW/results/heldout/human_ruizorera/released_${M}/dropin
echo "model=$M variant=$V host=$(hostname) start=$(date -Iseconds)"
"$RPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$D/pred_profiles.npz" \
  --annot "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/human_ribocode_annot_primary" --variant "$V" --out "$D" --min_aa 5 --pval 0.05
echo "done=$(date -Iseconds): $(($(wc -l < "$D/${V}_collapsed.txt") - 1)) calls"
