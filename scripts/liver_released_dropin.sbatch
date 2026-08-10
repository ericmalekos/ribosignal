#!/usr/bin/env bash
#SBATCH --job-name=liver_dropin
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=110G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/liver_dropin_%A_%a.out
set -euo pipefail
# RiboCode drop-in on the six RELEASED-model mouse-liver dumps (results/liver_released/).
#
# Replaces every model row in the mouse-liver analysis. The previous ones came from
# `orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes` -- the pre-union Fibroblast-universe
# checkpoint -- with no mamba4 arm at all, despite mamba4 being the primary release.
#
# 3 datasets x 2 models x 3 density variants = 18 tasks:
#   real            RiboCode on the pack's REAL counts -- the per-dataset reference
#   pred_obsdepth   predicted SHAPE at real per-tx depth -- isolates shape from depth
#   pred_preddepth  predicted shape AND predicted depth -- the standalone, no-Ribo-seq arm (theta=1)
#
# The calibrated (Poisson) arm is a separate job: theta must come from the Ribo-seq-FREE pgx anchor
# (annotated-CDS relative recall), not from matching an observed between-experiment ceiling, or the
# calibrated numbers would not be reachable on a tissue without Ribo-seq.
#
#   sbatch --array=0-17%4 scripts/liver_released_dropin.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
RPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python
ANNOT="/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/mouse_ribocode_annot"

DUMPS=(attn_mouse_wang_liver mamba4_mouse_wang_liver
       attn_mouse_janich_liver_decon mamba4_mouse_janich_liver_decon
       attn_mouse_gse243134_liver mamba4_mouse_gse243134_liver)
VARS=(real pred_obsdepth pred_preddepth)
D_I=$(( SLURM_ARRAY_TASK_ID / 3 )); V_I=$(( SLURM_ARRAY_TASK_ID % 3 ))
DS=${DUMPS[$D_I]}; V=${VARS[$V_I]}
DIR=$NEW/results/liver_released/$DS
[[ -s "$DIR/pred_profiles.npz" ]] || { echo "missing $DIR/pred_profiles.npz" >&2; exit 1; }
echo "dump=$DS variant=$V host=$(hostname) start=$(date -Iseconds)"
"$RPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$DIR/pred_profiles.npz" \
  --annot "$ANNOT" --variant "$V" --out "$DIR" --min_aa 5 --pval 0.05
echo "done=$(date -Iseconds): $(($(wc -l < "$DIR/${V}_collapsed.txt") - 1)) calls"
