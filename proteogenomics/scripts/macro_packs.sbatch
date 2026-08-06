#!/usr/bin/env bash
#SBATCH --job-name=macro_pack
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --array=0-11
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/macro_pack_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/macro_pack_%A_%a.err
set -euo pipefail
# Per-population coverage-only pack + no-Kozak ORF track for the 12 macrophage populations. One population per
# array task. Feeds dump_pred_profiles (the model prediction) -> enumerate_score_orfs (three-frame DBs).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PD=$NEW/proteogenomics/data/macrophage_tissue
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
mkdir -p "$NEW/proteogenomics/logs"

mapfile -t POPS < <(ls "$PD"/universes/*_universe_tx.txt | xargs -n1 basename | sed 's/_universe_tx.txt//' | sort)
POP="${POPS[$SLURM_ARRAY_TASK_ID]:?no population at index $SLURM_ARRAY_TASK_ID}"
DS="macro_$POP"
PACK="$NEW/data/packed_heldout_$DS"
echo "population=$POP dataset=$DS host=$(hostname) start=$(date -Iseconds)"

# 1. coverage-only pack (pools this population's 3 RNA coverage reps onto its expressed universe)
"$PY" "$NEW/scripts/heldout/build_heldout_pack.py" "$DS"

# 2. no-Kozak ORF track on this population's universe (matches the deployed model's training input)
"$PY" "$NEW/scripts/build_orf_track.py" --mode ext --kozak none --pack "$PACK" \
  --fasta "$PD/universes/${POP}_universe.fa" --out "$PACK/orf_track_v2_nokozak.npy"

echo "DONE $POP: pack=$(wc -l < "$PACK/tx_order.txt") tx, orf_track=$(ls -la "$PACK/orf_track_v2_nokozak.npy" | awk '{print $5}') bytes $(date -Iseconds)"
