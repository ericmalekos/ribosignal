#!/usr/bin/env bash
#SBATCH --job-name=heldout_pack_mm10
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=8:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout_mm10/pack_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout_mm10/pack_%A_%a.err
set -euo pipefail
# Swap ONLY the RNA coverage channel of a HELD-OUT eval pack to mm10; the Ribo target stays mm1.
# Same contract as the training-side rebuild: --reuse-target copies universe + target_counts
# verbatim, and the job FAILS if target_counts.npy is not byte-identical afterwards.
#
# DS is set per array index. Add datasets here as their mm10 coverage lands.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/python
DS_LIST=(mouse_janich_liver)
DS=${DS_LIST[$SLURM_ARRAY_TASK_ID]}
SRC=$NEW/data/packed_heldout_${DS}
COVD=$NEW/data/heldout_rnaseq_coverage_mm10/${DS}
OUT=$NEW/data/packed_heldout_${DS}_mm10
[ -d "$SRC" ] || { echo "no source pack $SRC" >&2; exit 1; }
[ -s "$OUT/coverage.npy" ] && { echo "SKIP $DS (pack exists)"; exit 0; }
mapfile -t COV < <(ls "$COVD"/*_coverage.hd5 2>/dev/null)
[ "${#COV[@]}" -gt 0 ] || { echo "no mm10 coverage in $COVD" >&2; exit 1; }
echo "[$(date +%T)] $DS: ${#COV[@]} mm10 RNA libs, target reused from packed_heldout_${DS}"
export PYTHONNOUSERSITE=1
"$PY" "$NEW/scripts/prepare/prepare_pack.py" \
  --reuse-target "$SRC" --rna-coverage "${COV[@]}" \
  --group "${DS}_mm10" --species mouse --out "$OUT"
[ -s "$OUT/coverage.npy" ] || { echo "no coverage.npy written for $DS" >&2; exit 2; }
if cmp -s "$SRC/target_counts.npy" "$OUT/target_counts.npy"; then
  echo "[$(date +%T)] DONE $DS -- target byte-identical, coverage is NEW (mm10)"
else
  echo "FATAL: target_counts.npy CHANGED for $DS; --reuse-target did not hold" >&2; exit 3
fi
