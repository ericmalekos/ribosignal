#!/usr/bin/env bash
#SBATCH --job-name=macro_mzml
#SBATCH --partition=medium
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --array=0-215%12
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/mzml_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/mzml_%A_%a.err
set -euo pipefail
# Convert the 216 IPX0001245001 macrophage .raw files to INDEXED mzML (ThermoRawFileParser 2.0, default
# centroiding) for the proteogenomics MS search. One file per array task, throttled to 12 concurrent.
# Idempotent: skips a file whose valid mzML already exists, so it can be re-run to mop up stragglers.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PD=$NEW/proteogenomics/data/macrophage_tissue
TRFP=/private/groups/carpenterlab/emalekos/conda_envs/proteomics/bin/ThermoRawFileParser
mkdir -p "$PD/mzml" "$NEW/proteogenomics/logs"

mapfile -t RAWS < <(ls "$PD"/proteomics_raw/*.raw | sort)
[ "$SLURM_ARRAY_TASK_ID" -lt "${#RAWS[@]}" ] || { echo "no file at index $SLURM_ARRAY_TASK_ID"; exit 0; }
RAW="${RAWS[$SLURM_ARRAY_TASK_ID]}"
B=$(basename "$RAW" .raw)
OUT="$PD/mzml/$B.mzML"

# skip if a non-trivial mzML already exists (mop-up re-runs)
if [ -s "$OUT" ] && [ "$(stat -c%s "$OUT")" -gt 1000000 ]; then
  echo "SKIP $B (mzML exists, $(stat -c%s "$OUT") bytes)"; exit 0
fi
echo "CONVERT $B host=$(hostname) start=$(date -Iseconds)"
# -f 2 = indexed mzML; default = centroid peak-picking (matches the HBL-1 processing)
"$TRFP" -i "$RAW" -b "$OUT" -f 2
echo "DONE $B ($(stat -c%s "$OUT" 2>/dev/null || echo 0) bytes) $(date -Iseconds)"
