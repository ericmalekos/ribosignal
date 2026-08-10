#!/usr/bin/env bash
#SBATCH --job-name=thp1_mzml
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/thp1_mzml_%A_%a.out
set -euo pipefail
# THP-1 BB7.2 immunopeptidome (PRIDE PXD015039) .raw -> INDEXED mzML, same converter and settings as
# the macrophage and HBL-1 arms (ThermoRawFileParser 2.0, -f 2, default centroiding) so search
# results stay comparable across pilots.
#
# This conversion IS the integrity check for these files. PRIDE publishes no per-file md5, so the
# download could only be size-checked; a truncated .raw fails to parse here rather than silently
# producing a short mzML. A converted file smaller than 1 MB is treated as a failure below.
#
#   sbatch --array=0-14%6 proteogenomics/scripts/thp1_raw2mzml.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PD=$NEW/proteogenomics/data/THP1_pilot
TRFP=/private/groups/carpenterlab/emalekos/conda_envs/proteomics/bin/ThermoRawFileParser
mkdir -p "$PD/mzml"

mapfile -t RAWS < <(ls "$PD"/immunopeptidome/*.raw | sort)
[ "$SLURM_ARRAY_TASK_ID" -lt "${#RAWS[@]}" ] || { echo "no file at index $SLURM_ARRAY_TASK_ID"; exit 0; }
RAW="${RAWS[$SLURM_ARRAY_TASK_ID]}"
B=$(basename "$RAW" .raw)
OUT="$PD/mzml/$B.mzML"

if [ -s "$OUT" ] && [ "$(stat -c%s "$OUT")" -gt 1000000 ]; then
  echo "SKIP $B (mzML exists, $(stat -c%s "$OUT") bytes)"; exit 0
fi
echo "CONVERT $B host=$(hostname) raw=$(stat -c%s "$RAW") start=$(date -Iseconds)"
"$TRFP" -i "$RAW" -b "$OUT" -f 2
SZ=$(stat -c%s "$OUT" 2>/dev/null || echo 0)
[ "$SZ" -gt 1000000 ] || { echo "FAIL $B: mzML is $SZ bytes -- raw likely truncated" >&2; exit 1; }
echo "DONE $B ($SZ bytes) $(date -Iseconds)"
