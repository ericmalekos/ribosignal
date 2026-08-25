#!/usr/bin/env bash
#SBATCH --job-name=pack_mm10
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=8:00:00
#SBATCH --array=0-7
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm10/pack_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/mm10/pack_%A_%a.err
set -euo pipefail
# Swap ONLY the RNA coverage channel to mm10; the Ribo target stays mm1, verbatim.
#
# prepare_pack.py --reuse-target reuses the source pack's universe AND target_counts unchanged and
# writes only new coverage -- its docstring names this the "swap RNA-seq" path (precedent: mm20->mm1).
# That is exactly right here: the label must stay mm1 (RiboCode does not drop multimappers and they
# destroy periodicity), while the INPUT is what was demonstrably corrupted -- GTF2I's CDS window
# reads 9% of the CDS mean at mm1 vs 138% at mm10 (results.md 2026-08-23).
#
# Output packs are NEW directories (packed_mm10cov_<tissue>); the existing packs are untouched, so
# the old and new can be compared and a retrain can be run both ways.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/python
COVD=$NEW/data/rnaseq_coverage/per_sample_mm10
TIS=(Fibroblast ES Fat HA_EC HCAEC Hepatocytes HUVEC VSMC)
SRC=(packed_union packed_union_ES packed_union_Fat packed_union_HA_EC packed_union_HCAEC \
     packed_union_Hepatocytes packed_union_HUVEC packed_union_VSMC)
T=${TIS[$SLURM_ARRAY_TASK_ID]}; S=${SRC[$SLURM_ARRAY_TASK_ID]}
# Fibroblast is the UNSUFFIXED member of a pack family -- dataset.tissue_pack_dir() maps
# RIBO_PACK_SUFFIX=mm10cov -> data/packed_mm10cov (Fibroblast) + data/packed_mm10cov_<T> (others).
if [ "$T" = "Fibroblast" ]; then OUT=$NEW/data/packed_mm10cov; else OUT=$NEW/data/packed_mm10cov_${T}; fi
[ -s "$OUT/coverage.npy" ] && { echo "SKIP $T (pack exists)"; exit 0; }
# this tissue's RNA runs -> their mm10 coverage hd5
if [ "$T" = "Fibroblast" ]; then
  mapfile -t RUNS < $NEW/data/fibroblast_rnaseq_srr.txt
else
  mapfile -t RUNS < <(awk -F'\t' -v t="$T" '$2==t{print $1}' $NEW/data/external/chothani_rna/srr_tissue.tsv)
fi
COV=(); for s in "${RUNS[@]}"; do
  f=$COVD/${s}_coverage.hd5
  [ -s "$f" ] && COV+=("$f") || echo "  WARN missing mm10 coverage for $s" >&2
done
[ "${#COV[@]}" -gt 0 ] || { echo "no mm10 coverage for $T" >&2; exit 1; }
echo "[$(date +%T)] $T: ${#COV[@]} of ${#RUNS[@]} RNA libs, target reused from $S"
export PYTHONNOUSERSITE=1
"$PY" $NEW/scripts/prepare/prepare_pack.py \
  --reuse-target "$NEW/data/$S" \
  --rna-coverage "${COV[@]}" \
  --group "$T" --species human --out "$OUT"
[ -s "$OUT/coverage.npy" ] || { echo "no coverage.npy written for $T" >&2; exit 2; }
# the target MUST be byte-identical to the source -- that is the whole contract of --reuse-target
if cmp -s "$NEW/data/$S/target_counts.npy" "$OUT/target_counts.npy"; then
  echo "[$(date +%T)] DONE $T -- target byte-identical to $S, coverage is NEW (mm10)"
else
  echo "FATAL: target_counts.npy CHANGED for $T; --reuse-target did not hold" >&2; exit 3
fi
