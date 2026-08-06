#!/usr/bin/env bash
#SBATCH --job-name=ms2r_line
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/ms2r_line_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/ms2r_line_%A_%a.err
set -uo pipefail
# Generalized per-fraction immunopeptidome MS2Rescore. Pass via --export: PILOT. Submit array 0-(3*Nreps-1).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
MS2R=/private/groups/carpenterlab/emalekos/conda_envs/ms2rescore/bin/ms2rescore
DBS=(canonical model null)
mapfile -t REPS < <(ls "$PILOT"/immunopeptidome/*.mzML | xargs -n1 basename | sed 's/.mzML//' | sort)
NR=${#REPS[@]}
D=${DBS[$((SLURM_ARRAY_TASK_ID / NR))]}
R=${REPS[$((SLURM_ARRAY_TASK_ID % NR))]}
OUT=$PILOT/rescore${DBSUB:-}${RTAG:-}/$D; mkdir -p "$OUT"    # DBSUB=_<tag> pairs with search_<tag>_<db>; RTAG=_seed3 for the stochasticity study
if [[ -s "$OUT/$R.psms.tsv" ]]; then echo "already done: $D/$R"; exit 0; fi
echo "PILOT=$(basename $PILOT) DB=$D rep=$R host=$(hostname) start=$(date -Iseconds)"
"$MS2R" -p "$PILOT/search${DBSUB:-}_$D/$R.pepXML" -s "$PILOT/immunopeptidome/$R.mzML" \
  -t pepxml -c "$NEW/proteogenomics/msfragger/ms2rescore_hbl1.json" \
  --output-path "$OUT/$R" 2>&1 | tail -8
echo "exit ${PIPESTATUS[0]} DB=$D R=$R done=$(date -Iseconds)"
[[ -s "$OUT/$R.psms.tsv" ]] && echo "  PSMs: $(($(wc -l < "$OUT/$R.psms.tsv")-1))" || echo "  NO psms.tsv"