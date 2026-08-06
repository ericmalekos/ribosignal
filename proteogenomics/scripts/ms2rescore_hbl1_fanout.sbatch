#!/usr/bin/env bash
#SBATCH --job-name=hbl1_ms2r
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --array=0-7
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/hbl1_ms2r_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/hbl1_ms2r_%A_%a.err
set -uo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
MS2R=/private/groups/carpenterlab/emalekos/conda_envs/ms2rescore/bin/ms2rescore
HB=$NEW/proteogenomics/data/HBL1_pilot
DBS=(canonical model null price)
REPS=(HBL1_DMSO_200M_050219_1 HBL1_DMSO_200M_050219_2)
NR=${#REPS[@]}                                    # 2
D=${DBS[$((SLURM_ARRAY_TASK_ID / NR))]}
R=${REPS[$((SLURM_ARRAY_TASK_ID % NR))]}
OUT=$HB/rescore/$D; mkdir -p "$OUT"
if [[ -s "$OUT/$R.psms.tsv" ]]; then echo "already done: $D/$R"; exit 0; fi
echo "host=$(hostname) DB=$D rep=$R start=$(date -Iseconds)"
# immunopeptidome = label-free high-res: basic + ms2pip immuno-HCD + DeepLC (config ms2rescore_hbl1.json)
"$MS2R" -p "$HB/search_$D/$R.pepXML" -s "$HB/immunopeptidome/$R.mzML" \
  -t pepxml -c "$NEW/proteogenomics/msfragger/ms2rescore_hbl1.json" \
  --output-path "$OUT/$R" 2>&1 | tail -12
rc=${PIPESTATUS[0]}
echo "exit $rc DB=$D rep=$R done=$(date -Iseconds)"
[[ -s "$OUT/$R.psms.tsv" ]] && echo "  rescored PSMs: $(($(wc -l < "$OUT/$R.psms.tsv")-1))" || echo "  NO psms.tsv (exit $rc)"