#!/usr/bin/env bash
#SBATCH --job-name=a549_ms2r
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --array=0-35
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/a549_ms2r_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/a549_ms2r_%A_%a.err
set -uo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
MS2R=/private/groups/carpenterlab/emalekos/conda_envs/ms2rescore/bin/ms2rescore
P=$NEW/proteogenomics/data/A549_pilot
DBS=(canonical model null)
FRACS=(d00615_Prot_12_01 d00616_Prot_12_02 d00617_Prot_12_03 d00618_Prot_12_04 d00619_Prot_12_05 d00620_Prot_12_06 \
       d00621_Prot_12_07 d00622_Prot_12_08 d00623_Prot_12_09 d00624_Prot_12_10 d00625_Prot_12_11 d00626_Prot_12_12)
NF=${#FRACS[@]}                                   # 12
D=${DBS[$((SLURM_ARRAY_TASK_ID / NF))]}
FR=${FRACS[$((SLURM_ARRAY_TASK_ID % NF))]}
OUT=$P/rescore/$D; mkdir -p "$OUT"
if [[ -s "$OUT/$FR.psms.tsv" ]]; then echo "already done: $D/$FR"; exit 0; fi
echo "host=$(hostname) DB=$D fraction=$FR start=$(date -Iseconds)"
# per-fraction: one pepXML + its specific mzML (the invocation that works); basic + ms2pip CID (no DeepLC/TMT)
"$MS2R" -p "$P/search_fresh_$D/$FR.pepXML" -s "$P/mzml/$FR.mzML" \
  -t pepxml -c "$NEW/proteogenomics/msfragger/ms2rescore_a549.json" \
  --output-path "$OUT/$FR" 2>&1 | tail -8
echo "exit ${PIPESTATUS[0]} DB=$D FR=$FR done=$(date -Iseconds)"
[[ -s "$OUT/$FR.psms.tsv" ]] && echo "  rescored PSMs: $(($(wc -l < "$OUT/$FR.psms.tsv")-1))" || echo "  NO psms.tsv"