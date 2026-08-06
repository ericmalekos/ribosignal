#!/usr/bin/env bash
#SBATCH --job-name=ms2r_test
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/ms2r_test_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/ms2r_test_%j.err
set -uo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
MS2R=/private/groups/carpenterlab/emalekos/conda_envs/ms2rescore/bin/ms2rescore
P=$NEW/proteogenomics/data/A549_pilot
OUT=$P/rescore_test; mkdir -p "$OUT"
echo "host=$(hostname) start=$(date -Iseconds)"
"$MS2R" -p "$P/search_fresh_model/d00615_Prot_12_01.pepXML" \
  -s "$P/mzml/d00615_Prot_12_01.mzML" \
  -t pepxml -c "$NEW/proteogenomics/msfragger/ms2rescore_a549.json" \
  --output-path "$OUT/d00615_model" 2>&1 | tail -30
echo "exit ${PIPESTATUS[0]} done=$(date -Iseconds)"
ls -la "$OUT"/d00615_model* 2>/dev/null | awk '{print $5,$9}'
