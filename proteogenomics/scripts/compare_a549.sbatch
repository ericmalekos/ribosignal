#!/usr/bin/env bash
#SBATCH --job-name=cmp_a549
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/cmp_a549_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/cmp_a549_%j.err
set -euo pipefail
# A549 class-aware compare (afterok:rescore). compare_rescored_classaware.py is A549-hardcoded and
# overwrites db_comparison_rescored_classaware.md -- so first preserve the mm20 result, then snapshot the
# fresh mm1 result to a _mm1 copy for a side-by-side (expected ~identical per Task 25).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
# compare_rescored_classaware.py imports mokapot (class-aware re-training) -> use the ms2rescore env, not cas12a
PY=/private/groups/carpenterlab/emalekos/conda_envs/ms2rescore/bin/python
D=$NEW/proteogenomics/data/A549_pilot
CA=$D/db_comparison_rescored_classaware.md
[ -f "$CA" ] && [ ! -f "${CA%.md}_mm20.md" ] && cp "$CA" "${CA%.md}_mm20.md" && echo "preserved mm20 -> ${CA%.md}_mm20.md"
$PY $NEW/proteogenomics/scripts/compare_rescored_classaware.py 2>&1 | tail -30
cp "$CA" "${CA%.md}_mm1.md" && echo "wrote mm1 -> ${CA%.md}_mm1.md"
