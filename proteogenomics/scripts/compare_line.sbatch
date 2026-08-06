#!/usr/bin/env bash
#SBATCH --job-name=cmp_resc
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/cmp_%x_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/cmp_%x_%j.err
set -euo pipefail
# Final class-FDR model/null/canonical comparison on the rescored PSMs (run afterok:rescore).
# LINE via --export. Writes the mm1 result to db_comparison_rescored_mm1.md alongside the mm20 one, so the
# two postures can be diffed (expected ~identical per results.md Task 25).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
OUT=$NEW/proteogenomics/data/${LINE}_pilot/db_comparison_rescored_mm1.md
echo "=== compare(mm1) $LINE $(date -Iseconds) ===" | tee "$OUT"
$PY $NEW/proteogenomics/scripts/compare_rescored.py --pilot ${LINE}_pilot --db_dir db --dbs canonical model null 2>&1 | tee -a "$OUT"
echo "wrote $OUT"
