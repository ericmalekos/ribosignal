#!/usr/bin/env bash
#SBATCH --job-name=ribo_dropin_cmp
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/compare_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/compare_%j.err
set -euo pipefail
# Task 14 step 3: compare drop-in ORF calls vs real (genomic locus, 0.5x enrichment, 90 nt floor --
# the headline operating point). Backend chosen by RUN. Submit AFTER the dropin array, e.g.:
#   sbatch --dependency=afterok:<dropin_jobid> --export=ALL,RUN=<run_dir> ribocode_dropin_compare.sbatch
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
RUN=${RUN:-results/loto/orf_v2_attn_orthrus_holdout_Hepatocytes}
OFFICIAL=${OFFICIAL:-$ECH/data/ribocode_per_tissue/Hepatocytes/Hepatocytes_collapsed.txt}
DROPIN=$NEW/$RUN/dropin
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
export PYTHONNOUSERSITE=1
echo "run=$RUN host=$(hostname) start=$(date -Iseconds)"
"$PY" "$NEW/scripts/compare_dropin_calls.py" \
  --dropin_dir "$DROPIN" --official "$OFFICIAL" --profiles "$DROPIN/pred_profiles.npz" \
  --key genomic --min_enrichment 0.5 --min_len 90 --out "$DROPIN"
echo "=== per-type recall (pred_obsdepth vs real) ==="
"$PY" -c "
import json
d=json.load(open('$DROPIN/dropin_metrics.json'))
for t,v in d.get('pred_obsdepth_vs_real',{}).get('recall_per_type',{}).items():
    print(f'  {t:<16} {v[\"recovered\"]:>6}/{v[\"n_real\"]:<6} recall {v[\"recall\"]:.3f}')
"
echo "done=$(date -Iseconds)"
