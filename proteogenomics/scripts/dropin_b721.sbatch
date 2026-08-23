#!/usr/bin/env bash
# B721 drop-in: RiboCode on the model's PREDICTED profile, scored against the MEASURED calls from
# the same cells. First ground-truth check in the project on an immunopeptidome line -- every prior
# MS comparison was model-vs-null, not model-vs-measured-translation.
#
# Two arms per model, per the standing two-arm rule:
#   pred_preddepth            theta = 1, deterministic  (standalone: no real Ribo-seq at all)
#   pred_preddepth + poisson  calibrated theta          (the deployed recipe)
#
#   sbatch --array=0-3 proteogenomics/scripts/dropin_b721.sbatch
#SBATCH --job-name=dropin_b721
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=110G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/dropin_b721_%A_%a.out
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
RPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python
ANNOT=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/human_ribocode_annot_primary
P=$NEW/proteogenomics/data/B721_pilot

MODELS=(mamba4 mamba4 attn attn)
ARMS=(theta1 poisson theta1 poisson)
M=${MODELS[$SLURM_ARRAY_TASK_ID]}; A=${ARMS[$SLURM_ARRAY_TASK_ID]}
# theta comes from that model's own pgx calibration, so the arm matches the DB that was built
THETA=$(/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 -c \
  "import json;print(json.load(open('$P/pgx_$M/calib/theta.json'))['theta'])")
OUT=$P/dropin_$M/$A; mkdir -p "$OUT"
echo "host=$(hostname) model=$M arm=$A theta=$THETA start=$(date -Iseconds)"

EXTRA=()
if [ "$A" = "poisson" ]; then EXTRA=(--pred_poisson --pred_scale "$THETA"); fi
"$RPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$P/pred_$M/pred_profiles.npz" \
  --annot "$ANNOT" --variant pred_preddepth --out "$OUT" --min_aa 5 --pval 0.05 "${EXTRA[@]}"
echo "done=$(date -Iseconds): $(wc -l < "$OUT"/pred_preddepth_collapsed.txt 2>/dev/null || echo 0) calls"
