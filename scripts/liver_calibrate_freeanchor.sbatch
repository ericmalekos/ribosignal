#!/usr/bin/env bash
#SBATCH --job-name=liver_calib
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=110G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/liver_calib_%A_%a.out
set -euo pipefail
# Ribo-seq-FREE Poisson calibration for the released mouse-liver + Ruiz-Orera dumps.
#
# THE POINT. pgx.calibrate anchors on recall of ANNOTATED CDS among expressed transcripts, in
# `relative` mode (hits vs hits at the saturating theta) -- both inputs available without any
# Ribo-seq: annotated CDS from the annotation, expression from the sample's own RNA-seq. That is the
# deployable anchor and the one every shipped pilot used (B721/THP-1/DoHH2 all chose theta=0.1).
#
# The OLD O2 liver analysis instead picked theta=0.02 by matching the observed Wang->Janich
# between-experiment CDS recall -- which needs two real Ribo-seq experiments. Quoting numbers from
# that theta as achievable on a Ribo-seq-free tissue would overstate the claim. This job produces the
# honest operating point, and the gap between the two thetas is itself a reportable result.
#
#   sbatch --array=0-7%3 scripts/liver_calibrate_freeanchor.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
export PYTHONPATH="$NEW/proteogenomics/scripts:${PYTHONPATH:-}"

NAMES=(attn_mouse_wang_liver mamba4_mouse_wang_liver
       attn_mouse_janich_liver_decon mamba4_mouse_janich_liver_decon
       attn_mouse_gse243134_liver mamba4_mouse_gse243134_liver
       released_attn released_mamba4)
SPECIES=(mouse mouse mouse mouse mouse mouse human human)
N=${NAMES[$SLURM_ARRAY_TASK_ID]}; SP=${SPECIES[$SLURM_ARRAY_TASK_ID]}
if [[ "$SP" == "human" ]]; then
  PROF=$NEW/results/heldout/human_ruizorera/$N/dropin/pred_profiles.npz
  ANNOT=$ECH/data/ribocode_annot
  OUT=$NEW/results/heldout/human_ruizorera/$N/calib
else
  PROF=$NEW/results/liver_released/$N/pred_profiles.npz
  ANNOT="/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/mouse_ribocode_annot"
  OUT=$NEW/results/liver_released/$N/calib
fi
[[ -s "$PROF" ]] || { echo "missing $PROF" >&2; exit 1; }
mkdir -p "$OUT"
echo "name=$N species=$SP host=$(hostname) start=$(date -Iseconds)"
cd "$NEW/proteogenomics/scripts"
"$PY" -m pgx.calibrate --profiles "$PROF" --annot "$ANNOT" --species "$SP" \
  --out "$OUT" --cds-recall 0.90 --recall-mode relative --min-aa 5 --pval 0.05
echo "done=$(date -Iseconds)"
"$PY" -c "
import json,pathlib
f=pathlib.Path('$OUT/theta.json')
print('  picked:', json.load(open(f)) if f.exists() else 'no theta.json')" 2>/dev/null || true
