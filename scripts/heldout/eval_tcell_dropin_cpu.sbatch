#!/usr/bin/env bash
#SBATCH --job-name=tcell_eval_cpu
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/tcell_eval_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/tcell_eval_%j.err
set -euo pipefail
# GSE155087 mouse CD4+ T-cell held-out eval, CPU (GPU partition jammed). Deployed no-Kozak mm1 model predicts
# the T-cell per-nt P-site density on the T-cell pack, then RiboCode detectORF (MOUSE annotation) on real /
# pred_obsdepth / pred_preddepth + the Poisson CDS-anchored arm -> orf_call_metrics (20-aa primary + sweep +
# two-arm) appended to the manuscript doc. Mirrors the union eval but for the mouse held-out pack.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
DATASET=mouse_gse155087_tcell
RUN=results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes    # deployed no-Kozak mm1 model
ORFTRACK=$NEW/data/packed_heldout_${DATASET}/orf_track_v2_nokozak.npy
FASTA=$NEW/data/heldout_refs/${DATASET}_universe.fa
ANNOT=$ECH/data/ribocode_annot_mouse                                           # MOUSE annotation for the drop-in
OUT=$NEW/results/heldout/$DATASET/onehot
RIBOPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
CASPY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
mkdir -p "$OUT/dropin" "$NEW/logs/heldout"
echo "T-cell eval DATASET=$DATASET run=$RUN host=$(hostname) CPU start=$(date -Iseconds)"

# 1. dump predicted profiles on the T-cell held-out pack (CPU, deployed model, mouse FASTA + T-cell ORF track)
# shellcheck disable=SC2086
singularity exec --no-home --cleanenv --env PYTHONNOUSERSITE=1 --env OMP_NUM_THREADS=16 \
  --env RIBO_ORF_TRACK="$ORFTRACK" --env RIBO_ONEHOT_FASTA="$FASTA" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" --heldout "$DATASET" \
    --out "$OUT/dropin" --fasta "$FASTA" --min_signal 50 --device cpu
echo "dump done=$(date -Iseconds)"

# guard: heavy skip = FASTA/pack mismatch -> abort before publishing wrong numbers
"$CASPY" -c "
import numpy as np,sys
d=np.load('$OUT/dropin/pred_profiles.npz',allow_pickle=True)
nt=len(d['tx_ids']); ns=int(str(d['meta'][3])); frac=ns/(nt+ns) if (nt+ns) else 0
print(f'dump: {nt} tx, {ns} skipped ({frac:.1%})'); sys.exit(1 if frac>0.05 else 0)" \
  || { echo 'ABORT: >5% tx skipped (FASTA/pack mismatch)'; exit 1; }

# 2. RiboCode drop-in (MOUSE annotation): real / pred_obsdepth / pred_preddepth
export PYTHONNOUSERSITE=1
for V in real pred_obsdepth pred_preddepth; do
  echo "=== dropin $V ==="
  "$RIBOPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$OUT/dropin/pred_profiles.npz" \
    --annot "$ANNOT" --variant "$V" --out "$OUT/dropin" --min_aa 5 --pval 0.05
done

# 3. two-arm Poisson CDS-anchored sweep
SW=$OUT/dropin_sweep_poisson
for TH in 0.02 0.05 0.1 0.2; do
  O=$SW/theta_$TH; rm -rf "$O"; mkdir -p "$O"
  "$RIBOPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$OUT/dropin/pred_profiles.npz" \
    --annot "$ANNOT" --variant pred_preddepth --pred_scale "$TH" --pred_poisson --seed 0 \
    --out "$O" --min_aa 5 --pval 0.05
done

# 4. standing-rule ORF-call metrics (20-aa primary + sweep + two-arm) -> manuscript doc
mkdir -p "$NEW/results/orf_call_metrics"
"$CASPY" "$NEW/scripts/heldout/orf_call_metrics.py" --run "$OUT" \
  --label "T-cell GSE155087 (mouse held-out)" --universe "mouse-tcell" \
  --poisson-sweep "$SW" --cds-recall 0.90 \
  --out-json "$NEW/results/orf_call_metrics/tcell_gse155087.json" \
  --append-md "$NEW/docs/manuscript_orf_call_eval.md"
echo "T-cell eval done=$(date -Iseconds)"
