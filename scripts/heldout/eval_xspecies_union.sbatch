#!/usr/bin/env bash
#SBATCH --job-name=xsp_union
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=0-4
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/xsp_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/xsp_%A_%a.err
set -euo pipefail
# Cross-species eval of a UNION-trained checkpoint on the 5 mouse held-out sets. MODEL env picks the
# checkpoint (attn_union|mamba4_union) + SIF + default device; the array index picks the dataset. Chain:
# dump -> RiboCode drop-in (mouse annot; real/pred_obsdepth/pred_preddepth) -> Poisson CDS-anchored sweep ->
# orf_call_metrics (20-aa primary + two-arm) -> xspecies doc. Mirrors eval_tcell_dropin_cpu.sbatch exactly,
# only the checkpoint + per-dataset paths are parameterized. min_signal 50 matches the flagship T-cell run
# so union results are directly comparable to the deployed-fibroblast cross-species numbers.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
REF=$NEW/data/heldout_refs
GSE=$NEW/proteogenomics/data/gse120762
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
ANNOT="/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/mouse_ribocode_annot"                       # MOUSE annotation for the drop-in
RIBOPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
CASPY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
MODEL=${MODEL:?set MODEL=attn_union|mamba4_union}

case "$MODEL" in
  attn_union)   RUN=results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes;   SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img;  DEVDEF=cpu ;;
  mamba4_union) RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes; SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img;  DEVDEF=cuda ;;
  *) echo "bad MODEL=$MODEL (want attn_union|mamba4_union)" >&2; exit 1 ;;
esac
DEV=${DEV:-$DEVDEF}; NV=""; [ "$DEV" = cuda ] && NV="--nv"

DATASETS=(mouse_gse155087_tcell mouse_wang_liver mouse_janich_liver mouse_gse120762_nt mouse_gse120762_lps)
DS=${DATASETS[$SLURM_ARRAY_TASK_ID]:?no dataset at index $SLURM_ARRAY_TASK_ID}
# PACKSUF appends a suffix to the PACK / heldout name so a variant eval pack can be scored with the
# identical chain (e.g. PACKSUF=_mm10 -> packed_heldout_<DS>_mm10). The FASTA, annotation, label and
# universe stay tied to the base DS, so the two runs are directly comparable. Empty = original.
PACKSUF=${PACKSUF:-}
DSPACK=${DS}${PACKSUF}
case "$DS" in
  mouse_gse155087_tcell) FA=$REF/mouse_gse155087_tcell_universe.fa; LBL="T-cell GSE155087" ;;
  mouse_wang_liver)      FA=$REF/mouse_wang_universe.fa;            LBL="Wang liver" ;;
  mouse_janich_liver)    FA=$REF/mouse_janich_universe.fa;          LBL="Janich liver" ;;
  mouse_gse120762_nt)    FA=$GSE/gse120762_nt_universe.fa;          LBL="GSE120762 NT" ;;
  mouse_gse120762_lps)   FA=$GSE/gse120762_lps_universe.fa;         LBL="GSE120762 LPS" ;;
esac
PACK=$NEW/data/packed_heldout_$DSPACK
TRACK=$PACK/orf_track_v2_nokozak.npy
OUT=$NEW/results/heldout_xspecies/$MODEL/$DSPACK
mkdir -p "$OUT/dropin" "$NEW/logs/heldout"
echo "xspecies MODEL=$MODEL DS=$DS pack=$DSPACK dev=$DEV host=$(hostname) start=$(date -Iseconds)"
[ -s "$NEW/$RUN/best.pt" ] || { echo "ABORT: missing checkpoint $RUN/best.pt"; exit 1; }
[ -s "$TRACK" ] || { echo "ABORT: missing nokozak track $TRACK"; exit 1; }

# 1. dump predicted profiles (mouse FASTA + dataset nokozak ORF track; min_signal 50 = tx with observed signal)
# shellcheck disable=SC2086
singularity exec $NV --no-home --cleanenv --env PYTHONNOUSERSITE=1 --env OMP_NUM_THREADS=16 \
  --env RIBO_ORF_TRACK="$TRACK" --env RIBO_ONEHOT_FASTA="$FA" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" --heldout "$DSPACK" \
    --out "$OUT/dropin" --fasta "$FA" --min_signal 50 --device "$DEV"

# guard: heavy skip = FASTA/pack mismatch
"$CASPY" -c "
import numpy as np,sys
d=np.load('$OUT/dropin/pred_profiles.npz',allow_pickle=True)
nt=len(d['tx_ids']); ns=int(str(d['meta'][3])); frac=ns/(nt+ns) if (nt+ns) else 0
print(f'dump: {nt} tx, {ns} skipped ({frac:.1%})'); sys.exit(1 if frac>0.05 else 0)" \
  || { echo 'ABORT: >5% tx skipped (FASTA/pack mismatch)'; exit 1; }

# 2. RiboCode drop-in (mouse annot): real / pred_obsdepth / pred_preddepth
export PYTHONNOUSERSITE=1
for V in real pred_obsdepth pred_preddepth; do
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

# 4. standing-rule ORF-call metrics (20-aa primary + two-arm) -> xspecies doc
mkdir -p "$NEW/results/orf_call_metrics"
"$CASPY" "$NEW/scripts/heldout/orf_call_metrics.py" --run "$OUT" \
  --label "$LBL ($MODEL xspecies${PACKSUF:+ $PACKSUF})" --universe "xsp-$MODEL${PACKSUF}" \
  --poisson-sweep "$SW" --cds-recall 0.90 \
  --out-json "$NEW/results/orf_call_metrics/xsp_${MODEL}_${DSPACK}.json" \
  --append-md "$NEW/docs/manuscript_xspecies_union.md"
echo "xspecies done MODEL=$MODEL DS=$DS $(date -Iseconds)"
