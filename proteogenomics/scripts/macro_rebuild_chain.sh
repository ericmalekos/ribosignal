#!/usr/bin/env bash
set -uo pipefail
# Drive the macrophage churn rebuild for ONE checkpoint, from finished dumps through the DB build.
#
#   macro_rebuild_chain.sh <TAG> <dump_jobid>
#
# Stages (each waits on the previous):
#   1. wait for the dump array to drain
#   2. macro_enumerate  -> orfs<TAG>/<pop>/candidates.{tsv,faa}
#   3. calibrate_f0_threshold.py -> CDS-anchored pred_frame0 cut (replaces the uncalibrated 0.5)
#   4. macro_dbs x2     -> db<TAG>/        at the legacy 0.5 cut
#                       -> db<TAG>_cal/    at the CDS-anchored cut
#
# TWO ARMS ON PURPOSE. Changing the checkpoint and the selection threshold in the same step would
# confound them: a shift in novel-peptide yield could be the new model or the new cut. The legacy-0.5 arm
# holds the threshold fixed so the checkpoint effect is isolated, and the calibrated arm holds the
# checkpoint fixed so the threshold effect is isolated. This also satisfies the standing two-arm rule.
#
# MSFragger is NOT launched here -- it is the 12 h/task pole and should be started deliberately once the
# DB sizes have been eyeballed.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PD=$NEW/proteogenomics/data/macrophage_tissue
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
S=$NEW/proteogenomics/scripts
TAG=${1:?usage: macro_rebuild_chain.sh <TAG> <dump_jobid>}
DUMPJOB=${2:?usage: macro_rebuild_chain.sh <TAG> <dump_jobid>}

echo "=== chain $TAG (dumps=$DUMPJOB) start $(date -Iseconds) ==="
while squeue -j "$DUMPJOB" -h 2>/dev/null | grep -q .; do sleep 120; done
NPRED=$(ls "$PD"/pred"$TAG"/*/pred_profiles.npz 2>/dev/null | wc -l)
echo "dumps drained: $NPRED/12 pred_profiles.npz"
[ "$NPRED" -eq 12 ] || { echo "ABORT: only $NPRED/12 dumps present"; exit 1; }

echo "--- enumerate ---"
JE=$(sbatch --parsable --job-name="enum${TAG}" --export=ALL,TAG="$TAG" "$S/macro_enumerate.sbatch")
echo "  job $JE"; while squeue -j "$JE" -h 2>/dev/null | grep -q .; do sleep 60; done
NC=$(ls "$PD"/orfs"$TAG"/*/candidates.faa 2>/dev/null | wc -l)
echo "  candidates: $NC/12"
[ "$NC" -eq 12 ] || { echo "ABORT: only $NC/12 enumerations"; exit 1; }

echo "--- calibrate CDS-anchored pred_frame0 threshold ---"
"$PY" "$S/calibrate_f0_threshold.py" --orfs-dir "$PD/orfs$TAG" --out "$PD/calib/${TAG#_}" \
  --cds-recall 0.90 --label "${TAG#_}" | tee "$PD/calib/${TAG#_}_calibration.txt"
TH=$("$PY" -c "import json;print(json.load(open('$PD/calib/${TAG#_}/f0_calibration.json'))['pooled_threshold'])")
echo "  CDS-anchored threshold = $TH   (legacy uncalibrated cut = 0.5)"

echo "--- DBs, arm 1: legacy 0.5 (isolates the checkpoint effect) ---"
J1=$(sbatch --parsable --job-name="db${TAG}" --export=ALL,TAG="$TAG" "$S/macro_dbs.sbatch")
echo "--- DBs, arm 2: CDS-anchored $TH (isolates the threshold effect) ---"
# macro_dbs passes --thresh 0.5 literally, so the calibrated arm overrides it via THRESH.
J2=$(sbatch --parsable --job-name="dbcal${TAG}" --export=ALL,TAG="${TAG}_cal",ORFS_TAG="$TAG",THRESH="$TH" "$S/macro_dbs.sbatch")
echo "  legacy job $J1   calibrated job $J2"
for J in $J1 $J2; do while squeue -j "$J" -h 2>/dev/null | grep -q .; do sleep 60; done; done

echo
echo "=== DB SIZES ($TAG) ==="
printf "  %-18s %12s %12s %12s\n" population "model@0.5" "model@$TH" "null"
for p in $(ls "$PD/db$TAG" 2>/dev/null); do
  m=$(grep -c '^>' "$PD/db$TAG/$p/db_model.fasta" 2>/dev/null || echo 0)
  c=$(grep -c '^>' "$PD/db${TAG}_cal/$p/db_model.fasta" 2>/dev/null || echo 0)
  n=$(grep -c '^>' "$PD/db$TAG/$p/db_null.fasta" 2>/dev/null || echo 0)
  printf "  %-18s %12s %12s %12s\n" "$p" "$m" "$c" "$n"
done
echo
echo "NEXT (deliberate, 12 h/task): submit MSFragger for both arms -- DBLIST=model ONLY."
echo "  The canonical and null searches are REUSED from the original untagged run:"
echo "    canonical  contains zero model content."
echo "    null       VERIFIED 2026-07-31: target-SEQUENCE sets are byte-identical across checkpoints"
echo "               (BMDM 244,269 both). Only the ACCESSIONS differ (76,300 of them), because"
echo "               build_a549_dbs dedups by sequence and the f0 tie-break picks a different"
echo "               representative ORF id. Class-specific FDR keys on the shared nuORF| prefix,"
echo "               so the peptide-level result is unchanged. Do NOT skip this check on a rebuild"
echo "               that alters enumeration (min_aa, NOVEL_CLASSES, or the universe FASTA)."
echo "  TAG=$TAG      DBLIST=model sbatch --array=0-11%6 $S/macro_msf.sbatch"
echo "  TAG=${TAG}_cal DBLIST=model sbatch --array=0-11%6 $S/macro_msf.sbatch"
echo "THEN:"
echo "  macro_churn_aggregate.py --tag $TAG      --reuse-tag \"\" "
echo "  macro_churn_aggregate.py --tag ${TAG}_cal --reuse-tag \"\" "
echo "=== chain $TAG done $(date -Iseconds) ==="
