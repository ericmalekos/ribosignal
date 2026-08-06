#!/usr/bin/env bash
# Re-run one immunopeptidome line on a RELEASED model with FROZEN search params.
#
# WHY. The original HBL-1 / DoHH2 / SU-DHL-4 results (2026-07-21) are wrong in two independent ways:
#   1. MODEL: they used results/loto/orf_v2_attn_onehot_holdout_Hepatocytes + the heuristic-Kozak
#      track -- five changes behind the released primary (architecture, union universe, Brain
#      dropped, --kozak none, mm1 coverage).
#   2. SEARCH PARAMS: fragger_hbl1_hla.params sets calibrate_mass = 2, which re-derives search
#      parameters PER DATABASE. On HBL-1 the null arm was scored with use_topN_peaks=100 /
#      intensity_transform=0 against the model arm's 150 / 1, cutoff 16.742 vs 18.483. The
#      model-vs-null comparison IS the claim, so it was never apples-to-apples.
#
# Both are fixed here in ONE search pass rather than two.
#
# SCOPE (deliberate): the DB-BUILDING METHOD IS HELD CONSTANT -- still build_a549_dbs.py with the
# pred_frame0 >= 0.5 cut, exactly as the original runs. So the delta is attributable to model +
# search params and nothing else. Switching to the pgx RiboCode-called builder (which the macrophage
# results use) is a separate, later comparison; doing it here would confound this one.
#
#   bash rerun_line_released.sh <LINE> [TAG]
#     LINE = HBL1 | DoHH2 | SUDHL4
#     TAG  = output suffix (default mamba4frozen) -> pred_<TAG> orfs_<TAG> db_<TAG> search_<TAG>_*
#
# Prereq: the predict must already have written <PILOT>/pred_<TAG>/pred_profiles.npz.
set -euo pipefail
LINE=${1:?usage: rerun_line_released.sh <LINE> [TAG]}
TAG=${2:-mamba4frozen}
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
P=$NEW/proteogenomics/data/${LINE}_pilot
# The predict dir is named by the MODEL, the output tag by model+params, so they differ. Pass it
# explicitly rather than deriving it by string surgery on TAG.
PRED=${PRED:-$P/pred_mamba4}
NPZ=$PRED/pred_profiles.npz

echo "=== rerun $LINE tag=$TAG $(date -Iseconds) ==="
[[ -s "$NPZ" ]] || { echo "ABORT: no $NPZ -- run the predict first" >&2; exit 1; }
echo "profiles: $NPZ ($(du -h "$NPZ" | cut -f1))"

echo "--- enumerate + score ORFs"
"$PY" "$NEW/proteogenomics/scripts/enumerate_score_orfs.py" \
  --profiles "$NPZ" --fasta "$P/${LINE}_universe.fa" --out_dir "$P/orfs_$TAG" 2>&1 | tail -3
[[ -s "$P/orfs_$TAG/candidates.faa" ]] || { echo "ABORT: no candidates.faa" >&2; exit 1; }

echo "--- build canonical / model / null DBs"
"$PY" "$NEW/proteogenomics/scripts/build_a549_dbs.py" \
  --candidates "$P/orfs_$TAG/candidates.faa" --out_dir "$P/db_$TAG" 2>&1 | tail -4
for d in canonical model null; do
  [[ -s "$P/db_$TAG/db_${d}.fasta" ]] || { echo "ABORT: missing db_${d}.fasta" >&2; exit 1; }
done

echo "--- submit HLA search (FROZEN params) + rescore"
SJ=$(sbatch --parsable --export=ALL,PILOT=$P,TEMPLATE=fragger_hbl1_hla_frozen.params,DBSUB=_$TAG \
      "$NEW/proteogenomics/scripts/msf_line_hla.sbatch")
NR=$(ls "$P"/immunopeptidome/*.mzML | wc -l); NA=$((3 * NR - 1))
RJ=$(sbatch --parsable --dependency=afterok:$SJ --array=0-$NA \
      --export=ALL,PILOT=$P,DBSUB=_$TAG "$NEW/proteogenomics/scripts/ms2rescore_line_fanout.sbatch")
echo "$LINE: search=$SJ  rescore=$RJ (afterok:$SJ, array 0-$NA over $NR mzML)"
echo "when both clear:  $PY $NEW/proteogenomics/scripts/compare_rescored.py --pilot ${LINE}_pilot --db_dir db_$TAG --dbs canonical model null"
