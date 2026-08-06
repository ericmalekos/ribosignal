#!/usr/bin/env bash
#SBATCH --job-name=postpred_a549
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/postpred_a549_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/postpred_a549_%j.err
set -euo pipefail
# A549 (tryptic, fresh universe) post-predict, run afterok:dump_a549_fresh. enumerate/score ORFs ->
# model+null DBs (db_fresh, with class maps) -> tryptic MSFragger search (msf_full, DBDIR=db_fresh
# TAG=_fresh -> search_fresh_*) -> MS2Rescore fanout (12 fractions x 3 DBs) -> class-aware compare.
# Fully afterok-chained; no polling.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
S=$NEW/proteogenomics/scripts
P=$NEW/proteogenomics/data/A549_pilot
echo "=== post_predict_a549(mm1) $(date -Iseconds) ==="
[ -s "$P/pred_fresh/pred_profiles.npz" ] || { echo "NO pred_fresh/pred_profiles.npz -- predict failed" >&2; exit 1; }
echo "enumerate/score ORFs (fresh universe) ..."
$PY $S/enumerate_score_orfs.py \
  --profiles "$P/pred_fresh/pred_profiles.npz" --fasta "$P/A549_universe.fa" --out_dir "$P/orfs_fresh" 2>&1 | tail -3
echo "build model+null DBs (db_fresh) ..."
$PY $S/build_a549_dbs.py \
  --candidates "$P/orfs_fresh/candidates.faa" --out_dir "$P/db_fresh" 2>&1 | tail -4
# clear stale rescore/search PSMs from the prior (mm20) run so the fanout's resume guard does NOT skip
# and the class-aware compare cannot read stale results.
rm -f "$P"/rescore/*/*.psms.tsv "$P"/search_fresh_*/*.pepXML "$P"/search_fresh_*/*.pin 2>/dev/null || true
# tryptic search over the 12 TMT fractions x 3 DBs (msf_full is array 0-2 over DBs, each does all 12 mzML)
SJ=$(sbatch --parsable --export=ALL,DBDIR=db_fresh,TAG=_fresh $S/msf_full.sbatch)
# rescore: array 0-35 (3 DBs x 12 fractions), reads search_fresh_$D/$FR.pepXML
RJ=$(sbatch --parsable --dependency=afterok:$SJ --array=0-35 $S/ms2rescore_fanout.sbatch)
CJ=$(sbatch --parsable --dependency=afterok:$RJ $S/compare_a549.sbatch)
echo "A549 mm1 chain: search=$SJ (0-2)  rescore=$RJ (afterok:$SJ, 0-35)  compare=$CJ (afterok:$RJ)"
echo "=== post_predict_a549(mm1) submitted downstream $(date -Iseconds) ==="
