#!/usr/bin/env bash
#SBATCH --job-name=postpred
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/postpred_%x_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/postpred_%x_%j.err
set -euo pipefail
# SLURM-native post-predict (run afterok:predict, so pred_profiles.npz is guaranteed present -- no polling,
# robust to arbitrary GPU-queue delays). LINE via --export. enumerate/score ORFs -> model+null DBs ->
# submit HLA search -> MS2Rescore (afterok) -> class-FDR compare (afterok). Then exit; the compare job
# writes the final result. Mirrors post_predict_line.sh but fully dependency-chained.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
S=$NEW/proteogenomics/scripts
P=$NEW/proteogenomics/data/${LINE}_pilot
echo "=== post_predict(mm1) $LINE $(date -Iseconds) ==="
[ -s "$P/pred/pred_profiles.npz" ] || { echo "NO pred_profiles.npz -- predict failed" >&2; exit 1; }
echo "enumerate/score ORFs ..."
$PY $S/enumerate_score_orfs.py \
  --profiles "$P/pred/pred_profiles.npz" --fasta "$P/${LINE}_universe.fa" --out_dir "$P/orfs" 2>&1 | tail -3
echo "build model+null DBs ..."
$PY $S/build_a549_dbs.py \
  --candidates "$P/orfs/candidates.faa" --out_dir "$P/db" 2>&1 | tail -4
# clear stale rescore/search PSMs from any prior (e.g. mm20) run so the fanout's resume guard does NOT
# skip and the compare cannot read stale results (this bug silently reused mm20 PSMs on the first mm1 pass).
rm -f "$P"/rescore/*/*.psms.tsv "$P"/search_*/*.pepXML "$P"/search_*/*.pin 2>/dev/null || true
SJ=$(sbatch --parsable --export=ALL,PILOT=$P $S/msf_line_hla.sbatch)
NR=$(ls "$P"/immunopeptidome/*.mzML | wc -l); NA=$((3*NR-1))
RJ=$(sbatch --parsable --dependency=afterok:$SJ --array=0-$NA --export=ALL,PILOT=$P $S/ms2rescore_line_fanout.sbatch)
CJ=$(sbatch --parsable --dependency=afterok:$RJ --export=ALL,LINE=$LINE $S/compare_line.sbatch)
echo "$LINE mm1 chain: search=$SJ  rescore=$RJ (afterok:$SJ, array 0-$NA)  compare=$CJ (afterok:$RJ)"
echo "=== post_predict(mm1) $LINE submitted downstream $(date -Iseconds) ==="
