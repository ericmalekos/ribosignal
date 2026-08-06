#!/usr/bin/env bash
# Turnkey post-predict chain for one MHC-I line. Args: LINE PREDICT_JOBID
# Waits for the GPU predict, then: enumerate/score ORFs -> model+null DBs -> HLA search -> MS2Rescore -> class-FDR.
LINE=$1; PJ=$2
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
P=$NEW/proteogenomics/data/${LINE}_pilot
echo "=== post_predict $LINE (waiting on predict $PJ) $(date -Iseconds) ==="
n=0; until ! squeue -j "$PJ" -h 2>/dev/null | grep -q .; do n=$((n+1)); [ $n -ge 90 ] && break; sleep 120; done
[ -s "$P/pred/pred_profiles.npz" ] || { echo "NO pred_profiles.npz -- predict failed"; tail -5 $NEW/proteogenomics/logs/dump_line_${PJ}.err 2>/dev/null; exit 1; }
echo "predict done. enumerate/score ORFs ..."
$PY $NEW/proteogenomics/scripts/enumerate_score_orfs.py \
  --profiles "$P/pred/pred_profiles.npz" --fasta "$P/${LINE}_universe.fa" --out_dir "$P/orfs" 2>&1 | tail -3
echo "build model+null DBs ..."
$PY $NEW/proteogenomics/scripts/build_a549_dbs.py \
  --candidates "$P/orfs/candidates.faa" --out_dir "$P/db" 2>&1 | tail -4
# submit HLA search (canonical/model/null) then MS2Rescore (afterok)
SJ=$(sbatch --parsable --export=ALL,PILOT=$P $NEW/proteogenomics/scripts/msf_line_hla.sbatch)
NR=$(ls "$P"/immunopeptidome/*.mzML | wc -l); NA=$((3*NR-1))
RJ=$(sbatch --parsable --dependency=afterok:$SJ --array=0-$NA --export=ALL,PILOT=$P $NEW/proteogenomics/scripts/ms2rescore_line_fanout.sbatch)
echo "search=$SJ  rescore=$RJ (afterok:$SJ, array 0-$NA)"
n=0; until ! squeue -j $SJ,$RJ -h 2>/dev/null | grep -q .; do n=$((n+1)); [ $n -ge 90 ] && break; sleep 120; done
echo "=== $LINE search+rescore cleared $(date -Iseconds) ==="
for d in canonical model null; do echo "  rescore/$d: $(ls $P/rescore/$d/*.psms.tsv 2>/dev/null|wc -l)/$NR reps"; done
echo "=== $LINE RESCORED model/null (raw baseline unknown until raw compare) ==="
$PY $NEW/proteogenomics/scripts/compare_rescored.py --pilot ${LINE}_pilot --db_dir db --dbs canonical model null 2>&1 | head -30
echo "=== post_predict $LINE DONE $(date -Iseconds) ==="
