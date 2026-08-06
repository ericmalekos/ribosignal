#!/usr/bin/env bash
#SBATCH --job-name=macro_db
#SBATCH --partition=medium
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=2:00:00
#SBATCH --array=0-11
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/macro_db_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/macro_db_%A_%a.err
set -euo pipefail
# Build the canonical/model/null search DBs per macrophage population from candidates.faa + the MOUSE vM38
# canonical proteome. All three DBs share the canonical base + inline REV_ decoys and differ ONLY in the
# novel-ORF set: db_model = pred_frame0>=0.5 selected, db_null = ALL candidates, db_canonical = baseline.
# build_a549_dbs.py is species-agnostic; --canonical points it at the mouse proteome. cas12a env (stdlib).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PD=$NEW/proteogenomics/data/macrophage_tissue
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
CANON=/private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations/gencode_proteins/gencode.vM38.pc_translations.fa.gz
# TAG namespaces artefacts per checkpoint (see macro_dumps.sbatch). Empty = original single-model layout.
# ORFS_TAG lets the DB dir differ from the candidates dir, so the SAME candidates can be cut at two
# thresholds into two DB trees (db<TAG> at the legacy 0.5, db<TAG>_cal at the CDS-anchored cut) without
# re-enumerating. THRESH is the pred_frame0 selection cut.
#
# The default 0.5 is UNCALIBRATED and far too permissive: on the deployed checkpoint it retains 99.7% of
# annotated CDS (i.e. it barely filters on the canonical axis) and admits ~7.9x more novel ORFs than a
# CDS-recall-anchored cut. Use calibrate_f0_threshold.py to derive a per-checkpoint anchored value.
TAG=${TAG:-}
ORFS_TAG=${ORFS_TAG:-$TAG}
THRESH=${THRESH:-0.5}
mkdir -p "$NEW/proteogenomics/logs"

mapfile -t POPS < <(ls "$PD"/universes/*_universe_tx.txt | xargs -n1 basename | sed 's/_universe_tx.txt//' | sort)
POP="${POPS[$SLURM_ARRAY_TASK_ID]:?no population at index $SLURM_ARRAY_TASK_ID}"
OUT="$PD/db${TAG}/$POP"; mkdir -p "$OUT"
echo "build DBs population=$POP tag=${TAG:-<none>} orfs_tag=${ORFS_TAG:-<none>} thresh=$THRESH host=$(hostname) start=$(date -Iseconds)"

"$PY" "$NEW/proteogenomics/scripts/build_a549_dbs.py" \
  --candidates "$PD/orfs${ORFS_TAG}/$POP/candidates.faa" \
  --canonical "$CANON" \
  --thresh "$THRESH" --min_aa 7 \
  --out_dir "$OUT"

# guard: model DB must be canonical + a non-trivial selected-novel set (else threshold/parse broke)
mseq=$(grep -c '^>' "$OUT/db_model.fasta"); cseq=$(grep -c '^>' "$OUT/db_canonical.fasta")
echo "DONE DBs $POP: canonical=$cseq model=$mseq null=$(grep -c '^>' "$OUT/db_null.fasta") seqs $(date -Iseconds)"
[ "$mseq" -gt "$cseq" ] || { echo "ABORT $POP: model DB not larger than canonical (no novel selected)"; exit 1; }
