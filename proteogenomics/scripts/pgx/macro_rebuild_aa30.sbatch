#!/usr/bin/env bash
#SBATCH --job-name=pgx_aa30_db
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --partition=medium
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/pgx_aa30_db_%A_%a.out
set -euo pipefail
# STAGE 1 of the tryptic 30-AA rebuild: databases. One array task per macrophage population.
#
# STANDING RULE (docs/PIPELINE_POLICY.md, memory feedback_orf_min_length_by_assay): the ORF-length
# floor is set BY ASSAY.
#   tryptic whole-cell lysate -> 30 aa   |   MHC / HLA immunopeptidomics -> 7 aa
# The 12 macrophage populations are TRYPTIC -- verified below from the search parameters themselves
# (search_enzyme_name_1 = trypsin, num_enzyme_termini = 2), never from a hardcoded name list, so a
# population that is actually MHC cannot be silently stripped of its microproteins.
#
# METHOD: LENGTH-FILTER THE EXISTING 7-AA DATABASES, do not rebuild from calls. In pgx.build_dbs,
# min_aa is a PURE LENGTH FILTER on final protein sequences, so `--min-aa 30` is exactly the >= 30 aa
# subset of `--min-aa 7`. Filtering reuses the ORIGINAL call sets, which matters because
# db_percall/ carries no db_summary.json -- a from-scratch rebuild could silently substitute a
# different call set and confound the floor change with a provenance change. Verified on BMDM: the
# filtered novel-target counts (1,533->1,195, 699->506, 2,173->1,265, 235,686->81,250) reconcile
# exactly with the independent sub-30-aa counts from orf_length_audit.py.
#
# WRITES TO A NEW TREE (pgx_xsubtype_aa30/), so the 7-aa results stay intact and comparable.
#
# STAGE 2 is the re-search and is NOT optional: shrinking the database changes decoy competition,
# the score cutoff and rank-1 winners, so post-hoc filtering of the old results is not a substitute.
#
#   sbatch --array=0-11 proteogenomics/scripts/pgx/macro_rebuild_aa30.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PD=$NEW/proteogenomics/data/macrophage_tissue
OLD=$PD/pgx_xsubtype
X=$PD/pgx_xsubtype_aa30
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3

mapfile -t POPS < <(awk -F'\t' 'NF{print $1}' "$PD/pop_mzml_map.tsv" | sort -u)
POP="${POPS[$SLURM_ARRAY_TASK_ID]:?no population at index $SLURM_ARRAY_TASK_ID}"

PAR=$(ls "$OLD/frozen/$POP"/*/fragger.params 2>/dev/null | head -1)
[[ -s "$PAR" ]] || { echo "ABORT: no fragger.params for $POP -- cannot verify the assay" >&2; exit 1; }
ENZ=$(awk -F'= *' '/^search_enzyme_name_1/{print $2}' "$PAR" | tr -d ' \r')
TERM=$(awk -F'= *' '/^num_enzyme_termini/{print $2}' "$PAR" | tr -d ' \r')
if [[ "$ENZ" != "trypsin" || "$TERM" != "2" ]]; then
  echo "ABORT: $POP is enzyme=$ENZ termini=$TERM, not tryptic. The 30-aa floor is for TRYPTIC" >&2
  echo "       whole-cell lysate only; MHC data stays at 7 aa." >&2
  exit 1
fi
echo "pop=$POP assay=tryptic (enzyme=$ENZ termini=$TERM) -> min_aa=30 host=$(hostname) start=$(date -Iseconds)"

cd "$NEW/proteogenomics/scripts"
for SUB in db_percall db_null; do
  SRC=$OLD/$SUB/$POP
  [[ -d "$SRC" ]] || { echo "  (no $SUB/$POP -- skipped, not zeroed)"; continue; }
  "$PY" -m pgx.filter_db_min_aa --in "$SRC" --out "$X/$SUB/$POP" --min-aa 30 --assay tryptic
done
echo "pop=$POP done=$(date -Iseconds)"
