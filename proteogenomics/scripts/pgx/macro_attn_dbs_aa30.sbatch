#!/usr/bin/env bash
#SBATCH --job-name=pgx_attn_db30
#SBATCH --cpus-per-task=4
#SBATCH --mem=110G
#SBATCH --time=3:00:00
#SBATCH --partition=medium
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/pgx_attn_db30_%A_%a.out
set -euo pipefail
# STAGE B of the attn macrophage sweep: build the ATTENTION model's per-population database at 30 AA
# and place it BESIDE the mamba one, so a single search can score both against identical spectra.
#
# STANDING RULE 5 (docs/PIPELINE_POLICY.md): the ORF-length floor is set by ASSAY. Macrophages are
# tryptic whole-cell lysate -- verified from their own fragger.params (search_enzyme_name_1 = trypsin,
# num_enzyme_termini = 2) -- so 30 AA. `--assay tryptic` SETS min_aa=30 and aborts if a conflicting
# --min-aa is passed, so the floor cannot drift back to 7 by an edit here.
#
# WHY BOTH MODELS AT 30 AND NOT attn-AT-30 ALONE. Comparing attn@30 against the existing mamba@7
# would change two variables at once and answer nothing. The mamba 30-AA databases already exist
# (job 36784385 length-filtered them into this same tree), so putting attn@30 next to them makes the
# model comparison valid at the cost of ONE extra database per population.
#
# DATABASE NAMING. mamba's arm is `db_model_predicted`; attn's is `db_model_predicted_attn`. Distinct
# names let both live in one DBDIR and be searched as separate arms of the same frozen run, which is
# what makes the comparison spectra-identical rather than merely parameter-matched.
#
# PER-POPULATION, NOT POOLED. Each subtype gets its own database from its OWN prediction -- 566 to
# 2,131 novel sequences across the 12 for mamba. That is the design P12 rests on. It also means
# cross-population "unique peptide" claims are confounded (a peptide can only be found where its
# database contains it); compare arms WITHIN a population, never populations to each other.
#
#   sbatch --array=0-11 proteogenomics/scripts/pgx/macro_attn_dbs_aa30.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PD=$NEW/proteogenomics/data/macrophage_tissue
SRC=$PD/pgx_macro_attn
X=$PD/pgx_xsubtype_aa30
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3

mapfile -t POPS < <(awk -F'\t' 'NF{print $1}' "$PD/pop_mzml_map.tsv" | sort -u)
POP="${POPS[$SLURM_ARRAY_TASK_ID]:?no population at index $SLURM_ARRAY_TASK_ID}"

CALLS=$SRC/$POP/calls/calls_poisson.tsv
EXTS=$SRC/$POP/extensions/poisson/extensions.tsv
for f in "$CALLS" "$EXTS"; do
  [[ -s "$f" ]] || { echo "ABORT: missing attn input $f (did stage A finish for $POP?)" >&2; exit 1; }
done
DST=$X/db_percall/$POP
[[ -d "$DST" ]] || { echo "ABORT: no 30-AA tree at $DST -- run macro_rebuild_aa30.sbatch first" >&2; exit 1; }
# The mamba 30-AA database must already be here, or the two arms would not be comparable.
[[ -s "$DST/db_model_predicted.fasta" ]] || { echo "ABORT: mamba 30-AA DB absent from $DST" >&2; exit 1; }

echo "pop=$POP arm=predicted_attn assay=tryptic(30aa) host=$(hostname) start=$(date -Iseconds)"
echo "  attn calls=$(( $(wc -l < "$CALLS") - 1 ))  extensions=$(( $(wc -l < "$EXTS") - 1 ))"

if [[ -s "$DST/db_model_predicted_attn.fasta" ]]; then
  echo "  skip (already built)"
else
  cd "$NEW/proteogenomics/scripts"
  "$PY" -m pgx.build_dbs --species mouse --assay tryptic --out "$DST" \
    --calls "$CALLS" --extensions "$EXTS" --arm predicted_attn --dbs model
fi

m=$(grep -c "^>nuORF" "$DST/db_model_predicted.fasta" 2>/dev/null || echo 0)
a=$(grep -c "^>nuORF" "$DST/db_model_predicted_attn.fasta" 2>/dev/null || echo 0)
echo "pop=$POP done=$(date -Iseconds)  novel@30aa: mamba=$m attn=$a"
