#!/usr/bin/env bash
#SBATCH --job-name=macro_msf
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --time=12:00:00
#SBATCH --array=0-35%6
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/macro_msf_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/macro_msf_%A_%a.err
set -euo pipefail
# Per-population MSFragger search. task/NDB -> population (from pop_mzml_map.tsv), task%NDB -> DB.
#
# DBLIST selects which DBs to search and therefore the required array size (12 * NDB - 1):
#   DBLIST=canonical,model,null  -> --array=0-35   (default, full rebuild)
#   DBLIST=model                 -> --array=0-11   (checkpoint swap: canonical is model-independent, and
#                                                   db_null's SEQUENCE set is too -- enumeration is pure
#                                                   over the universe FASTA, only the class map's f0
#                                                   column moves. VERIFY the null FASTA is unchanged
#                                                   before reusing its search, do not assume it.)
# TAG namespaces DB + search dirs per checkpoint (see macro_dumps.sbatch). Each task searches ONE population's valid mzML
# against ONE DB with the label-free tryptic params (fragger_macro_lfq.params: TMT stripped, Cys-carbamidomethyl
# fixed, ox-M + acetyl-protein-Nterm variable, 20 ppm Orbitrap precursor, 0.6 Da ion-trap fragment -- matches
# the Orbitrap-Fusion-Lumos ITMS/HCD acquisition). The 2 TRFP-unparseable mzML stubs are skipped by the size gate.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PD=$NEW/proteogenomics/data/macrophage_tissue
JAR=/private/groups/carpenterlab/emalekos/conda_envs/msfragger/share/msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar
JAVA=/private/groups/carpenterlab/emalekos/conda_envs/msfragger/bin/java
PARAMS=$NEW/proteogenomics/msfragger/fragger_macro_lfq.params
IFS=',' read -ra DBS <<< "${DBLIST:-canonical,model,null}"
NDB=${#DBS[@]}
TAG=${TAG:-}
mkdir -p "$NEW/proteogenomics/logs"

mapfile -t POPS < <(cut -f1 "$PD/pop_mzml_map.tsv" | sort)
POP="${POPS[$((SLURM_ARRAY_TASK_ID / NDB))]:?no population for task $SLURM_ARRAY_TASK_ID}"
D="${DBS[$((SLURM_ARRAY_TASK_ID % NDB))]}"
PREF=$(awk -v s="$POP" '$1==s{print $2}' "$PD/pop_mzml_map.tsv")
DB="$PD/db${TAG}/$POP/db_${D}.fasta"
OUT="$PD/search${TAG}/$POP/$D"; mkdir -p "$OUT"
SCR=/data/tmp/emalekos/macromsf_${POP}_${D}_${SLURM_JOB_ID}; mkdir -p "$SCR"
echo "pop=$POP db=$D tag=${TAG:-<none>} ndb=$NDB prefix=$PREF host=$(hostname) start=$(date -Iseconds)"
[ "$(ls "$OUT"/*.tsv 2>/dev/null | wc -l)" -ge 17 ] && { echo "SKIP $POP/$D (already searched)"; exit 0; }
[ -s "$DB" ] || { echo "ABORT: missing DB $DB (DB-build step incomplete?)"; exit 1; }

# stage COMPLETE mzML only: >1MB (excludes the <1MB TRFP vendor-bug stubs) AND a proper closing tag
# (excludes OOM-truncated files -- MSFragger aborts the whole batch on any corrupt file, so a size gate
# alone is not enough; a truncated >1MB file must be caught here or it fails the entire population).
n=0
for m in "$PD"/mzml/${PREF}_Rep*_F*.mzML; do
  [ "$(stat -c%s "$m" 2>/dev/null || echo 0)" -gt 1000000 ] || { echo "skip stub $(basename "$m")"; continue; }
  tail -c 400 "$m" 2>/dev/null | grep -qE "</indexedmzML>|</mzML>" || { echo "skip truncated $(basename "$m")"; continue; }
  cp "$m" "$SCR/"; n=$((n + 1))
done
echo "staged $n mzML for $POP"
[ "$n" -gt 0 ] || { echo "ABORT: no mzML staged for $POP (prefix $PREF)"; exit 1; }

awk -v db="$DB" '/^database_name[[:space:]]*=/{print "database_name = " db; next}{print}' "$PARAMS" > "$SCR/fragger.params"
cd "$SCR"
$JAVA -Xmx140g -jar "$JAR" "$SCR/fragger.params" "$SCR"/*.mzML
cp "$SCR"/*.pin "$SCR"/*.tsv "$SCR"/*.pepXML "$OUT"/ 2>/dev/null || true   # pepXML kept for MS2Rescore
rm -rf "$SCR"
echo "DONE pop=$POP db=$D $(date -Iseconds): $(ls "$OUT"/*.pin 2>/dev/null | wc -l) fraction pin in $OUT"
