#!/usr/bin/env bash
#SBATCH --job-name=msf_line
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --array=0-2
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/msf_line_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/msf_line_%A_%a.err
set -euo pipefail
# Generalized immunopeptidome HLA search. Pass via --export: PILOT (the <LINE>_pilot dir).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
JAR=/private/groups/carpenterlab/emalekos/conda_envs/msfragger/share/msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar
JAVA=/private/groups/carpenterlab/emalekos/conda_envs/msfragger/bin/java
DBS=(canonical model null); D=${DBS[$SLURM_ARRAY_TASK_ID]}
# TEMPLATE: default reproduces the original (calibrate_mass=2) runs. Use the FROZEN template for any
# cross-database comparison -- calibrate_mass=2 re-derives search params PER DATABASE, and a
# divergence test on HBL-1 showed the null arm scored with use_topN_peaks=100 / intensity_transform=0
# against the model arm's 150 / 1, cutoff 16.742 vs 18.483. See proteogenomics/methods.md.
#   TEMPLATE=fragger_hbl1_hla_frozen.params DBSUB=_frozen sbatch ...
TEMPLATE=${TEMPLATE:-fragger_hbl1_hla.params}
DBSUB=${DBSUB:-}
DB=$PILOT/db${DBSUB}/db_${D}.fasta
OUT=$PILOT/search${DBSUB}_${D}; mkdir -p "$OUT"
echo "TEMPLATE=$TEMPLATE  DBDIR=db${DBSUB}  OUT=$(basename "$OUT")"
SCR=/data/tmp/emalekos/msfline_$(basename $PILOT)_${D}_${SLURM_JOB_ID}; mkdir -p "$SCR"
echo "PILOT=$(basename $PILOT) DB=$D host=$(hostname) start=$(date -Iseconds)"
[[ -s "$DB" ]] || { echo "missing DB $DB" >&2; exit 1; }
awk -v db="$DB" '/^database_name[[:space:]]*=/{print "database_name = " db; next}{print}' \
  $NEW/proteogenomics/msfragger/$TEMPLATE > "$SCR/fragger.params"
cp "$PILOT"/immunopeptidome/*.mzML "$SCR/"
cd "$SCR"
$JAVA -Xmx110g -jar "$JAR" "$SCR/fragger.params" "$SCR"/*.mzML
cp "$SCR"/*.pin "$SCR"/*.tsv "$SCR"/*.pepXML "$OUT"/ 2>/dev/null || true
rm -rf "$SCR"; echo "DB=$D done=$(date -Iseconds); $(ls $OUT/*.pin 2>/dev/null|wc -l) mzML searched"