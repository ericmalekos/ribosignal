#!/usr/bin/env bash
#SBATCH --job-name=msf_hbl1
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --array=0-2
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/msf_hbl1_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/msf_hbl1_%A_%a.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
JAR=/private/groups/carpenterlab/emalekos/conda_envs/msfragger/share/msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar
JAVA=/private/groups/carpenterlab/emalekos/conda_envs/msfragger/bin/java
DBS=(canonical model null price ribocode); D=${DBS[$SLURM_ARRAY_TASK_ID]}
HB=$NEW/proteogenomics/data/HBL1_pilot
DB=$HB/db/db_${D}.fasta
OUT=$HB/search_${D}; mkdir -p "$OUT"
SCR=/data/tmp/emalekos/msfhbl1_${D}_${SLURM_JOB_ID}; mkdir -p "$SCR"
echo "DB=$D host=$(hostname) start=$(date -Iseconds)"
awk -v db="$DB" '/^database_name[[:space:]]*=/{print "database_name = " db; next}{print}' \
  $NEW/proteogenomics/msfragger/fragger_hbl1_hla.params > "$SCR/fragger.params"
cp "$HB"/immunopeptidome/*.mzML "$SCR/"   # ThermoRawFileParser-converted (msfragger conda build lacks the Thermo .raw reader)
cd "$SCR"
$JAVA -Xmx110g -jar "$JAR" "$SCR/fragger.params" "$SCR"/*.mzML
cp "$SCR"/*.pin "$SCR"/*.tsv "$SCR"/*.pepXML "$OUT"/ 2>/dev/null || true
rm -rf "$SCR"; echo "DB=$D done=$(date -Iseconds); $(ls $OUT/*.pin 2>/dev/null|wc -l) raw searched"
