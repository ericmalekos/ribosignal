#!/usr/bin/env bash
#SBATCH --job-name=msf_a549
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --time=12:00:00
#SBATCH --array=0-2
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/msf_full_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/msf_full_%A_%a.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
JAR=/private/groups/carpenterlab/emalekos/conda_envs/msfragger/share/msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar
JAVA=/private/groups/carpenterlab/emalekos/conda_envs/msfragger/bin/java
DBS=(canonical model null)
D=${DBS[$SLURM_ARRAY_TASK_ID]}
DB=$NEW/proteogenomics/data/A549_pilot/${DBDIR:-db}/db_${D}.fasta
OUT=$NEW/proteogenomics/data/A549_pilot/search${TAG:-}_${D}; mkdir -p "$OUT"
SCR=/data/tmp/emalekos/msf_${D}_${SLURM_JOB_ID}; mkdir -p "$SCR"
echo "DB=$D host=$(hostname) start=$(date -Iseconds)"
awk -v db="$DB" '/^database_name[[:space:]]*=/{print "database_name = " db; next}{print}' \
  $NEW/proteogenomics/msfragger/fragger_a549_tmt.params > "$SCR/fragger.params"
cp $NEW/proteogenomics/data/A549_pilot/mzml/*.mzML "$SCR/"
cd "$SCR"
$JAVA -Xmx140g -jar "$JAR" "$SCR/fragger.params" "$SCR"/*.mzML
cp "$SCR"/*.pin "$SCR"/*.tsv "$SCR"/*.pepXML "$OUT"/ 2>/dev/null || true   # keep pepXML for MS2Rescore
rm -rf "$SCR"
echo "DB=$D done=$(date -Iseconds); result tsv/pin in $OUT ($(ls $OUT/*.pin 2>/dev/null|wc -l) fractions)"
