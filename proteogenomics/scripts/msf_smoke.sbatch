#!/usr/bin/env bash
#SBATCH --job-name=msf_a549_smoke
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/msf_smoke_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/msf_smoke_%j.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
JAR=/private/groups/carpenterlab/emalekos/conda_envs/msfragger/share/msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar
JAVA=/private/groups/carpenterlab/emalekos/conda_envs/msfragger/bin/java
DB=$NEW/proteogenomics/data/A549_pilot/db/db_canonical.fasta
SCR=/data/tmp/emalekos/msf_smoke_${SLURM_JOB_ID}; mkdir -p "$SCR"
awk -v db="$DB" '/^database_name[[:space:]]*=/{print "database_name = " db; next}{print}' \
  $NEW/proteogenomics/msfragger/fragger_a549_tmt.params > "$SCR/fragger.params"
cp $NEW/proteogenomics/data/A549_pilot/mzml/d00615_Prot_12_01.mzML $NEW/proteogenomics/data/A549_pilot/mzml/d00616_Prot_12_02.mzML "$SCR/"
echo "host=$(hostname) start=$(date -Iseconds)"
cd "$SCR"
$JAVA -Xmx110g -jar "$JAR" "$SCR/fragger.params" "$SCR"/*.mzML
echo "== PSM counts (target vs REV_ decoy) per fraction =="
for t in "$SCR"/*.tsv; do
  [ -f "$t" ] && echo "$(basename $t): $(tail -n +2 $t | wc -l) PSMs, decoy $(tail -n +2 $t | grep -c REV_ || true)"
done
# keep the pin/tsv results, drop the staged mzML
mkdir -p $NEW/proteogenomics/data/A549_pilot/search_smoke
cp "$SCR"/*.pin "$SCR"/*.tsv $NEW/proteogenomics/data/A549_pilot/search_smoke/ 2>/dev/null || true
rm -f "$SCR"/*.mzML
echo "done=$(date -Iseconds)"
