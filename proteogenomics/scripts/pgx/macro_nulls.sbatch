#!/usr/bin/env bash
#SBATCH --job-name=pgx_nulls
#SBATCH --cpus-per-task=4
#SBATCH --mem=110G
#SBATCH --time=2:00:00
#SBATCH --partition=medium
#SBATCH --output=%x_%A_%a.out
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PD=$NEW/proteogenomics/data/macrophage_tissue
X=$PD/pgx_xsubtype
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
mapfile -t POPS < <(awk -F'\t' '{print $1}' "$PD/pop_mzml_map.tsv" | sort)
POP="${POPS[$SLURM_ARRAY_TASK_ID]:?}"
[ -s "$X/db_null/$POP/db_null_nc.fasta" ] && { echo "SKIP $POP (already built)"; exit 0; }
echo "building nulls for $POP host=$(hostname) start=$(date -Iseconds)"
cd "$NEW/proteogenomics/scripts"
"$PY" -m pgx.build_dbs --species mouse --out "$X/db_null/$POP" \
  --universe-fa "$PD/universes/${POP}_universe.fa" --dbs null_atg,null_nc
echo "done $POP $(date -Iseconds)"
