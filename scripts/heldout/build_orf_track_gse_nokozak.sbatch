#!/usr/bin/env bash
#SBATCH --job-name=orf_gse_nk
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=0-1
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orf/orf_gse_nk_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orf/orf_gse_nk_%A_%a.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
mkdir -p $NEW/logs/orf
CONDS=(nt lps); COND=${CONDS[$SLURM_ARRAY_TASK_ID]}
PACK=$NEW/data/packed_heldout_mouse_gse120762_${COND}
FASTA=$NEW/proteogenomics/data/gse120762/gse120762_${COND}_universe.fa
# NO-KOZAK track for the GSE universe, to match the mm1 no-kozak model's training input.
$PY $NEW/scripts/build_orf_track.py --mode ext --kozak none --pack "$PACK" --fasta "$FASTA" \
  --out "$PACK/orf_track_v2_nokozak.npy"
echo "no-kozak orf track $COND -> $PACK/orf_track_v2_nokozak.npy"
