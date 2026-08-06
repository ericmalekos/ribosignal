#!/usr/bin/env bash
#SBATCH --job-name=orf_gse
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orf/orf_gse_%x_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orf/orf_gse_%x_%j.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
mkdir -p $NEW/logs/orf
PACK=$NEW/data/packed_heldout_mouse_gse120762_${COND:?}
FASTA=$NEW/proteogenomics/data/gse120762/gse120762_${COND}_universe.fa
# same track as the deployed model: --mode ext --kozak heuristic, aligned to this pack's tx
$PY $NEW/scripts/build_orf_track.py --mode ext --kozak heuristic --pack "$PACK" --fasta "$FASTA" --out "$PACK/orf_track_v2.npy"
echo "orf track $COND -> $PACK/orf_track_v2.npy"
