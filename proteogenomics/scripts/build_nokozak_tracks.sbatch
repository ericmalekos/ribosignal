#!/usr/bin/env bash
# Build --kozak none ORF tracks for the MS pilot universes.
#
# WHY. The released models (orf_v2_{mamba4,attn}_onehot_union_noBrain_nokozak_mm1) are trained on a
# no-Kozak ORF track. The pilot packs only carry kozak=heuristic tracks, built for the OLD deployed
# model. Feeding a nokozak-trained model the heuristic track is a train/inference mismatch, so the
# re-run needs these first.
#
#   sbatch --array=0-3 proteogenomics/scripts/build_nokozak_tracks.sbatch
#SBATCH --job-name=nokozak_track
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=6:00:00
#SBATCH --partition=medium
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/nokozak_track_%A_%a.out
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
PACKS=(human_hbl1 human_dohh2 human_sudhl4 human_a549_fresh)
FA=(HBL1_pilot/HBL1_universe.fa DoHH2_pilot/DoHH2_universe.fa SUDHL4_pilot/SUDHL4_universe.fa A549_pilot/A549_universe.fa)
P=${PACKS[$SLURM_ARRAY_TASK_ID]}; F=${FA[$SLURM_ARRAY_TASK_ID]}
PACK=$NEW/data/packed_heldout_$P
FASTA=$NEW/proteogenomics/data/$F
OUT=$PACK/orf_track_v2_nokozak.npy
echo "pack=$P fasta=$(basename $FASTA) host=$(hostname) start=$(date -Iseconds)"
[[ -d "$PACK" ]]  || { echo "ABORT: no pack $PACK" >&2; exit 1; }
[[ -s "$FASTA" ]] || { echo "ABORT: no fasta $FASTA" >&2; exit 1; }
if [[ "${FORCE:-0}" != "1" && -s "$OUT" ]]; then echo "SKIP (exists)"; exit 0; fi
"$PY" "$NEW/scripts/build_orf_track.py" --mode ext --kozak none --pack "$PACK" --fasta "$FASTA" --out "$OUT"
# name must match content: assert the meta really says kozak=none
META=${OUT%.npy}_meta.json
[[ -s "$META" ]] && grep -q '"kozak": *"none"' "$META" || { echo "ABORT: $META does not say kozak=none" >&2; exit 1; }
echo "DONE $P: $(du -h "$OUT" | cut -f1)  $(date -Iseconds)"
