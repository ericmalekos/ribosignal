#!/usr/bin/env bash
# Dump predicted profiles over a gene_type universe (one model per array task).
#
# The gene_type universe adds NMD / retained_intron / TEC isoforms of coding genes, which the
# transcript_type dumps do NOT cover, so the model arms need fresh predictions before any pgx
# database can be built on it. Nulls do not: they depend only on the universe FASTA.
#
#   PACK=<...>/data/packed_heldout_macro_BMDM_gt \
#   FASTA=<...>/universes_genetype/BMDM_universe.fa \
#   TXLIST=<...>/universes_genetype/BMDM_universe_tx.txt \
#   OUTROOT=<...>/macrophage_tissue  POP=BMDM \
#     sbatch proteogenomics/scripts/pgx/dump_genetype.sbatch
#SBATCH --job-name=pgx_dump_gt
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --gres=gpu:1
#SBATCH --time=6:00:00
#SBATCH --partition=gpu
#SBATCH --output=%x_%A_%a.out
set -euo pipefail

: "${PACK:?}" "${FASTA:?}" "${TXLIST:?}" "${OUTROOT:?}" "${POP:?}"
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache

# mamba4 needs mamba_ssm (orthrus SIF); attn runs in either. One model per array index.
RUNS=(orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
      orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes)
TAGS=(mamba4 attn)
SIFS=("$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img"
      "$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img")

RUN=results/loto/${RUNS[$SLURM_ARRAY_TASK_ID]}
TAG=${TAGS[$SLURM_ARRAY_TASK_ID]}
SIF=${SIFS[$SLURM_ARRAY_TASK_ID]}
OUT=$OUTROOT/pred_genetype_${TAG}/$POP
mkdir -p "$OUT"

[ -s "$NEW/$RUN/best.pt" ] || { echo "ABORT: no checkpoint at $NEW/$RUN/best.pt" >&2; exit 1; }
[ -s "$SIF" ] || { echo "ABORT: no SIF at $SIF" >&2; exit 1; }
[ -s "$OUT/pred_profiles.npz" ] && { echo "SKIP $TAG (already dumped)"; exit 0; }
echo "model=$TAG pop=$POP host=$(hostname) sif=$(basename "$SIF") start=$(date -Iseconds)"

singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 --env OMP_NUM_THREADS=8 \
  --env RIBO_ORF_TRACK="$PACK/orf_track_v2_nokozak.npy" \
  --env RIBO_ONEHOT_FASTA="$FASTA" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" --heldout "macro_${POP}_gt" \
    --out "$OUT" --fasta "$FASTA" --tx_list "$TXLIST" --device cuda

# A heavy skip rate means the FASTA and the pack disagree -- fail loudly rather than silently
# reporting on a biased subset (see feedback_dump_silent_skip_onehot_fasta).
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 -c "
import numpy as np, sys
d = np.load('$OUT/pred_profiles.npz', allow_pickle=True)
nt = len(d['tx_ids']); ns = int(str(d['meta'][3])); frac = ns / (nt + ns) if (nt + ns) else 0
print(f'$TAG dump: {nt} tx, {ns} skipped ({frac:.1%})')
sys.exit(1 if frac > 0.05 else 0)" || { echo "ABORT $TAG: >5% tx skipped" >&2; exit 1; }
echo "done=$(date -Iseconds)"
