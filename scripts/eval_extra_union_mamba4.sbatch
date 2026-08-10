#!/usr/bin/env bash
#SBATCH --job-name=eval_extra_m4
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/eval_extra_m4_%j.out
set -euo pipefail
# eval_extra on the released MAMBA4 model. Companion to eval_extra_union.sbatch (which does attn).
#
# WHY THIS IS NEEDED. Locked decision D1b (2026-07-31) makes `orf_v2_mamba4` the MAIN Figure 1 model,
# with attn moving to supplemental but staying shipped as the CPU inference path. Every Fig 1 panel
# was nevertheless still built on attn. Repointing them to mamba4 needs its per-stratum metrics, and
# the mamba4 released run has `test_metrics.json` and `localization_metrics.json` but no
# `extra_metrics.json` -- which is exactly what figure B4 reads.
#
# GPU-only, unavoidably: mamba_ssm dispatches to causal_conv1d_cuda with no CPU fallback in the
# installed build, so unlike the attn arms this cannot be diverted to a CPU partition when the gpu
# queue is busy.
#
# Env matches train_loto_union.sbatch exactly (pack suffix + union ORF track + union FASTA move
# together, all three naming the same universe). The union track is kozak=none despite the plain
# filename; both facts are asserted below rather than assumed, because the Fibroblast-universe track
# is ALSO kozak=none and would pass a kozak-only check while mis-aligning every ORF.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img
RUN=$NEW/results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
ORFTRACK=$NEW/data/packed_union/orf_track_v2.npy
FASTA=$NEW/data/union_universe.fa

for f in "$ORFTRACK" "$FASTA" "$RUN/best.pt"; do
  [[ -s "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
META=${ORFTRACK%.npy}_meta.json
[[ -s "$META" ]] && { grep -q '"kozak": *"none"' "$META" || \
  { echo "ABORT: $META is not kozak=none" >&2; exit 1; }; }
WANT_BYTES=2126471848
GOT_BYTES=$(stat -c%s "$ORFTRACK")
[[ "$GOT_BYTES" -eq "$WANT_BYTES" ]] || { echo "ABORT: $ORFTRACK is $GOT_BYTES bytes, expected \
$WANT_BYTES (the union track)." >&2; exit 1; }
echo "run=$RUN track=$GOT_BYTES B host=$(hostname) start=$(date -Iseconds)"

singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_PACK_SUFFIX=union \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  --env RIBO_ONEHOT_FASTA="$FASTA" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/eval_extra.py" --run "$RUN" --eval_cap 0 --device cuda

OUT=$RUN/extra_metrics.json
[[ -s "$OUT" ]] || { echo "ABORT: no $OUT written" >&2; exit 1; }
echo "done=$(date -Iseconds) -> $OUT"
