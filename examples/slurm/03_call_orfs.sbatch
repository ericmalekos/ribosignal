#!/usr/bin/env bash
#SBATCH --job-name=ribosignal_calls
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail
# Both callers, both arms. `real` is the reference; `pred_preddepth` calls from predicted shape
# and depth, so no observed P-site enters the call.
: "${RIBO_IMAGE:?}"; : "${RIBO_WORK:?}"; : "${RIBO_REF:=$RIBO_WORK/ref}"
ARCH=${ARCH:-mamba4}
T=${SLURM_CPUS_PER_TASK:-8}
S() { singularity exec --cleanenv -B "$RIBO_WORK" -B "$RIBO_REF" "$RIBO_IMAGE" "$@"; }
P=$RIBO_WORK/pred/$ARCH/pred_profiles.npz

for V in real pred_preddepth; do
  S python /opt/riboseq-model/scripts/ribocode_dropin.py --profiles "$P" \
    --annot "$RIBO_WORK/calls/annot" --variant "$V" --out "$RIBO_WORK/calls/$ARCH" \
    --min_aa 30 --pval 0.05
done
S python /opt/riboseq-model/scripts/ribocode_dropin.py --profiles "$P" \
  --annot "$RIBO_WORK/calls/annot" --variant pred_preddepth --pred_poisson --pred_scale 0.05 \
  --out "$RIBO_WORK/calls/${ARCH}_poisson" --min_aa 30 --pval 0.05

# Ribo-TISH needs a GTF restricted to the predicted transcripts (whole genes only);
# --emit-gtf writes it. Building it from pack_meta.tsv is wrong: the pack holds every transcript.
S python /opt/riboseq-model/scripts/orfcallers/ribotish_dropin.py --profiles "$P" \
  --gtf "$RIBO_REF/annotation.gtf" --emit-gtf "$RIBO_WORK/calls/scored_$ARCH.gtf" \
  --genome "$RIBO_REF/genome.fa" \
  --variant pred_preddepth --out "$RIBO_WORK/calls/ribotish_$ARCH" \
  --longest --minaalen 5 --fpth 0.05 --numproc "$T"
