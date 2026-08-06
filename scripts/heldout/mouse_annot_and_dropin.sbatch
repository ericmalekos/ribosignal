#!/usr/bin/env bash
#SBATCH --job-name=ho_mouse_annot
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/mouse_annot_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/mouse_annot_%j.err
set -euo pipefail
# Build the vM38 RiboCode annotation DB that the mouse cross-species drop-in needs (never persisted),
# then run the drop-in (3 density variants) + compare against the official mouse RiboCode calls on the
# mouse GPU dump's pred_profiles.npz. The mouse_wang_liver_collapsed.txt official calls were produced
# by RiboCode from vM38 + GRCm39, so rebuilding the annot from the same GTF/genome reproduces the same
# transcript and genomic coordinates the --key genomic compare relies on.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
PYR=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
PREP=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/prepare_transcripts
GTF=$BASE/annotations/gencode.vM38.annotation.gtf
GENOME=/private/groups/carpenterlab/emalekos/genomes/GRCm39.primary_assembly.genome.fa
ANNOT=$NEW/data/mouse_ribocode_annot
DATASET=mouse_wang_liver
DROPIN=$NEW/results/heldout/$DATASET/onehot/dropin
PROFILES=$DROPIN/pred_profiles.npz
OFFICIAL=$NEW/data/heldout_psites/$DATASET/${DATASET}_collapsed.txt
export PYTHONNOUSERSITE=1
echo "host=$(hostname) start=$(date -Iseconds)"

# 0) Ensure a samtools/pysam .fai sidecar for the genome (prepare_transcripts reads per-tx via faidx).
if [ ! -f "$GENOME.fai" ]; then
  echo "-- faidx $GENOME --"
  "$PYR" -c "import pysam,sys; pysam.faidx(sys.argv[1])" "$GENOME"
fi

# 1) Build the vM38 RiboCode annotation DB (idempotent: skip if transcripts.pickle already present).
if [ -f "$ANNOT/transcripts.pickle" ]; then
  echo "annot exists, skip prepare_transcripts: $ANNOT"
else
  echo "-- prepare_transcripts (vM38 + GRCm39) -> $ANNOT --"
  "$PREP" -g "$GTF" -f "$GENOME" -o "$ANNOT"
fi
[ -f "$ANNOT/transcripts.pickle" ] || { echo "ERROR: annot build failed, no $ANNOT/transcripts.pickle"; exit 1; }

# 2) The mouse GPU dump (separate fast job) writes pred_profiles.npz. Wait bounded for it (prepare
#    takes ~20 min, the dump takes a few, so it is normally already present by the time we get here).
for _ in $(seq 1 30); do [ -f "$PROFILES" ] && break; echo "waiting for $PROFILES ..."; sleep 30; done
[ -f "$PROFILES" ] || { echo "ERROR: no $PROFILES (mouse dump 35335047 not done); rerun this job after it lands"; exit 1; }

# 3) RiboCode drop-in (real / predicted-shape x obs-depth / x pred-depth) + compare vs official.
for V in real pred_obsdepth pred_preddepth; do
  echo "-- variant $V --"
  "$PYR" "$NEW/scripts/ribocode_dropin.py" \
    --profiles "$PROFILES" --annot "$ANNOT" --variant "$V" --out "$DROPIN" --min_aa 5 --pval 0.05
done
echo "-- compare --"
"$PYR" "$NEW/scripts/compare_dropin_calls.py" \
  --dropin_dir "$DROPIN" --official "$OFFICIAL" --profiles "$PROFILES" \
  --key genomic --min_len 90 --min_enrichment 0.5 --out "$DROPIN" \
  --tx2gene "$NEW/data/mouse_tx2biotype.tsv"
echo "done $(date -Iseconds)"
