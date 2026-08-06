#!/usr/bin/env bash
#SBATCH --job-name=ho_gpu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/gpu_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/heldout/gpu_%j.err
set -euo pipefail
# Held-out predicted-profile dump + RiboCode drop-in, all in ONE self-contained GPU job.
# Replaces the fragile short-partition sharded-CPU chain (dump array -> afterok merge -> afterok
# finalize): a single timed-out shard poisoned the afterok dependency and killed the whole finalize.
# On a GPU the O(L^2) transformer attention that blew the 55-min short wall runs in minutes, so we
# dump the full test set un-sharded (nshards=1 -> dump_pred_profiles.py writes pred_profiles.npz
# directly, no merge), then run the 3 density variants + compare in-place. DATASET/BACKEND via --export;
# mouse overrides ORFTRACK/FASTA_OVERRIDE/ANNOT.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
PYR=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
DATASET=${DATASET:?set DATASET=<name>}
BACKEND=${BACKEND:?set BACKEND=onehot|orthrus|rinalmo}
ORFTRACK=${ORFTRACK:-$NEW/data/packed/orf_track_v2.npy}
FASTA_OVERRIDE=${FASTA_OVERRIDE:-}
ANNOT=${ANNOT:-$ECH/data/ribocode_annot}
OFFICIAL=$NEW/data/heldout_psites/$DATASET/${DATASET}_collapsed.txt
case "$BACKEND" in
  onehot)  RUN=results/loto/orf_v2_attn_onehot_holdout_Hepatocytes ;;
  orthrus) RUN=results/loto/orf_v2_attn_orthrus_holdout_Hepatocytes ;;
  rinalmo) RUN=results/loto/orf_v2_attn_holdout_Hepatocytes ;;
  *) echo "bad BACKEND $BACKEND"; exit 1 ;;
esac
OUT=$NEW/results/heldout/$DATASET/$BACKEND
DROPIN=$OUT/dropin
PROFILES=$DROPIN/pred_profiles.npz
mkdir -p "$DROPIN" "$NEW/logs/heldout"
export PYTHONNOUSERSITE=1
echo "GPU dump+finalize dataset=$DATASET backend=$BACKEND host=$(hostname) gpu=${CUDA_VISIBLE_DEVICES:-?} start=$(date -Iseconds)"

FASTA_ARG=(); ENV_FASTA=()
if [ -n "$FASTA_OVERRIDE" ]; then FASTA_ARG=(--fasta "$FASTA_OVERRIDE"); ENV_FASTA=(--env RIBO_ONEHOT_FASTA="$FASTA_OVERRIDE"); fi

# 1) GPU dump -- nshards defaults to 1, so this writes pred_profiles.npz directly (no shard/merge).
# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 --env OMP_NUM_THREADS=8 \
  --env RIBO_ORF_TRACK="$ORFTRACK" "${ENV_FASTA[@]}" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" \
    --heldout "$DATASET" --out "$DROPIN" "${FASTA_ARG[@]}" --device cuda
[ -f "$PROFILES" ] || { echo "dump failed, no $PROFILES"; exit 1; }
echo "dump OK: $PROFILES ($(stat -c%s "$PROFILES" 2>/dev/null) bytes) -- profile Pearson computable from here."

# 2) RiboCode drop-in: 3 density variants (real / predicted-shape x obs-depth / x pred-depth) + compare.
# Requires a prepared RiboCode annotation DB (transcripts.pickle). Human reuses the v49 fibroblast
# annot; mouse has no persisted annot yet, so the drop-in auto-skips (the dump/profile Pearson above
# still lands) until a vM38 annot is built and passed via ANNOT=.
if [ -f "$ANNOT/transcripts.pickle" ]; then
  for V in real pred_obsdepth pred_preddepth; do
    echo "-- variant $V --"
    "$PYR" "$NEW/scripts/ribocode_dropin.py" \
      --profiles "$PROFILES" --annot "$ANNOT" --variant "$V" --out "$DROPIN" --min_aa 5 --pval 0.05
  done
  echo "-- compare --"
  "$PYR" "$NEW/scripts/compare_dropin_calls.py" \
    --dropin_dir "$DROPIN" --official "$OFFICIAL" --profiles "$PROFILES" \
    --key genomic --min_len 90 --min_enrichment 0.5 --out "$DROPIN"
else
  echo "SKIP drop-in: no $ANNOT/transcripts.pickle (dump/profile-Pearson only for this dataset)"
fi
echo "done dataset=$DATASET backend=$BACKEND $(date -Iseconds)"
