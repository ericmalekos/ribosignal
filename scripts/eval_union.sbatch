#!/usr/bin/env bash
#SBATCH --job-name=eval_union
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=140G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/eval_union_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/eval_union_%j.err
set -euo pipefail
# Predict on Hepatocytes holdout for an FM-embedding variant: dump per-nt profiles (GPU; the trained model
# reads the pre-computed FM embeddings via tx_index, so the rinalmo SIF/torch works for every backend) then
# RiboCode detectORF on real / pred_obsdepth / pred_preddepth (ribocode env).
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
ORFTRACK=$NEW/data/packed_union/orf_track_v2.npy
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
BACKEND=${BACKEND:?set BACKEND=onehot|rinalmo|orthrus|hydrarna}
MIXER=${MIXER:-transformer}                           # match the training: transformer or mamba (mamba dump needs the orthrus SIF)
N_ATTN=${N_ATTN:-2}                                   # match the training: 2 (default) or 4 (attn4)
D_STATE=${D_STATE:-16}                                # match the training: mamba SSM state size (16 default, or 32/64)
case "$MIXER" in
  transformer) NAME=orf_v2_attn; [ "$N_ATTN" = 2 ] || NAME=orf_v2_attn${N_ATTN}; ARCHLBL=""; [ "$N_ATTN" = 2 ] || ARCHLBL="attn${N_ATTN} "; SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba)       NAME=orf_v2_mamba; ARCHLBL="mamba"; [ "$N_ATTN" = 2 ] || { NAME=orf_v2_mamba${N_ATTN}; ARCHLBL="mamba${N_ATTN}"; }; [ "$D_STATE" = 16 ] || { NAME=${NAME}_ds${D_STATE}; ARCHLBL="${ARCHLBL} ds${D_STATE}"; }; ARCHLBL="${ARCHLBL} "; SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img ;;
  *) echo "bad MIXER=$MIXER (want transformer|mamba)" >&2; exit 1 ;;
esac
RUN=results/loto/${NAME}_${BACKEND}_union_noBrain_nokozak_mm1_holdout_Hepatocytes
RIBOPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
ANNOT=$ECH/data/ribocode_annot
mkdir -p "$NEW/logs/dropin"
echo "eval BACKEND=$BACKEND run=$RUN host=$(hostname) start=$(date -Iseconds)"

# DUMP_DEVICE=cpu runs the (only) GPU step on idle CPUs -- the 5M-param model is CPU-deployable, so this
# bypasses a jammed GPU queue (submit with --partition=gpu --gres=none). --nv is dropped for CPU.
DUMP_DEVICE=${DUMP_DEVICE:-cuda}
[ "$DUMP_DEVICE" = cpu ] && NVFLAG="" || NVFLAG="--nv"
# shellcheck disable=SC2086
singularity exec $NVFLAG --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True --env RIBO_PACK_SUFFIX=union --env RIBO_ORF_TRACK="$ORFTRACK" --env RIBO_ONEHOT_FASTA="$NEW/data/union_universe.fa" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" \
    --fasta "$NEW/data/union_universe.fa" --device "$DUMP_DEVICE"

# guard: a heavily-skipped dump = wrong onehot FASTA / pack (silently biases the eval to a subset).
# Fail loudly BEFORE the drop-in + doc append, rather than publishing wrong numbers.
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 -c "
import numpy as np,sys
d=np.load('$NEW/$RUN/dropin/pred_profiles.npz',allow_pickle=True)
nt=len(d['tx_ids']); ns=int(str(d['meta'][3])); frac=ns/(nt+ns) if (nt+ns) else 0
print(f'dump: {nt} tx, {ns} skipped ({frac:.1%})')
sys.exit(1 if frac>0.05 else 0)" || { echo 'ABORT: >5% tx skipped in dump (FASTA/pack mismatch)'; exit 1; }

export PYTHONNOUSERSITE=1
for V in real pred_obsdepth pred_preddepth; do
  echo "=== dropin $V ==="
  "$RIBOPY" "$NEW/scripts/ribocode_dropin.py" \
    --profiles "$NEW/$RUN/dropin/pred_profiles.npz" --annot "$ANNOT" \
    --variant "$V" --out "$NEW/$RUN/dropin" --min_aa 5 --pval 0.05
done

# --- two-arm Poisson CDS-anchored sweep + full standing-rule ORF-call metrics -> manuscript doc ---
SW=$NEW/$RUN/dropin_sweep_poisson
for TH in 0.02 0.05 0.1 0.2; do
  O=$SW/theta_$TH; rm -rf "$O"; mkdir -p "$O"
  "$RIBOPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$NEW/$RUN/dropin/pred_profiles.npz" \
    --annot "$ANNOT" --variant pred_preddepth --pred_scale "$TH" --pred_poisson --seed 0 \
    --out "$O" --min_aa 5 --pval 0.05
done
CASPY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
mkdir -p "$NEW/results/orf_call_metrics"
"$CASPY" "$NEW/scripts/heldout/orf_call_metrics.py" --run "$NEW/$RUN" --label "${ARCHLBL}${BACKEND} union" \
  --universe union --poisson-sweep "$SW" --cds-recall 0.90 \
  --out-json "$NEW/results/orf_call_metrics/${NAME}_${BACKEND}_union.json" \
  --append-md "$NEW/docs/manuscript_orf_call_eval.md"
echo "eval BACKEND=$BACKEND done=$(date -Iseconds)"
