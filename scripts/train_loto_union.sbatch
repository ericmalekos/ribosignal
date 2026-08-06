#!/usr/bin/env bash
#SBATCH --job-name=loto_union
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --time=24:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_union_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/loto_union_%j.err
set -euo pipefail
# Broader UNION-universe (84,472 tx) training. Same recipe as the Fibroblast-universe sweep (nokozak + mm1 +
# noBrain, holdout Hepatocytes, one-hot/rinalmo/orthrus/hydrarna) -- ONLY the universe differs. RIBO_PACK_SUFFIX
# routes tissue_pack_dir to data/packed_union*; the ORF track + one-hot FASTA are the union ones. FM backends
# need the union FM embeddings (extracted separately) before they can run; onehot needs only the union FASTA.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ORFTRACK=$NEW/data/packed_union/orf_track_v2.npy
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
BACKEND=${BACKEND:?set BACKEND=onehot|rinalmo|orthrus|hydrarna}
MIXER=${MIXER:-transformer}                           # transformer (default) or mamba (needs the orthrus SIF for mamba_ssm)
N_ATTN=${N_ATTN:-2}                                   # transformer layers: 2 (default) or 4 (attn4, the arch-sweep winner)
D_STATE=${D_STATE:-16}                                # mamba SSM state size (16 default); 32/64 add long-range capacity WITHOUT depth
SEED=${SEED:-0}                                       # 0 = base run (no suffix); 1/2 = seed replicates for the seed-confirm sweep
case "$MIXER" in
  transformer) ARCHARGS="--n_attn_layers $N_ATTN"; NAME=orf_v2_attn; [ "$N_ATTN" = 2 ] || NAME=orf_v2_attn${N_ATTN}; SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img ;;
  mamba)       ARCHARGS="--n_attn_layers $N_ATTN --mixer mamba --mamba_d_state $D_STATE"; NAME=orf_v2_mamba; [ "$N_ATTN" = 2 ] || NAME=orf_v2_mamba${N_ATTN}; [ "$D_STATE" = 16 ] || NAME=${NAME}_ds${D_STATE}; SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img ;;
  *) echo "bad MIXER=$MIXER (want transformer|mamba)" >&2; exit 1 ;;
esac
OUT=$NEW/results/loto/${NAME}_${BACKEND}_union_noBrain_nokozak_mm1_holdout_Hepatocytes
[ "$SEED" = 0 ] || OUT=${OUT}_seed${SEED}             # seed0 keeps the base name; seed1/2 get a suffix
mkdir -p "$NEW/logs/train" "$OUT"
echo "UNION BACKEND=$BACKEND MIXER=$MIXER name=$NAME ($ARCHARGS) SIF=$(basename "$SIF") holdout=Hepatocytes host=$(hostname) start=$(date -Iseconds)"

# shellcheck disable=SC2086
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_PACK_SUFFIX=union \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  --env RIBO_ONEHOT_FASTA="$NEW/data/union_universe.fa" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train_loto.py" \
    --out "$OUT" \
    --holdout Hepatocytes \
    --train_tissues Fibroblast,VSMC,ES,Fat,HA_EC,HCAEC,HUVEC \
    --emb_backend "$BACKEND" --cov_norm global_mean \
    --min_train_signal 50 --val_fold 0 \
    --max_tx_per_tissue 6000 \
    --resume \
    --use_orf_track $ARCHARGS \
    --epochs 30 --patience 6 --warmup 2 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget 16000 --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed "$SEED"
echo "train_done=$(date -Iseconds) -> $OUT"
