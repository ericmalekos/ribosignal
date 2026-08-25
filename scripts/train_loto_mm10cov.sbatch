#!/usr/bin/env bash
#SBATCH --job-name=loto_mm10cov
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --time=24:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/mm10cov_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/train/mm10cov_%j.err
set -euo pipefail
# mm10-RNA-coverage retrain, union universe, mamba4. Two arms, run SEPARATELY:
#
#   ARM=A  mm10 RNA coverage, FULL training set.
#          Fixes the model INPUT. mm1 RNA coverage is artifactually empty over multimapped
#          regions -- GTF2I ENST00000443166.5 reads 9% of CDS mean across a 1,617 nt window at
#          mm1 and 138% at mm10. A vs the mm1 baseline isolates the input effect.
#
#   ARM=B  mm10 RNA coverage, training set with multimap-corrupted tx REMOVED.
#          Fixes the training LABEL. The Ribo target stays mm1 everywhere (RiboCode does not drop
#          multimappers and they destroy periodicity), so tx whose mm1 Ribo label is missing most
#          of its reads teach an artifact. Those 7,640 tx (frac_discarded > 0.5) are dropped from
#          TRAIN and VAL and KEPT IN TEST. B vs A isolates the exclusion effect.
#
# A and B differ from each other in exactly one thing (the training set) and both differ from the
# mm1 baseline in exactly one other (the coverage channel). Test split is identical in all three.
#
# Baseline to compare against:
#   results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes
#   (best val pearson_median 0.6190 @ epoch 23; seeds 1,2 also on disk)
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img     # mamba_ssm lives here
# mm10cov packs reuse the union tx_order verbatim (prepare_pack.py --reuse-target), so the union
# ORF track and one-hot FASTA align 1:1 and are correct unchanged.
ORFTRACK=$NEW/data/packed_union/orf_track_v2.npy              # kozak:"none" (verified in its meta json)
EXCL=$NEW/data/mm25_diagnostic/train_exclude_mmcost_gt0.5.txt
ARM=${ARM:?set ARM=A or ARM=B}
SEED=${SEED:-0}
N_ATTN=${N_ATTN:-4}; D_STATE=${D_STATE:-16}                   # mamba4 = the standing mamba arm (FIG_MODEL default)
NAME=orf_v2_mamba${N_ATTN}_onehot_union_noBrain_nokozak_mm10cov
case "$ARM" in
  A) EXCLARGS=(); NAME=${NAME} ;;
  B) [ -s "$EXCL" ] || { echo "missing exclusion list $EXCL" >&2; exit 1; }
     EXCLARGS=(--env RIBO_TRAIN_EXCLUDE_TX="$EXCL"); NAME=${NAME}_exclmm50 ;;
  *) echo "bad ARM=$ARM (want A|B)" >&2; exit 1 ;;
esac
OUT=$NEW/results/loto/${NAME}_holdout_Hepatocytes
[ "$SEED" = 0 ] || OUT=${OUT}_seed${SEED}
mkdir -p "$NEW/logs/train" "$OUT"
# record what this arm actually changed, so the run is self-describing later
cat > "$OUT/arm.json" <<JSON
{"arm":"$ARM","model":"mamba$N_ATTN","rna_coverage":"mm10","ribo_target":"mm1",
 "pack_family":"packed_mm10cov","universe":"union_84472",
 "train_exclusion":$([ "$ARM" = B ] && echo "\"frac_discarded>0.5, 7640 tx, train+val only, test kept\"" || echo null),
 "baseline":"results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes",
 "seed":$SEED}
JSON
echo "ARM=$ARM name=$NAME seed=$SEED host=$(hostname) gpu=${CUDA_VISIBLE_DEVICES:-?} start=$(date -Iseconds)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
  --env RIBO_PACK_SUFFIX=mm10cov \
  --env RIBO_ORF_TRACK="$ORFTRACK" \
  --env RIBO_ONEHOT_FASTA="$NEW/data/union_universe.fa" \
  "${EXCLARGS[@]}" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/train_loto.py" \
    --out "$OUT" \
    --holdout Hepatocytes \
    --train_tissues Fibroblast,VSMC,ES,Fat,HA_EC,HCAEC,HUVEC \
    --emb_backend onehot --cov_norm global_mean \
    --min_train_signal 50 --val_fold 0 \
    --max_tx_per_tissue 6000 \
    --resume \
    --use_orf_track --n_attn_layers "$N_ATTN" --mixer mamba --mamba_d_state "$D_STATE" \
    --epochs 30 --patience 6 --warmup 2 \
    --channels 256 --n_blocks 10 --dropout 0.1 \
    --lr 3e-4 --weight_decay 1e-2 \
    --budget 16000 --count_weight 0.1 \
    --eval_cap 2500 --num_workers 8 \
    --device cuda --seed "$SEED"
echo "train_done=$(date -Iseconds) -> $OUT"
