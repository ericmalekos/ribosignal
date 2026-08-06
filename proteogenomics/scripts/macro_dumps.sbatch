#!/usr/bin/env bash
#SBATCH --job-name=macro_dump
#SBATCH --partition=gpu
#SBATCH --gres=none
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=8:00:00
#SBATCH --array=0-11%6
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/macro_dump_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/macro_dump_%A_%a.err
set -euo pipefail
# Model prediction on each macrophage population (CPU, idle GPU-node cores). Produces the per-nt
# predicted P-site profiles that enumerate_score_orfs uses to SELECT translated ORFs for the MODEL DB.
# --fasta is REQUIRED (dump_pred_profiles hardcodes the fibroblast FASTA otherwise -> silent tx skips).
#
# CHECKPOINT SELECTION: the architecture is read from <run>/args.json (mixer, n_attn_layers, channels,
# emb_backend, use_orf_track), so switching models means changing MODEL_RUN and nothing else. Both the
# attn and mamba4 union checkpoints are emb_backend=onehot + use_orf_track=true, so the existing
# RIBO_ONEHOT_FASTA / orf_track_v2_nokozak.npy wiring carries over unchanged. An FM-backed checkpoint
# (rinalmo/orthrus/hydrarna) would NOT work here -- it needs per-tx embeddings over the mouse universe
# that do not exist.
#
# TAG namespaces every downstream artefact so two checkpoints can coexist. Defaults reproduce the
# original single-model layout exactly ($PD/pred/$POP, etc).
#   MODEL_RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes \
#   TAG=_mamba4_union sbatch macro_dumps.sbatch
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
PD=$NEW/proteogenomics/data/macrophage_tissue
# CONTAINER + DEVICE are mixer-dependent:
#   transformer (attn) -> rinalmo image, runs on CPU (--gres=none, idle GPU-node cores)
#   mamba       (mamba4) -> ORTHRUS image (the only one carrying mamba_ssm 1.2.0.post1, torch 2.2.0+cu118)
#                           and it is GPU-ONLY: mamba_ssm ships CUDA kernels with no CPU fallback, so a
#                           --device cpu run dies at import with ModuleNotFoundError/no-CUDA. Submit with
#                           --gres=gpu:1 and DEVICE=cuda.
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
SIF=${SIF:-$SIFDIR/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img}
DEVICE=${DEVICE:-cpu}
# default = the pre-union deployed checkpoint: 7-tissue LOTO (Fibroblast,VSMC,ES,Fat,HA_EC,HCAEC,HUVEC),
# hold out Hepatocytes, applied cross-species to mouse. NOT a "fibroblast model".
RUN=${MODEL_RUN:-results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes}
TAG=${TAG:-}
mkdir -p "$NEW/proteogenomics/logs"
[ -s "$NEW/$RUN/best.pt" ]   || { echo "ABORT: no checkpoint at $NEW/$RUN/best.pt"; exit 1; }
[ -s "$NEW/$RUN/args.json" ] || { echo "ABORT: no args.json at $NEW/$RUN/args.json"; exit 1; }
echo "MODEL_RUN=$RUN"; echo "TAG=${TAG:-<none>}"
grep -oE '"(mixer|n_attn_layers|channels|emb_backend|use_orf_track)": *[^,}]*' "$NEW/$RUN/args.json" | sed 's/^/  /'

mapfile -t POPS < <(ls "$PD"/universes/*_universe_tx.txt | xargs -n1 basename | sed 's/_universe_tx.txt//' | sort)
POP="${POPS[$SLURM_ARRAY_TASK_ID]:?no population at index $SLURM_ARRAY_TASK_ID}"
DS="macro_$POP"
PACK="$NEW/data/packed_heldout_$DS"
FASTA="$PD/universes/${POP}_universe.fa"
TXLIST="$PD/universes/${POP}_universe_tx.txt"   # predict the FULL expressed universe
OUT="$PD/pred${TAG}/$POP"
mkdir -p "$OUT"
echo "dump population=$POP host=$(hostname) start=$(date -Iseconds)"
[ -s "$OUT/pred_profiles.npz" ] && { echo "SKIP $POP (pred_profiles.npz exists)"; exit 0; }

NVFLAG=""; [ "$DEVICE" = cuda ] && NVFLAG="--nv"
echo "SIF=$(basename "$SIF")  DEVICE=$DEVICE"
# shellcheck disable=SC2086
singularity exec $NVFLAG --no-home --cleanenv --env PYTHONNOUSERSITE=1 --env OMP_NUM_THREADS=16 \
  --env RIBO_ORF_TRACK="$PACK/orf_track_v2_nokozak.npy" \
  --env RIBO_ONEHOT_FASTA="$FASTA" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$NEW/$RUN" --heldout "$DS" \
    --out "$OUT" --fasta "$FASTA" --tx_list "$TXLIST" --device "$DEVICE"

# guard: heavy skip = wrong FASTA/pack
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 -c "
import numpy as np,sys
d=np.load('$OUT/pred_profiles.npz',allow_pickle=True)
nt=len(d['tx_ids']); ns=int(str(d['meta'][3])); frac=ns/(nt+ns) if (nt+ns) else 0
print(f'$POP dump: {nt} tx, {ns} skipped ({frac:.1%})'); sys.exit(1 if frac>0.05 else 0)" \
  || { echo "ABORT $POP: >5% tx skipped"; exit 1; }
echo "DONE dump $POP $(date -Iseconds)"
