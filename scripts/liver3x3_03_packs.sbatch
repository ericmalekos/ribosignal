#!/usr/bin/env bash
#SBATCH --job-name=l3x3_pack
#SBATCH --partition=medium
#SBATCH --cpus-per-task=6
#SBATCH --mem=180G
#SBATCH --time=8:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/l3x3_pack_%A_%a.out
set -euo pipefail
# STAGE 3 of the mouse-liver 3x3: assemble all NINE cells from the stage-1 cached hd5.
#
#            RNA input ->     janich      gse243134    wang
#   Ribo   janich                0            1          2
#   source gse243134             3            4          5
#          wang                  6            7          8
#
# EVERY cell is assembled here, including the three diagonals that stage 1 already wrote as a
# by-product. That is deliberate: if six cells came from prepare_pack and three from
# prepare_from_bams, the diagonal would differ from the off-diagonal by code path rather than by
# data, and the diagonal is exactly where the factorial's "matched Ribo and RNA" baseline lives.
# Nine identical assemblies, one comparison.
#
# The three expensive BAM passes happened once in stage 1; this stage is pure hd5 pooling.
#
#   sbatch --array=0-8%5 scripts/liver3x3_03_packs.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
L=$NEW/data/liver3x3
UNI=$NEW/data/heldout_refs/mouse_liver_union3_tx.txt
TRACK=$L/orf_track_union3_nokozak.npy

DS=(janich gse243134 wang)
RIBO_I=$((SLURM_ARRAY_TASK_ID / 3))
RNA_I=$((SLURM_ARRAY_TASK_ID % 3))
RB=${DS[$RIBO_I]}
RN=${DS[$RNA_I]}
OUT=$L/pack_ribo-${RB}_rna-${RN}

mapfile -t PS < <(ls $L/pool_$RB/_work/psites/*_psites.hd5 | sort)
mapfile -t CV < <(ls $L/pool_$RN/_work/coverage/*_coverage.hd5 | sort)
[[ ${#PS[@]} -ge 1 && ${#CV[@]} -ge 1 ]] || { echo "no cached hd5 for ribo=$RB rna=$RN" >&2; exit 1; }
[[ -s "$TRACK" ]] || { echo "missing ORF track $TRACK -- run stage 2 first" >&2; exit 1; }
echo "cell=$SLURM_ARRAY_TASK_ID ribo=$RB (${#PS[@]} hd5) rna=$RN (${#CV[@]} hd5) host=$(hostname) start=$(date -Iseconds)"

"$PY" "$NEW/scripts/prepare/prepare_pack.py" \
  --ribo-psites "${PS[@]}" \
  --rna-coverage "${CV[@]}" \
  --universe-tx "$UNI" \
  --group "l3x3_ribo-${RB}_rna-${RN}" --species mouse --out "$OUT" \
  --tx2biotype "$NEW/data/tx2biotype_mouse.tsv"

# Gate the contract, not the exit status. Every cell must share the universe exactly, or the nine
# numbers are not comparable and the whole factorial is void.
"$PY" -c "
import numpy as np, sys, json, hashlib
o='$OUT'
want=[l.strip() for l in open('$UNI') if l.strip()]
got=[l.strip() for l in open(o+'/tx_order.txt')]
if got!=want:
    print(f'ABORT: tx_order != universe ({len(got)} vs {len(want)})'); sys.exit(1)
t=np.load(o+'/target_counts.npy'); c=np.load(o+'/coverage.npy'); off=np.load(o+'/offsets.npy')
tr=np.load('$TRACK', mmap_mode='r')
assert len(t)==len(c)==int(off[-1])==tr.shape[0], 'axis mismatch across target/coverage/offsets/track'
print(f'  {len(got):,} tx  {t.sum():,} P-sites  {c.sum():,} coverage  sum_L {len(t):,}')
print('  OK: universe + axes agree with the shared ORF track')
"
echo "done=$(date -Iseconds) -> $OUT"
