#!/usr/bin/env bash
#SBATCH --job-name=l3x3_orf
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=8:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/l3x3_orf_%j.out
set -euo pipefail
# STAGE 2 of the mouse-liver 3x3: build the ORF track ONCE on the shared UNION universe.
#
# All nine cells share one universe -- tx_order.txt is byte-identical across the three stage-1 pools
# (md5 34e432db66e0e30659a80b5b6873ad32) -- so one track is correct for every cell and is referenced
# through RIBO_ORF_TRACK rather than copied into nine packs.
#
# --kozak none, ALWAYS. The Kozak ablation (task 20) found the heuristic redundant at best and harmful
# cross-tissue: four arms tied in-distribution (0.642-0.646) and V1/no-Kozak generalised best, while a
# learned kernel rediscovered Kozak unsupervised (r=0.807) making the explicit gate superfluous. This
# is a standing project rule (feedback_kozak_never_set), not a per-run choice. The released models
# were trained no-Kozak, so anything else would also break train/inference matching.
#
#   sbatch scripts/liver3x3_02_orftrack.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
L=$NEW/data/liver3x3
PACK=$L/pool_janich                     # any of the three: identical tx_order/offsets/lengths
FASTA=$NEW/data/heldout_refs/mouse_liver_union3_universe.fa
OUT=$L/orf_track_union3_nokozak.npy

for f in "$PACK/tx_order.txt" "$PACK/offsets.npy" "$FASTA"; do
  [[ -s "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
[[ -s "$OUT" ]] && { echo "SKIP: $OUT exists"; exit 0; }
echo "host=$(hostname) start=$(date -Iseconds)"
echo "  pack=$PACK  tx=$(wc -l < "$PACK/tx_order.txt")  fasta=$(grep -c '^>' "$FASTA")"

"$PY" "$NEW/scripts/build_orf_track.py" --pack "$PACK" --fasta "$FASTA" \
  --kozak none --mode ext --out "$OUT"

[[ -s "$OUT" ]] || { echo "ABORT: no track written" >&2; exit 1; }
# Name must match content -- a track called nokozak that was built with the heuristic is the exact
# silent mismatch the byte-size assertion elsewhere in this project exists to catch.
META=${OUT%.npy}_meta.json
if [[ -s "$META" ]]; then
  grep -q '"kozak": *"none"' "$META" || { echo "ABORT: $META is not kozak=none" >&2; exit 1; }
  echo "  meta confirms kozak=none"
fi
"$PY" -c "
import numpy as np
t=np.load('$OUT', mmap_mode='r'); o=np.load('$PACK/offsets.npy')
print(f'  track shape {t.shape} dtype {t.dtype}; pack sum_L {int(o[-1]):,}')
assert t.shape[0]==int(o[-1]), 'track rows != pack sum_L'
print('  OK: track rows match pack sum_L')
"
echo "done=$(date -Iseconds) -> $OUT"
