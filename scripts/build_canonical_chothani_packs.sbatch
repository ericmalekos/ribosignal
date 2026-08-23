#!/usr/bin/env bash
#SBATCH --job-name=canon_pack
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/canon_pack_%A_%a.out
set -euo pipefail
# Build CANONICAL Chothani training packs: data/packed_canon (Fibroblast) + data/packed_canon_<T>.
#
# These are the packs a retrain reads. They differ from the deployed training packs
# (data/packed_union*) in exactly one thing: the P-sites come from alignments produced by
# scripts/riboseq_align.sbatch (EndToEnd + mm1 + ncRNA/cross-gene filter) instead of the 2026-07
# recipe. Measured divergence of that swap is 3.04-11.66% of signal depending on tissue
# (results/chothani_regeneration/divergence_by_tissue.md).
#
# WHAT IS REUSED, AND WHY THAT IS SOUND:
#   coverage.npy      RNA-seq is NOT re-aligned. It stays Local (only Ribo needed EndToEnd), and the
#                     universe is byte-identical, so the stored coverage is still exactly right.
#                     Re-deriving it would cost ~2.5 h/sample and attest to nothing.
#   tx_order/lengths  identical universe by construction; asserted below, not assumed.
#   ORF track         a function of sequence + annotation, not reads. data/packed_union/orf_track_v2.npy
#                     (meta says kozak=none, shape 212647172x5) stays valid and is passed to training
#                     via RIBO_ORF_TRACK -- no copy needed.
#
# WHAT IS REBUILT: target_counts.npy, and pack_meta.tsv with it. pack_meta carries per-transcript
# total_psites and the trainer filters on it (--min_train_signal 50), so a stale meta would silently
# select a different training set than the counts justify. write_pack regenerates both together.
#
# The pooled psites hd5 are a free by-product of the divergence runs
# (results/chothani_regeneration/_canonical_<T>/), so this stage is pooling only -- no BAM passes.
#
#   sbatch --array=0-7 scripts/build_canonical_chothani_packs.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3

TISSUES=(ES Fat Fibroblast HA_EC HCAEC Hepatocytes HUVEC VSMC)
NBAM=(6 6 32 6 5 5 3 11)          # expected sample count per tissue, from the manifest
T=${TISSUES[$SLURM_ARRAY_TASK_ID]}
WANT=${NBAM[$SLURM_ARRAY_TASK_ID]}

if [[ "$T" == "Fibroblast" ]]; then REF=$NEW/data/packed_union;    OUT=$NEW/data/packed_canon
else                                REF=$NEW/data/packed_union_$T; OUT=$NEW/data/packed_canon_$T; fi
W=$NEW/results/chothani_regeneration/_canonical_$T

for p in "$REF/tx_order.txt" "$REF/coverage.npy" "$REF/pack_meta.tsv" "$W"; do
  [[ -e "$p" ]] || { echo "missing $p" >&2; exit 1; }
done
shopt -s nullglob
HD5=("$W"/*_psites.hd5)
shopt -u nullglob
# Gate on the FULL sample set. A tissue pooled from a subset would look like a valid pack and train
# without complaint, while quietly being a different experiment from its deployed counterpart.
[[ ${#HD5[@]} -eq "$WANT" ]] || {
  echo "ABORT: $T has ${#HD5[@]} psites hd5, expected $WANT (divergence run not finished?)" >&2; exit 1; }
echo "tissue=$T ref=$(basename "$REF") out=$(basename "$OUT") hd5=${#HD5[@]} host=$(hostname) start=$(date -Iseconds)"

"$PY" - "$T" "$REF" "$OUT" "${HD5[@]}" <<'PYEOF'
import sys, pathlib, numpy as np
sys.path.insert(0, "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/scripts")
from prepare.packlib import pool_per_nt, write_pack

tissue, ref, out = sys.argv[1], pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
hd5 = sys.argv[4:]

order = [l.strip() for l in open(ref / "tx_order.txt") if l.strip()]
lengths = np.load(ref / "lengths.npy").astype(np.int64)
offs = np.concatenate([[0], np.cumsum(lengths)])
len_of = dict(zip(order, lengths))

# biotype from the reference pack_meta so the new pack keeps the 6-column form the eval code expects
biotype = {}
with open(ref / "pack_meta.tsv") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ci = {c: i for i, c in enumerate(hdr)}
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        biotype[f[ci["tx_id"]]] = (f[ci["chrom"]], f[ci["transcript_type"]])

cov_flat = np.load(ref / "coverage.npy")
cov = {t: cov_flat[offs[i]:offs[i + 1]].astype(np.int64) for i, t in enumerate(order)}
del cov_flat

tgt = pool_per_nt(hd5, "p_sites", set(order))
for t in order:                      # pool_per_nt omits transcripts with no reads
    if t not in tgt:
        tgt[t] = np.zeros(len_of[t], dtype=np.int64)

old = np.load(ref / "target_counts.npy").astype(np.int64)
new_tot = int(sum(int(v.sum()) for v in tgt.values()))
print(f"  P-sites {old.sum():,} -> {new_tot:,} ({100*(new_tot-old.sum())/old.sum():+.2f}%)")

write_pack(out, order, tgt, cov, group=f"{tissue}_canon", species="human",
           biotype=biotype, ref_len_of=len_of, n_ribo=len(hd5), n_rna=0)

# The ORF track is shared and row-indexed, so any drift in row order silently misaligns every
# feature. Assert rather than trust.
new_order = [l.strip() for l in open(out / "tx_order.txt")]
assert new_order == order, "tx_order changed"
assert np.array_equal(np.load(out / "lengths.npy").astype(np.int64), lengths), "lengths changed"
assert np.array_equal(np.load(out / "coverage.npy").astype(np.int64), np.load(ref / "coverage.npy").astype(np.int64)), \
    "coverage was altered; it must be reused verbatim"
print(f"  OK: {len(order):,} tx, order+lengths+coverage identical to {ref.name}")
PYEOF
echo "done=$(date -Iseconds) -> $OUT"
