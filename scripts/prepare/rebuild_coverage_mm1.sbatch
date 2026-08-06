#!/usr/bin/env bash
#SBATCH --job-name=rebuild_mm1
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --array=1-8
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/prepare/rebuild_mm1_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/prepare/rebuild_mm1_%A_%a.err
set -euo pipefail
# Rebuild each Chothani tissue pack's COVERAGE with mm1 (unique-mapper) RNA-seq, via
# prepare_pack.py --reuse-target (Ribo target unchanged -> reused byte-for-byte; only coverage
# is re-pooled). Builds to a TEMP dir, asserts the target stayed byte-identical, backs up the
# mm20 coverage, then swaps ONLY coverage.npy/coverage_norm.json/pack_meta.tsv into the real
# pack -- the real pack's target/offsets/lengths/tx_order are NEVER overwritten. Brain dropped.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
PER=$NEW/data/rnaseq_coverage/per_sample
MAP=$NEW/data/loto_rnaseq_srr_tissue.tsv
mkdir -p "$NEW/logs/prepare"

TISSUES=(Fibroblast VSMC ES Fat HA_EC HCAEC HUVEC Hepatocytes)
T=${TISSUES[$((SLURM_ARRAY_TASK_ID - 1))]}
PACK=$([ "$T" = Fibroblast ] && echo "$NEW/data/packed" || echo "$NEW/data/packed_$T")
TMP=/data/tmp/emalekos/mm1_rebuild/$T
rm -rf "$TMP"; mkdir -p "$TMP"

# tissue's mm1 RNA coverage SRRs: Fibroblast = per-sample files NOT in the tissue map; else map
if [ "$T" = Fibroblast ]; then
  SRRS=$(comm -23 \
    <(ls "$PER"/SRR*_coverage.hd5 | grep -oE 'SRR[0-9]+' | sort -u) \
    <(cut -f1 "$MAP" | sort -u))
else
  SRRS=$(awk -F'\t' -v t="$T" '$2==t{print $1}' "$MAP")
fi
COV_ARGS=""
for s in $SRRS; do
  f="$PER/${s}_coverage.hd5"
  [ -s "$f" ] || { echo "MISSING mm1 coverage $f" >&2; exit 1; }
  COV_ARGS="$COV_ARGS $f"
done
N=$(echo $SRRS | wc -w)
echo "[$T] pack=$PACK  mm1_coverage_samples=$N  host=$(hostname)  start=$(date -Iseconds)"

# 1. build the mm1 pack to TMP, reusing the existing pack's universe + Ribo target
# shellcheck disable=SC2086
"$PY" "$NEW/scripts/prepare/prepare_pack.py" \
  --reuse-target "$PACK" \
  --rna-coverage $COV_ARGS \
  --tx2biotype "$NEW/data/tx2biotype.tsv" \
  --group "$T" --species human --out "$TMP"

# 2. safety: target must be byte-identical (reuse path re-derives it)
"$PY" - "$TMP" "$PACK" <<'EOF'
import sys, numpy as np
tmp, pack = sys.argv[1], sys.argv[2]
a = np.load(f"{tmp}/target_counts.npy"); b = np.load(f"{pack}/target_counts.npy")
assert np.array_equal(a, b), "TARGET CHANGED -- aborting swap"
ca = np.load(f"{tmp}/coverage.npy");     cb = np.load(f"{pack}/coverage.npy")
print(f"  target identical=True  coverage changed={not np.array_equal(ca,cb)}  "
      f"mm1_total={int(ca.sum())}  mm20_total={int(cb.sum())}")
EOF

# 3. back up mm20 (once) then swap ONLY the coverage-derived files into the real pack
[ -f "$PACK/coverage_mm20.npy" ]       || cp "$PACK/coverage.npy"        "$PACK/coverage_mm20.npy"
[ -f "$PACK/coverage_norm_mm20.json" ] || cp "$PACK/coverage_norm.json"  "$PACK/coverage_norm_mm20.json"
[ -f "$PACK/pack_meta_mm20.tsv" ]      || cp "$PACK/pack_meta.tsv"       "$PACK/pack_meta_mm20.tsv"
cp "$TMP/coverage.npy"       "$PACK/coverage.npy"
cp "$TMP/coverage_norm.json" "$PACK/coverage_norm.json"
cp "$TMP/pack_meta.tsv"      "$PACK/pack_meta.tsv"
rm -rf "$TMP"
echo "[$T] mm1 coverage swapped into $PACK (mm20 backed up as *_mm20.*)  done=$(date -Iseconds)"
