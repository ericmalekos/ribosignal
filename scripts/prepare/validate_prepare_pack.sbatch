#!/usr/bin/env bash
#SBATCH --job-name=validate_prepare_pack
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/prepare/validate_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/prepare/validate_%j.err
set -euo pipefail
# Validate prepare_pack.py end-to-end via the FULL future-user path: pool HUVEC's Ribo target
# straight from the raw RiboCode per-sample psites hd5 + the mm1 RNA coverage, reusing the
# shared v49 universe, and assert the resulting target_counts is BYTE-IDENTICAL to the existing
# data/packed_HUVEC pack (Ribo unchanged) while coverage differs (mm1 vs mm20).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
OUT=/data/tmp/emalekos/prepare_pack_validate/huvec_mm1
mkdir -p "$NEW/logs/prepare" "$(dirname "$OUT")"

HUVEC_SRR="SRR15513270 SRR15513271 SRR15513272"   # mm1 RNA per-sample coverage
COV_ARGS=""
for s in $HUVEC_SRR; do COV_ARGS="$COV_ARGS $NEW/data/rnaseq_coverage/per_sample/${s}_coverage.hd5"; done

# shellcheck disable=SC2086
"$PY" "$NEW/scripts/prepare/prepare_pack.py" \
  --ribo-psites "$ECH/data/ribocode_per_tissue/HUVEC" \
  --rna-coverage $COV_ARGS \
  --ref-pack "$NEW/data/packed" \
  --tx2biotype "$NEW/data/tx2biotype.tsv" \
  --group HUVEC --species human --out "$OUT"

echo "=== VALIDATION ==="
"$PY" - "$OUT" "$NEW/data/packed_HUVEC" <<'EOF'
import sys, numpy as np
new, old = sys.argv[1], sys.argv[2]
ta = np.load(f"{new}/target_counts.npy"); tb = np.load(f"{old}/target_counts.npy")
ca = np.load(f"{new}/coverage.npy");      cb = np.load(f"{old}/coverage.npy")
la = np.load(f"{new}/lengths.npy");       lb = np.load(f"{old}/lengths.npy")
print("lengths identical   :", np.array_equal(la, lb))
print("target byte-identical:", np.array_equal(ta, tb), "(Ribo unchanged; MUST be True)")
print("coverage differs (mm1):", not np.array_equal(ca, cb),
      f"  mm1_total={int(ca.sum())}  mm20_total={int(cb.sum())}")
ok = np.array_equal(la, lb) and np.array_equal(ta, tb) and not np.array_equal(ca, cb)
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
EOF
echo "done=$(date -Iseconds)"
