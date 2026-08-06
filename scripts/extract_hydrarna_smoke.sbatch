#!/usr/bin/env bash
#SBATCH --job-name=hydra_smoke
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH --time=00:30:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/hydrarna/smoke_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/hydrarna/smoke_%j.err
set -euo pipefail
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-hydrarna-latest.img
OUT=$NEW/data/hydrarna_token_emb/smoke
mkdir -p "$OUT" "$NEW/logs/hydrarna"
# Triton (Mamba's GPU-kernel compiler) writes a cache to $HOME/.triton; --no-home leaves $HOME on the
# read-only group fs, so redirect HOME + TRITON_CACHE_DIR to node-local writable scratch.
SCRATCH=/data/tmp/emalekos/hydra_${SLURM_JOB_ID}
mkdir -p "$SCRATCH/.triton"
echo "host=$(hostname) scratch=$SCRATCH start=$(date -Iseconds)"
# extract first 3 transcripts
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env HOME="$SCRATCH" --env TRITON_CACHE_DIR="$SCRATCH/.triton" \
  -B "$BASE" -B /data/tmp/emalekos "$SIF" \
  python3 "$NEW/scripts/extract_hydrarna_pertoken.py" \
    -i "$NEW/data/fibroblast_universe.fa" -o "$OUT" --start 0 --end 3 --log-every 1
echo "=== verify: emb shape (L,1024) matches each transcript's length ==="
/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 - "$NEW" "$OUT" <<'EOF'
import sys, numpy as np
NEW, OUT = sys.argv[1], sys.argv[2]
seqs = []
h = None; p = []
for line in open(f"{NEW}/data/fibroblast_universe.fa"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if h: seqs.append((h, "".join(p)))
        h = line[1:].split()[0]; p = []
    elif line: p.append(line.strip())
    if len(seqs) >= 3 and h and "".join(p) == "": break
seqs = seqs[:3]
import os
ok = True
for txid, seq in seqs:
    f = f"{OUT}/{txid}_tokens.npy"
    if not os.path.exists(f): print(f"  MISSING {txid}"); ok = False; continue
    a = np.load(f)
    match = a.shape == (len(seq), 1024)
    print(f"  {txid}: emb {a.shape} vs seq {len(seq)} dtype={a.dtype}  {'OK' if match else 'MISMATCH'}")
    ok = ok and match
print("SMOKE PASS" if ok else "SMOKE FAIL")
EOF
echo "done=$(date -Iseconds)"
