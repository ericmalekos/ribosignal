#!/usr/bin/env bash
#SBATCH --job-name=orthrus_smoke
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orthrus/smoke_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orthrus/smoke_%j.err
set -euo pipefail
REPO=/private/groups/carpenterlab/emalekos/RNAZoo_meta
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img
WEIGHTS_HOST=/private/groups/carpenterlab/emalekos/orthrus6t_weights
NEW=$REPO/RNAZoo/experiments/riboseq_signal_model
mkdir -p "$NEW/logs/orthrus"
echo "=== node: $(hostname) ==="; nvidia-smi -L

# 1) Equivalence check: model.forward (per-token, L x 512) mean-pools to model.representation.
#    Confirms forward() is the pre-pool tensor and one-hot has no CLS/EOS to strip.
singularity exec --nv --no-home --cleanenv \
  --env PYTHONNOUSERSITE=1 --env HF_HUB_OFFLINE=1 --env TRANSFORMERS_OFFLINE=1 \
  --bind ${REPO}:${REPO} --bind ${WEIGHTS_HOST}:/weights \
  "$SIF" /opt/conda/envs/orthrus/bin/python - <<'PY'
import numpy as np, torch
from transformers import AutoModel
dev = torch.device("cuda")
model = AutoModel.from_pretrained("/weights/orthrus-large-4-track",
                                  trust_remote_code=True).to(dev).eval()
def oh(seq):
    L=len(seq); o=np.zeros((L,4),np.float32); idx={"A":0,"C":1,"G":2,"T":3,"U":3}
    for i,b in enumerate(seq.upper()):
        j=idx.get(b)
        if j is not None: o[i,j]=1
    return o
for seq in ["ACGUACGUACGUACGGCUAGCUAGCUAGCUAGCUAGCUAGC",
            "AUGGCCUAG"*30,
            "GC"*500]:
    L=len(seq)
    x=torch.from_numpy(oh(seq)).to(dev).unsqueeze(0)          # (1,L,4)
    lengths=torch.tensor([L],device=dev,dtype=torch.float32)
    with torch.no_grad():
        h=model.forward(x, channel_last=True)                 # (1,L,512) per-token
        r=model.representation(x, lengths, channel_last=True) # (1,512) pooled
    per=h.squeeze(0).float().cpu().numpy()
    pooled=r.squeeze(0).float().cpu().numpy()
    mean_of_per=per.mean(axis=0)
    dmax=float(np.abs(mean_of_per-pooled).max())
    assert per.shape==(L,512), f"bad shape {per.shape} for L={L}"
    print(f"L={L:5d}  per-token{per.shape}  |mean(forward)-representation|max={dmax:.2e}")
    assert dmax < 1e-3, "forward-mean does not match representation -> pooling mismatch"
print("EQUIVALENCE OK: forward() is per-token (L,512), no CLS/EOS, mean==representation")
PY

# 2) Real code-path check: run the actual extractor on the first 24 records of chunk_1.
SMOKE_OUT="$NEW/data/orthrus4t_token_emb_smoke/chunk_1"
rm -rf "$NEW/data/orthrus4t_token_emb_smoke"
singularity exec --nv --no-home --cleanenv \
  --env PYTHONNOUSERSITE=1 --env HF_HUB_OFFLINE=1 --env TRANSFORMERS_OFFLINE=1 \
  --bind ${REPO}:${REPO} --bind ${WEIGHTS_HOST}:/weights \
  "$SIF" /opt/conda/envs/orthrus/bin/python \
    $NEW/scripts/extract_orthrus4t_pertoken.py \
    --input  $NEW/data/rinalmo_chunks/chunk_1.fa \
    --output $SMOKE_OUT \
    --max-length 10000 --limit 24 --dtype float16

# 3) Verify a produced file aligns to its RiNALMo counterpart (same tx, same L).
singularity exec --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --bind ${REPO}:${REPO} \
  "$SIF" /opt/conda/envs/orthrus/bin/python - "$SMOKE_OUT" "$NEW" <<'PY'
import sys, glob, os, numpy as np
smoke, NEW = sys.argv[1], sys.argv[2]
files = sorted(glob.glob(os.path.join(smoke, "*_tokens.npy")))
print(f"produced {len(files)} per-token files")
idx = {}
with open(os.path.join(NEW, "data/rinalmo_token_emb/tx_index.tsv")) as fh:
    fh.readline()
    for ln in fh:
        t, p = ln.rstrip("\n").split("\t"); idx[t] = p
checked = 0
for f in files[:8]:
    tx = os.path.basename(f)[:-len("_tokens.npy")]
    o = np.load(f)
    assert o.ndim == 2 and o.shape[1] == 512, f"{tx}: bad orthrus shape {o.shape}"
    assert o.dtype == np.float16, f"{tx}: dtype {o.dtype}"
    if tx in idx:
        r = np.load(idx[tx])
        assert r.shape[0] == o.shape[0], f"{tx}: L mismatch orthrus {o.shape[0]} vs rinalmo {r.shape[0]}"
        print(f"  {tx}: orthrus{o.shape} rinalmo{r.shape} L-aligned OK")
        checked += 1
assert checked > 0, "no overlap with rinalmo index to cross-check"
print(f"ALIGNMENT OK: {checked} tx cross-checked, orthrus L == rinalmo L (1:1 nt axis)")
PY
echo "=== smoke done ==="
