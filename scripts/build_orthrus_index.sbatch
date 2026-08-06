#!/usr/bin/env bash
#SBATCH --job-name=orthrus_idx
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orthrus/index_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/orthrus/index_%j.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3

$PY $NEW/scripts/build_embedding_index.py --emb-dir $NEW/data/orthrus4t_token_emb

# Compare the Orthrus and RiNALMo tx sets: they should match 1:1 (same universe, same
# 10000-nt cap). Report the intersection so the 3-way comparison split is on record.
$PY - <<'PY'
NEW="/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model"
def load(p):
    s=set()
    with open(p) as fh:
        fh.readline()
        for ln in fh: s.add(ln.split("\t",1)[0])
    return s
r=load(f"{NEW}/data/rinalmo_token_emb/tx_index.tsv")
o=load(f"{NEW}/data/orthrus4t_token_emb/tx_index.tsv")
print(f"rinalmo tx: {len(r):,}")
print(f"orthrus tx: {len(o):,}")
print(f"intersection: {len(r & o):,}")
print(f"rinalmo-only: {len(r - o):,}   orthrus-only: {len(o - r):,}")
PY
echo "index build done=$(date -Iseconds)"
