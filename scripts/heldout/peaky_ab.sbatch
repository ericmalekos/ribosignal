#!/usr/bin/env bash
#SBATCH --job-name=peaky_ab
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=01:30:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/peaky_ab_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/peaky_ab_%j.err
set -euo pipefail
# O3 track-2 A/B: run the standard O3 Phase-1 flow on the PK=0.2 anti-smoothing model, apples-to-apples vs the
# standard nokozak mm1 model (feature-AUC 0.669, and its over-calling). CPU on the gpu partition idle cores.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
PK=$NEW/results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_peaky0p2_holdout_Hepatocytes
STD=$NEW/results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes
PROF=$PK/dropin/pred_profiles.npz
ANNOT=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/human_ribocode_annot_primary
RIBOPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
CASPY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
export PYTHONNOUSERSITE=1
mkdir -p "$NEW/logs/dropin"

# 1. drop-in: real (truth) + pred_preddepth (standalone over-calling at theta=1)
for V in real pred_preddepth; do
  echo "=== dropin $V ==="
  "$RIBOPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$PROF" --annot "$ANNOT" \
    --variant "$V" --out "$PK/dropin" --min_aa 5 --pval 0.05
done

# 2. Poisson theta=2.5 candidate pool (same candidate set as the standard O3 Phase-1 features)
OUT25=$PK/dropin_sweep_poisson/theta_2.5; rm -rf "$OUT25"; mkdir -p "$OUT25"
echo "=== poisson theta=2.5 candidates ==="
"$RIBOPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$PROF" --annot "$ANNOT" \
  --variant pred_preddepth --pred_scale 2.5 --pred_poisson --seed 0 --out "$OUT25" --min_aa 5 --pval 0.05

# 3. shape features on the PEAKY predicted density (theta=2.5 candidates, real truth)
echo "=== orf_features (peaky) ==="
"$CASPY" "$NEW/scripts/heldout/orf_features.py" --profiles "$PROF" \
  --candidates "$OUT25/pred_preddepth_collapsed.txt" --truth "$PK/dropin/real_collapsed.txt" \
  --out "$PK/o3_features_peaky.tsv"

# 4. fp_filter train -> feature-separation AUC (vs the standard model's 0.669)
echo "=== fp_filter train (peaky) ==="
"$CASPY" "$NEW/scripts/heldout/fp_filter.py" train --features "$PK/o3_features_peaky.tsv" \
  --out "$PK/fp_filter_peaky.joblib"

# 5. head-to-head summary: over-calling (call composition) peaky vs standard, both pred_preddepth theta=1
echo ""
echo "==================== A/B SUMMARY: peaky (PK=0.2) vs standard ===================="
"$CASPY" - "$PK" "$STD" <<'EOF'
import sys
from collections import Counter
PK, STD = sys.argv[1], sys.argv[2]
def comp(path):
    n=Counter()
    with open(path) as fh:
        h=fh.readline().rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
        for ln in fh:
            n[ln.split("\t")[ci["ORF_type"]]]+=1
    return n
for lab,base in [("STANDARD",STD),("PEAKY(PK=0.2)",PK)]:
    c=comp(f"{base}/dropin/pred_preddepth_collapsed.txt")
    tot=sum(c.values())
    print(f"  {lab:15s} total={tot:6d}  annotated={c['annotated']:6d}  uORF={c['uORF']:5d}  "
          f"novel={c['novel']:5d}  dORF={c['dORF']:5d}")
print("  (fewer non-canonical at similar annotated = less over-calling = anti-smoothing helped)")
EOF
echo "peaky A/B done $(date -Iseconds)"
