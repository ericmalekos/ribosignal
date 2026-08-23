#!/usr/bin/env bash
#SBATCH --job-name=huvec_canon
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH --time=11:45:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/huvec_canon_%j.out
set -euo pipefail
# THE reproducibility test, on the final recipe. Supersedes arms A/B/C, all three of which were
# flawed: A and B omitted the filter entirely, and C applied it to COORDINATE-SORTED bams where
# cross-gene detection silently fails (it reported 0.30% instead of the true ~9.5%).
#
# These bams come from scripts/riboseq_align.sbatch: EndToEnd + mm1 + filter applied to raw
# query-grouped STAR output, in that order. Same pooling code as every other pack.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
RBIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
ANNOT=$NEW/data/human_ribocode_annot_primary
PACK=$NEW/data/packed_union_HUVEC
BAMDIR=$NEW/data/ribo_bam_e2e_human/HUVEC
W=$NEW/results/chothani_regeneration/_canonical_HUVEC
mkdir -p "$W"
echo "host=$(hostname) start=$(date -Iseconds)"
if ! ls "$W"/*_psites.hd5 >/dev/null 2>&1; then
  ls "$BAMDIR"/*.Aligned.toTranscriptome.out.bam > "$W/bam_list.txt"
  echo "  bams: $(wc -l < "$W/bam_list.txt")"
  ( cd "$W" && "$RBIN/metaplots" -a "$ANNOT" -i "$W/bam_list.txt" -o canon_HUVEC )
  CFG=$(ls "$W"/canon_HUVEC_pre_config.txt | head -1)
  ( cd "$W" && "$RBIN/RiboCode" -a "$ANNOT" -c "$CFG" -l no -g -o canon_HUVEC )
fi
"$PY" - <<PYEOF
import numpy as np, glob, json, pathlib, sys
sys.path.insert(0, "$NEW/scripts")
from prepare.packlib import pool_per_nt
pack = pathlib.Path("$PACK")
order = [l.strip() for l in open(pack / "tx_order.txt")]
lens = dict(zip(order, np.load(pack / "lengths.npy").astype(np.int64)))
hd5 = sorted(glob.glob("$W/*_psites.hd5"))
pooled = pool_per_nt(hd5, "p_sites", set(order))
new = np.concatenate([np.asarray(pooled[t], dtype=np.int64) if t in pooled
                      else np.zeros(lens[t], dtype=np.int64) for t in order])
old = np.load(pack / "target_counts.npy").astype(np.int64)
d = np.abs(new - old)
res = dict(tissue="HUVEC", recipe="CANONICAL: EndToEnd + mm1 + filter on query-grouped bams",
           n_bams=len(hd5), old_psites=int(old.sum()), new_psites=int(new.sum()),
           pct_of_signal=float(100*d.sum()/max(old.sum(),1)),
           pearson=float(np.corrcoef(new, old)[0,1]))
pathlib.Path("$W/divergence_canonical.json").write_text(json.dumps(res, indent=2))
print(f"  CANONICAL: old {res['old_psites']:,} | new {res['new_psites']:,} | "
      f"{res['pct_of_signal']:.2f}% of signal | per-nt r={res['pearson']:.4f}")
print("  prior (all flawed): A my recon 11.69% r=0.935 | B minimal 12.82% r=0.583 | C broken-filter n/a")
print("  reference: Janich gate via current pipeline = 1.69%")
PYEOF
echo "done=$(date -Iseconds)"
