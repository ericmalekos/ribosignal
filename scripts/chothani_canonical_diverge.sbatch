#!/usr/bin/env bash
#SBATCH --job-name=cho_canon
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH --time=11:45:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/cho_canon_%A_%a.out
set -euo pipefail
# THE reproducibility test, per tissue, on the final recipe.
#
# Generalises scripts/huvec_canonical_diverge.sbatch from one tissue to all eight. HUVEC alone
# answered "can the training data be regenerated?" (3.04%, r=0.957). Eight tissues answer the
# question a reviewer actually asks: is that number a property of the pipeline, or of HUVEC?
#
# Bams come from scripts/riboseq_align.sbatch: EndToEnd + mm1 + filter applied to raw query-grouped
# STAR output, in that order. Pooling is the same packlib code path as every pack in the project.
#
# The comparison is per-nt P-sites, new vs the training pack's target_counts.npy, over the pack's
# OWN universe and row order. Two figures are reported because they answer different questions:
#   pct_of_signal  sum|new-old| / sum(old)   -- how much of the signal moved
#   pearson        per-nt correlation        -- whether it moved coherently or was reshaped
# A recipe can score well on one and badly on the other; the soft-clipped arm did exactly that
# (12.82% of signal but r=0.583, i.e. the signal was displaced, not merely rescaled).
#
#   sbatch --array=0-7 scripts/chothani_canonical_diverge.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
RBIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
ANNOT=$NEW/data/human_ribocode_annot_primary

TISSUES=(ES Fat Fibroblast HA_EC HCAEC Hepatocytes HUVEC VSMC)
T=${TISSUES[$SLURM_ARRAY_TASK_ID]}
# Fibroblast's union pack is `data/packed_union` (no tissue suffix -- it was the base tissue when
# the union universe was built); every other tissue got `packed_union_<T>`.
#
# NOT `data/packed`. That is the ORIGINAL single-tissue Fibroblast pack on the pre-union universe
# (36,668 tx vs 84,472) and its coverage_norm records `n_ribo_samples: 0` -- it is coverage-only, so
# diffing P-sites against it would compare this tissue's reads to a pack that has none. Verified via
# coverage_norm.json: packed = Fibroblast/0 ribo/36,668 tx, packed_union = Fibroblast/32 ribo/84,472
# tx. All eight tissues therefore share the 84,472-tx union universe, which is what makes the eight
# divergence numbers comparable to each other.
if [[ "$T" == "Fibroblast" ]]; then PACK=$NEW/data/packed_union; else PACK=$NEW/data/packed_union_$T; fi
BAMDIR=$NEW/data/ribo_bam_e2e_human/$T
W=$NEW/results/chothani_regeneration/_canonical_$T
mkdir -p "$W"

for p in "$PACK/tx_order.txt" "$PACK/target_counts.npy" "$BAMDIR"; do
  [[ -e "$p" ]] || { echo "missing $p" >&2; exit 1; }
done
shopt -s nullglob
BAMS=("$BAMDIR"/*.Aligned.toTranscriptome.out.bam)
shopt -u nullglob
[[ ${#BAMS[@]} -ge 1 ]] || { echo "no canonical bams in $BAMDIR" >&2; exit 1; }
echo "tissue=$T pack=$(basename "$PACK") bams=${#BAMS[@]} host=$(hostname) start=$(date -Iseconds)"

# Resumable: metaplots + RiboCode is the expensive half, and it writes hd5 before the compare runs.
shopt -s nullglob
_h=("$W"/*_psites.hd5)
shopt -u nullglob
if [[ ${#_h[@]} -eq 0 ]]; then
  printf '%s\n' "${BAMS[@]}" > "$W/bam_list.txt"
  ( cd "$W" && "$RBIN/metaplots" -a "$ANNOT" -i "$W/bam_list.txt" -o "canon_$T" )
  CFG=$(ls "$W/canon_${T}_pre_config.txt" | head -1)
  ( cd "$W" && "$RBIN/RiboCode" -a "$ANNOT" -c "$CFG" -l no -g -o "canon_$T" )
fi

"$PY" - <<PYEOF
import numpy as np, glob, json, pathlib, sys
sys.path.insert(0, "$NEW/scripts")
from prepare.packlib import pool_per_nt
pack = pathlib.Path("$PACK")
order = [l.strip() for l in open(pack / "tx_order.txt")]
lens = dict(zip(order, np.load(pack / "lengths.npy").astype(np.int64)))
hd5 = sorted(glob.glob("$W/*_psites.hd5"))
if not hd5:
    sys.exit("no psites hd5 written -- RiboCode stage produced nothing")
pooled = pool_per_nt(hd5, "p_sites", set(order))
new = np.concatenate([np.asarray(pooled[t], dtype=np.int64) if t in pooled
                      else np.zeros(lens[t], dtype=np.int64) for t in order])
old = np.load(pack / "target_counts.npy").astype(np.int64)
if new.shape != old.shape:
    sys.exit(f"axis mismatch: pooled {new.shape} vs pack {old.shape}")
d = np.abs(new - old)
res = dict(tissue="$T", pack=pack.name,
           recipe="CANONICAL: EndToEnd + mm1 + filter on query-grouped bams",
           n_bams=len(hd5), n_tx=len(order),
           old_psites=int(old.sum()), new_psites=int(new.sum()),
           pct_of_signal=float(100*d.sum()/max(old.sum(),1)),
           pct_total_change=float(100*(new.sum()-old.sum())/max(old.sum(),1)),
           pearson=float(np.corrcoef(new, old)[0,1]))
pathlib.Path("$W/divergence_canonical.json").write_text(json.dumps(res, indent=2))
print(f"  {res['tissue']}: old {res['old_psites']:,} | new {res['new_psites']:,} "
      f"({res['pct_total_change']:+.2f}%) | {res['pct_of_signal']:.2f}% of signal moved | "
      f"per-nt r={res['pearson']:.4f}")
PYEOF
echo "done=$(date -Iseconds) -> $W/divergence_canonical.json"
