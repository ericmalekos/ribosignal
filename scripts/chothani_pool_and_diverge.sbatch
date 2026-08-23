#!/usr/bin/env bash
#SBATCH --job-name=cho_pool
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH --time=11:45:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/cho_pool_%A_%a.out
set -euo pipefail
# Phase 1.2/1.3: re-derive each Chothani tissue's P-site target from the re-fetched, re-aligned BAMs,
# and MEASURE how far it lands from the pack the released models were actually trained on.
#
# TARGETS ONLY, NOT A FULL PACK. The divergence question is entirely about `target_counts.npy`. The
# RNA coverage hd5 for these tissues were deleted along with the P-sites and are not being
# regenerated here, so forcing a full pack build would require inventing an RNA arm that has nothing
# to do with the measurement. This runs the Ribo stage and pools onto the existing pack's transcript
# axis, which is exactly what a comparison needs and nothing more.
#
# EXPECT DIVERGENCE. RiboCode metaplots re-derives which read lengths carry periodicity PER SAMPLE,
# and that alone cost 1.69% on the Janich rebuild with byte-identical inputs. The original recipe is
# unrecoverable -- data/packed has no provenance.json at all. So the deliverable is a number and a
# reproducible path, not a match.
#
#   sbatch --array=0-7%3 scripts/chothani_pool_and_diverge.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
RBIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
ANNOT=$NEW/data/human_ribocode_annot_primary
BAM=$NEW/data/chothani_bam
MAN=$NEW/data/external/chothani_ribo_refetch/manifest.tsv
OUT=$NEW/results/chothani_regeneration
mkdir -p "$OUT"

TISSUES=(ES Fat Fibroblast HA_EC HCAEC Hepatocytes HUVEC VSMC)
T=${TISSUES[$SLURM_ARRAY_TASK_ID]}
# Fibroblast is the original data/packed_union; the rest are packed_union_<tissue>.
if [[ "$T" == "Fibroblast" ]]; then PACK=$NEW/data/packed_union; else PACK=$NEW/data/packed_union_$T; fi
[[ -s "$PACK/target_counts.npy" ]] || { echo "no reference pack at $PACK" >&2; exit 1; }

mapfile -t RUNS < <(awk -F'\t' -v t="$T" 'NR>1{s=$2; sub(/_[0-9]+ .*/,"",s); if (s==t) print $1}' "$MAN")
[[ ${#RUNS[@]} -ge 1 ]] || { echo "no runs for tissue $T" >&2; exit 1; }
BAMS=(); for r in "${RUNS[@]}"; do
  f=$BAM/${r}.Aligned.toTranscriptome.out.bam
  [[ -s "$f" ]] || { echo "missing BAM $f" >&2; exit 1; }
  BAMS+=("$f")
done
W=$OUT/_work_$T
mkdir -p "$W"
echo "tissue=$T runs=${#BAMS[@]} pack=$(basename "$PACK") host=$(hostname) start=$(date -Iseconds)"

# ---- metaplots + RiboCode: per-sample P-site offsets, exactly as the standard path does ----
if ! ls "$W"/*_psites.hd5 >/dev/null 2>&1; then
  printf '%s\n' "${BAMS[@]}" > "$W/bam_list.txt"
  ( cd "$W" && "$RBIN/metaplots" -a "$ANNOT" -i "$W/bam_list.txt" -o "cho_$T" )
  CFG=$(ls "$W"/cho_${T}_pre_config.txt 2>/dev/null | head -1)
  [[ -s "$CFG" ]] || { echo "ABORT: metaplots wrote no config for $T" >&2; exit 1; }
  grep -cE "^[^#]" "$CFG" | sed 's/^/  read-length rows selected: /'
  ( cd "$W" && "$RBIN/RiboCode" -a "$ANNOT" -c "$CFG" -l no -g -o "cho_$T" )
fi

# ---- pool onto the reference pack's transcript axis and measure the divergence ----
"$PY" - <<PYEOF
import numpy as np, glob, json, pathlib, sys
sys.path.insert(0, "$NEW/scripts")
from prepare.packlib import pool_per_nt          # same pooling the packs were built with
pack = pathlib.Path("$PACK")
order = [l.strip() for l in open(pack / "tx_order.txt")]
lens = np.load(pack / "lengths.npy")
hd5 = sorted(glob.glob("$W/*_psites.hd5"))
print(f"  pooling {len(hd5)} psites hd5 onto {len(order):,} tx")
if not hd5:
    print("ABORT: RiboCode produced no psites hd5", file=sys.stderr); raise SystemExit(1)
pooled = pool_per_nt(hd5, "p_sites", set(order))   # key verified vs prepare_pack.py:122
lens_by_tx = dict(zip(order, np.asarray(lens, dtype=np.int64)))
missing = [t for t in order if t not in pooled]
if missing:
    print(f"  {len(missing):,} tx had no pooled reads; filled with zeros")
new = np.concatenate([np.asarray(pooled[t], dtype=np.int64) if t in pooled
                      else np.zeros(lens_by_tx[t], dtype=np.int64) for t in order])
old = np.load(pack / "target_counts.npy").astype(np.int64)
if new.shape != old.shape:
    print(f"ABORT: axis mismatch new={new.shape} old={old.shape}", file=sys.stderr); raise SystemExit(1)
d = np.abs(new - old)
rel = d.sum() / max(old.sum(), 1)
res = dict(tissue="$T", pack=str(pack), n_tx=len(order), n_bams=len(hd5),
           old_psites=int(old.sum()), new_psites=int(new.sum()),
           nt_differing=int(np.count_nonzero(d)), abs_diff=int(d.sum()),
           pct_of_signal=float(100 * rel),
           pearson=float(np.corrcoef(new, old)[0, 1]))
pathlib.Path("$OUT/${T}_divergence.json").write_text(json.dumps(res, indent=2))
print(f"  OLD {res['old_psites']:,} P-sites | NEW {res['new_psites']:,} | "
      f"{res['nt_differing']:,} nt differ | {res['pct_of_signal']:.4f}% of signal | "
      f"per-nt r={res['pearson']:.4f}")
PYEOF
echo "done=$(date -Iseconds)"
