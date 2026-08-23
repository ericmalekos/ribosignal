#!/usr/bin/env bash
#SBATCH --job-name=huvec_filt
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH --time=11:45:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/huvec_filt_%j.out
set -euo pipefail
# Third arm of the HUVEC reproducibility ablation: the ORIGINAL recipe, filter step included.
#
# The original recipe was never lost -- it is documented in build_psite_target.py's own docstring:
#   "Each input hd5 already has per-sample per-read-length P-site offsets + sense-only strand baked
#    in (RiboCode metaplots/pre_config), and rRNA/tRNA/miRNA transcripts + cross-gene multimappers
#    removed; within-gene isoform multimappers are retained and each counted (posture A)."
# and implemented in filter_tx_heldout.py, which is "byte-for-byte the same logic" as the
# expression_context_human filter that produced the training inputs.
#
# My Chothani reconstruction skipped that filter entirely, which is why every tissue came out with
# ~5% MORE P-sites than the pack. The three arms isolate each deviation:
#   A  my reconstruction  (bowtie2 contam + snoRNA + EndToEnd, NO filter)  -> 11.69%, r=0.9350
#   B  minimal            (no bowtie2, no EndToEnd, NO filter)             -> job 36750912
#   C  minimal + FILTER   (this job)                                        -> the real test
#
# If C lands near the 1.69% the Janich gate achieved, the training data IS reproducible and the
# earlier "cannot regenerate to better than 13%" was my recipe drift, not a property of the data.
#
#   sbatch scripts/huvec_filtered_test.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
RPY=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/python3
RBIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
ANNOT=$NEW/data/human_ribocode_annot_primary
PACK=$NEW/data/packed_union_HUVEC
SRCB=$NEW/results/chothani_regeneration/_minimal_HUVEC/bam
W=$NEW/results/chothani_regeneration/_filtered_HUVEC
BAMDIR=$W/bam
NCRNA=$NEW/data/ncrna_tx_human_v49.txt
T2G=$NEW/data/tx_to_gene_human_v49.tsv
T=${SLURM_CPUS_PER_TASK:-12}
mkdir -p "$BAMDIR"
for f in "$NCRNA" "$T2G" "$PACK/target_counts.npy"; do
  [[ -s "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
# Depends on arm B having produced the minimal BAMs; they are COPIED, never filtered in place,
# because filter_tx_heldout.py overwrites its input and arm B's numbers must stay reproducible.
ls "$SRCB"/*.Aligned.toTranscriptome.out.bam >/dev/null 2>&1 \
  || { echo "arm B BAMs not ready at $SRCB -- run huvec_pipeline_alignment_test.sbatch first" >&2; exit 1; }
echo "host=$(hostname) start=$(date -Iseconds)"

for src in "$SRCB"/*.Aligned.toTranscriptome.out.bam; do
  b=$(basename "$src"); dst=$BAMDIR/$b
  [[ -s "$dst" ]] && { echo "  skip $b"; continue; }
  cp "$src" "$dst"
  before=$("$ENV/bin/samtools" view -@ "$T" -c "$dst")
  "$RPY" "$NEW/scripts/heldout/filter_tx_heldout.py" "$dst" "$NCRNA" "$T2G"
  after=$("$ENV/bin/samtools" view -@ "$T" -c "$dst")
  echo "  $b: $before -> $after records ($(awk -v a=$before -v b=$after 'BEGIN{printf "%.2f", 100*(a-b)/a}')% dropped)"
done

if ! ls "$W"/*_psites.hd5 >/dev/null 2>&1; then
  ls "$BAMDIR"/*.Aligned.toTranscriptome.out.bam > "$W/bam_list.txt"
  ( cd "$W" && "$RBIN/metaplots" -a "$ANNOT" -i "$W/bam_list.txt" -o filt_HUVEC )
  CFG=$(ls "$W"/filt_HUVEC_pre_config.txt | head -1)
  ( cd "$W" && "$RBIN/RiboCode" -a "$ANNOT" -c "$CFG" -l no -g -o filt_HUVEC )
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
res = dict(tissue="HUVEC", recipe="minimal + filter_tx_heldout (ORIGINAL recipe)", n_bams=len(hd5),
           old_psites=int(old.sum()), new_psites=int(new.sum()),
           pct_of_signal=float(100 * d.sum() / max(old.sum(), 1)),
           pearson=float(np.corrcoef(new, old)[0, 1]))
pathlib.Path("$W/divergence_filtered.json").write_text(json.dumps(res, indent=2))
print(f"  ARM C (original recipe): old {res['old_psites']:,} | new {res['new_psites']:,} | "
      f"{res['pct_of_signal']:.2f}% of signal | per-nt r={res['pearson']:.4f}")
print(f"  ARM A (my reconstruction): 11.69%, r=0.9350")
PYEOF
echo "done=$(date -Iseconds)"
