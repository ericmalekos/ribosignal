#!/usr/bin/env bash
#SBATCH --job-name=l3x3_valid
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=110G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/l3x3_valid_%j.out
set -euo pipefail
# STEP 0 of the mouse-liver 3x3 factorial, and it BLOCKS everything else.
#
# The factorial needs Janich and GSE243134 Ribo counts on universes their existing packs do not cover,
# so both must be re-pooled from archived BAMs. Wang has no BAMs and is sliced from its
# full-annotation pack instead. Two code paths, one table -- unless they are shown to agree.
#
# This re-pools Janich from its 7 BAMs onto its OWN universe, where an authoritative pack already
# exists, and requires target_counts.npy to match. If they disagree, the factorial would be measuring
# the difference between two P-site pipelines rather than between three experiments, and every number
# downstream would be void. Nothing else runs until this passes.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
ANNOT="/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/mouse_ribocode_annot"
REF=$NEW/data/heldout_refs
GOLD=$NEW/data/packed_heldout_mouse_janich_liver_decon
OUT=$NEW/data/liver3x3/_validate_janich_own
mkdir -p "$(dirname "$OUT")"
echo "host=$(hostname) start=$(date -Iseconds)"

# mouse_janich_liver/ holds the RNA BAMs (SRR19302xx) -- the RIBO BAMs are the re-fetched
# SRR19301xx in the arm-suffixed dir. Using the wrong one would pool RNA through the Ribo path.
mapfile -t BAMS < <(ls $NEW/data/heldout_bam/mouse_janich_liver_ribo/*.bam 2>/dev/null | sort)
[[ ${#BAMS[@]} -ge 1 ]] || { echo "no Janich BAMs" >&2; exit 1; }
echo "Janich Ribo BAMs: ${#BAMS[@]}"
printf '  %s\n' "${BAMS[@]##*/}"

"$PY" "$NEW/scripts/prepare/prepare_from_bams.py" \
  --ribo-bam "${BAMS[@]}" \
  --rna-coverage "$NEW/data/heldout_rnaseq_coverage/mouse_janich_liver" \
  --universe-tx "$REF/mouse_janich_decon_universe_tx.txt" \
  --annot "$ANNOT" --tx2biotype "$NEW/data/tx2biotype_mouse.tsv" \
  --ribocode-bin /private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin \
  --group validate_janich_own --species mouse --out "$OUT" 2>&1 | tail -12

"$PY" - <<PYEOF
import numpy as np, sys
g="$GOLD"; o="$OUT"
go=[l.strip() for l in open(g+"/tx_order.txt")]; oo=[l.strip() for l in open(o+"/tx_order.txt")]
if go!=oo:
    print(f"FAIL tx_order differs ({len(go)} vs {len(oo)})"); sys.exit(1)
a=np.load(g+"/target_counts.npy"); b=np.load(o+"/target_counts.npy")
print(f"gold  target_counts: {a.sum():,} P-sites over {len(a):,} nt")
print(f"rebuilt              {b.sum():,} P-sites over {len(b):,} nt")
if np.array_equal(a,b):
    print("PASS: target_counts byte-identical -- the BAM path reproduces the historical pack"); sys.exit(0)
d=np.abs(a.astype(np.int64)-b.astype(np.int64))
rel=d.sum()/max(a.sum(),1)
print(f"DIFFER: {np.count_nonzero(d):,} nt differ, total abs diff {d.sum():,} ({100*rel:.4f}% of signal)")
# A tiny difference can come from RiboCode offset ties; a large one means a different pipeline.
if rel < 0.001:
    print("WARN: <0.1% divergence -- usable, but record it and do not call the paths identical")
    sys.exit(0)
print("FAIL: the two P-site paths disagree materially; do not build the factorial on them")
sys.exit(1)
PYEOF
echo "done=$(date -Iseconds)"
