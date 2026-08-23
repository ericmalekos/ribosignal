#!/usr/bin/env bash
#SBATCH --job-name=quant_wang
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/quant_wang_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/quant_wang_%j.err
set -euo pipefail
# Quantify the Wang mouse-liver RNA-seq so the O2 arms can be scored on a REAL expressed universe.
#
# Wang previously had NO salmon quant: its RNA was processed into per-nt coverage and the fastqs deleted,
# so its "universe" was the whole 221,835-tx vM38 annotation, unfiltered. That is the one place in the
# pipeline where the expression gate was skipped, and it is why the O2 model arms predicted across ~205k
# transcripts the reference experiments never covered.
#
# Adapter status (checked before running, per the Janich lesson): SRR5262874 is 50 nt reads with 0.1%
# TruSeq adapter -> NO trimming needed. Contrast Janich, which was 95% adapter and required cutadapt.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
SALMON=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/salmon
IDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_vM38   # decoy-aware (rule)
FQ=$NEW/data/external/wang2021_mouse/fastq/mouse_liver
QD=$NEW/data/tpm/salmon_quant_wang_liver
REF=$NEW/data/heldout_refs
mkdir -p "$QD" "$NEW/logs/janich"
echo "quant_wang start=$(date -Iseconds) host=$(hostname)"

for f in "$FQ"/SRR*.fastq.gz; do
  s=$(basename "$f" .fastq.gz)
  # re-assert the adapter check per run so a surprise library cannot slip through silently
  A=$(zcat "$f" | head -400000 | awk 'NR%4==2' | grep -c AGATCGGAAGAGC || true)
  echo "  $s adapter-positive reads: $A / 100000"
  [ "$A" -lt 50000 ] || { echo "  ABORT: $s is >50% adapter, needs cutadapt (see rebuild_janich.sbatch)"; exit 1; }
  if [ -s "$QD/$s/quant.sf" ]; then echo "  SKIP $s (quant exists)"; continue; fi
  "$SALMON" quant -i "$IDX" -l A -r "$f" -p 16 --seqBias --gcBias -o "$QD/$s" > "$QD/$s.log" 2>&1 \
    || { echo "salmon failed on $s"; exit 1; }
  echo "    $(grep -oE 'mapped [0-9]+ \([0-9.]+%\)' "$QD/$s.log" | tail -1)"
done

# mean-TPM merge across the two replicates
"$PY" - "$QD" <<'PYEOF'
import glob, sys
qd = sys.argv[1]
qs = sorted(glob.glob(f"{qd}/SRR*/quant.sf"))
data, n = {}, 0
for q in qs:
    n += 1
    with open(q) as f:
        f.readline()
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            d = data.setdefault(p[0], [p[1], p[2], 0.0, 0.0])
            d[2] += float(p[3]); d[3] += float(p[4])
out = f"{qd}/wang_merged_quant.sf"
with open(out, "w") as w:
    w.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
    for name, (L, eff, t, nr) in data.items():
        w.write(f"{name}\t{L}\t{eff}\t{t/max(n,1):.6f}\t{nr:.3f}\n")
det = sum(1 for v in data.values() if v[2]/max(n,1) > 0)
ge1 = sum(1 for v in data.values() if v[2]/max(n,1) >= 1)
print(f"  merged {n} runs: {len(data):,} tx; detected {det:,}; TPM>=1 {ge1:,} -> {out}")
PYEOF

# expressed pc+lncRNA universe, same construction as every other dataset
"$PY" "$NEW/proteogenomics/scripts/build_line_universe.py" \
  --salmon "$QD/wang_merged_quant.sf" --species mouse \
  --out_prefix "$REF/mouse_wang_expressed" --min_tpm 1.0 --max_len 10000
echo "  Wang expressed universe: $(wc -l < "$REF/mouse_wang_expressed_universe_tx.txt") tx"
echo "  (the old Wang 'universe' was 221,835 = the entire vM38 annotation, unfiltered)"
echo "quant_wang done=$(date -Iseconds)"
