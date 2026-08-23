#!/usr/bin/env bash
#SBATCH --job-name=rebuild_janich
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/rebuild_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/rebuild_%j.err
set -euo pipefail
# Rebuild the Janich mouse-liver UNIVERSE, which was built from an untrimmed salmon quant.
#
# ROOT CAUSE (diagnosed 2026-07-30): the fastqs are NOT truncated or corrupt -- salmon processed all
# 45,344,157 fragments. They are ARTseq small-insert total-RNA libraries (GEO GSM1644100: "ZT0, total
# RNA, rep1"; both the RNA and the ribosome-profiling libraries were built with the ARTseq kit), so the
# median insert is 34 nt inside a 101 nt read and 95% of reads carry the 3' TruSeq adapter.
# build_janich_universe.sbatch ran salmon on the RAW fastqs, so it mapped 11,811 / 45,344,157 = 0.03%,
# TPM collapsed onto a handful of transcripts, and the universe came out at 2,434 tx.
#
# What was NEVER affected: align_janich_rna.sbatch always ran cutadapt, so the BAMs are clean (0%
# soft-clip, 28-35M CIGARs) and the coverage hd5 used as the model input are correct. Only the universe
# (and everything derived from it) is wrong. Those are quarantined under
# data/_quarantine/janich_stale_universe_2026-07-30/ (see its README).
#
# Steps: cutadapt -> salmon (decoy-aware, all 7 runs) -> mean-TPM merge -> expressed universe -> pack +
# no-Kozak ORF track, reusing the existing (good) coverage. No GPU needed; predict/eval is a separate job.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
SALMON=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/salmon
CUTADAPT=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/cutadapt
IDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_vM38   # decoy-aware (rule)
FQ=$NEW/data/external/janich_liver_rnaseq
QD=$NEW/data/tpm/salmon_quant_janich_liver
REF=$NEW/data/heldout_refs
mkdir -p "$QD" "$NEW/logs/janich"
echo "rebuild_janich start=$(date -Iseconds) host=$(hostname)"

# ---- 1. cutadapt THEN salmon, single-end, decoy-aware, all 7 runs ----
# THE BUG THIS FIXES: build_janich_universe.sbatch ran salmon on the RAW fastqs. These are ARTseq
# small-insert libraries (median insert 34 nt in a 101 nt read), so 95% of reads carry the 3' TruSeq
# adapter and salmon mapped 0.03% (11,811 / 45,344,157). STAR was always fine because
# align_janich_rna.sbatch DOES trim (BAM CIGARs are clean 28-35M, 0% soft-clip), which is why the
# coverage/model input was never affected -- only the universe. Trim flags mirror that script exactly.
TRIM=$NEW/data/external/janich_liver_rnaseq_trimmed
mkdir -p "$TRIM"
for f in "$FQ"/SRR*.fastq.gz; do
  s=$(basename "$f" .fastq.gz)
  if [ ! -s "$TRIM/$s.fastq.gz" ]; then
    echo "  cutadapt $s"
    "$CUTADAPT" -j 16 -a AGATCGGAAGAGC -a "A{18}" --times 2 --minimum-length 20 \
      -o "$TRIM/$s.fastq.gz" "$f" > "$TRIM/$s.cutadapt.log" 2>&1 \
      || { echo "cutadapt failed on $s"; exit 1; }
    grep -E "Reads written|Total reads processed" "$TRIM/$s.cutadapt.log" | sed "s/^/    /"
  fi
  if [ -s "$QD/$s/quant.sf" ]; then echo "  SKIP salmon $s (quant exists)"; continue; fi
  echo "  salmon $s"
  "$SALMON" quant -i "$IDX" -l A -r "$TRIM/$s.fastq.gz" -p 16 \
    --seqBias --gcBias -o "$QD/$s" > "$QD/$s.log" 2>&1 || { echo "salmon failed on $s"; exit 1; }
  MR=$(grep -oE "mapped [0-9]+ \([0-9.]+%\)" "$QD/$s.log" | tail -1)
  echo "    $s $MR"
  # guard: the whole point of the rebuild is a sane mapping rate
  pct=$(echo "$MR" | grep -oE "\([0-9.]+%" | tr -d '(%')
  awk -v p="${pct:-0}" 'BEGIN{ if (p+0 < 20) { print "    WARNING: mapping rate still low ("p"%)"; } }'
done

# ---- 2. mean-TPM merge across the runs (mirrors build_macro_universes.py) ----
"$PY" - "$QD" <<'PYEOF'
import glob, os, sys
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
out = f"{qd}/janich_merged_quant.sf"
with open(out, "w") as w:
    w.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
    for name, (L, eff, tsum, nr) in data.items():
        w.write(f"{name}\t{L}\t{eff}\t{tsum / max(n,1):.6f}\t{nr:.3f}\n")
det = sum(1 for v in data.values() if v[2] / max(n,1) > 0)
ge1 = sum(1 for v in data.values() if v[2] / max(n,1) >= 1)
print(f"  merged {n} runs: {len(data):,} tx; detected {det:,}; TPM>=1 {ge1:,} -> {out}")
PYEOF

# ---- 3. expressed universe (pc + lncRNA, TPM>=1, <=10kb, no chrM) ----
"$PY" "$NEW/proteogenomics/scripts/build_line_universe.py" \
  --salmon "$QD/janich_merged_quant.sf" --species mouse \
  --out_prefix "$REF/mouse_janich" --min_tpm 1.0 --max_len 10000
NTX=$(wc -l < "$REF/mouse_janich_universe_tx.txt")
echo "  new Janich universe: $NTX tx (was 2,434 from the broken quant)"
[ "$NTX" -gt 8000 ] || { echo "ABORT: universe still small ($NTX); investigate before proceeding"; exit 1; }

# ---- 4. pack + no-Kozak ORF track on the NEW universe ----
"$PY" "$NEW/scripts/prepare/prepare_pack.py" \
  --ribo-psites "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/results/merged_liver_ribocode/merged_liver" \
  --rna-coverage "$NEW/data/heldout_rnaseq_coverage/mouse_janich_liver" \
  --universe-tx "$REF/mouse_janich_universe_tx.txt" \
  --fasta "$REF/mouse_janich_universe.fa" \
  --tx2biotype "$NEW/data/tx2biotype_mouse.tsv" \
  --build-orf-track --kozak none \
  --group janich_liver --species mouse \
  --out "$NEW/data/packed_heldout_mouse_janich_liver"
echo "  pack: $(wc -l < "$NEW/data/packed_heldout_mouse_janich_liver/tx_order.txt") tx"
ls -la "$NEW/data/packed_heldout_mouse_janich_liver/orf_track_v2_nokozak.npy"   # named by kozak setting
echo "rebuild_janich done=$(date -Iseconds)"
