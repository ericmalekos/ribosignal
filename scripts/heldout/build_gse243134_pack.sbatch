#!/usr/bin/env bash
#SBATCH --job-name=pack_243134
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=140G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/gse243134/pack_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/gse243134/pack_%j.err
set -euo pipefail
# THIRD LIVER ARM, stage 2 (rebuilt 2026-07-31).
#
#   Ribo-seq  : GSE243134, 21 WT RPF libraries  (stage 1 output, unchanged, good)
#   RNA-seq   : GSE302188, 4 "Liver, 3 Months" poly(A) replicates  <-- REPLACEMENT
#
# GSE243134's own totalRNA (Ribo-Zero) arm was REJECTED: salmon 17.7% median, 93.6% of TPM absorbed by
# ten tRNA/7SL structural transcripts (typed `lncRNA`, so the biotype filter cannot drop them), TPM>=1
# = 3,534, universe 3,321 tx vs Wang 13,041. Alignment and trimming were correct; the library chemistry
# was wrong. See data/external/mouse_liver_polya_rnaseq/README.md and methods.md 7.2.
#
# CAVEAT: this RNA is NOT donor-matched to the Ribo-seq (different study), unlike the Wang and Janich
# arms. Liver is strongly circadian and the RPF is ZT5/ZT12 while GSE302188 records no time of day. The
# arm's primary contribution is an independent set of OBSERVED RiboCode calls, which needs the universe
# but not donor-matched per-nt RNA; the model-prediction sub-arm inherits the mismatch.
#
# The GSE243134 _340_/_347_ flowcell pairs are pooled here: prepare_pack sums every psites/coverage file
# in the given directory, so per-run files from one library add together correctly.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
RIBOPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
SALMON=$ENV/bin/salmon
SIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_vM38
STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
D=$NEW/data/external/gse243134_mouse_liver
RD=$NEW/data/external/mouse_liver_polya_rnaseq
B=$NEW/data/heldout_bam/mouse_gse243134_liver
QD=$NEW/data/tpm/salmon_quant_gse302188_liver
COVDIR=$NEW/data/heldout_rnaseq_coverage/mouse_gse302188_liver
REF=$NEW/data/heldout_refs
mkdir -p "$QD" "$COVDIR"
echo "start=$(date -Iseconds) host=$(hostname)"

RNARUNS="SRR34449471 SRR34449553 SRR34449454 SRR34449569"

# ---- 1. salmon on the poly(A) RNA runs (PE) ----
for r in $RNARUNS; do
  [ -s "$QD/$r/quant.sf" ] && { echo "  SKIP salmon $r"; continue; }
  T=/data/tmp/emalekos/gse302188_q/$r; mkdir -p "$T"
  "$ENV/bin/cutadapt" -j 16 -a AGATCGGAAGAGC -A AGATCGGAAGAGC --minimum-length 20 \
    -o "$T/1.fq.gz" -p "$T/2.fq.gz" "$RD/fastq/${r}_1.fastq.gz" "$RD/fastq/${r}_2.fastq.gz" > "$T/cut.log" 2>&1
  "$SALMON" quant -i "$SIDX" -l A -1 "$T/1.fq.gz" -2 "$T/2.fq.gz" -p 16 --seqBias --gcBias \
    -o "$QD/$r" > "$QD/$r.log" 2>&1 || { echo "salmon failed $r"; exit 1; }
  rm -rf "$T"
done
echo "--- salmon mapping rates (poly(A); Wang 91-93%, REJECTED GSE243134 totalRNA 9.75-23.64%) ---"
"$PY" - "$QD" <<'PYEOF'
import glob, json, sys
rs = []
for m in sorted(glob.glob(f"{sys.argv[1]}/SRR*/aux_info/meta_info.json")):
    rs.append(json.load(open(m))["percent_mapped"])
print(f"    min={min(rs):.2f}% median={sorted(rs)[len(rs)//2]:.2f}% max={max(rs):.2f}% n={len(rs)}")
PYEOF

# ---- 2. mean-TPM merge + expressed universe ----
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
out = f"{qd}/gse302188_merged_quant.sf"
with open(out, "w") as w:
    w.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
    for name, (L, eff, t, nr) in data.items():
        w.write(f"{name}\t{L}\t{eff}\t{t/max(n,1):.6f}\t{nr:.3f}\n")
tpm1 = sum(1 for v in data.values() if v[2]/max(n,1) >= 1)
print(f"  merged {n} runs: {len(data):,} tx; TPM>=1 {tpm1:,} -> {out}")
PYEOF
"$PY" "$NEW/proteogenomics/scripts/build_line_universe.py" \
  --salmon "$QD/gse302188_merged_quant.sf" --species mouse \
  --out_prefix "$REF/mouse_gse243134" --min_tpm 1.0 --max_len 10000
NTX=$(wc -l < "$REF/mouse_gse243134_universe_tx.txt")
echo "  universe: $NTX tx   (Wang 13,041 / Janich 10,431 / rejected totalRNA 3,321)"
# Guard that caught the totalRNA failure before it reached the eval. Keep it on every new arm.
[ "$NTX" -gt 6000 ] || { echo "ABORT: universe implausibly small ($NTX)"; exit 1; }

# ---- 3. STAR the poly(A) RNA runs -> per-nt coverage hd5 ----
# --outFilterMultimapNmax 1 matches every other RNA-coverage build in this project.
for r in $RNARUNS; do
  OUTCOV=$COVDIR/${r}_coverage.hd5
  [ -s "$OUTCOV" ] && { echo "  SKIP coverage $r"; continue; }
  S=/data/tmp/emalekos/gse302188_star/$r; rm -rf "$S"; mkdir -p "$S"
  "$ENV/bin/cutadapt" -j 16 -a AGATCGGAAGAGC -A AGATCGGAAGAGC --minimum-length 20 \
    -o "$S/1.fq.gz" -p "$S/2.fq.gz" "$RD/fastq/${r}_1.fastq.gz" "$RD/fastq/${r}_2.fastq.gz" > "$S/cut.log" 2>&1
  "$ENV/bin/STAR" --genomeDir "$STARIDX" --readFilesIn "$S/1.fq.gz" "$S/2.fq.gz" --readFilesCommand zcat \
    --runThreadN 16 --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM \
    --outFilterMultimapNmax 1 --outSAMattributes NH HI AS nM --outFileNamePrefix "$S/${r}."
  grep -E "Uniquely mapped reads %" "$S/${r}.Log.final.out" | sed "s/^/    $r /"
  "$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$S/${r}.Aligned.toTranscriptome.out.bam" "$OUTCOV" "$r"
  rm -rf "$S"
done

# ---- 4+5. psites + pack via the BYO-BAM front-end (prepare_from_bams chains:
#      filter_tx_heldout -> samtools index -> RiboCode metaplots -> RiboCode -g -> *_psites.hd5
#      -> prepare_pack). Ribo BAMs are GSE243134 stage 1 output, unchanged. ----
export RIBOCODE_BIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
export RIBO_PYSAM_PYTHON=$RIBOPY
export RIBO_SAMTOOLS=$ENV/bin/samtools
RIBOBAMS=$(awk -F'\t' 'NR>1 && $2 ~ /RPF/{print $1}' "$D/manifest_wt.tsv" \
           | while read r; do [ -s "$B/${r}.Aligned.toTranscriptome.out.bam" ] && echo "$B/${r}.Aligned.toTranscriptome.out.bam"; done)
COVS=$(ls "$COVDIR"/*_coverage.hd5)
echo "  ribo BAMs (GSE243134): $(echo "$RIBOBAMS" | wc -l)   rna coverage (GSE302188): $(echo "$COVS" | wc -l)"

# shellcheck disable=SC2086
"$PY" "$NEW/scripts/prepare/prepare_from_bams.py" \
  --ribo-bam $RIBOBAMS \
  --rna-coverage $COVS \
  --universe-tx "$REF/mouse_gse243134_universe_tx.txt" \
  --fasta "$REF/mouse_gse243134_universe.fa" \
  --annot "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/mouse_ribocode_annot" \
  --ncrna-tx "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/ncrna_tx_mouse_vM38.txt" \
  --tx-to-gene "$NEW/data/heldout_refs/mouse_tx_to_gene.tsv" \
  --tx2biotype "$NEW/data/tx2biotype_mouse.tsv" \
  --build-orf-track --kozak none \
  --group gse243134_liver --species mouse --threads 16 \
  --workdir /data/tmp/emalekos/gse243134_pack \
  --out "$NEW/data/packed_heldout_mouse_gse243134_liver"

P=$NEW/data/packed_heldout_mouse_gse243134_liver
echo "  pack: $(wc -l < "$P/tx_order.txt") tx"
ls -la "$P"/orf_track_v2_nokozak.npy
echo
echo "===== THREE LIVER ARMS SIDE BY SIDE ====="
for u in mouse_wang_expressed mouse_janich mouse_gse243134; do
  [ -s "$REF/${u}_universe_tx.txt" ] && printf "  %-22s universe %6d tx\n" "$u" "$(wc -l < "$REF/${u}_universe_tx.txt")"
done
echo "done=$(date -Iseconds)"
