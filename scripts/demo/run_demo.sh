#!/usr/bin/env bash
# End-to-end RiboSignal demo on public data: RNA-seq in, ORF calls out, scored against
# real Ribo-seq from the same donor.
#
#   bash scripts/demo/run_demo.sh
#
# LINUX. Do not run this on macOS: bioconda's osx-64 STAR 2.7.11b aborts on any
# --quantMode, and 2.7.10b accepts --readFilesCommand but silently yields ZERO reads and
# still exits 0. Both were reproduced; the second is the dangerous one.
#
# Resumable: every stage skips if its output exists. Delete the output to redo a stage.
# Every stage ends in an assertion, because the failure mode throughout this pipeline is
# a plausible-looking empty or shifted file rather than a crash.
#
# ---------------------------------------------------------------------------------------
# The data. Chothani et al., GEO GSE182371 (Ribo-seq) + GSE182372 (RNA-seq), Hepatocytes_1
# -- one donor, both assays. Hepatocytes is the tissue HELD OUT of training for both
# released checkpoints, so this is a genuine held-out test and not a memorised one.
#
#   SRR15513269   RNA-seq,  paired, 22,062,488 pairs, 76 nt   -> the model's input
#   SRR15513208   Ribo-seq, single, 182,236,906 reads         -> the ground truth
#
# MEASURED, not assumed: SRR15513269 is ALREADY adapter-trimmed (TruSeq in 4 reads per
# 400,000; lengths ragged 74-76). It gets no -a. SRR15513208 has NOT been measured -- the
# download did not finish -- so step 3 measures it and step 5 acts on the answer. Do not
# skip that.
# ---------------------------------------------------------------------------------------
set -euo pipefail

REPO=${RIBO_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}
WORK=${WORK:-$PWD/ribosignal_demo}
THREADS=${THREADS:-16}
# chr22 is the fast path: 1,747 genes / 11,615 transcripts, index ~400 MB, ~1 GB RAM.
# CHROM=all uses the whole primary assembly and needs ~32 GB RAM for genomeGenerate.
CHROM=${CHROM:-chr22}
GENCODE=${GENCODE:-49}
PY=${PY:-python}

R=$WORK/ref; F=$WORK/fastq; L=$WORK/logs
mkdir -p "$R" "$F" "$L" "$WORK"/{rna,ribo,pack,pred,calls}
cd "$WORK"
say() { printf '\n\033[1m=== %s\033[0m\n' "$*"; }
die() { printf '\nFAILED: %s\n' "$*" >&2; exit 1; }

# =========================================================================================
say "0. tools"
# The project's own image has STAR/samtools/cutadapt/RiboCode/torch/mamba_ssm at the exact
# pinned versions AND the weights baked in, so prefer it where Docker exists:
#   docker run --rm -it -v "$WORK:$WORK" -w "$WORK" ghcr.io/ericmalekos/riboseq-model:latest
# It has no gffread and no ribotish; install those alongside.
for t in STAR samtools cutadapt gffread; do command -v $t >/dev/null || die "$t not on PATH"; done
$PY -c "import torch, numpy, pysam, h5py" || die "need torch numpy pysam h5py"
STAR --version; samtools --version | head -1; cutadapt --version

# =========================================================================================
say "1. reference (GENCODE v$GENCODE, $CHROM)"
GTF=$R/$CHROM.gtf; FA=$R/$CHROM.fa; TX=$R/transcripts.fa; IDX=$R/star_index
if [ ! -s "$GTF" ]; then
  [ -s "$R/gencode.gtf.gz" ] || curl -sSL -o "$R/gencode.gtf.gz" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$GENCODE/gencode.v$GENCODE.annotation.gtf.gz"
  [ -s "$R/genome.fa.gz" ] || curl -sSL -o "$R/genome.fa.gz" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$GENCODE/GRCh38.primary_assembly.genome.fa.gz"
  if [ "$CHROM" = all ]; then
    gunzip -c "$R/genome.fa.gz" > "$FA"; gunzip -c "$R/gencode.gtf.gz" > "$GTF"
  else
    gunzip -c "$R/genome.fa.gz" | awk -v c=">$CHROM" '/^>/{p=($1==c)} p' > "$FA"
    gunzip -c "$R/gencode.gtf.gz" | awk -F'\t' -v c="$CHROM" '$1==c' > "$GTF"
  fi
  samtools faidx "$FA"
fi
[ -s "$TX" ] || { gffread -w "$R/tx.raw" -g "$FA" "$GTF"
                  awk '/^>/{print $1; next}{print}' "$R/tx.raw" > "$TX"; rm -f "$R/tx.raw"; }
NTX=$(grep -c '^>' "$TX"); echo "  transcripts: $NTX"
[ "$NTX" -gt 100 ] || die "only $NTX transcripts in $TX"

if [ ! -s "$IDX/SA" ]; then
  GLEN=$(awk '{s+=$2} END{print s}' "$FA.fai")
  NB=$($PY -c "import math;print(min(14,int(math.log2($GLEN)/2-1)))")
  echo "  genome $GLEN bp, --genomeSAindexNbases $NB, --sjdbOverhang 75 (76 nt reads)"
  STAR --runMode genomeGenerate --runThreadN "$THREADS" --genomeDir "$IDX" \
       --genomeFastaFiles "$FA" --sjdbGTFfile "$GTF" --sjdbOverhang 75 \
       --genomeSAindexNbases "$NB" --outFileNamePrefix "$IDX/" > "$L/index.log" 2>&1
fi
[ -s "$IDX/geneInfo.tab" ] || die "index has no transcriptome tables"

# =========================================================================================
say "2. fastq"
ena() { # ena <SRR> <path-suffix...>
  local s=$1; shift
  for f in "$@"; do
    [ -s "$F/$f" ] && continue
    curl -sSL -o "$F/$f" "https://ftp.sra.ebi.ac.uk/vol1/fastq/${s:0:6}/0${s: -2}/$s/$f"
  done
}
ena SRR15513269 SRR15513269_1.fastq.gz SRR15513269_2.fastq.gz
ena SRR15513208 SRR15513208.fastq.gz
ls -la "$F"

# =========================================================================================
say "3. MEASURE the adapters (never default them)"
$PY "$REPO/scripts/demo/measure_adapter.py" "$F/SRR15513269_1.fastq.gz" | tee "$L/adapter_rna.txt"
$PY "$REPO/scripts/demo/measure_adapter.py" "$F/SRR15513208.fastq.gz"   | tee "$L/adapter_ribo.txt"
# RNA verdict is known and asserted; the Ribo verdict is read out of the file.
grep -q "already trimmed" "$L/adapter_rna.txt" || \
  die "SRR15513269 now shows an adapter; it did not before. Read $L/adapter_rna.txt and set RNA_ADAPTER."
RIBO_ADAPTER=$(awk '/Use `cutadapt -a /{print $4}' "$L/adapter_ribo.txt" | head -1)
echo "  ribo adapter -> ${RIBO_ADAPTER:-<none detected>}"

# =========================================================================================
say "4. RNA-seq -> coverage"
# -m 20 only. NOT --maximum-length and NOT --discard-untrimmed: both are right for
# footprints and wrong for RNA (the latter kept 16.4% of reads on a library where the
# correct call kept 99.0%).
if [ ! -s rna/t1.fq.gz ]; then
  cutadapt --trim-n -m 20 -j "$THREADS" -o rna/t1.fq.gz -p rna/t2.fq.gz \
    "$F/SRR15513269_1.fastq.gz" "$F/SRR15513269_2.fastq.gz" > "$L/cutadapt_rna.log" 2>&1
fi
if [ ! -s rna/SRR15513269.Aligned.toTranscriptome.out.bam ]; then
  STAR --genomeDir "$IDX" --readFilesIn rna/t1.fq.gz rna/t2.fq.gz --readFilesCommand zcat \
       --runThreadN "$THREADS" --outSAMtype None --quantMode TranscriptomeSAM \
       --outFilterMultimapNmax 10 --outSAMattributes NH HI AS nM \
       --outFileNamePrefix rna/SRR15513269.
fi
NIN=$(awk -F'\t' '/Number of input reads/{gsub(/ /,"",$2);print $2}' rna/SRR15513269.Log.final.out)
echo "  input reads: $NIN"; grep -E "Uniquely mapped reads %" rna/SRR15513269.Log.final.out
[ "${NIN:-0}" -gt 0 ] || die "STAR read 0 input reads -- the silent-empty-BAM failure. Check readFilesCommand."
samtools quickcheck rna/SRR15513269.Aligned.toTranscriptome.out.bam || die "bad RNA bam"

# A chromosome subset has far fewer transcripts than the whole annotation, so the guard's
# whole-transcriptome default floor of 5,000 would fire spuriously. Scale it, do not remove it.
FLOOR=$([ "$CHROM" = all ] && echo 5000 || echo $((NTX / 10)))
[ -s pack/coverage.hd5 ] || $PY "$REPO/scripts/rnaseq_coverage.py" \
    --bam rna/SRR15513269.Aligned.toTranscriptome.out.bam \
    --out pack/coverage.hd5 --sample SRR15513269 --min-covered-tx "$FLOOR"

# =========================================================================================
say "5. Ribo-seq -> P-sites (the ground truth)"
# Footprints ARE length-selected, so here -M and --discard-untrimmed are correct -- the
# opposite of step 4. EndToEnd matters: soft-clipping shifts the inferred P-site.
if [ ! -s ribo/t.fq.gz ]; then
  if [ -n "$RIBO_ADAPTER" ]; then
    cutadapt -a "$RIBO_ADAPTER" -m 20 -M 40 --discard-untrimmed -j "$THREADS" \
      -o ribo/t.fq.gz "$F/SRR15513208.fastq.gz" > "$L/cutadapt_ribo.log" 2>&1
  else
    cutadapt --trim-n -m 20 -M 40 -j "$THREADS" \
      -o ribo/t.fq.gz "$F/SRR15513208.fastq.gz" > "$L/cutadapt_ribo.log" 2>&1
  fi
fi
if [ ! -s ribo/SRR15513208.Aligned.toTranscriptome.out.bam ]; then
  STAR --genomeDir "$IDX" --readFilesIn ribo/t.fq.gz --readFilesCommand zcat \
       --runThreadN "$THREADS" --outSAMtype None --quantMode TranscriptomeSAM \
       --outFilterMultimapNmax 1 --alignEndsType EndToEnd \
       --outSAMattributes NH HI AS nM --outFileNamePrefix ribo/SRR15513208.
fi
grep -E "Number of input reads|Uniquely mapped reads %" ribo/SRR15513208.Log.final.out
[ -s pack/psites.hd5 ] || $PY "$REPO/scripts/ribo_psites.py" \
    --bam ribo/SRR15513208.Aligned.toTranscriptome.out.bam --gtf "$GTF" \
    --out pack/psites.hd5 --sample SRR15513208 --lengths 25:35 2>&1 | tee "$L/psites.log"
# The frame-0 fraction is the whole QC. 1/3 is noise; a real library reaches 0.5-0.7.
# ribo_psites.py already exits non-zero below --min-frame0; this echoes it for the record.
grep -E "frame-0 fraction" "$L/psites.log" || true

# =========================================================================================
say "6. pack + ORF track"
[ -s pack/demo/tx_order.txt ] || $PY "$REPO/scripts/build_pack.py" \
    --fasta "$TX" --coverage pack/coverage.hd5 --psites pack/psites.hd5 --out pack/demo
[ -s pack/demo/orf_track_v2_nokozak.npy ] || $PY "$REPO/scripts/build_orf_track.py" \
    --pack pack/demo --fasta "$TX" --mode ext --kozak none
export RIBO_ONEHOT_FASTA=$TX
export RIBO_ORF_TRACK=$WORK/pack/demo/orf_track_v2_nokozak.npy

# =========================================================================================
say "7. weights"
if [ ! -s weights/attn_best.pt ]; then
  $PY -c "from huggingface_hub import snapshot_download; snapshot_download('emalek/RiboSignal', local_dir='weights')"
fi
(cd weights && sha256sum -c SHA256SUMS)

# =========================================================================================
say "8. predict, both checkpoints"
# On a GPU box set RIBO_MAMBA_IMPL=cuda to REFUSE the pure-PyTorch fallback, so a missing
# mamba_ssm is an error instead of a silent 100x slowdown.
for ARCH in attn mamba4; do
  [ -s "pred/$ARCH/pred_profiles.npz" ] && continue
  $PY "$REPO/scripts/dump_pred_profiles.py" --run weights/ --arch "$ARCH" \
      --pack pack/demo --out "pred/$ARCH" \
      --device "$([ -n "${CUDA_VISIBLE_DEVICES:-}" ] && echo cuda || echo cpu)"
done

# =========================================================================================
say "9. RiboCode annotation + ORF calls"
[ -s calls/annot/transcripts.pickle ] || \
  prepare_transcripts -g "$GTF" -f "$FA" -o calls/annot
for ARCH in attn mamba4; do
  P=pred/$ARCH/pred_profiles.npz
  # real  = the caller on OBSERVED P-sites, the reference the others are judged against
  # pred_preddepth = fully de novo: predicted shape AND predicted depth, no Ribo-seq used
  # + the Poisson dial, which trades ~6 points of CDS recall for uORF/non-canonical precision
  [ -s "calls/$ARCH/real_collapsed.txt" ] || $PY "$REPO/scripts/ribocode_dropin.py" \
      --profiles "$P" --annot calls/annot --variant real --out "calls/$ARCH" --min_aa 30 --pval 0.05
  [ -s "calls/$ARCH/pred_preddepth_collapsed.txt" ] || $PY "$REPO/scripts/ribocode_dropin.py" \
      --profiles "$P" --annot calls/annot --variant pred_preddepth --out "calls/$ARCH" --min_aa 30 --pval 0.05
  [ -s "calls/${ARCH}_poisson/pred_preddepth_collapsed.txt" ] || $PY "$REPO/scripts/ribocode_dropin.py" \
      --profiles "$P" --annot calls/annot --variant pred_preddepth --pred_poisson --pred_scale 0.05 \
      --out "calls/${ARCH}_poisson" --min_aa 30 --pval 0.05
done

# =========================================================================================
say "10. Ribo-TISH, as a second caller"
# Ribo-TISH only skips its bam path when the profile covers every transcript of a gene, so
# the GTF must be restricted to exactly the scored transcripts.
if command -v ribotish >/dev/null; then
  cut -f1 pack/demo/pack_meta.tsv | tail -n +2 > calls/tx_ids.txt
  awk -v FL=calls/tx_ids.txt 'BEGIN{while((getline l < FL)>0) k[l]=1}
       !/^#/ { if (match($0,/transcript_id "[^"]+"/)) {
         t=substr($0,RSTART+15,RLENGTH-16); if (t in k) print } }' "$GTF" > calls/scored.gtf
  for ARCH in attn mamba4; do
    [ -s "calls/ribotish_$ARCH/pred_preddepth.txt" ] || $PY "$REPO/scripts/orfcallers/ribotish_dropin.py" \
        --profiles "pred/$ARCH/pred_profiles.npz" --gtf calls/scored.gtf --genome "$FA" \
        --variant pred_preddepth --out "calls/ribotish_$ARCH" \
        --longest --minaalen 5 --fpth 0.05 --numproc "$THREADS"
  done
else
  echo "  ribotish not on PATH -- skipping (pip install ribotish, or its own env via --ribotish)"
fi

# =========================================================================================
say "11. score: predicted vs observed"
$PY "$REPO/scripts/demo/score_demo.py" --pred-dir pred --calls-dir calls \
    --out "$WORK/RESULTS.json" | tee "$WORK/RESULTS.txt"

say "done -- $WORK/RESULTS.txt"
