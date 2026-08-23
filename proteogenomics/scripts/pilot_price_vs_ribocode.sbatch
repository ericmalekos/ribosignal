#!/usr/bin/env bash
#SBATCH --job-name=hbl1_pricevrc
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/hbl1_pricevrc_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/hbl1_pricevrc_%j.err
set -euo pipefail
RS=/private/groups/carpenterlab/emalekos/conda_envs/riboseq;  export PATH="$RS/bin:$PATH"
RC=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
GEDI=/private/groups/carpenterlab/emalekos/conda_envs/gedi/bin
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ANNOT=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/human_ribocode_annot_primary
PRICE_IDX=$ECH/data/orfcaller_smoke/price/gedi_genome/gencode.v49.annotation
NCRNA=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/ncrna_tx_human_v49.txt
TX2GENE=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/tx_to_gene_human_v49.tsv
RB=$NEW/proteogenomics/data/HBL1_pilot/riboseq
OUT=$NEW/proteogenomics/data/HBL1_pilot/pricevrc; mkdir -p "$OUT"
SCR=/data/tmp/emalekos/pricevrc_$SLURM_JOB_ID; mkdir -p "$SCR"; cd "$SCR"

echo "[$(date +%T)] subset 100M reads + trim (paper-matched: -u 2, NNNNNAGATCGGAAGAGC, 26-34)"
set +o pipefail; zcat $RB/SRR12285192.fastq.gz | head -400000000 > raw.fq; set -o pipefail
cutadapt -j12 -u 2 -a NNNNNAGATCGGAAGAGC --trim-n -m 26 -M 34 -o trim.fq raw.fq > "$OUT/cutadapt.log" 2>&1

echo "[$(date +%T)] STAR (paper params; genome SortedByCoord BAM + transcriptome BAM)"
STAR --runThreadN 12 --genomeDir "$IDX" --readFilesIn trim.fq \
  --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
  --alignEndsType EndToEnd --seedSearchStartLmax 15 \
  --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 25 --outFilterMultimapNmax 20 \
  --outSAMattributes NH HI AS nM MD --limitBAMsortRAM 40000000000 \
  --outFileNamePrefix "$SCR/hbl1." --outTmpDir "$SCR/_star" > "$OUT/star.log" 2>&1
grep -E "Uniquely mapped reads number|Uniquely mapped reads %" "$SCR/hbl1.Log.final.out" | sed 's/^/  /'
cp "$SCR/hbl1.Log.final.out" "$OUT/"

# ---------- PRICE (genome BAM; single-mapper NH==1 = MAPQ 255, matching the training recipe) ----------
echo "[$(date +%T)] PRICE: filter unique genome reads + index"
samtools view -@12 -q 255 -b "$SCR/hbl1.Aligned.sortedByCoord.out.bam" > "$SCR/hbl1.genome.uniq.bam"
samtools index -@12 "$SCR/hbl1.genome.uniq.bam"
echo "  unique genome reads: $(samtools view -c -@12 $SCR/hbl1.genome.uniq.bam)"
mkdir -p "$OUT/price"
cp "$SCR/hbl1.genome.uniq.bam" "$SCR/hbl1.genome.uniq.bam.bai" "$OUT/price/"   # persist for PRICE iteration
echo "[$(date +%T)] PRICE run (-genomic grch38_v49 = REGISTERED name, not a file path -- the bug that silently killed prior runs)"
set +e
$GEDI/gedi -e Price -reads "$SCR/hbl1.genome.uniq.bam" -genomic grch38_v49 \
    -nthreads 12 -D -progress -prefix "$OUT/price/hbl1" > "$OUT/price/price.log" 2>&1
echo "  PRICE exit=$?"
echo "  PRICE output files:"; ls -la "$OUT/price/" | awk '{print $5,$9}'

# ---------- RiboCode (transcriptome BAM, cleaned) ----------
echo "[$(date +%T)] RiboCode: clean transcriptome BAM"
cp "$SCR/hbl1.Aligned.toTranscriptome.out.bam" "$SCR/hbl1.tx.bam"
$RS/bin/python3 $NEW/proteogenomics/scripts/filter_tx_bam_generic.py "$SCR/hbl1.tx.bam" "$NCRNA" "$TX2GENE" >> "$OUT/star.log" 2>&1
mkdir -p "$OUT/ribocode"; echo "$SCR/hbl1.tx.bam" > "$OUT/ribocode/bamlist.txt"
echo "[$(date +%T)] RiboCode metaplots + ORF call"
"$RC/metaplots" -a "$ANNOT" -i "$OUT/ribocode/bamlist.txt" -o "$OUT/ribocode/hbl1" > "$OUT/ribocode/metaplots.log" 2>&1
CONFIG="$OUT/ribocode/hbl1_pre_config.txt"
grep -E "^# [0-9]+\s" "$CONFIG" 2>/dev/null | awk -F'\t' '{printf "  periodicity len=%s offset=%s f0=%s\n",$1,$3,$7}' | head
"$RC/RiboCode" -a "$ANNOT" -c "$CONFIG" -l no -g -o "$OUT/ribocode/hbl1" > "$OUT/ribocode/ribocode.log" 2>&1

# ---------- compare ----------
echo; echo "########## PRICE vs RiboCode (same 100M-read subset alignment) ##########"
PT=$(ls "$OUT/price/hbl1"*.orfs.tsv 2>/dev/null | head -1)
echo "PRICE orfs.tsv: ${PT:-NONE PRODUCED}"
[ -n "${PT:-}" ] && { echo "  PRICE ORFs: $(($(wc -l < "$PT")-1))"; head -1 "$PT" | tr '\t' '\n' | grep -in "type\|orf\|codon" | head; }
RCT="$OUT/ribocode/hbl1_collapsed.txt"
echo "RiboCode collapsed: $([ -s "$RCT" ] && echo "$(($(wc -l < "$RCT")-1)) ORFs" || echo NONE)"
[ -s "$RCT" ] && awk -F'\t' 'NR>1{print $2}' "$RCT" | sort | uniq -c | sort -rn | head
echo "[$(date +%T)] DONE"
rm -rf "$SCR"