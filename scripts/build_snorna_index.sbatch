#!/usr/bin/env bash
#SBATCH --job-name=snorna_idx
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/snorna_idx_%j.out
set -euo pipefail
# Build a snoRNA/scaRNA bowtie2 depletion index for the human Ribo-seq + RNA-seq path.
#
# WHY THIS EXISTS. build_contaminant_index.sbatch deliberately routed snoRNA to a POST-alignment BED
# filter rather than to up-front depletion. That leaves a hole: the salmon quant that DEFINES the
# transcript universe runs on raw FASTQ, so it never sees a post-alignment filter. Measured on
# GSE304796 (CAR-T, Ribo-Zero total RNA): 43.1% of TPM sat in `ncRNA_host` genes -- snoRNA host
# transcripts, which GENCODE types as lncRNA so no biotype filter can drop them -- against 15.3% for
# the poly(A) control GSE208041. One gene, ENSG00000282885, carried ~18-20% of the library alone.
# Removing the reads is cleaner and more consistent than excluding host genes downstream.
#
# SEPARATE INDEX, NOT A REBUILD OF contam_grch38_v49. Two reasons. (1) Jobs were actively reading the
# contaminant index when this was written, and overwriting a file another job holds open is how you
# get an output of exactly the right size that is wrong. (2) Keeping the stages separate keeps the
# rRNA/tRNA depletion rate and the snoRNA depletion rate independently measurable in the logs, which
# is what makes a regression visible later.
#
# snoRNA + scaRNA ONLY. scaRNAs are a snoRNA subclass and behave identically here. snRNA is
# deliberately EXCLUDED, preserving the original decision in build_contaminant_index.sbatch: snRNA
# loci are long enough that up-front depletion also strips legitimately overlapping mRNA reads.
# snoRNAs are 60-300 nt and sit in host introns, so they do not carry that risk. misc_RNA, vault_RNA
# and the rest stay in the post-alignment BED as before -- this changes one class, not the policy.
#
#   sbatch scripts/build_snorna_index.sbatch
GTF=/private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations/gencode.v49.annotation.gtf
GEN=/private/groups/carpenterlab/emalekos/genomes/GRCh38.primary_assembly.genome.fa
IDXDIR=/private/groups/carpenterlab/emalekos/Bowtie2_indexes
OUT=$IDXDIR/snorna_grch38_v49
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
BED=$NEW/data/ncrna_intervals
BT=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin
BEDTOOLS=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/bedtools
T=${SLURM_CPUS_PER_TASK:-12}
mkdir -p "$IDXDIR" "$BED"
for f in "$GTF" "$GEN"; do [[ -s "$f" ]] || { echo "missing $f" >&2; exit 1; }; done
[[ -s "$OUT.1.bt2" ]] && { echo "SKIP: $OUT already built"; exit 0; }
echo "host=$(hostname) start=$(date -Iseconds)"

SNO='snoRNA|scaRNA'
# gene_type, not transcript_type -- a snoRNA gene's transcripts are all snoRNA, and gene_type is the
# field every other contaminant rule in this project keys on.
awk -v pat="$SNO" -F'\t' '$3=="exon" && $9 ~ ("gene_type \"("pat")\"") {
    match($9, /gene_id "[^"]+"/); g=substr($9, RSTART+9, RLENGTH-10)
    printf "%s\t%d\t%d\t%s\t0\t%s\n", $1, $4-1, $5, g, $7 }' "$GTF" \
  | sort -k1,1 -k2,2n | "$BEDTOOLS" merge -s -c 4,6 -o distinct,distinct -i - \
  > "$BED/human_v49_snorna.bed"
N=$(wc -l < "$BED/human_v49_snorna.bed")
echo "  snoRNA/scaRNA merged intervals: $N"
[[ "$N" -gt 500 ]] || { echo "ABORT: only $N snoRNA intervals -- expected >500 in GENCODE v49" >&2; exit 1; }

"$BEDTOOLS" getfasta -s -fi "$GEN" -bed "$BED/human_v49_snorna.bed" -fo "$OUT.fa"
NF=$(grep -c '^>' "$OUT.fa")
echo "  sequences written: $NF  ($(du -h "$OUT.fa" | cut -f1))"
[[ "$NF" -eq "$N" ]] || { echo "ABORT: $NF sequences from $N intervals" >&2; exit 1; }

"$BT/bowtie2-build" --threads "$T" "$OUT.fa" "$OUT" > /dev/null
[[ -s "$OUT.1.bt2" ]] || { echo "ABORT: bowtie2-build produced no index" >&2; exit 1; }

cat > "$OUT.README" <<EOF
snorna_grch38_v49 -- snoRNA + scaRNA bowtie2 depletion index
built $(date -Iseconds) by scripts/build_snorna_index.sbatch

source GTF:  $GTF
genome:      $GEN
gene_types:  snoRNA, scaRNA   (snRNA deliberately excluded -- up-front snRNA depletion also strips
                               legitimately overlapping mRNA reads; see build_contaminant_index.sbatch)
intervals:   $N merged, strand-aware
sequences:   $NF

Used as a SECOND bowtie2 depletion stage after contam_grch38_v49 (rRNA/tRNA/miRNA/rDNA), so the two
depletion rates stay separately visible in the per-sample logs. Motivated by 43.1% of GSE304796 TPM
landing in snoRNA host genes typed as lncRNA.
EOF
echo "done=$(date -Iseconds) -> $OUT"
