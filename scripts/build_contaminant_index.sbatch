#!/usr/bin/env bash
#SBATCH --job-name=contam_idx
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/contam_idx_%j.out
set -euo pipefail
# Build the rRNA/tRNA contaminant bowtie2 index for the two-stage Ribo-seq alignment, plus the
# ncRNA-locus BED the post-STAR filter uses. Neither existed in this project before.
#
# WHY BOTH. The project's standing Ribo-seq rule allows either (a) two-stage align -- bowtie2 against
# an rRNA/tRNA/miRNA FASTA first, then STAR on the unmapped -- or (b) a post-STAR interval filter.
# The CAR-T recipe calls for (a), and (b) is still wanted for the snoRNA/canonical-chrom step after
# alignment, so this builds the inputs for both from ONE definition of "contaminant" -- otherwise the
# two stages can disagree about what they are removing.
#
# SHARED, NOT PER-PROJECT (CLAUDE.md): the index lands in a dedicated top-level Bowtie2_indexes/ dir
# named for (tool, species, annotation), so later projects reuse it instead of rebuilding.
#
# BIOTYPES. rRNA + Mt_rRNA + rRNA_pseudogene + cytoplasmic tRNA + Mt_tRNA are depleted by bowtie2 (they dominate
# a Ribo-seq library and carry no ORF information). miRNA joins them per the project rule. snoRNA,
# snRNA, scaRNA, misc_RNA and vault_RNA go into the BED for the POST-alignment filter rather than the
# depletion index -- the CAR-T recipe removes snoRNA after alignment, and depleting snRNA up front
# would also strip legitimately overlapping mRNA reads.
GTF=/private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations/gencode.v49.annotation.gtf
# CYTOPLASMIC tRNAs live in a SEPARATE GENCODE file. The main annotation has zero `gene_type "tRNA"`
# genes -- only 22 Mt_tRNA -- so a depletion index built from it alone leaves cytoplasmic tRNA
# fragments in, and those are one of the largest contaminants in a Ribo-seq library. Verified: 0 tRNA
# genes in the main GTF, 929 in the tRNAs file.
TRNA=/private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations/gencode.v49.tRNAs.gtf.gz
# THE 45S rDNA REPEAT IS NOT IN THE GENOME ASSEMBLY AS ANNOTATED rRNA. GENCODE labels only 47 genes
# `rRNA`, and the tandem rDNA arrays that produce nearly all cellular rRNA are absent from GRCh38
# primary. A GENCODE-derived index therefore misses the single largest Ribo-seq contaminant.
# Measured on GSE208041 (a library with no submitter pre-depletion): the GENCODE-only index caught
# 6.27% of reads, U13369 + 5S caught 61.68%. Without these two files that library aligned at only
# 22.5% unique with 75.8% of reads discarded as multi-mapping rDNA.
RDNA="/private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations/rdna_U13369.fa \
      /private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations/rdna_5S.fa"
GEN=/private/groups/carpenterlab/emalekos/genomes/GRCh38.primary_assembly.genome.fa
IDXDIR=/private/groups/carpenterlab/emalekos/Bowtie2_indexes
OUT=$IDXDIR/contam_grch38_v49
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
BED=$NEW/data/ncrna_intervals
BT=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin
BEDTOOLS=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/bedtools
T=${SLURM_CPUS_PER_TASK:-12}
mkdir -p "$IDXDIR" "$BED"
echo "host=$(hostname) start=$(date -Iseconds)"

DEPLETE='rRNA|Mt_rRNA|rRNA_pseudogene|Mt_tRNA|miRNA'
POSTFILT='snoRNA|snRNA|scaRNA|misc_RNA|vault_RNA|sRNA|scRNA|ribozyme'

# ---- 1. exon intervals per class, straight from the GTF ----
# gene_type, not transcript_type: a contaminant gene's transcripts are all contaminant, and gene_type
# is the field the project's existing chrM/ncRNA rules key on.
{ awk -v pat="$DEPLETE" -F'\t' '$3=="exon" && $9 ~ ("gene_type \"("pat")\"") {
    match($9, /gene_id "[^"]+"/); g=substr($9, RSTART+9, RLENGTH-10)
    printf "%s\t%d\t%d\t%s\t0\t%s\n", $1, $4-1, $5, g, $7 }' "$GTF"
  # tRNA file: feature column is "tRNA"/"exon" and every record is by definition a contaminant here
  zcat "$TRNA" | awk -F'\t' '!/^#/ && ($3=="exon" || $3=="tRNA") {
    match($9, /gene_id "[^"]+"/); g=(RSTART?substr($9, RSTART+9, RLENGTH-10):"tRNA")
    printf "%s\t%d\t%d\t%s\t0\t%s\n", $1, $4-1, $5, g, $7 }'
} | sort -k1,1 -k2,2n | "$BEDTOOLS" merge -s -c 4,6 -o distinct,distinct -i - \
  > "$BED/human_v49_deplete.bed"
awk -v pat="$POSTFILT" -F'\t' '$3=="exon" && $9 ~ ("gene_type \"("pat")\"") {
    match($9, /gene_id "[^"]+"/); g=substr($9, RSTART+9, RLENGTH-10)
    printf "%s\t%d\t%d\t%s\t0\t%s\n", $1, $4-1, $5, g, $7 }' "$GTF" \
  | sort -k1,1 -k2,2n | "$BEDTOOLS" merge -s -c 4,6 -o distinct,distinct -i - \
  > "$BED/human_v49_postfilter.bed"
echo "deplete intervals:    $(wc -l < "$BED/human_v49_deplete.bed")"
echo "postfilter intervals: $(wc -l < "$BED/human_v49_postfilter.bed")"
[[ -s "$BED/human_v49_deplete.bed" ]] || { echo "ABORT: empty deplete BED" >&2; exit 1; }

# ---- 2. contaminant FASTA ----
# -s keeps strand. Sequence identity is what bowtie2 needs; exon-level intervals are fine because a
# footprint sitting inside an rRNA exon aligns to it regardless of splicing.
"$BEDTOOLS" getfasta -fi "$GEN" -bed "$BED/human_v49_deplete.bed" -s -fo "$OUT.genome.fa"
# rDNA goes in as whole sequences, not genome intervals -- it has no reliable primary-assembly locus.
cat "$OUT.genome.fa" $RDNA > "$OUT.fa"
rm -f "$OUT.genome.fa"
NSEQ=$(grep -c '^>' "$OUT.fa"); BP=$(grep -v '^>' "$OUT.fa" | tr -d '\n' | wc -c)
echo "contaminant FASTA: $NSEQ sequences, $BP bp"
[[ "$NSEQ" -gt 100 ]] || { echo "ABORT: only $NSEQ contaminant sequences -- GTF parse likely wrong" >&2; exit 1; }

# ---- 3. bowtie2 index ----
"$BT/bowtie2-build" --threads "$T" -q "$OUT.fa" "$OUT"
ls -la "$OUT"*.bt2 | awk '{printf "  %10s  %s\n",$5,$9}'
cat > "$OUT.README" <<TXT
Contaminant depletion index for human Ribo-seq (two-stage alignment, stage 1).
Built $(date -Iseconds) from:
  GTF    $GTF
  genome $GEN
Depleted gene_types: $DEPLETE  + cytoplasmic tRNA + 45S rDNA repeat (U13369) + 5S
  -> $NSEQ merged exon intervals, $BP bp
Post-alignment filter BED (NOT in this index): $POSTFILT
  -> $NEW/data/ncrna_intervals/human_v49_postfilter.bed
Usage: bowtie2 --very-fast-local -x $OUT -U reads.fq --un-gz unaligned.fq.gz
TXT
echo "done=$(date -Iseconds) -> $OUT"
