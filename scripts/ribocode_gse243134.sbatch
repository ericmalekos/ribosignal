#!/usr/bin/env bash
#SBATCH --job-name=rc_gse243134
#SBATCH --partition=long
#SBATCH --cpus-per-task=8
#SBATCH --mem=140G
#SBATCH --time=48:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rc_gse243134_%j.out
set -euo pipefail
# RiboCode ORF calls on the GSE243134 mouse-liver Ribo-seq -- the MISSING third arm of the
# between-dataset concordance analysis.
#
# WHY THIS DID NOT EXIST. GSE243134 entered the project as a third liver arm, but its matched
# totalRNA was REJECTED on library chemistry (Ribo-Zero; salmon 17.7%, universe crushed to 3,321 tx),
# and the arm was shelved with it. That conflated two independent things: the RNA side is unusable,
# the Ribo-seq side is excellent (21 WT runs, 273 M unique RPF, 133.4 M P-sites in its pack -- more
# than Wang's 97.5 M or Janich's 88.5 M). **Observed-vs-observed ORF concordance needs no RNA-seq at
# all**, so the Ribo-seq can carry a third point in the O2 ceiling on its own.
#
# 21 RPF runs (manifest `strategy` = OTHER; the other 19 rows are the rejected RNA-Seq). Same
# annotation, caller and defaults as the Wang and Janich collapsed call sets it will be compared
# against, so the three are directly comparable: prepared vM38 transcripts, start ATG, min_AA 5,
# pval 0.05, longest_orf, fdr_bh.
#
# Depth footnote for the writeup: 133.4 M vs 97.5 M (Wang) vs 88.5 M (Janich) P-sites, i.e. within
# 1.5x, so the three-way F1s are not dominated by a depth mismatch. Depth still matters (see the
# depth-crossover result) and is reported alongside.
#
#   sbatch scripts/ribocode_gse243134.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
RC=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
SAMTOOLS=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/samtools
ANNOT="/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/mouse_ribocode_annot"
BAMDIR=$NEW/data/heldout_bam/mouse_gse243134_liver
MANIFEST=$NEW/data/external/gse243134_mouse_liver/manifest_wt.tsv
OUT=$NEW/data/heldout_psites/mouse_gse243134_liver; mkdir -p "$OUT"
[[ -d "$ANNOT" ]] || { echo "no mouse annot at $ANNOT" >&2; exit 1; }

# The archived transcriptome BAMs are coordinate-sorted + indexed IN PLACE by
# scripts/sort_archived_tx_bams.sbatch, which this job is chained after. STAR emits
# --quantMode TranscriptomeSAM output in READ order (no @HD SO: line), so metaplots' pysam fetch
# cannot use them raw -- the first attempt at this job died on exactly that.
T=${SLURM_CPUS_PER_TASK:-8}
BAMLIST=$OUT/rpf_bams.txt; : > "$BAMLIST"
n=0
while read -r srr; do
  b=$BAMDIR/$srr.Aligned.toTranscriptome.out.bam
  [[ -s "$b" ]] || { echo "missing $b" >&2; exit 1; }
  # sed -n 1p, NOT `head -1 | grep -q`: under `set -o pipefail` the early-exiting head/grep SIGPIPE
  # the upstream samtools and the test reads FALSE even when the tag is present. That exact bug made
  # this job abort on a BAM that had just been correctly sorted.
  # quickcheck as well as the header: a BAM truncated mid-write keeps a VALID header and would
  # be waved through. That is how job 36834132 died after 92 min of work (2026-08-21).
  "$SAMTOOLS" quickcheck "$b" 2>/dev/null \
    || { echo "$b is TRUNCATED (quickcheck failed)" >&2; exit 1; }
  [[ "$("$SAMTOOLS" view -H "$b" | sed -n 1p)" == *SO:coordinate* ]] \
    || { echo "$b is NOT coordinate-sorted; run sort_archived_tx_bams.sbatch first" >&2; exit 1; }
  [[ -s "$b.bai" ]] || "$SAMTOOLS" index -@ "$T" "$b"
  echo "$b" >> "$BAMLIST"; n=$((n+1))
done < <(awk -F'\t' 'NR>1 && $3!="RNA-Seq"{print $1}' "$MANIFEST")
echo "[$(date +%T)] host=$(hostname)  $n RPF BAMs pooled"
[[ $n -eq 21 ]] || { echo "expected 21 RPF BAMs, got $n" >&2; exit 1; }

cd "$OUT"
echo "[$(date +%T)] metaplots (P-site offsets per BAM)"
"$RC/metaplots" -a "$ANNOT" -i "$BAMLIST" -o mouse_gse243134_liver > "$OUT/metaplots.log" 2>&1
CONFIG=$OUT/mouse_gse243134_liver_pre_config.txt
[[ -f "$CONFIG" ]] || { echo "no pre_config produced" >&2; tail -20 "$OUT/metaplots.log"; exit 2; }
grep -c "^SRR" "$CONFIG" | xargs echo "  libraries with a P-site offset:"

echo "[$(date +%T)] RiboCode joint ORF call"
"$RC/RiboCode" -a "$ANNOT" -c "$CONFIG" -l no -g -o mouse_gse243134_liver > "$OUT/ribocode.log" 2>&1
echo "[$(date +%T)] DONE: $(($(wc -l < "$OUT"/mouse_gse243134_liver_collapsed.txt) - 1)) collapsed calls"
