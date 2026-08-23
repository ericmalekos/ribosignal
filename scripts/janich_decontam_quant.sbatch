#!/usr/bin/env bash
#SBATCH --job-name=janich_decon
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich_decon_%A_%a.out
set -euo pipefail
# Remove the structural-RNA reads that dominate Janich liver TPM, then re-quantify.
#
# THE PROBLEM. All 7 Janich liver RNA replicates put 52-70% of salmon TPM into just two transcripts:
# Gm59647 (ENSMUST00000292750.1, a 7SL/SRP duplicate, 302 nt) and n-TKctt14 (ENSMUST00000359371.1,
# tRNA-derived, 100 nt). Both are typed `lncRNA`, so no biotype filter removes them. TPM is
# compositional, so those two absorb over half the pool and depress every other transcript's TPM by
# ~2x. Measured consequence: the TPM >= 1 universe is 11,532 tx where renormalisation without them
# gives 16,469 -- about 4,900 real transcripts per replicate suppressed below threshold.
#
# NOTE the per-nt coverage the MODEL consumes is already clean (these two carry 0.000% of the pack's
# coverage mass) because STAR runs with --outFilterMultimapNmax 1 plus the ncRNA-locus filter, which
# strips multi-copy structural RNA. Salmon has no such filter and uses length-normalised EM, which is
# why the two quantifications disagree by orders of magnitude here. This job fixes the SALMON side,
# i.e. universe selection; no trained model is affected.
#
# WHY FILTER READS RATHER THAN EDIT THE INDEX. CLAUDE.md keeps ONE canonical index per
# (tool, species, annotation) and reuses it across projects. Building a Janich-specific index minus
# two transcripts would make its TPMs non-comparable with Wang and GSE302188, which are the other two
# liver arms in the three-way ORF-call comparison. Removing the READS keeps the shared
# decoy-aware vM38 index untouched, and is also what the user asked for.
#
# Equivalence note: with the reads removed, TPM for the survivors reduces to the same arithmetic as
# renormalising the existing quant.sf, EXCEPT for reads that ambiguously map between a contaminant
# and a real transcript -- those vanish here but would be reassigned under renormalisation. The
# comparison against the renormalisation estimate is therefore a check on how much ambiguity exists.
#
#   sbatch --array=0-6 scripts/janich_decontam_quant.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
BT=/private/groups/carpenterlab/emalekos/conda_envs/crackling/bin       # bowtie2 lives here
SAL=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/salmon
SALIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_vM38
CONTAM=$NEW/data/external/janich_contaminants.fa
FQDIR=$NEW/data/external/janich_liver_rnaseq_trimmed                     # adapter-fixed (task 35)
OUTROOT=$NEW/data/tpm/salmon_quant_janich_liver_decontam
RUNS=(SRR1930216 SRR1930217 SRR1930218 SRR1930221 SRR1930222 SRR1930225 SRR1930226)
SRR=${RUNS[$SLURM_ARRAY_TASK_ID]}
T=${SLURM_CPUS_PER_TASK:-8}
mkdir -p "$OUTROOT"

# The bowtie2 index is two short sequences; build it once, guarded so concurrent array tasks do not
# race on the same path (one writer per path, project rule).
IDXDIR=$NEW/data/external/janich_contam_bt2
if [[ ! -s "$IDXDIR/contam.1.bt2" ]]; then
  mkdir -p "$IDXDIR"
  LOCK="$IDXDIR/.building"
  if mkdir "$LOCK" 2>/dev/null; then
    "$BT/bowtie2-build" --quiet "$CONTAM" "$IDXDIR/contam"
    rmdir "$LOCK"
  else
    while [[ -d "$LOCK" ]]; do sleep 5; done
  fi
fi
[[ -s "$IDXDIR/contam.1.bt2" ]] || { echo "no bowtie2 index at $IDXDIR" >&2; exit 1; }

SCR=/data/tmp/emalekos/janich_decon_${SRR}_${SLURM_JOB_ID}; mkdir -p "$SCR"
trap 'rm -rf "$SCR"' EXIT
echo "[$(date +%T)] $SRR host=$(hostname)"

# --very-sensitive-local: these are short structural RNAs and a footprint may only partially overlap;
# end-to-end would leave partial-overlap reads behind, which is the opposite of what decontamination
# wants. Reads that FAIL to align are the ones kept (--un-gz).
"$BT/bowtie2" -p "$T" --very-sensitive-local -x "$IDXDIR/contam" \
   -U "$FQDIR/$SRR.fastq.gz" --un-gz "$SCR/$SRR.clean.fastq.gz" -S /dev/null \
   2> "$OUTROOT/$SRR.bowtie2.log"
grep -E "reads; of these|aligned exactly 1 time|aligned >1 times|overall alignment rate" \
   "$OUTROOT/$SRR.bowtie2.log" | sed 's/^/  /'

"$SAL" quant -i "$SALIDX" -l A -r "$SCR/$SRR.clean.fastq.gz" \
   --validateMappings --gcBias --seqBias -p "$T" -o "$OUTROOT/$SRR" \
   > "$OUTROOT/$SRR.salmon.log" 2>&1
echo "[$(date +%T)] $SRR done; tx TPM>=1: $(awk 'NR>1 && $4>=1' "$OUTROOT/$SRR/quant.sf" | wc -l)"
