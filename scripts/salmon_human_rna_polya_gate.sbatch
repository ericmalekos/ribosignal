#!/usr/bin/env bash
#SBATCH --job-name=hs_salmon
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/hs_salmon_%A_%a.out
set -euo pipefail
# Salmon quant of the human RNA arms, to run the poly(A) gate BEFORE any universe is defined.
#
# WHY THIS GATE. A universe built from Ribo-Zero total RNA is corrupted: the library keeps tRNA and
# 7SL, both of which are annotated as ordinary transcripts, so a biotype filter cannot remove them and
# they soak up the TPM mass. That is what crushed the GSE243134 liver universe to 3,321 transcripts.
# SRA `library_selection` does NOT settle it -- it says `cDNA` for poly(A) and total-RNA libraries
# alike. The composition of the quant does settle it.
#
# GSE304796 (CAR-T) is the live risk: its RNA is rRNA-depleted TOTAL RNA per the submitters, not
# poly(A)-selected. GSE208041 (THP-1) is included as the comparison arm -- a dataset expected to pass
# makes a failure interpretable rather than ambiguous.
#
# Decoy-aware index, always (feedback_salmon_decoy_aware): transcriptome-only quant misattributes
# intronic and repeat-element reads to lncRNAs and biases exactly the TPMs this gate reads.
#
# NOTE this gate judges the UNIVERSE, not the pack. A TPM-contaminated library can still be a clean
# model input -- salmon EM + length normalisation and STAR mm1 + ncRNA filtering disagree strongly on
# multi-copy loci (feedback_salmon_tpm_vs_pack_coverage). A FAIL here blocks universe definition from
# that RNA, not use of the RNA as coverage.
#
#   sbatch --array=0-6 scripts/salmon_human_rna_polya_gate.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
SALMON=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/salmon
IDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_v49
OUT=$NEW/data/tpm/salmon_human_rna
T=${SLURM_CPUS_PER_TASK:-12}

RUNS=(SRR21228005 SRR21228006 SRR21228015 SRR21228016 SRR34901749 SRR34901750 SRR34901751)
DIRS=(GSE208041_thp1 GSE208041_thp1 GSE208041_thp1 GSE208041_thp1 GSE304796_cart GSE304796_cart GSE304796_cart)
SRR=${RUNS[$SLURM_ARRAY_TASK_ID]}
SUB=${DIRS[$SLURM_ARRAY_TASK_ID]}
D=$NEW/data/external/$SUB
O=$OUT/$SRR

[[ -x "$SALMON" ]] || { echo "no salmon at $SALMON" >&2; exit 1; }
[[ -d "$IDX" ]] || { echo "no index at $IDX" >&2; exit 1; }
[[ -s "$O/quant.sf" ]] && { echo "SKIP $SRR (quant.sf exists)"; exit 0; }
mkdir -p "$O"
echo "run=$SRR dir=$SUB host=$(hostname) start=$(date -Iseconds)"

# -l A lets salmon infer the library type; paired and single arms both appear in this array.
if [[ -s "$D/${SRR}_1.fastq.gz" && -s "$D/${SRR}_2.fastq.gz" ]]; then
  echo "  paired-end"
  "$SALMON" quant -i "$IDX" -l A -1 "$D/${SRR}_1.fastq.gz" -2 "$D/${SRR}_2.fastq.gz" \
    -p "$T" --validateMappings -o "$O"
elif [[ -s "$D/${SRR}.fastq.gz" ]]; then
  echo "  single-end"
  "$SALMON" quant -i "$IDX" -l A -r "$D/${SRR}.fastq.gz" \
    -p "$T" --validateMappings -o "$O"
else
  echo "no fastq for $SRR in $D" >&2; exit 1
fi

[[ -s "$O/quant.sf" ]] || { echo "ABORT: salmon produced no quant.sf" >&2; exit 1; }
echo "  mapping rate: $(grep -o '"percent_mapped": *[0-9.]*' "$O/aux_info/meta_info.json" 2>/dev/null | head -1)"
echo "done=$(date -Iseconds) -> $O"
