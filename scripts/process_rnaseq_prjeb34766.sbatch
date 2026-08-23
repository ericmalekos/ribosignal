#!/usr/bin/env bash
#SBATCH --job-name=rna_34766
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rna_34766_%j.out
set -euo pipefail
# PRJEB34766 mouse-liver RNA -> the two model inputs, for the RNA-quality factorial (#68).
#   (1) salmon decoy-aware vM38 isoform TPM  -> the alternative UNIVERSE
#   (2) STAR --quantMode TranscriptomeSAM -> rnaseq_coverage.py -> the alternative COVERAGE track
#
# WHY THIS ARM. Every existing mouse-liver RNA arm is SINGLE-end (Janich 50-51 nt, GSE243134
# 75-80 nt), so RNA quality has never been varied independently of anything else. This is PAIRED,
# library_selection=PolyA, 7 CreNeg control runs / 223 M pairs, and -- confirmed from the ENA sample
# `strain` field, not inferred from a filename -- C57BL/6 male adult liver, matching the Ribo-seq.
# Reverba flox/AlbCreERT2 circadian liver (PNAS 2020, doi 10.1073/pnas.2005330117). CrePos/KO/WT runs
# are excluded: Reverba deletion perturbs the circadian liver transcriptome.
#
# The 7 runs are pooled into ONE quant and ONE alignment: they are matched controls of one tissue and
# the model consumes a single per-nt coverage track per arm. NOTE this is a circadian study, so the
# controls may span timepoints; pooling averages over a rhythmic transcriptome, which is acceptable
# for a coverage track but is stated rather than hidden.
#
# Same recipe as every other coverage arm (align_rnaseq_sample.sh): STAR --outFilterMultimapNmax 1,
# transcriptome axis, rnaseq_coverage.py auto-strand. Deviating would confound the RNA comparison
# with a pipeline difference.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grcm39_vM38
SALIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_vM38
FQ=$NEW/data/external/PRJEB34766_liver_rna
OUTS=$NEW/data/tpm/salmon_quant_prjeb34766_liver
OUTC=$NEW/data/heldout_rnaseq_coverage/mouse_prjeb34766_liver
T=${SLURM_CPUS_PER_TASK:-12}
SCR=/data/tmp/emalekos/rna34766_${SLURM_JOB_ID}; mkdir -p "$SCR" "$OUTS" "$OUTC"
trap 'rm -rf "$SCR"' EXIT
echo "host=$(hostname) start=$(date -Iseconds)"
[[ -d "$STARIDX" ]] || { echo "no STAR index $STARIDX" >&2; exit 1; }

mapfile -t RUNS < <(awk -F'\t' 'NR>1{print $1}' "$FQ/manifest.tsv")
R1=(); R2=()
for s in "${RUNS[@]}"; do
  a=$FQ/${s}_1.fastq.gz; b=$FQ/${s}_2.fastq.gz
  [[ -s "$a" && -s "$b" ]] || { echo "missing pair for $s" >&2; exit 1; }
  R1+=("$a"); R2+=("$b")
done
echo "runs: ${#RUNS[@]} paired"

# --- 1. salmon: per-run (so the factorial can mean-merge like the other arms) ---
for i in "${!RUNS[@]}"; do
  s=${RUNS[$i]}
  [[ -s "$OUTS/$s/quant.sf" ]] && { echo "skip salmon $s"; continue; }
  "$ENV/bin/salmon" quant -i "$SALIDX" -l A -1 "${R1[$i]}" -2 "${R2[$i]}" \
    --validateMappings --gcBias --seqBias -p "$T" -o "$OUTS/$s" > "$OUTS/$s.log" 2>&1
  echo "  salmon $s: $(awk 'NR>1 && $4>=1' "$OUTS/$s/quant.sf" | wc -l) tx TPM>=1"
done

# --- gate: library chemistry, before anything expensive downstream (feedback_rnaseq_must_be_polya) ---
# polya_check takes no --tx2biotype: the biotype approach was retired because the decoy-aware index
# holds pc+lncRNA only, so no library can ever show rRNA TPM. Distinguish a FLAGGED library (exit 1)
# from a BROKEN invocation (exit 2+), or a crashed gate reads as a passed gate.
set +e; "$PY" "$NEW/proteogenomics/scripts/polya_check.py" --quant "$OUTS/${RUNS[0]}/quant.sf"; gate=$?; set -e
case $gate in
  0) echo "  polya gate: PASS" ;;
  1) echo "  polya gate: FLAGGED -- inspect before use" ;;
  *) echo "  polya gate: ERRORED (exit $gate) -- gate did not run" >&2; exit 1 ;;
esac

# --- 2. STAR over ALL runs pooled -> one transcriptome BAM -> per-nt coverage ---
IFS=,; RR1="${R1[*]}"; RR2="${R2[*]}"; unset IFS
"$ENV/bin/STAR" --genomeDir "$STARIDX" \
  --readFilesIn "$RR1" "$RR2" --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 1 --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCR/PRJEB34766."
TBAM=$SCR/PRJEB34766.Aligned.toTranscriptome.out.bam
[[ -s "$TBAM" ]] || { echo "STAR produced no transcriptome BAM" >&2; exit 1; }
cp "$SCR/PRJEB34766.Log.final.out" "$OUTC/PRJEB34766.Log.final.out"
grep -E "Uniquely mapped reads %|Number of input reads" "$OUTC/PRJEB34766.Log.final.out"

"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$TBAM" \
  "$OUTC/PRJEB34766_coverage.hd5" PRJEB34766
echo "done=$(date -Iseconds) -> $OUTC/PRJEB34766_coverage.hd5"
