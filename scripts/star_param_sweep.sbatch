#!/usr/bin/env bash
#SBATCH --job-name=star_sweep
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/star_sweep_%A_%a.out
set -euo pipefail
# STAR parameter sweep for Ribo-seq, one-at-a-time around a fixed base.
#
# BASE (held constant in every arm, per the standing rule): --alignEndsType EndToEnd,
# --outFilterMultimapNmax 1, --quantMode TranscriptomeSAM. Arm 0 leaves the four swept parameters at
# STAR's defaults, so every other arm reads directly as "what this one constraint costs".
#
# One-at-a-time, not a grid: 4 params x 3-4 levels is 13 arms instead of ~150, and the question is
# the marginal cost of each constraint, not their interaction. Interactions can be checked later on
# whichever constraints turn out to matter.
#
# THE READOUT IS PERIODICITY, NOT RETENTION. Losing reads is not itself bad -- the point of a
# mismatch filter is to lose the wrong ones. A constraint is worth paying for if it RAISES the
# in-frame fraction (f0%) measured against annotated CDS starts. A constraint that drops reads
# without raising f0 is pure loss.
#
# Single run (SRR15513213, HUVEC_1, ~100 M reads) for speed; this measures alignment behaviour, not
# biology, so one deep library is sufficient and all arms see identical input.
#
#   sbatch --array=0-12 scripts/star_param_sweep.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
CDS=$NEW/data/human_ribocode_annot_primary/transcripts_cds.txt
FQ=$NEW/data/external/chothani_ribo_refetch/SRR15513213.fastq.gz
OUT=$NEW/results/star_param_sweep
mkdir -p "$OUT"
T=${SLURM_CPUS_PER_TASK:-8}

#            label                      extra STAR flags
ARMS=(
  "base_defaults|"
  "mismatchNoverL_0.20|--outFilterMismatchNoverLmax 0.20"
  "mismatchNoverL_0.10|--outFilterMismatchNoverLmax 0.10"
  "mismatchNoverL_0.04|--outFilterMismatchNoverLmax 0.04"
  "matchNminOverL_0.75|--outFilterMatchNminOverLread 0.75"
  "matchNminOverL_0.85|--outFilterMatchNminOverLread 0.85"
  "matchNminOverL_0.95|--outFilterMatchNminOverLread 0.95"
  "scoreMinOverL_0.75|--outFilterScoreMinOverLread 0.75"
  "scoreMinOverL_0.85|--outFilterScoreMinOverLread 0.85"
  "scoreMinOverL_0.95|--outFilterScoreMinOverLread 0.95"
  "seedStartLmax_30|--seedSearchStartLmax 30"
  "seedStartLmax_20|--seedSearchStartLmax 20"
  "seedStartLmax_15|--seedSearchStartLmax 15"
)
IFS='|' read -r LAB EXTRA <<< "${ARMS[$SLURM_ARRAY_TASK_ID]}"
ROW=$OUT/${LAB}.tsv
[[ -s "$ROW" ]] && { echo "SKIP $LAB (done)"; exit 0; }
SCR=/data/tmp/emalekos/sweep_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT
echo "arm=$LAB extra='${EXTRA}' host=$(hostname) start=$(date -Iseconds)"

# Same trim for every arm: N-padded reads, verified (no ligated adapter present).
"$ENV/bin/cutadapt" --trim-n -m 20 -M 40 -j "$T" -o "$SCR/t.fq.gz" "$FQ" > "$SCR/cutadapt.log" 2>&1

# shellcheck disable=SC2086  # EXTRA is intentionally word-split into flags
"$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$SCR/t.fq.gz" --readFilesCommand zcat \
  --runThreadN "$T" --alignEndsType EndToEnd --outFilterMultimapNmax 1 \
  --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM \
  --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCR/s." $EXTRA
L=$SCR/s.Log.final.out
IN=$(awk -F'\t' '/Number of input reads/{gsub(/ /,"",$2); print $2}' "$L")
UNIQ=$(awk -F'\t' '/Uniquely mapped reads number/{gsub(/ /,"",$2); print $2}' "$L")
UPCT=$(awk -F'\t' '/Uniquely mapped reads %/{gsub(/[ %]/,"",$2); print $2}' "$L")
MMR=$(awk -F'\t' '/Mismatch rate per base/{gsub(/[ %]/,"",$2); print $2}' "$L")
SHORT=$(awk -F'\t' '/too short/{gsub(/[ %]/,"",$2); print $2}' "$L" | head -1)

# Periodicity: in-frame fraction of 5' ends relative to annotated CDS starts. This is the metric
# that decides whether a constraint bought anything.
"$PY" - <<PYEOF > "$SCR/frame.txt"
import subprocess, collections
cds={}
for line in open("$CDS"):
    f=line.rstrip("\n").split("\t")
    if len(f)>=3:
        try: cds[f[0]]=(int(f[1]),int(f[2]))
        except ValueError: pass
S="$ENV/bin/samtools"
p=subprocess.Popen([S,"view","$SCR/s.Aligned.toTranscriptome.out.bam"],stdout=subprocess.PIPE,text=True,bufsize=1<<20)
fr=[0,0,0]; n=0
for line in p.stdout:
    f=line.split("\t",10)
    if len(f)<10: continue
    c=cds.get(f[2])
    if c is None: continue
    L=len(f[9])
    if not (25<=L<=34): continue
    fr[(int(f[3])-1-c[0])%3]+=1; n+=1
    if n>=3_000_000: break
p.terminate()
tot=sum(fr) or 1
print(f"{100*max(fr)/tot:.2f}\t{n}")
PYEOF
read -r F0 NFR < "$SCR/frame.txt"
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$LAB" "${EXTRA:-<defaults>}" "$IN" "$UNIQ" "$UPCT" "$MMR" "$F0" "$NFR" > "$ROW"
echo "  in=$IN uniq=$UNIQ (${UPCT}%) mismatch_rate=${MMR}% max_frame=${F0}% over $NFR reads"
echo "done=$(date -Iseconds)"
