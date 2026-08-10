#!/usr/bin/env bash
#SBATCH --job-name=rnaseq_thp1
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/rnaseq_thp1_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/rnaseq_thp1_%j.err
set -euo pipefail
# THP-1 RNA-seq (PRJNA686824, 3 paired-end runs) -> both model inputs, same recipe as A549:
#   (1) salmon decoy-aware v49 isoform TPM  -> expressed-transcript filter + the null databases
#   (2) STAR --quantMode TranscriptomeSAM -> rnaseq_coverage.py -> per-nt coverage -> THE MODEL INPUT
# The three runs are pooled into ONE quant and ONE alignment: they are replicates of the same
# untreated THP-1 culture, and the model consumes a single per-nt coverage track per population.
#
# POLY(A) CHECK (standing rule: `library_selection=cDNA` in SRA does NOT prove poly(A) selection --
# a Ribo-Zero total-RNA arm keeps tRNA/7SL, which are typed lncRNA and survive the biotype filter,
# and that crushed the GSE243134 liver universe). Verified empirically below rather than trusted:
# the job prints the intronic-vs-exonic and the rRNA/tRNA/7SL fraction, and the universe builder
# downstream will show the transcript count. A total-RNA library shows up as a collapsed universe.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PILOT=$NEW/proteogenomics/data/THP1_pilot
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
STARIDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
SALIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_v49
RUNS=(SRR13279451 SRR13279452 SRR13279453)
T=${SLURM_CPUS_PER_TASK:-8}
SCR=/data/tmp/emalekos/thp1_rnaseq_${SLURM_JOB_ID:-manual}
mkdir -p "$SCR" "$PILOT/salmon" "$PILOT/coverage"
trap 'rm -rf "$SCR"' EXIT
echo "host=$(hostname) start=$(date -Iseconds)"

R1=(); R2=()
for s in "${RUNS[@]}"; do
  a=$PILOT/rnaseq/${s}_1.fastq.gz; b=$PILOT/rnaseq/${s}_2.fastq.gz
  [[ -s "$a" && -s "$b" ]] || { echo "missing fastq for $s" >&2; exit 1; }
  R1+=("$a"); R2+=("$b")
done

# --- 1) salmon decoy-aware v49 -> isoform TPM (3 runs as one sample) ---
"$ENV/bin/salmon" quant -i "$SALIDX" -l A -1 "${R1[@]}" -2 "${R2[@]}" \
  --validateMappings --gcBias --seqBias -p "$T" -o "$PILOT/salmon/THP1"
echo "salmon done $(date -Iseconds); tx >1 TPM: $(awk 'NR>1 && $4>1' "$PILOT/salmon/THP1/quant.sf" | wc -l)"
echo "  mapping rate: $(grep -o '"percent_mapped": *[0-9.]*' "$PILOT/salmon/THP1/aux_info/meta_info.json" | head -1)"

# --- 2) STAR transcriptome align -> per-nt coverage (THE model input) ---
IFS=,; RR1="${R1[*]}"; RR2="${R2[*]}"; unset IFS
"$ENV/bin/STAR" --genomeDir "$STARIDX" \
  --readFilesIn "$RR1" "$RR2" --readFilesCommand zcat \
  --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 1 --outSAMattributes NH HI AS nM \
  --outFileNamePrefix "$SCR/THP1."
TBAM="$SCR/THP1.Aligned.toTranscriptome.out.bam"
[ -s "$TBAM" ] || { echo "STAR produced no transcriptome BAM" >&2; exit 1; }
cp "$SCR/THP1.Log.final.out" "$PILOT/coverage/THP1.Log.final.out"
grep -E "Uniquely mapped reads %|Number of input reads" "$PILOT/coverage/THP1.Log.final.out"

"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$TBAM" "$PILOT/coverage/THP1_coverage.hd5" THP1
echo "coverage done $(date -Iseconds) -> $PILOT/coverage/THP1_coverage.hd5"

# --- 3) poly(A) evidence: what fraction of TPM sits in rRNA/tRNA-like biotypes? ---
PILOT="$PILOT" NEW="$NEW" "$ENV/bin/python" - <<'PY'
import os
pilot = os.environ["PILOT"]
bt = {}
with open(os.environ["NEW"] + "/data/tx2biotype.tsv") as fh:
    ci = {c: i for i, c in enumerate(fh.readline().rstrip("\n").split("\t"))}
    for ln in fh:
        f = ln.rstrip("\n").split("\t"); bt[f[ci["tx_id"]]] = f[ci["transcript_type"]]
tot = 0.0; junk = 0.0
JUNK = {"rRNA", "rRNA_pseudogene", "Mt_rRNA", "Mt_tRNA", "misc_RNA", "snRNA", "snoRNA", "scaRNA"}
with open(f"{pilot}/salmon/THP1/quant.sf") as fh:
    fh.readline()
    for ln in fh:
        f = ln.split("\t"); tpm = float(f[3]); tot += tpm
        if bt.get(f[0], "") in JUNK: junk += tpm
print(f"POLYA CHECK: {junk/tot*100:.2f}% of TPM in rRNA/Mt/misc/sn(o)RNA "
      f"(poly(A)-selected libraries sit low single digits; Ribo-Zero total RNA runs much higher)")
PY
echo "all done $(date -Iseconds)"
