#!/usr/bin/env bash
#SBATCH --job-name=hs_sort
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/hs_sort_%A_%a.out
set -euo pipefail
# Coordinate-sort + index the human *.toTranscriptome.bam files.
#
# `STAR --quantMode TranscriptomeSAM` ALWAYS emits a read-ordered BAM with no `@HD SO:` tag, whatever
# --outSAMtype says, and the processing scripts copy it out untouched (only the GENOME bam gets
# sorted, for umi_tools dedup). prepare_from_bams calls `samtools index`, which then fails. This is
# the same gap that cost two runs on the Janich rebuild; it is fixed here for all three human
# datasets at once rather than per-dataset.
#
# Sort -> VERIFY -> replace -> index. The verify step checks the header tag, quickcheck, AND that the
# record count survived: a truncated sort yields a valid-looking BAM with fewer reads, and replacing
# the original with that is silent data loss.
#
# Idempotent: an already-sorted BAM is skipped, so this can be re-run after the last CAR-T sample
# lands without redoing the rest.
#
#   sbatch --array=0-17%6 scripts/sort_human_tx_bams.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
S=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/samtools
T=${SLURM_CPUS_PER_TASK:-8}
mapfile -t BAMS < <(ls $NEW/data/processed/*_snodep/*.toTranscriptome.bam 2>/dev/null | sort)
[[ ${#BAMS[@]} -gt 0 ]] || { echo "no BAMs found" >&2; exit 1; }
i=${SLURM_ARRAY_TASK_ID:-0}
[[ "$i" -lt "${#BAMS[@]}" ]] || { echo "task $i beyond ${#BAMS[@]} BAMs; nothing to do"; exit 0; }
f=${BAMS[$i]}
n=$(basename "$f")
SCR=/data/tmp/emalekos/hssort_${SLURM_JOB_ID}_${i}
mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT
echo "file=$n host=$(hostname) start=$(date -Iseconds)"

# `head -1 | grep -q` under `set -o pipefail` gives a SIGPIPE false negative; use sed -n 1p.
if [[ "$("$S" view -H "$f" 2>/dev/null | sed -n 1p)" == *SO:coordinate* ]]; then
  [[ -s "$f.bai" ]] || "$S" index -@ "$T" "$f"
  echo "already coordinate-sorted; ensured index"; exit 0
fi

before=$("$S" view -@ "$T" -c "$f")
"$S" sort -@ "$T" -m 3G -o "$SCR/sorted.bam" "$f"
after=$("$S" view -@ "$T" -c "$SCR/sorted.bam")
echo "  records before=$before after=$after"
[[ "$before" -eq "$after" ]] || { echo "ABORT: record count changed ($before -> $after)" >&2; exit 1; }
"$S" view -H "$SCR/sorted.bam" | sed -n 1p | grep -q SO:coordinate \
  || { echo "ABORT: sorted BAM lacks SO:coordinate" >&2; exit 1; }
"$S" quickcheck "$SCR/sorted.bam" || { echo "ABORT: quickcheck failed" >&2; exit 1; }

mv "$SCR/sorted.bam" "$f"
"$S" index -@ "$T" "$f"
echo "done=$(date -Iseconds): $n sorted + indexed"
