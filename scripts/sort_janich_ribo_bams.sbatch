#!/usr/bin/env bash
#SBATCH --job-name=janich_sort
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich_sort_%A_%a.out
set -euo pipefail
# Coordinate-sort the re-aligned Janich Ribo transcriptome BAMs.
#
# `STAR --quantMode TranscriptomeSAM` ALWAYS emits a read-ordered BAM with no `@HD SO:` tag, whatever
# --outSAMtype says. `samtools index` then fails, and prepare_from_bams calls index. The archived
# GSE243134 BAMs carry `SO:coordinate` because a separate sorting pass was run over them earlier;
# align_gse243134.sbatch on its own does not produce sorted output, so copying that recipe verbatim
# reproduced the gap.
#
# Sort -> VERIFY -> replace -> index, in that order. The verify step checks the header tag, quickcheck,
# AND that the record count survived: a truncated sort yields a valid-looking BAM with fewer reads,
# and replacing the original with that would be silent data loss.
#
#   sbatch --array=0-4 scripts/sort_janich_ribo_bams.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
S=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/samtools
B=$NEW/data/heldout_bam/mouse_janich_liver_ribo
T=${SLURM_CPUS_PER_TASK:-12}
mapfile -t BAMS < <(ls "$B"/*.bam | sort)
f=${BAMS[$SLURM_ARRAY_TASK_ID]}
n=$(basename "$f")
SCR=/data/tmp/emalekos/jsort_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$SCR"; trap 'rm -rf "$SCR"' EXIT
echo "file=$n host=$(hostname) start=$(date -Iseconds)"

# `head -1 | grep -q` under `set -o pipefail` gives a SIGPIPE false negative; use sed -n 1p.
if [[ "$("$S" view -H "$f" 2>/dev/null | sed -n 1p)" == *SO:coordinate* ]]; then
  echo "already coordinate-sorted; nothing to do"; exit 0
fi

before=$("$S" view -@ "$T" -c "$f")
"$S" sort -@ "$T" -m 3G -o "$SCR/sorted.bam" "$f"
after=$("$S" view -@ "$T" -c "$SCR/sorted.bam")
hdr=$("$S" view -H "$SCR/sorted.bam" 2>/dev/null | sed -n 1p)
echo "  records before=$before after=$after"
[[ "$before" -eq "$after" ]] || { echo "ABORT: record count changed ($before -> $after)" >&2; exit 1; }
[[ "$hdr" == *SO:coordinate* ]] || { echo "ABORT: sorted BAM lacks SO:coordinate" >&2; exit 1; }
"$S" quickcheck "$SCR/sorted.bam" || { echo "ABORT: quickcheck failed" >&2; exit 1; }

mv "$SCR/sorted.bam" "$f"
"$S" index -@ "$T" "$f"
echo "done=$(date -Iseconds): $n sorted + indexed"
