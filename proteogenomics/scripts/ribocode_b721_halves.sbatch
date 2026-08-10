#!/usr/bin/env bash
# B721 split-half ceiling: RiboCode called independently on two disjoint halves of the SAME
# Ribo-seq, so the two call sets differ only by sampling noise and depth.
#
# WHY. The drop-in F1 (predicted vs measured) is uninterpretable without this. An F1 of 0.67 could
# mean the model is mediocre, or it could mean RiboCode itself only reproduces its own calls at 0.7
# when handed an independent measurement of the same cells. Every prior calibration in this project
# (O2) was anchored the same way: compare the model to the reproducibility ceiling of the assay, not
# to a perfect oracle. Same script, same annotation, same pval/min_len as the drop-in scoring.
#
# The split is by HLA allele so each half is a set of complete biological libraries, never two
# technical re-sequencings of one library split across halves (that would understate the noise):
#   half A: A0101 x2, A3303 x1   (SRR10846525,26,27)
#   half B: B1501 x2, B4402 x2   (SRR10846528,29,30,31)
#
#   sbatch --array=0-1 proteogenomics/scripts/ribocode_b721_halves.sbatch
#SBATCH --job-name=rc_b721_half
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=110G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/rc_b721_half_%A_%a.out
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
RC=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
ANNOT=$ECH/data/ribocode_annot
RB=$NEW/proteogenomics/data/B721_pilot/riboseq

case "$SLURM_ARRAY_TASK_ID" in
  0) H=A; SRRS=(SRR10846525 SRR10846526 SRR10846527) ;;
  1) H=B; SRRS=(SRR10846528 SRR10846529 SRR10846530 SRR10846531) ;;
  *) echo "bad array index" >&2; exit 1 ;;
esac
OUT=$NEW/proteogenomics/data/B721_pilot/ribocode_half$H; mkdir -p "$OUT"

BAMLIST=$OUT/bams.txt; : > "$BAMLIST"
for s in "${SRRS[@]}"; do
  b=$RB/$s.Aligned.toTranscriptome.out.bam
  [[ -f "$b" && -f "$b.bai" ]] || { echo "missing/unindexed $b" >&2; exit 1; }
  echo "$b" >> "$BAMLIST"
done
echo "[$(date +%T)] host=$(hostname) half=$H  $(wc -l < "$BAMLIST") BAMs"

cd "$OUT"
"$RC/metaplots" -a "$ANNOT" -i "$BAMLIST" -o half$H > "$OUT/metaplots.log" 2>&1
CONFIG=$OUT/half${H}_pre_config.txt
[[ -f "$CONFIG" ]] || { echo "no pre_config" >&2; tail -20 "$OUT/metaplots.log"; exit 2; }
cat "$CONFIG"
"$RC/RiboCode" -a "$ANNOT" -c "$CONFIG" -l no -g -o half$H > "$OUT/ribocode.log" 2>&1
echo "[$(date +%T)] DONE half=$H: $(($(wc -l < "$OUT"/half${H}_collapsed.txt) - 1)) calls"
