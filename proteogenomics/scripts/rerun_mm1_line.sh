#!/usr/bin/env bash
# Strict-uniform mm1 re-run of one Fig 2 immunopeptidome line, fully SLURM-native (afterok chain, no
# head-node babysitting): align RNA with --outFilterMultimapNmax 1 (SKIP_SALMON=1, universe unchanged)
# -> rebuild the model-input pack from the new mm1 coverage -> GPU predict -> enumerate/score ORFs ->
# model+null DBs -> HLA search -> MS2Rescore -> class-FDR compare. Result: db_comparison_rescored_mm1.md.
#
# Usage: rerun_mm1_line.sh <LINE> <HELDOUT> <SRR>
#   e.g. rerun_mm1_line.sh DoHH2 human_dohh2 SRR12285182
set -euo pipefail
LINE=$1; HELDOUT=$2; SRR=$3
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
S=$NEW/proteogenomics/scripts
P=$NEW/proteogenomics/data/${LINE}_pilot
R1=$P/rnaseq/${SRR}_1.fastq.gz; R2=$P/rnaseq/${SRR}_2.fastq.gz
[ -f "$R1" ] && [ -f "$R2" ] || { echo "missing FASTQ: $R1 / $R2" >&2; exit 1; }
[ -f "$P/${LINE}_universe.fa" ] || { echo "missing universe $P/${LINE}_universe.fa (run build_line_universe first)" >&2; exit 1; }

AJ=$(sbatch --parsable --export=ALL,LINE=$LINE,R1=$R1,R2=$R2,SKIP_SALMON=1 $S/align_line.sbatch)
PKJ=$(sbatch --parsable --dependency=afterok:$AJ --export=ALL,DS=$HELDOUT $NEW/scripts/heldout/build_heldout_pack.sbatch)
DJ=$(sbatch --parsable --dependency=afterok:$PKJ --export=ALL,LINE=$LINE,HELDOUT=$HELDOUT $S/dump_line.sbatch)
XJ=$(sbatch --parsable --dependency=afterok:$DJ --export=ALL,LINE=$LINE $S/post_predict_line.sbatch)
echo "$LINE ($HELDOUT, $SRR) mm1 chain:  align=$AJ -> pack=$PKJ -> predict=$DJ -> postpred=$XJ"
