#!/usr/bin/env bash
# QC ONLY: run RiboCode metaplots (P-site offset + 3-nt periodicity) on the cleaned HBL-1 DMSO transcriptome
# BAMs. Does NOT call ORFs -- ORF calling is held pending QC review.
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
RC=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
ANNOT=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/human_ribocode_annot_primary
RB=$NEW/proteogenomics/data/HBL1_pilot/riboseq
OUT=$NEW/proteogenomics/data/HBL1_pilot/ribocode; mkdir -p "$OUT"; cd "$OUT"
D1=$RB/SRR12285192.Aligned.toTranscriptome.out.bam
D2=$RB/SRR12285191.Aligned.toTranscriptome.out.bam
for b in "$D1" "$D2"; do [[ -f "$b" && -f "$b.bai" ]] || { echo "missing/unindexed $b" >&2; exit 1; }; done
printf "%s\n%s\n" "$D1" "$D2" > "$OUT/hbl1_dmso_bams.txt"
"$RC/metaplots" -a "$ANNOT" -i "$OUT/hbl1_dmso_bams.txt" -o HBL1_DMSO_QC > "$OUT/metaplots_qc.log" 2>&1
echo "metaplots QC done -> $OUT/HBL1_DMSO_QC_pre_config.txt"
