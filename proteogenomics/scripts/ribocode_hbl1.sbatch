#!/usr/bin/env bash
#SBATCH --job-name=hbl1_ribocode
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/hbl1_ribocode_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/hbl1_ribocode_%j.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
RC=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
ANNOT=$ECH/data/ribocode_annot                      # existing v49 prepared-transcripts (same as training)
RB=$NEW/proteogenomics/data/HBL1_pilot/riboseq
OUT=$NEW/proteogenomics/data/HBL1_pilot/ribocode; mkdir -p "$OUT"

D1=$RB/SRR12285192.Aligned.toTranscriptome.out.bam
D2=$RB/SRR12285191.Aligned.toTranscriptome.out.bam
for b in "$D1" "$D2"; do [[ -f "$b" && -f "$b.bai" ]] || { echo "missing/unindexed $b" >&2; exit 1; }; done
BAMLIST=$OUT/hbl1_dmso_bams.txt; printf "%s\n%s\n" "$D1" "$D2" > "$BAMLIST"
cd "$OUT"
echo "[$(date +%T)] metaplots (pool D1+D2) using annot=$ANNOT"
"$RC/metaplots" -a "$ANNOT" -i "$BAMLIST" -o HBL1_DMSO > "$OUT/metaplots.log" 2>&1
CONFIG=$OUT/HBL1_DMSO_pre_config.txt
[[ -f "$CONFIG" ]] || { echo "no pre_config produced" >&2; tail -20 "$OUT/metaplots.log"; exit 2; }
echo "[$(date +%T)] pre_config rows: $(wc -l < "$CONFIG")"; cat "$CONFIG"
echo "[$(date +%T)] RiboCode joint ORF call"
"$RC/RiboCode" -a "$ANNOT" -c "$CONFIG" -l no -g -o HBL1_DMSO > "$OUT/ribocode.log" 2>&1
echo "[$(date +%T)] DONE. outputs:"; ls -lh "$OUT" | grep -iE "HBL1_DMSO" | head
echo "-- collapsed ORF count --"; wc -l "$OUT"/HBL1_DMSO*collapsed.txt 2>/dev/null || true