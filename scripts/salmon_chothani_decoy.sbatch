#!/usr/bin/env bash
#SBATCH --job-name=salmon_dec
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --array=1-25%6
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/salmon_chothani/dec_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/salmon_chothani/dec_%A_%a.err
set -euo pipefail
# Restore per-tissue Chothani TPM by the ORIGINAL method: decoy-aware selective alignment FROM
# FASTQ. This supersedes data/tpm/chothani_alnmode/, which was an interim recovery from STAR
# transcriptome BAMs and therefore lost 18.6% of reads to --outFilterMultimapNmax 1 (see methods.md).
#
# Reads are 74-76 nt, comfortably above the index k=31. Index carries 194 decoys (primary-assembly
# chromosomes), per the project's standing "always quant against the decoy-aware index" rule --
# transcriptome-only quant misattributes intronic and repeat reads to lncRNAs.
#
# Runs ONLY on md5-verified FASTQ pairs; a run whose pair is absent or short is skipped loudly
# rather than silently quantified on a truncated file.
#
#   sbatch scripts/salmon_chothani_decoy.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
SALMON=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/salmon
IDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_v49
D=$NEW/data/external/chothani_rna
MAN=$D/manifest_25.tsv
row=$(awk 'NR>1' "$MAN" | sed -n "${SLURM_ARRAY_TASK_ID}p")
SRR=$(echo "$row" | cut -f1); TIS=$(echo "$row" | cut -f2)
BYTES=$(echo "$row" | cut -f5)
R1=$D/fastq/${SRR}_1.fastq.gz; R2=$D/fastq/${SRR}_2.fastq.gz
B1=$(echo "$BYTES" | cut -d';' -f1); B2=$(echo "$BYTES" | cut -d';' -f2)
for pair in "$R1:$B1" "$R2:$B2"; do
  f=${pair%%:*}; want=${pair##*:}
  [ -s "$f" ] || { echo "SKIP $SRR: missing $(basename "$f")" >&2; exit 0; }
  got=$(stat -c %s "$f")
  [ "$got" -ge "$want" ] || { echo "SKIP $SRR: $(basename "$f") is $got, expected $want (truncated)" >&2; exit 0; }
done
OUT=$NEW/data/tpm/chothani_decoy/${TIS}/${SRR}
mkdir -p "$OUT" "$NEW/logs/salmon_chothani"
echo "[$(date +%T)] $SRR ($TIS) decoy-aware -> $OUT"
"$SALMON" quant -i "$IDX" -l A -1 "$R1" -2 "$R2" -o "$OUT" -p "${SLURM_CPUS_PER_TASK:-12}"
[ -s "$OUT/quant.sf" ] || { echo "FATAL no quant.sf for $SRR" >&2; exit 3; }
pm=$(python3 -c "import json;print(json.load(open('$OUT/aux_info/meta_info.json')).get('percent_mapped'))" 2>/dev/null || echo "?")
echo "[$(date +%T)] $SRR done -- percent_mapped=$pm"
