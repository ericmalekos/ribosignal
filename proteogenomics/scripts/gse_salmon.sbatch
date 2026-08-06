#!/usr/bin/env bash
#SBATCH --job-name=gse_salmon
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/gse_salmon_%x_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/gse_salmon_%x_%j.err
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
SALIDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_vM38
DEST=$NEW/proteogenomics/data/gse120762/raw
case "${COND:?set COND=nt|lps}" in
  nt)  SRRS="SRR7956038 SRR7956039 SRR7956040" ;;   # BMDM WT/NT rep1-3
  lps) SRRS="SRR7956041 SRR7956042" ;;              # BMDM LPS rep1-2
  *) echo "bad COND"; exit 1 ;;
esac
R1=""; R2=""; for s in $SRRS; do R1+="$DEST/${s}_1.fastq.gz "; R2+="$DEST/${s}_2.fastq.gz "; done
OUT=$NEW/proteogenomics/data/gse120762/salmon_$COND
echo "COND=$COND reps=$SRRS start=$(date -Iseconds)"
"$ENV/bin/salmon" quant -i "$SALIDX" -l A -1 $R1 -2 $R2 --gcBias --seqBias -p 8 -o "$OUT"
echo "salmon $COND done: tx>1TPM $(awk 'NR>1 && $4>1' "$OUT/quant.sf" | wc -l)"
