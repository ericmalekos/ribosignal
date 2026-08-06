#!/usr/bin/env bash
#SBATCH --job-name=hbl1_price
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=80G
#SBATCH --time=08:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/hbl1_price_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/hbl1_price_%j.err
set -euo pipefail
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq; export PATH="$ENV/bin:$PATH"
GEDI=/private/groups/carpenterlab/emalekos/conda_envs/gedi/bin
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
RB=$NEW/proteogenomics/data/HBL1_pilot/riboseq
OUT=$NEW/proteogenomics/data/HBL1_pilot/price; mkdir -p "$OUT"
D1=$RB/SRR12285192.genome.uniq.bam
D2=$RB/SRR12285191.genome.uniq.bam
for b in "$D1" "$D2"; do [[ -f "$b" && -f "$b.bai" ]] || { echo "missing $b" >&2; exit 1; }; done

echo "[$(date +%T)] pool D1+D2 unique genome BAMs"
POOL=$RB/HBL1_DMSO.genome.uniq.bam
samtools merge -f -@12 "$POOL" "$D1" "$D2"
samtools index -@12 "$POOL"
echo "  pooled unique genome reads: $(samtools view -c -@12 $POOL)"

echo "[$(date +%T)] PRICE (-genomic grch38_v49 registered name)"
$GEDI/gedi -e Price -reads "$POOL" -genomic grch38_v49 -nthreads 12 -progress \
    -prefix "$OUT/HBL1_DMSO" > "$OUT/price.log" 2>&1
echo "  PRICE exit=$?"
grep -iE "Remaining after multiple testing|ORFs" "$OUT/price.log" | tail -5

PT="$OUT/HBL1_DMSO.orfs.tsv"
[[ -s "$PT" ]] || { echo "no PRICE orfs.tsv!" >&2; tail -20 "$OUT/price.log"; exit 3; }
echo "[$(date +%T)] PRICE ORFs: $(($(wc -l < "$PT")-1))"
echo "-- Type distribution --"
awk -F'\t' 'NR>1{print $6}' "$PT" | sort | uniq -c | sort -rn
echo "[$(date +%T)] DONE -> $PT"