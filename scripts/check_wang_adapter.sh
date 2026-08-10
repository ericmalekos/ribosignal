#!/usr/bin/env bash
set -euo pipefail
# Adapter / insert-length check for the Wang Ribo runs, on the FULL downloaded files.
#
# Run this BEFORE align_wang_ribo.sbatch and read the output. Two candidate 3' adapters are trimmed
# in-place on a large read sample and the resulting insert-length distribution is compared. The
# correct adapter is the one that yields a mode near 30 nt with most inserts in 26-34 nt -- that is
# ribosome-footprint biology, and it is the check that distinguished TruSeq from a spurious poly-A hit
# on the Janich runs (poly-A "fired" on 93% of reads but produced mode 94 nt and 0.0% in range).
#
#   bash scripts/check_wang_adapter.sh
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
D=$NEW/data/external/wang2021_ribo_refetch
N=${N:-400000}          # reads sampled (4x lines)

for FQ in "$D"/*.fastq.gz; do
  [[ -s "$FQ" ]] || continue
  SRR=$(basename "$FQ" .fastq.gz)
  echo "=== $SRR ==="
  echo -n "  raw read length mode: "
  zcat "$FQ" | head -$((N*4)) | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c | sort -rn | head -1
  for A in AGATCGGAAGAGC CTGTAGGCACCATCAAT; do
    hits=$(zcat "$FQ" | head -$((N*4)) | awk 'NR%4==2' | grep -c "$A" || true)
    tot=$((N))
    pct=$(awk -v h="$hits" -v t="$tot" 'BEGIN{printf "%.1f", 100*h/t}')
    stats=$(zcat "$FQ" | head -$((N*4)) \
      | "$ENV/bin/cutadapt" -a "$A" --minimum-length 1 -o - - 2>/dev/null \
      | awk 'NR%4==2{L=length($0); n++; c[L]++; if(L>=26&&L<=34) inr++}
             END{ best=0; bl=0; for(l in c) if(c[l]>best){best=c[l]; bl=l}
                  printf "mode %s nt, 26-34 = %.1f%%", bl, 100*inr/n }')
    printf "  %-18s present in %5s%% of reads -> %s\n" "$A" "$pct" "$stats"
  done
done
echo
echo "Pick the adapter whose trimmed distribution peaks near 30 nt with a high 26-34 fraction."
echo "Then: ADAPTER=<chosen> sbatch --array=0-1 scripts/align_wang_ribo.sbatch"
