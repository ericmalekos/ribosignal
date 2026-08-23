#!/usr/bin/env bash
#SBATCH --job-name=repack_union_fib
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=6
#SBATCH --mem=200G
#SBATCH --time=10:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/prepare/repack_union_fib_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/prepare/repack_union_fib_%j.err
set -euo pipefail
# Broader-universe repack, REFERENCE tissue (Fibroblast): build the fresh 84,472-tx UNION pack from the
# per-sample Ribo psites (pooled target) + mm1 RNA coverage, and the shared no-Kozak ORF track. The other 7
# tissues reuse this pack's tx_order via --ref-pack (repack_union_others.sbatch). New universe -> no reuse-target.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
PER=$NEW/data/rnaseq_coverage/per_sample
MAP=$NEW/data/loto_rnaseq_srr_tissue.tsv
OUT=$NEW/data/packed_union
mkdir -p "$NEW/logs/prepare"
T=Fibroblast
# Fibroblast RNA coverage SRRs = per-sample files NOT in the tissue map
SRRS=$(comm -23 <(ls "$PER"/SRR*_coverage.hd5 | grep -oE 'SRR[0-9]+' | sort -u) <(cut -f1 "$MAP" | sort -u))
COV_ARGS=""
for s in $SRRS; do f="$PER/${s}_coverage.hd5"; [ -s "$f" ] || { echo "MISSING cov $f" >&2; exit 1; }; COV_ARGS="$COV_ARGS $f"; done
PSITES="/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/ribocode_per_tissue/$T"/*_psites.hd5
echo "[$T] union repack: cov_samples=$(echo $SRRS|wc -w) psites_samples=$(ls $PSITES|wc -l) host=$(hostname) start=$(date -Iseconds)"
# shellcheck disable=SC2086
"$PY" "$NEW/scripts/prepare/prepare_pack.py" \
  --ribo-psites $PSITES \
  --rna-coverage $COV_ARGS \
  --universe-tx "$NEW/data/union_tx_order.txt" \
  --tx2biotype "$NEW/data/tx2biotype.tsv" \
  --build-orf-track --fasta "$NEW/data/union_universe.fa" --kozak none \
  --group "$T" --species human --out "$OUT"
echo "[$T] done=$(date -Iseconds) -> $OUT"
ls -la "$OUT" | grep -E "tx_order|target|coverage|orf_track" | awk '{print "  "$5" "$NF}'
