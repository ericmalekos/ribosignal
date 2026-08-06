#!/usr/bin/env bash
#SBATCH --job-name=repack_union_oth
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=6
#SBATCH --mem=200G
#SBATCH --time=10:00:00
#SBATCH --array=1-7
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/prepare/repack_union_oth_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/prepare/repack_union_oth_%A_%a.err
set -euo pipefail
# Broader-universe repack, the 7 non-Fibroblast tissues: build each UNION pack by pooling its per-sample Ribo
# psites (target) + mm1 RNA coverage, reusing Fibroblast's union tx_order (--ref-pack) so all packs are row-
# aligned to the same 84,472-tx axis + the shared ORF track. Run after repack_union_fib (dependency).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
PER=$NEW/data/rnaseq_coverage/per_sample
MAP=$NEW/data/loto_rnaseq_srr_tissue.tsv
REF=$NEW/data/packed_union
TISSUES=(VSMC ES Fat HA_EC HCAEC HUVEC Hepatocytes)
T=${TISSUES[$((SLURM_ARRAY_TASK_ID - 1))]}
OUT=$NEW/data/packed_union_$T
SRRS=$(awk -F'\t' -v t="$T" '$2==t{print $1}' "$MAP")
COV_ARGS=""
for s in $SRRS; do f="$PER/${s}_coverage.hd5"; [ -s "$f" ] || { echo "MISSING cov $f" >&2; exit 1; }; COV_ARGS="$COV_ARGS $f"; done
PSITES="$ECH/data/ribocode_per_tissue/$T"/*_psites.hd5
echo "[$T] union repack: cov_samples=$(echo $SRRS|wc -w) psites_samples=$(ls $PSITES|wc -l) host=$(hostname) start=$(date -Iseconds)"
# shellcheck disable=SC2086
"$PY" "$NEW/scripts/prepare/prepare_pack.py" \
  --ribo-psites $PSITES \
  --rna-coverage $COV_ARGS \
  --ref-pack "$REF" \
  --tx2biotype "$NEW/data/tx2biotype.tsv" \
  --group "$T" --species human --out "$OUT"
echo "[$T] done=$(date -Iseconds) -> $OUT ($(wc -l < $OUT/tx_order.txt) tx)"
