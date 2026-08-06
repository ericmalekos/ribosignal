#!/usr/bin/env bash
# Characterize MS2Rescore/mokapot run-to-run variance in the novel-peptide discovery counts (Fig 2b error
# bars). For each immunopeptidome, re-run the mokapot rescore N times on the SAME (fixed mm1) search, each
# into its own rescore_seed<K>/ dir, then an aggregator computes mean +/- SD of the class-FDR novel-peptide
# counts per DB. The search is deterministic; only mokapot is re-rolled, so the spread is pure rescoring
# stochasticity. Usage: run_mokapot_stochastic.sh [N=8]
set -euo pipefail
N=${1:-8}
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
S=$NEW/proteogenomics/scripts
DATASETS=(HBL1 DoHH2 SUDHL4)
deps=""
for L in "${DATASETS[@]}"; do
  P=$NEW/proteogenomics/data/${L}_pilot
  NR=$(ls "$P"/immunopeptidome/*.mzML | wc -l); NA=$((3*NR-1))
  for K in $(seq 1 "$N"); do
    rm -rf "$P/rescore_seed${K}"          # fresh mokapot roll each seed
    J=$(sbatch --parsable --array=0-$NA --export=ALL,PILOT=$P,RTAG=_seed${K} $S/ms2rescore_line_fanout.sbatch)
    deps="${deps}:${J}"
    echo "  $L seed $K -> rescore array $J (0-$NA) -> rescore_seed${K}/"
  done
done
AGG=$(sbatch --parsable --dependency=afterok${deps} --export=ALL,N=$N $S/aggregate_mokapot_stochastic.sbatch)
echo "aggregator=$AGG (afterok all ${N}x${#DATASETS[@]} rescore arrays)"
