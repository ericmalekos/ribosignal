#!/usr/bin/env bash
# THP-1 end-to-end, SLURM-native afterok chain (no head-node babysitting):
#   universe -> pack (+ nokozak ORF track) -> GPU predict (BOTH released models) -> pgx (both models)
#
# Prerequisite: process_rnaseq_thp1.sbatch has finished, leaving
#   THP1_pilot/salmon/THP1/quant.sf        (isoform TPM -> universe + the null databases)
#   THP1_pilot/coverage/THP1_coverage.hd5  (per-nt RNA depth -> THE model input)
#
# WHY NOT rerun_mm1_line.sh / post_predict_line.sbatch. That chain is the RETIRED f0-threshold route
# (enumerate_score_orfs -> build_a549_dbs -> MS2Rescore -> compare_line). Everything current goes
# through `pgx`, which selects ORFs by RiboCode call rather than by an f0 cut and builds both
# threshold arms. Using the old chain here would produce a THP-1 row that cannot sit in the same
# table as A549 / HBL-1 / SU-DHL-4 / DoHH2.
#
# ENZYME IS NONSPECIFIC (HLA-I peptides are not tryptic), which means:
#   * the frozen HLA template is used (fragger_hbl1_hla_frozen.params), and
#   * `null_nc` is dropped automatically by pgx.run -- a near-cognate null under nonspecific digestion
#     enumerates ~2.86e9 peptides against MSFragger's ~2e9 cap. `null_atg` is the null. The ncStart
#     columns survive because they are a substring scan over the DB, not a search.
#
# Usage: thp1_chain.sh [<dependency jobid>]
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
S=$NEW/proteogenomics/scripts
P=$NEW/proteogenomics/data/THP1_pilot
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
SIFDIR=/private/groups/carpenterlab/emalekos/singularity_cache
DEP=${1:-}
dep() { [[ -n "$1" ]] && echo "--dependency=afterok:$1" || echo ""; }

# ---- 1. universe + pack + nokozak ORF track (light; one CPU job) ----
# A real script file, NOT `sbatch --wrap`: --wrap bodies run under `sh` (dash), which rejects
# `set -o pipefail`. The first attempt died on that in under a second and took the whole afterok
# chain with it (every downstream job went DependencyNeverSatisfied).
PREP=$(sbatch --parsable $(dep "$DEP") $S/thp1_prep.sbatch)
echo "prep=$PREP"

# ---- 2. predict, one GPU job per released model ----
D_M=$(sbatch --parsable --dependency=afterok:$PREP --job-name=thp1_pred_mamba4 \
  --export=ALL,LINE=THP1,HELDOUT=human_thp1,OUTSUB=pred_mamba4,KOZAK=nokozak,\
MODEL_RUN=results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes,\
SIF=$SIFDIR/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img \
  $S/dump_line.sbatch)
D_A=$(sbatch --parsable --dependency=afterok:$PREP --job-name=thp1_pred_attn \
  --export=ALL,LINE=THP1,HELDOUT=human_thp1,OUTSUB=pred_attn,KOZAK=nokozak,\
MODEL_RUN=results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes \
  $S/dump_line.sbatch)
echo "predict: mamba4=$D_M attn=$D_A"

# ---- 3. pgx per model, each with its OWN search root ----
#
# NOT using --shared-search-root here, deliberately, even though gencode and null_atg are
# byte-identical between the two models. `pipeline.sbatch` submits its searches and exits, so
# `afterok` on it means "the other pipeline SUBMITTED", not "its searches finished". Both runs would
# then submit searches for the same arm into the same directory: two writers on one path, which the
# project's concurrency rule forbids. The `.pgx_hash` guard skips a COMPLETED search, not a running
# one, so it does not close that window.
#
# The cost of not sharing is one duplicate gencode + null_atg search. Measured on the other
# nonspecific HLA pilots that is ~5 min per arm, i.e. ~10-20 min total. Cheap insurance against a
# corrupt search directory of exactly the right size.
G_M=$(SPECIES=human LABEL=THP1 PROFILES=$P/pred_mamba4/pred_profiles.npz \
  UNIVERSE_FA=$P/THP1_universe.fa SALMON=$P/salmon/THP1/quant.sf \
  MZML="$P/mzml/*.mzML" ENZYME=nonspecific OUT=$P/pgx_mamba4 PGX_SCRIPTS=$S \
  sbatch --parsable --dependency=afterok:$D_M --job-name=thp1_pgx_mamba4 $S/pgx/pipeline.sbatch)
G_A=$(SPECIES=human LABEL=THP1 PROFILES=$P/pred_attn/pred_profiles.npz \
  UNIVERSE_FA=$P/THP1_universe.fa SALMON=$P/salmon/THP1/quant.sf \
  MZML="$P/mzml/*.mzML" ENZYME=nonspecific OUT=$P/pgx_attn PGX_SCRIPTS=$S \
  sbatch --parsable --dependency=afterok:$D_A --job-name=thp1_pgx_attn $S/pgx/pipeline.sbatch)
echo "pgx: mamba4=$G_M attn=$G_A"
echo "THP-1 chain: prep=$PREP -> predict($D_M,$D_A) -> pgx($G_M,$G_A)"
