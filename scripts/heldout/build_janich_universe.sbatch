#!/usr/bin/env bash
# ============================================================================
# DO NOT RUN -- RETIRED 2026-08-19.
#
# This points at expression_context_human/.../salmon_quant_mouse_liver_decoy/ZT0_rep1/quant.sf,
# which is the PRE-adapter-fix quant (0.03% salmon mapping -- the corruption fixed by task 35).
# That tree was deleted 2026-08-15 and the path is dead. It was deliberately NOT repointed at the
# live data/tpm/salmon_quant_janich_liver/, because doing so would make a script that must not run
# look runnable and would silently produce a DIFFERENT universe.
#
# The live 10,431-tx Janich universe was rebuilt post-fix from data/tpm/salmon_quant_janich_liver/.
# See scripts/janich_decon_universe_pack.sbatch (lines 24-28) for the original note.
# ============================================================================
#SBATCH --job-name=janich_univ
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/universe_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/universe_%j.err
set -euo pipefail
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
$PY /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/scripts/build_line_universe.py   --salmon /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human/data/salmon_quant_mouse_liver_decoy/ZT0_rep1/quant.sf   --species mouse --out_prefix /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/heldout_refs/mouse_janich
echo "janich universe -> /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/heldout_refs/mouse_janich_universe_tx.txt"
