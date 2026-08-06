#!/usr/bin/env bash
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
