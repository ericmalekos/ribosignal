#!/usr/bin/env bash
#SBATCH --job-name=janich_pack
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/pack_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/pack_%j.err
set -euo pipefail
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
REF=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/heldout_refs
$PY /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/scripts/prepare/prepare_pack.py   --ribo-psites /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human/data/ribocode_mouse_liver/Liver_5samp   --rna-coverage /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/heldout_rnaseq_coverage/mouse_janich_liver   --universe-tx $REF/mouse_janich_universe_tx.txt   --fasta $REF/mouse_janich_universe.fa --build-orf-track --kozak none   --tx2biotype /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/tx2biotype_mouse.tsv   --group janich_liver --species mouse --out /private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/packed_heldout_mouse_janich_liver
