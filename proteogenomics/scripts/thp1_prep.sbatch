#!/usr/bin/env bash
#SBATCH --job-name=thp1_prep
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/thp1_prep_%j.out
set -euo pipefail
# THP-1 universe + coverage-only pack + no-Kozak ORF track.
#
# THIS IS A SCRIPT FILE, NOT `sbatch --wrap`, ON PURPOSE. --wrap runs the body under `sh` (dash on
# this cluster), where `set -o pipefail` is an illegal option; the first attempt died instantly with
# "Illegal option -o pipefail" and took the whole afterok chain down with it
# (DependencyNeverSatisfied). Anything needing bash syntax gets a real script with a shebang.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
S=$NEW/proteogenomics/scripts
P=$NEW/proteogenomics/data/THP1_pilot
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
echo "host=$(hostname) start=$(date -Iseconds)"

# --- 0. library chemistry gate (project rule feedback_rnaseq_must_be_polya) ---
# Runs BEFORE anything expensive: a total-RNA library must not reach the model.
"$PY" "$S/polya_check.py" --quant "$P/salmon/THP1/quant.sf"

# --- 1. expressed universe from the THP-1 salmon quant ---
"$PY" "$S/build_line_universe.py" --salmon "$P/salmon/THP1/quant.sf" \
    --out_prefix "$P/THP1" --species human
echo "universe: $(wc -l < "$P/THP1_universe_tx.txt") tx"

# --- 2. coverage-only pack + no-Kozak ORF track ---
# --coverage-only: there is no THP-1 Ribo-seq, so target_counts stays zero and the pack carries only
# the RNA coverage the model consumes.
# --kozak none per feedback_kozak_never_set: the released models trained on the no-Kozak track, and
# feeding them the heuristic track is a train/inference mismatch. prepare_pack names the track by its
# content, so it lands as orf_track_v2_nokozak.npy.
"$PY" "$NEW/scripts/prepare/prepare_pack.py" \
    --rna-coverage "$P/coverage/THP1_coverage.hd5" \
    --universe-tx "$P/THP1_universe_tx.txt" --coverage-only \
    --build-orf-track --fasta "$P/THP1_universe.fa" --kozak none \
    --tx2biotype "$NEW/data/tx2biotype.tsv" \
    --group THP1 --species human --out "$NEW/data/packed_heldout_human_thp1"
echo "done=$(date -Iseconds)"
ls -la "$NEW/data/packed_heldout_human_thp1" | tail -n +2
