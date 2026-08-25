#!/usr/bin/env bash
#SBATCH --job-name=prep_line_full
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=8:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/prep_full_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/prep_full_%j.err
set -euo pipefail
# Expressed universe + coverage-only pack + no-Kozak ORF track, for ONE line at ONE multimap posture,
# built from the FULL-DEPTH ALL-RUN salmon quant and coverage.
#
# Generalises thp1_prep.sbatch. Two corrections baked in:
#   * PY is `riboseq`, not `cas12a`. prepare_pack.py -> packlib -> `import h5py`, and cas12a has
#     torch but NO h5py, so the cas12a path dies instantly on import.
#   * The universe is rebuilt PER RUN-SET. build_line_universe.py derives it from the salmon quant,
#     so adding runs changes the universe -- the existing packs/ORF tracks are NOT reusable and the
#     mm1 arm must be rebuilt too, not just mm10.
#
# The ORF track is built per pack (--build-orf-track --kozak none) because tx_order changes with the
# universe; a track from another pack would be silently misaligned.
#
# Env: LINE, MM (1|10), SPECIES (default human).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
S=$NEW/proteogenomics/scripts
PY=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/python
LINE=${LINE:?set LINE}; MM=${MM:?set MM=1 or 10}; SPECIES=${SPECIES:-human}
P=$NEW/proteogenomics/data/${LINE}_pilot
SAL=$P/salmon_full/${LINE}/quant.sf
COV=$P/coverage_mm${MM}_full/${LINE}_coverage.hd5
OUT=$NEW/data/packed_heldout_${LINE,,}_mm${MM}_full
export PYTHONNOUSERSITE=1
echo "host=$(hostname) LINE=$LINE MM=$MM start=$(date -Iseconds)"
[ -s "$SAL" ] || { echo "missing salmon quant $SAL" >&2; exit 1; }
[ -s "$COV" ] || { echo "missing coverage $COV" >&2; exit 1; }
# 0. library chemistry gate -- a total-RNA library must never reach the model (poly-A rule)
"$PY" "$S/polya_check.py" --quant "$SAL"
# 1. expressed universe from THIS run-set's quant
"$PY" "$S/build_line_universe.py" --salmon "$SAL" --out_prefix "$P/${LINE}_full" --species "$SPECIES"
echo "  universe: $(wc -l < "$P/${LINE}_full_universe_tx.txt") tx"
# 2. coverage-only pack + its own no-Kozak ORF track
"$PY" "$NEW/scripts/prepare/prepare_pack.py" \
    --rna-coverage "$COV" \
    --universe-tx "$P/${LINE}_full_universe_tx.txt" --coverage-only \
    --build-orf-track --fasta "$P/${LINE}_full_universe.fa" --kozak none \
    --group "${LINE}_mm${MM}_full" --species "$SPECIES" --out "$OUT"
[ -s "$OUT/coverage.npy" ] || { echo "no coverage.npy in $OUT" >&2; exit 2; }
ls -1 "$OUT" | sed 's/^/  /'
echo "done=$(date -Iseconds) -> $OUT"
