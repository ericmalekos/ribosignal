#!/usr/bin/env bash
#SBATCH --job-name=l3x3_depth
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=10:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/l3x3_depth_%A_%a.out
set -euo pipefail
# Task #75: the MATCHED-DEPTH arm of the mouse-liver 3x3 factorial.
#
# WHY IT IS NOW NECESSARY, NOT OPTIONAL. The task deferred this pending the natural-depth result,
# with the criterion "if the depth ordering reappears there, this becomes necessary". It reappears
# unambiguously. Natural depths are wang 78.3M < janich 110.0M < gse243134 150.5M, and the
# observed-vs-observed matrix tracks that ordering rather than lab or protocol:
#
#   novel recall   wang AS THE SCORED SET:      0.612 (vs janich), 0.575 (vs gse243134)   <- low
#                  wang AS THE REFERENCE:       0.805 (janich),    0.834 (gse243134)      <- high
#
# A shallow library finds fewer ORFs, so it misses a deep reference's calls (low recall as the scored
# set) while a deep experiment recovers most of ITS few calls (high recall as the reference). That is
# the signature of depth, and E2 uses this matrix as "the between-experiment CEILING" in a headline
# figure -- so the ceiling currently confounds protocol difference with sequencing depth.
#
# DESIGN: all three thinned to 75M, INCLUDING WANG. Thinning only the two deeper sets would leave
# wang as the single arm carrying no thinning noise, which is a subtler version of the same
# asymmetry. 75M sits just below wang's 78.3M so every arm is thinned and none is special.
#
#   gse243134  f=0.4982     janich  f=0.6818     wang  f=0.9584
#
# THREE SEEDS, because binomial thinning is stochastic and a single draw can move a recall by more
# than the effect being measured. The spread across seeds IS the error bar; report it, never one seed.
#
# Uses ribocode_dropin.py --subsample, which bin-thins the REAL per-nt P-sites (each footprint kept
# i.i.d. with probability f). No new machinery; the same code path as the depth_crossover work.
#
#   sbatch --array=0-8 scripts/liver3x3_matched_depth.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
RIBOPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
ANNOT=$NEW/data/mouse_ribocode_annot
[[ -d "$ANNOT" ]] || ANNOT=$(ls -d $NEW/data/*mouse*ribocode*annot* 2>/dev/null | head -1)
RES=$NEW/results/mouse_liver_3x3_canon
OUT=$NEW/results/liver3x3_matched_depth
mkdir -p "$OUT" "$NEW/logs"

DS=(gse243134 janich wang)
FRAC=(0.4982 0.6818 0.9584)
SEEDS=(0 1 2)
T=${SLURM_ARRAY_TASK_ID:-0}
DI=$((T / 3)); SI=$((T % 3))
D=${DS[$DI]}; F=${FRAC[$DI]}; S=${SEEDS[$SI]}

PROF=$RES/mamba4_ribo-${D}_rna-${D}/pred_profiles.npz
for f in "$PROF" "$ANNOT/transcripts_cds.txt"; do
  [[ -s "$f" ]] || { echo "ABORT: missing $f" >&2; exit 1; }
done
O=$OUT/${D}_seed${S}
mkdir -p "$O"
echo "task=$T dataset=$D frac=$F seed=$S annot=$(basename $ANNOT) host=$(hostname) start=$(date -Iseconds)"

if [[ -s "$O/real_depth${F}_seed${S}_collapsed.txt" ]]; then
  echo "  exists, skipping"
else
  export PYTHONNOUSERSITE=1
  "$RIBOPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$PROF" --annot "$ANNOT" \
    --variant real --subsample "$F" --seed "$S" --out "$O" --min_aa 5 --pval 0.05
fi
C=$(ls "$O"/real_depth*_seed${S}_collapsed.txt 2>/dev/null | head -1)
[[ -s "$C" ]] || { echo "ABORT: no collapsed output in $O" >&2; exit 1; }
echo "  calls: $(($(wc -l < "$C") - 1))"
echo "done=$(date -Iseconds) -> $O"
