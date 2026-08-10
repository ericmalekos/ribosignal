#!/usr/bin/env bash
#SBATCH --job-name=ribo_dropin_ctg
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --array=0-2
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/ctg_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/ctg_%A_%a.err
set -euo pipefail
# Non-AUG (ATG + CTG) RiboCode drop-in: re-run detectORF with ALTERNATIVE_START_CODON_LIST=["CTG"]
# on the SAME predicted/real per-nt density (reuses the existing pred_profiles.npz -- no re-dump), so
# real-vs-predicted ORF calling can be compared on CTG-initiated ORFs (which the ATG-only default
# never enumerated). CTG is fallback-only in RiboCode (an ORF is called at a CTG only where its
# in-frame stop has no ATG), so CTG calls are a DISJOINT new set on top of the unchanged ATG calls.
# One array task per density variant; writes to dropin_ctg/ (separate dir; filenames are standard so
# compare_dropin_calls.py runs verbatim on it).
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
RUN=${RUN:-results/loto/orf_v2_attn_onehot_holdout_Hepatocytes}
ALT=${ALT:-CTG}                            # e.g. "CTG,GTG,TTG,ACG,ATA,ATT,ATC,AAG,AGG" for all near-cognates
PROFILES=$NEW/$RUN/dropin/pred_profiles.npz
ANNOT=${ANNOT:-$ECH/data/ribocode_annot}   # mouse: "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/data/mouse_ribocode_annot"
OUT=$NEW/$RUN/${OUTSUB:-dropin_ctg}        # OUTSUB=dropin_allalt for the all-near-cognate run
PY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
mkdir -p "$NEW/logs/dropin" "$OUT"

VARIANTS=(real pred_obsdepth pred_preddepth)
V=${VARIANTS[$SLURM_ARRAY_TASK_ID]}
echo "variant=$V alt=$ALT run=$RUN host=$(hostname) start=$(date -Iseconds)"
[ -f "$PROFILES" ] || { echo "ERROR: missing $PROFILES (run the ATG drop-in dump first)"; exit 1; }
export PYTHONNOUSERSITE=1
"$PY" "$NEW/scripts/ribocode_dropin.py" \
  --profiles "$PROFILES" --annot "$ANNOT" --variant "$V" \
  --alt_start_codons "$ALT" --out "$OUT" --min_aa 5 --pval 0.05
echo "done variant=$V $(date -Iseconds)"
