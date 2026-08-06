#!/usr/bin/env bash
#SBATCH --job-name=ctg_compare
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/ctg_compare_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/dropin/ctg_compare_%j.err
set -euo pipefail
# Runs AFTER ribocode_dropin_ctg.sbatch. (1) compare_dropin_calls.py on dropin_ctg/ = the combined
# ATG+CTG drop-in numbers (real-with-CTG as reference); (2) compare_dropin_ctg.py = the start-codon-
# stratified ATG-vs-CTG breakdown (the actual answer). Submit with:
#   sbatch --dependency=afterok:<ctg_array_jobid> --export=ALL,RUN=<run_dir> ctg_compare.sbatch
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
RUN=${RUN:-results/loto/orf_v2_attn_onehot_holdout_Hepatocytes}
CTGDIR=$NEW/$RUN/dropin_ctg
# --official row is secondary here; default to this run's OWN ATG-only real calls, so real_vs_official
# = "how many calls does enabling CTG add over the ATG baseline". TX2GENE/FASTA empty = human defaults;
# set both for a non-human (mouse vM38) held-out.
OFFICIAL=${OFFICIAL:-$NEW/$RUN/dropin/real_collapsed.txt}
PROFILES=$NEW/$RUN/dropin/pred_profiles.npz
TX2GENE=${TX2GENE:-}
FASTA=${FASTA:-}
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
export PYTHONNOUSERSITE=1
echo "run=$RUN tx2gene=${TX2GENE:-<human default>} host=$(hostname) start=$(date -Iseconds)"

T2G_ARG=(); [ -n "$TX2GENE" ] && T2G_ARG=(--tx2gene "$TX2GENE")
FASTA_ARG=(); [ -n "$FASTA" ] && FASTA_ARG=(--fasta "$FASTA")

if [ -f "$OFFICIAL" ]; then
  echo "############ (1) COMBINED ATG+CTG drop-in (genomic key, 0.5x enrich, 90 nt) ############"
  "$PY" "$NEW/scripts/compare_dropin_calls.py" \
    --dropin_dir "$CTGDIR" --official "$OFFICIAL" --profiles "$PROFILES" \
    --key genomic --min_enrichment 0.5 --min_len 90 --out "$CTGDIR" "${T2G_ARG[@]}"
else
  echo "############ (1) SKIPPED: no ATG-only reference at $OFFICIAL ############"
  echo "(expected for the Kozak arms, which run only the CTG drop-in; the stratified step (2) below"
  echo " already reports ATG-vs-CTG, so the combined step is not needed here.)"
fi

echo
echo "############ (2) START-CODON-STRATIFIED ATG vs CTG ############"
"$PY" "$NEW/scripts/compare_dropin_ctg.py" \
  --dropin_dir "$CTGDIR" --profiles "$PROFILES" --min_enrichment 0.5 --min_len 90 --out "$CTGDIR" \
  "${T2G_ARG[@]}" "${FASTA_ARG[@]}"
echo "done=$(date -Iseconds)"
