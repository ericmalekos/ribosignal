#!/usr/bin/env bash
# Submit the standing-rule ORF-call evaluation for BOTH mm10-coverage retrain arms.
#
# Reuses scripts/eval_union.sbatch unchanged via its PSUF/RUN/TAG/LABEL overrides, so both arms get
# the identical eval the union mm1 baseline got: the >5%-skip dump guard, RiboCode drop-in on
# real / pred_obsdepth / pred_preddepth, and the two-arm CDS-anchored Poisson sweep
# (theta 0.02/0.05/0.1/0.2) required by the two-arm ORF-calling rule.
#
# Run only after BOTH trainings finish -- each needs its own best.pt.
#   bash scripts/eval_mm10cov_arms.sh            # submit both
#   ARMS=A bash scripts/eval_mm10cov_arms.sh     # just one
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
B=orf_v2_mamba4_onehot_union_noBrain_nokozak_mm10cov
for arm in ${ARMS:-A B}; do
  case "$arm" in
    A) RUN=results/loto/${B}_holdout_Hepatocytes;           LBL="mamba4 onehot mm10cov armA" ;;
    B) RUN=results/loto/${B}_exclmm50_holdout_Hepatocytes;  LBL="mamba4 onehot mm10cov armB excl" ;;
    *) echo "bad arm $arm" >&2; exit 1 ;;
  esac
  # Gate on test_metrics.json, NOT best.pt: train_loto.py writes best.pt after epoch 0, so best.pt
  # is present within ~28 min of launch and would submit an eval of a 1-epoch model -- which also
  # appends a bogus row to docs/manuscript_orf_call_eval.md. test_metrics.json is written only after
  # the training loop exits, so it is the correct "training finished" sentinel.
  if [ ! -s "$NEW/$RUN/test_metrics.json" ]; then
    ep=$(python3 -c "import json;print(len(json.load(open('$NEW/$RUN/history.json'))))" 2>/dev/null || echo 0)
    echo "SKIP arm $arm: training not finished (no test_metrics.json; $ep epochs so far)" >&2
    continue
  fi
  jid=$(BACKEND=onehot MIXER=mamba N_ATTN=4 PSUF=mm10cov RUN="$RUN" \
        TAG="${B}_arm${arm}" LABEL="$LBL" \
        sbatch --parsable "$NEW/scripts/eval_union.sbatch")
  echo "  arm $arm -> job $jid   ($RUN)"
done
