#!/usr/bin/env bash
# Orchestrate the Kozak-arms eval chain on SLURM. Per arm, three dependent stages:
#   (1) eval_kozak_arms.sbatch     [gpu]    dump pred_profiles.npz + localization AUROC
#   (2) ribocode_dropin_ctg.sbatch [medium] ATG+CTG RiboCode drop-in on the predicted profiles
#   (3) ctg_compare.sbatch         [short]  start-codon-stratified ATG-vs-CTG compare -> dropin_ctg_metrics.json
# Stage 2 for arm i depends on the arm-i element of the stage-1 array (afterok:<jid>_<i>); stage 3 on
# the whole stage-2 array. After ALL arms' stage 3 land, run:  aggregate_kozak.py + plot_kozak_weights.py
#
# Usage:  run_kozak_eval_chain.sh [arm indices ...]     (default: 0 1 2 3)
#   e.g.  run_kozak_eval_chain.sh 0 2 3     # v0/v2/v3 now (final); then  run_kozak_eval_chain.sh 1  when v1 lands
# Arm index -> slug: 0=v0_heuristic 1=v1_nokozak 2=v2_pwm 3=v3_learned.
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
SLUGS=(v0_heuristic v1_nokozak v2_pwm v3_learned)
ARMS=("$@"); [ ${#ARMS[@]} -eq 0 ] && ARMS=(0 1 2 3)

# Guard: every requested arm must have a FINAL best.pt (do not eval a mid-training checkpoint).
for i in "${ARMS[@]}"; do
  run=$NEW/results/kozak/${SLUGS[$i]}_onehot_fib2hep
  [ -f "$run/best.pt" ] || { echo "ERROR: no best.pt for arm $i (${SLUGS[$i]}) at $run -- training not done?"; exit 1; }
done

alist=$(IFS=,; echo "${ARMS[*]}")
echo "== stage 1: eval (dump + localization) for arms [$alist] =="
EVJ=$(sbatch --parsable --array="$alist" "$NEW/scripts/eval_kozak_arms.sbatch")
echo "  eval array job: $EVJ (elements: $alist)"

declare -a CMPJOBS=()
for i in "${ARMS[@]}"; do
  slug=${SLUGS[$i]}
  RUN=results/kozak/${slug}_onehot_fib2hep
  echo "== arm $i ($slug) =="
  DJ=$(sbatch --parsable --dependency="afterok:${EVJ}_${i}" \
       --export=ALL,RUN="$RUN" "$NEW/scripts/ribocode_dropin_ctg.sbatch")
  echo "  stage 2 drop-in (dep eval ${EVJ}_${i}): $DJ"
  CJ=$(sbatch --parsable --dependency="afterok:${DJ}" \
       --export=ALL,RUN="$RUN" "$NEW/scripts/ctg_compare.sbatch")
  echo "  stage 3 compare (dep drop-in ${DJ}): $CJ"
  CMPJOBS+=("$CJ")
done

echo
echo "submitted. stage-3 (compare) job ids: ${CMPJOBS[*]}"
echo "when ALL of those finish, run on the head node:"
echo "  PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3"
echo "  \$PY $NEW/scripts/plot_kozak_weights.py     # refresh the V3 learned-weights figure"
echo "  \$PY $NEW/scripts/aggregate_kozak.py         # -> results/kozak/kozak_summary.md"
