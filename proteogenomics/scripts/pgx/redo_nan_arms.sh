#!/usr/bin/env bash
# Re-run the extension arms that the bh() NaN bug silently zeroed, then rebuild their
# model_standard databases.
#
# The bug: numpy sorts NaN last, so bh()'s reverse cumulative minimum started on a NaN p-value and
# propagated it through every q-value. One degenerate Wilcoxon (flat extension region) nulled the
# whole arm -- reported as "passing: 0" beside "would_pass_whole_orf_f0_0.5: 4434", which reads as a
# finding rather than arithmetic. Fixed in pgx.seqtools.bh; pinned by test_bh_is_nan_safe.
#
# FULL BLAST RADIUS (2026-08-07). Scanned EVERY extension_summary.json under proteogenomics/data --
# 52 arms -- for the signature `mtime < fix && tested > 0 && passing == 0`. Five arms were affected,
# all of them `standard`; every `poisson` arm was clean, so the headline results do not move.
#
#   1-3. DoHH2/attn, SUDHL4/attn, B721/mamba4        -- fixed by the loop below
#   4.   macrophage pgx_attn_union/BMDM/standard      -- fixed 2026-08-06 (found on the second pass)
#   5.   macrophage pgx_genetype/attn/BMDM/standard   -- fixed 2026-08-07 (found on the THIRD pass,
#        by scanning the whole tree instead of the pilots). 0 -> 4,455 passing of 9,882 tested; its
#        db_model_standard.fasta went 126,700 -> 135,560 entries (novel 6,124 -> 10,554).
#        This one mattered beyond bookkeeping: the genetype table compared attn against mamba4, and
#        attn's model_standard database was missing its ENTIRE N-terminal-extension class while
#        mamba4's was not. That comparison was not a model difference.
#
# LESSON, recorded because it cost three passes: scoping a blast-radius scan to "the datasets I was
# just working on" finds the arms you already suspected. Enumerate the artifact type across the
# whole tree and let the signature decide.
#
# The two macrophage arms are NOT re-run by this script; they were fixed in place because they need
# per-population profile/universe paths that this pilot-shaped loop does not carry.
set -uo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PD=$NEW/proteogenomics/data
RPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
cd "$NEW/proteogenomics/scripts"
export PYTHONPATH="$NEW/proteogenomics/scripts:${PYTHONPATH:-}"

for spec in "DoHH2:attn:pred_attn" "SUDHL4:attn:pred_attn" "B721:mamba4:pred_mamba4"; do
  L=${spec%%:*}; M=$(echo "$spec" | cut -d: -f2); PRED=${spec##*:}
  P=$PD/${L}_pilot; G=$P/pgx_$M
  UNI=$P/${L}_universe.fa
  CALLS=$G/calls/standard/theta_1/pred_preddepth_collapsed.txt
  OUT=$G/extensions/standard
  echo "=== $L / $M / standard ==="
  [[ -s "$CALLS" ]] || { echo "  ABORT: no $CALLS"; continue; }
  before=$( [ -s "$OUT/extension_summary.json" ] && $PY -c "import json;print(json.load(open('$OUT/extension_summary.json'))['passing'])" || echo "?" )
  "$RPY" -m pgx.extensions --profiles "$P/$PRED/pred_profiles.npz" --calls-collapsed "$CALLS" \
      --fasta "$UNI" --species human --theta 1.0 --start-codons near_cognate --out "$OUT" \
      --ext-qvalue 0.05 --orf-qvalue 0.05 --min-ext-aa 5 --min-nonzero 5 2>&1 | grep -iE "tested|passing|wrote" | sed 's/^/  /'
  after=$($PY -c "import json;print(json.load(open('$OUT/extension_summary.json'))['passing'])")
  echo "  passing: $before -> $after"
  "$PY" -m pgx.build_dbs --species human --out "$G/db" --calls "$G/calls/calls_standard.tsv" \
      --extensions "$OUT/extensions.tsv" --arm standard --dbs model --min-aa 7 2>&1 | grep -iE "db_model" | sed 's/^/  /'
done
