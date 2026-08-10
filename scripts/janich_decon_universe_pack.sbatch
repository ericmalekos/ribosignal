#!/usr/bin/env bash
#SBATCH --job-name=janich_decon_pack
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich_decon_pack_%j.out
set -euo pipefail
# Rebuild the Janich liver universe + pack on the DECONTAMINATED salmon quants.
#
# WHY. All 7 Janich replicates put 52-70% of salmon TPM into two structural transcripts (Gm59647
# 7SL/SRP dup, n-TKctt14 tRNA-derived), both typed `lncRNA` so no biotype filter drops them. TPM is
# compositional, so they suppressed ~4,900 real transcripts per replicate below the TPM >= 1 cut:
# 11,532 mean before vs 16,379 after removing the reads.
#
# SCOPE: this corrects UNIVERSE SELECTION only. The per-nt coverage the model consumes was always
# clean -- those two transcripts carry 0.000% of the pack's coverage mass, because STAR runs with
# --outFilterMultimapNmax 1 plus the ncRNA-locus filter while salmon has no multimapper filter and
# uses length-normalised EM. No trained model or prediction was ever affected.
#
# NON-DESTRUCTIVE: the existing data/packed_heldout_mouse_janich_liver and
# data/heldout_refs/mouse_janich_universe_tx.txt are left untouched; everything lands in *_decon, so
# the liver ORF-call comparison can be run both ways and the difference attributed.
#
# NOTE scripts/heldout/build_janich_universe.sbatch is STALE -- it points at
# expression_context_human/.../salmon_quant_mouse_liver_decoy/ZT0_rep1/quant.sf, which is the
# PRE-adapter-fix quant (0.03% salmon mapping, the corruption task 35 fixed). The live 10,431-tx
# universe was rebuilt post-fix from data/salmon_quant_janich_liver/. Do not re-run that script.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
S=$NEW/proteogenomics/scripts
Q=$NEW/data/salmon_quant_janich_liver_decontam
REF=$NEW/data/heldout_refs
COV=$NEW/data/heldout_rnaseq_coverage/mouse_janich_liver
OLDPACK=$NEW/data/packed_heldout_mouse_janich_liver
OUTPACK=${OUTPACK:-$NEW/data/packed_heldout_mouse_janich_liver_decon}
# The Ribo P-sites live in the SIBLING project (Liver_5samp), not under data/heldout_psites/. Rather
# than hard-code that path and risk drifting from what the live pack actually consumed, read the
# exact input list out of the existing pack's provenance.json -- same 5 P-site hd5, so the ONLY thing
# that changes between old and new pack is the universe.
mapfile -t PS < <("$PY" -c "
import json;print('\n'.join(d['path'] for d in json.load(open('$OLDPACK/provenance.json'))['ribo_psites_inputs']))")
[[ ${#PS[@]} -ge 1 ]] || { echo 'could not read ribo_psites_inputs from the old pack' >&2; exit 1; }
for p in "${PS[@]}"; do [[ -s "$p" ]] || { echo "missing P-site input $p" >&2; exit 1; }; done
echo "reusing ${#PS[@]} Ribo P-site inputs from the existing pack's provenance"
export PYTHONPATH="$S:${PYTHONPATH:-}"
echo "host=$(hostname) start=$(date -Iseconds)"

mapfile -t QUANTS < <(ls "$Q"/*/quant.sf | sort)
[[ ${#QUANTS[@]} -eq 7 ]] || { echo "expected 7 decontaminated quants, got ${#QUANTS[@]}" >&2; exit 1; }
echo "quants: ${#QUANTS[@]}"

# ---- 1. universe (mean TPM across the 7 replicates; pgx.universe supports multi-quant merge) ----
cd "$S" && "$PY" -m pgx.universe --salmon "${QUANTS[@]}" --species mouse \
    --out-prefix "$REF/mouse_janich_decon" 2>&1 | tail -5
cd "$NEW"
NTX=$(wc -l < "$REF/mouse_janich_decon_universe_tx.txt")
echo "decontaminated universe: $NTX tx  (was $(wc -l < "$REF/mouse_janich_universe_tx.txt"))"
[[ $NTX -gt 6000 ]] || { echo "ABORT: universe $NTX looks crushed" >&2; exit 1; }

# ---- 2. pack on the new universe, reusing the UNCHANGED Ribo P-sites and RNA coverage ----
# --kozak none per feedback_kozak_never_set.
"$PY" "$NEW/scripts/prepare/prepare_pack.py" \
    --ribo-psites "${PS[@]}" --rna-coverage "$COV" \
    --universe-tx "$REF/mouse_janich_decon_universe_tx.txt" \
    --build-orf-track --fasta "$REF/mouse_janich_decon_universe.fa" --kozak none \
    --tx2biotype "$NEW/data/tx2biotype_mouse.tsv" \
    --group janich_liver_decon --species mouse --out "$OUTPACK" 2>&1 | tail -8
echo "done=$(date -Iseconds)"
ls -la "$OUTPACK" | tail -n +2
