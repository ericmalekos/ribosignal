#!/usr/bin/env bash
set -uo pipefail
# Phase 4.1/4.2: per-codon occupancy across every dump that has a known CDS table and universe FASTA.
#
# The dumps span two species and several universes, so each entry pairs an npz with the CDS table and
# FASTA it was actually built against. Getting that pairing wrong yields codon assignments shifted
# relative to the true reading frame -- which would look like a real (and wrong) result rather than an
# error, so the pairing is explicit here rather than inferred from a path.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
MCDS=$NEW/data/mouse_ribocode_annot/transcripts_cds.txt
HCDS=$NEW/data/human_ribocode_annot_primary/transcripts_cds.txt
REF=$NEW/data/heldout_refs
OUT=$NEW/results/codon_occupancy
mkdir -p "$OUT"

run () {  # label npz cds fasta
  local lab=$1 npz=$2 cds=$3 fa=$4
  [[ -s "$npz" && -s "$cds" && -s "$fa" ]] || { echo "  SKIP $lab (missing input)"; return; }
  [[ -s "$OUT/${lab}.tsv" ]] && { echo "  skip $lab (done)"; return; }
  "$PY" "$NEW/scripts/codon_occupancy.py" --npz "$npz" --cds "$cds" --fasta "$fa" \
    --label "$lab" --out "$OUT/${lab}.tsv"
}

# ---- mouse liver 3x3: the matched (diagonal) cells, both models ----
for M in attn mamba4; do
  for D in janich gse243134 wang; do
    run "mouse3x3_${M}_${D}" \
      "$NEW/results/mouse_liver_3x3/${M}_ribo-${D}_rna-${D}/pred_profiles.npz" \
      "$MCDS" "$REF/mouse_liver_union3_universe.fa"
  done
done

# ---- human ORF calls: the two datasets with a real observed arm ----
run "human_gse208041_attn"   "$NEW/results/human_orf_calls/attn_gse208041/pred_profiles.npz"   "$HCDS" "$REF/human_thp1_gse208041_universe.fa"
run "human_gse208041_mamba4" "$NEW/results/human_orf_calls/mamba4_gse208041/pred_profiles.npz" "$HCDS" "$REF/human_thp1_gse208041_universe.fa"
run "human_cart_attn"        "$NEW/results/human_orf_calls/attn_cart/pred_profiles.npz"        "$HCDS" "$REF/human_cart_gse304796_universe.fa"
run "human_cart_mamba4"      "$NEW/results/human_orf_calls/mamba4_cart/pred_profiles.npz"      "$HCDS" "$REF/human_cart_gse304796_universe.fa"

# ---- released-model mouse liver runs on their own universes ----
run "released_wang_mamba4"      "$NEW/results/liver_released/mamba4_mouse_wang_liver/pred_profiles.npz"        "$MCDS" "$REF/mouse_wang_universe.fa"
run "released_gse243134_mamba4" "$NEW/results/liver_released/mamba4_mouse_gse243134_liver/pred_profiles.npz"   "$MCDS" "$REF/mouse_gse243134_universe.fa"
run "released_janich_mamba4"    "$NEW/results/liver_released/mamba4_mouse_janich_liver_decon/pred_profiles.npz" "$MCDS" "$REF/mouse_janich_decon_universe.fa"

echo "done -> $OUT ($(ls "$OUT"/*.tsv 2>/dev/null | wc -l) tables)"
