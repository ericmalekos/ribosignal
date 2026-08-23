#!/usr/bin/env bash
#SBATCH --job-name=rnaq_packs
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=8:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/rnaq_packs_%j.out
set -euo pipefail
# Task 68, stage 1: build the two NEW packs for the RNA-quality factorial.
#
# THE QUESTION. Every mouse-liver arm so far varies RNA-seq and Ribo-seq together, so "better RNA"
# has never been separated from "different experiment". Here the Ribo-seq is HELD FIXED -- the same 5
# Janich P-site hd5 in all three arms, read out of the existing pack's provenance rather than
# re-listed -- and only the RNA changes.
#
# RNA-seq enters the model twice, and the two paths are separable:
#   (1) UNIVERSE  -- salmon TPM >= 1 decides which transcripts exist at all
#   (2) COVERAGE  -- the per-nt RNA depth track the model consumes as an input channel
#
#   arm a  universe Janich-decon (14,887 tx)     coverage Janich        ALREADY BUILT
#   arm b  universe Janich-decon (14,887 tx)     coverage PRJEB34766    <- built here
#   arm c  universe PRJEB34766                   coverage PRJEB34766    <- built here
#
#   a vs b  = pure COVERAGE effect (universe held identical, so the transcript space and the ORF
#             track are byte-identical and the comparison is not confounded by what got called)
#   b vs c  = pure UNIVERSE effect
#   a vs c  = the full "swap the RNA-seq experiment" effect
#
# WHY PRJEB34766 IS THE COMPARATOR. Both Janich (50-51 nt) and GSE243134 (75-80 nt) RNA are SINGLE-end.
# PRJEB34766 is PAIRED, poly(A)-selected per ENA `library_selection`, 7 CreNeg control runs / 223 M
# pairs, and C57BL/6 male adult liver per the ENA sample `strain` field -- matching the Ribo-seq strain.
# Verified on arrival: 90.6-92.4% salmon mapping, 81.1% STAR unique, top single transcript 4.3-5.0% of
# TPM (Janich pre-decontamination was 52-70%). PRJEB86747 was rejected as an arm: C3H/HeNRj strain plus
# a mutation-accumulation design, so "better RNA" would have been inseparable from "different mouse".
#
# arm b uses --ref-pack, so tx_order and lengths are copied verbatim from arm a and arm a's
# orf_track_v2_nokozak.npy applies unchanged -- no second track build, and no risk of a track that
# disagrees by a row. arm c has its own universe and therefore needs its own track (--kozak none per
# feedback_kozak_never_set).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
S=$NEW/proteogenomics/scripts
REF=$NEW/data/heldout_refs
ARM_A=$NEW/data/packed_heldout_mouse_janich_liver_decon
# NOTE the packed_heldout_ prefix is REQUIRED: dataset.heldout_pack_dir() resolves
# --heldout <name> to data/packed_heldout_<name> with no env override.
ARM_B=$NEW/data/packed_heldout_rnaq_b_uniJanich_cov34766
ARM_C=$NEW/data/packed_heldout_rnaq_c_uni34766_cov34766
COV34766=$NEW/data/heldout_rnaseq_coverage/mouse_prjeb34766_liver/PRJEB34766_coverage.hd5
Q34766=$NEW/data/tpm/salmon_quant_prjeb34766_liver
echo "host=$(hostname) start=$(date -Iseconds)"

[[ -s "$COV34766" ]] || { echo "no PRJEB34766 coverage at $COV34766" >&2; exit 1; }
[[ -s "$ARM_A/provenance.json" ]] || { echo "arm a pack missing: $ARM_A" >&2; exit 1; }

# Same 5 Ribo P-site inputs as arm a, read from arm a's provenance so the arms cannot drift.
mapfile -t PS < <("$PY" -c "
import json;print('\n'.join(d['path'] for d in json.load(open('$ARM_A/provenance.json'))['ribo_psites_inputs']))")
[[ ${#PS[@]} -ge 1 ]] || { echo 'could not read ribo_psites_inputs from arm a' >&2; exit 1; }
for p in "${PS[@]}"; do [[ -s "$p" ]] || { echo "missing P-site input $p" >&2; exit 1; }; done
echo "Ribo held fixed: ${#PS[@]} P-site inputs from arm a's provenance"

# ---- arm b: universe HELD (--ref-pack arm a), coverage swapped ----
if [[ -s "$ARM_B/coverage.npy" ]]; then
  echo "skip arm b (exists)"
else
  "$PY" "$NEW/scripts/prepare/prepare_pack.py" \
      --ribo-psites "${PS[@]}" --rna-coverage "$COV34766" \
      --ref-pack "$ARM_A" --tx2biotype "$NEW/data/tx2biotype_mouse.tsv" \
      --group rnaq_b_uniJanich_cov34766 --species mouse --out "$ARM_B" 2>&1 | tail -8
fi

# The whole point of --ref-pack is that arm a's ORF track stays valid. Assert it rather than assume:
# a row-order difference here would silently mis-align every ORF the caller sees.
"$PY" - <<PYEOF
import numpy as np, sys
a=[l.strip() for l in open("$ARM_A/tx_order.txt")]
b=[l.strip() for l in open("$ARM_B/tx_order.txt")]
assert a==b, f"arm b tx_order differs from arm a ({len(a)} vs {len(b)}) -- ORF track is NOT reusable"
la,lb=np.load("$ARM_A/lengths.npy"),np.load("$ARM_B/lengths.npy")
assert np.array_equal(la,lb), "arm b lengths differ from arm a"
ta,tb=np.load("$ARM_A/target_counts.npy"),np.load("$ARM_B/target_counts.npy")
assert np.array_equal(ta,tb), "arm b Ribo target differs from arm a -- the Ribo arm did NOT stay fixed"
ca,cb=np.load("$ARM_A/coverage.npy"),np.load("$ARM_B/coverage.npy")
assert not np.array_equal(ca,cb), "arm b coverage is IDENTICAL to arm a -- the RNA swap did not happen"
print(f"arm b OK: {len(a):,} tx, target byte-identical to arm a, coverage differs "
      f"(mean {ca.mean():.2f} -> {cb.mean():.2f})")
PYEOF
ln -sf "$ARM_A/orf_track_v2_nokozak.npy" "$ARM_B/orf_track_v2_nokozak.npy"

# ---- arm c: universe from PRJEB34766 salmon, coverage PRJEB34766 ----
if [[ ! -s "$REF/mouse_prjeb34766_universe_tx.txt" ]]; then
  mapfile -t QUANTS < <(ls "$Q34766"/*/quant.sf | sort)
  [[ ${#QUANTS[@]} -ge 5 ]] || { echo "expected >=5 PRJEB34766 quants, got ${#QUANTS[@]}" >&2; exit 1; }
  echo "building arm c universe from ${#QUANTS[@]} quants (mean TPM merge)"
  export PYTHONPATH="$S:${PYTHONPATH:-}"
  cd "$S" && "$PY" -m pgx.universe --salmon "${QUANTS[@]}" --species mouse \
      --out-prefix "$REF/mouse_prjeb34766" 2>&1 | tail -5
  cd "$NEW"
fi
NTX_C=$(wc -l < "$REF/mouse_prjeb34766_universe_tx.txt")
echo "arm c universe: $NTX_C tx  (arm a/b: $(wc -l < "$ARM_A/tx_order.txt"))"
[[ $NTX_C -gt 6000 ]] || { echo "ABORT: arm c universe $NTX_C looks crushed" >&2; exit 1; }

if [[ -s "$ARM_C/coverage.npy" ]]; then
  echo "skip arm c (exists)"
else
  "$PY" "$NEW/scripts/prepare/prepare_pack.py" \
      --ribo-psites "${PS[@]}" --rna-coverage "$COV34766" \
      --universe-tx "$REF/mouse_prjeb34766_universe_tx.txt" \
      --build-orf-track --fasta "$REF/mouse_prjeb34766_universe.fa" --kozak none \
      --tx2biotype "$NEW/data/tx2biotype_mouse.tsv" \
      --group rnaq_c_uni34766_cov34766 --species mouse --out "$ARM_C" 2>&1 | tail -8
fi

"$PY" - <<PYEOF
import numpy as np
for tag,d in [("a","$ARM_A"),("b","$ARM_B"),("c","$ARM_C")]:
    o=[l.strip() for l in open(f"{d}/tx_order.txt")]
    cov=np.load(f"{d}/coverage.npy"); tgt=np.load(f"{d}/target_counts.npy")
    print(f"arm {tag}: {len(o):,} tx  sum_L={len(cov):,}  cov_mean={cov.mean():.2f}  "
          f"psites={tgt.sum():,}")
PYEOF
echo "done=$(date -Iseconds)"
