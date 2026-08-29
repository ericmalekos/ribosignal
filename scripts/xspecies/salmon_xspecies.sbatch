#!/usr/bin/env bash
#SBATCH --job-name=xsp_salmon
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/xspecies/salmon_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/xspecies/salmon_%A_%a.err
#
# salmon quant for the cross-species RNA arms, feeding the per-species expressed-transcript
# universe that Phase 5 builds on.
#
#   sbatch --array=0-N --export=ALL,SHEET=<rna_<sp>.tsv> scripts/xspecies/salmon_xspecies.sbatch
#
# TWO THINGS THIS ADDS OVER scripts/salmon_chothani_decoy.sbatch:
#   * a SINGLE-END branch. That script is paired-end only (`-1`/`-2` with no `-r`), and
#     worm, fly and yeast RNA in this expansion are all single-end.
#   * per-species index resolution through scripts/xspecies/species.tsv, so the
#     annotation version cannot drift from the one the STAR index was built against.
#
# DECOY-AWARE INDEX ONLY. A transcriptome-only index misattributes intronic and repeat
# reads to lncRNAs and biases TPM, which then biases the universe this feeds.
# `-l A` lets salmon infer strandedness rather than asserting it per dataset.
set -euo pipefail

NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
SHEET=${SHEET:?set SHEET=<rna samplesheet.tsv>}
T=${SLURM_CPUS_PER_TASK:-12}

mapfile -t ROWS < <(awk -F'\t' 'NR>1' "$SHEET")
IFS=$'\t' read -r DS RUN R1 R2 ADAPTER SPECIES LAYOUT <<< "${ROWS[$SLURM_ARRAY_TASK_ID]}"
[[ "$R2" == "-" ]] && R2=""
[[ "$LAYOUT" == "PAIRED" || "$LAYOUT" == "SINGLE" ]] || {
  echo "samplesheet row parsed wrong: layout='$LAYOUT' species='$SPECIES'" >&2; exit 1; }

SPTSV=$NEW/scripts/xspecies/species.tsv
ANNVER=$(awk -F'\t' -v s="$SPECIES" '!/^#/ && $1==s {print $2; exit}' "$SPTSV")
[[ -n "$ANNVER" ]] || { echo "species '$SPECIES' has no row in $SPTSV" >&2; exit 1; }
IDX=/private/groups/carpenterlab/emalekos/Salmon_indexes/salmon_index_decoy_${SPECIES}_${ANNVER}
[[ -s "$IDX/info.json" ]] || { echo "missing salmon index $IDX" >&2; exit 1; }

OUT=$NEW/data/xspecies_tpm/$DS/$RUN
[[ -s "$OUT/quant.sf" ]] && { echo "SKIP $RUN (quant.sf exists)"; exit 0; }
mkdir -p "$OUT"

echo "dataset=$DS run=$RUN species=$SPECIES layout=$LAYOUT host=$(hostname) start=$(date -Iseconds)"
if [[ "$LAYOUT" == "PAIRED" ]]; then
  "$ENV/bin/salmon" quant -i "$IDX" -l A -1 "$R1" -2 "$R2" -o "$OUT" -p "$T"
else
  "$ENV/bin/salmon" quant -i "$IDX" -l A -r "$R1" -o "$OUT" -p "$T"
fi

[[ -s "$OUT/quant.sf" ]] || { echo "ABORT: no quant.sf written" >&2; exit 3; }
"$ENV/bin/python" - "$OUT" <<'PY'
import json, sys, csv
from pathlib import Path
out = Path(sys.argv[1])
meta = json.loads((out / "aux_info" / "meta_info.json").read_text())
n = sum(1 for _ in (out / "quant.sf").open()) - 1
expressed = 0
with (out / "quant.sf").open() as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if float(r["TPM"]) >= 1:
            expressed += 1
print(f"  mapped={meta.get('percent_mapped'):.1f}%  libtype={meta.get('library_types')}  "
      f"tx={n:,}  TPM>=1: {expressed:,}")
# A collapse in percent_mapped is the signal that the index or the library is wrong.
if meta.get("percent_mapped", 0) < 20:
    raise SystemExit(f"ABORT: only {meta['percent_mapped']:.1f}% mapped")
PY
echo "=== $RUN DONE=$(date -Iseconds)"
