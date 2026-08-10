#!/usr/bin/env bash
#SBATCH --job-name=hs_annot
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/hs_annot_%j.out
set -euo pipefail
# Build the HUMAN RiboCode prepared annotation (GENCODE v49 / GRCh38).
#
# WHY THIS EXISTS. There is no human RiboCode annotation anywhere under RNAZoo_meta -- the tree only
# has the mouse one. RiboCode's `metaplots` and `RiboCode` both refuse to run without a prepared
# annotation dir ("pleas run prepare_transcripts.py first"), so the THP-1 (GSE208041, GSE39561) and
# CAR-T (GSE304796) ORF calls cannot proceed until this lands. Mouse work is unaffected.
#
# The mouse annotation lives at data/mouse_ribocode_annot/ (it was relocated out of
# expression_context_human at some point; 17 scripts still pointed at the old path and were
# repointed). This one is written alongside it, in-project, under the matching name -- NOT back into
# another experiment's data dir, which is how the mouse one got lost track of in the first place.
#
# Annotation choice: GENCODE v49 primary assembly, matching the STAR index the human Ribo-seq BAMs
# were aligned against (STAR_indexes/star_index_grch38_v49). RiboCode joins on versioned transcript
# IDs, so a release mismatch here silently drops transcripts rather than erroring.
#
#   sbatch scripts/build_human_ribocode_annot.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
RBIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
GTF=/private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations/gencode.v49.annotation.gtf
FA=/private/groups/carpenterlab/emalekos/genomes/GRCh38.primary_assembly.genome.fa
OUT=$NEW/data/human_ribocode_annot

for f in "$GTF" "$FA" "$RBIN/prepare_transcripts"; do
  [[ -s "$f" ]] || { echo "missing $f" >&2; exit 1; }
done
if [[ -s "$OUT/transcripts.pickle" && -s "$OUT/transcripts_cds.txt" ]]; then
  echo "SKIP: $OUT already prepared"; exit 0
fi
mkdir -p "$OUT"
echo "host=$(hostname) start=$(date -Iseconds)"
echo "gtf=$GTF"
echo "fa=$FA"

"$RBIN/prepare_transcripts" -g "$GTF" -f "$FA" -o "$OUT"

# Gate on the three files RiboCode actually loads, not on exit status: a partial prepare leaves a
# directory that looks populated and then fails inside metaplots hours later, inside another job.
for f in transcripts_cds.txt transcripts_sequence.fa transcripts.pickle; do
  [[ -s "$OUT/$f" ]] || { echo "ABORT: $OUT/$f missing or empty after prepare" >&2; exit 1; }
done
NCDS=$(wc -l < "$OUT/transcripts_cds.txt")
echo "  transcripts_cds.txt: $NCDS coding transcripts"
[[ "$NCDS" -gt 50000 ]] || { echo "ABORT: only $NCDS CDS entries -- expected >50k for GENCODE v49" >&2; exit 1; }
head -1 "$OUT/transcripts_cds.txt" | grep -q '^ENST' || { echo "ABORT: cds table is not ENST-keyed" >&2; exit 1; }

cat > "$OUT/PROVENANCE.md" <<EOF
# human_ribocode_annot

RiboCode prepared annotation, built $(date -Iseconds) by scripts/build_human_ribocode_annot.sbatch.

- GTF:    $GTF
- genome: $FA
- tool:   $RBIN/prepare_transcripts
- coding transcripts: $NCDS

Matches STAR_indexes/star_index_grch38_v49, which is what the human Ribo-seq BAMs
(THP-1 GSE208041, THP-1 GSE39561, CAR-T GSE304796) are aligned against. RiboCode joins on
VERSIONED transcript IDs, so do not pair this annotation with BAMs from a different GENCODE release.
EOF
echo "done=$(date -Iseconds) -> $OUT"
