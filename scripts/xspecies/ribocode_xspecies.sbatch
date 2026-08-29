#!/usr/bin/env bash
#SBATCH --job-name=xsp_ribocode
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/xspecies/ribocode_%A.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/xspecies/ribocode_%A.err
#
# Pooled RiboCode P-site build for one cross-species Ribo arm:
#   metaplots (per-read-length P-site offsets, aggregated over the arm's BAMs)
#   -> RiboCode -l no -g  (joint ORF call + per-nt <run>_psites.hd5)
#
#   sbatch --export=ALL,SHEET=<ribo_<sp>.tsv> scripts/xspecies/ribocode_xspecies.sbatch
#
# This mirrors scripts/heldout/ribocode_heldout.sh but resolves its paths from
# scripts/xspecies/species.tsv and the samplesheet, rather than from
# scripts/heldout/heldout_config.sh. That config hardcodes
# BAM_DIR="$NEW/data/heldout_bam/$DATASET" and carries a case block per dataset; adding
# nine arms to it would mean editing a file five other datasets already depend on, for
# no gain. The RiboCode invocation itself is identical.
set -euo pipefail

NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
RIBOCODE_BIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
SHEET=${SHEET:?set SHEET=<ribo samplesheet.tsv>}
BAMROOT=${BAMROOT:-$NEW/data/xspecies_ribo_bam_mm1}
# P-site output root. MUST move whenever BAMROOT moves: this script writes
# <root>/<dataset>/ and the mm1 P-sites, the packs and every ORF call set already live
# under the default. Pointing a second posture at the same root would silently overwrite
# the call sets the model was scored against.
PSITEROOT=${PSITEROOT:-$NEW/data/xspecies_psites}

DS=$(awk -F'\t' 'NR==2{print $1}' "$SHEET")
SPECIES=$(awk -F'\t' 'NR==2{print $5}' "$SHEET")
[[ -n "$DS" && -n "$SPECIES" ]] || { echo "could not read dataset/species from $SHEET" >&2; exit 1; }

ANNOT=$NEW/data/ribocode_annot_$SPECIES
[[ -s "$ANNOT/transcripts_sequence.fa" ]] || { echo "missing RiboCode annot $ANNOT" >&2; exit 1; }

BAM_DIR=$BAMROOT/$DS
PSITE_DIR=$PSITEROOT/$DS
mkdir -p "$PSITE_DIR"

echo "dataset=$DS species=$SPECIES host=$(hostname) start=$(date -Iseconds)"
echo "  annot=$ANNOT"
echo "  bams =$BAM_DIR"
echo "  out  =$PSITE_DIR"

BAM_LIST="$PSITE_DIR/${DS}_bam_list.txt"; : > "$BAM_LIST"
n=0
while IFS=$'\t' read -r _ RUN _rest; do
  bam="$BAM_DIR/${RUN}.Aligned.toTranscriptome.out.bam"
  [[ -f "$bam" ]] || { echo "MISSING bam: $bam" >&2; exit 2; }
  # quickcheck every BAM we did not just write: one truncated mid-write keeps a VALID
  # header, and RiboCode then dies 90 minutes in on an opaque pysam error.
  "$ENV/bin/samtools" quickcheck "$bam" || { echo "ABORT: $bam failed quickcheck" >&2; exit 2; }
  [[ -f "$bam.bai" ]] || "$ENV/bin/samtools" index "$bam"
  echo "$bam" >> "$BAM_LIST"; n=$((n+1))
done < <(awk -F'\t' 'NR>1' "$SHEET")
echo "  $n pooled BAMs, all passing quickcheck"

cd "$PSITE_DIR"
echo "[$(date +%T)] metaplots"
"$RIBOCODE_BIN/metaplots" -a "$ANNOT" -i "$BAM_LIST" -o "$DS"
CONFIG="$PSITE_DIR/${DS}_pre_config.txt"
[[ -f "$CONFIG" ]] || { echo "metaplots produced no config: $CONFIG" >&2; exit 4; }
echo "  pre_config rows: $(( $(wc -l < "$CONFIG") - 1 ))"
# The config lists the read lengths whose 3-nt periodicity was good enough to assign a
# P-site offset. An empty one means no periodicity was found at all, which for a Ribo
# library means the trimming or the annotation is wrong -- fail here, not downstream.
[[ $(( $(wc -l < "$CONFIG") - 1 )) -ge 1 ]] || { echo "ABORT: no read length passed the periodicity check" >&2; exit 5; }
sed -n '1,12p' "$CONFIG" | sed 's/^/    /'

echo "[$(date +%T)] RiboCode"
"$RIBOCODE_BIN/RiboCode" -a "$ANNOT" -c "$CONFIG" -l no -g -o "$DS"

nh=$(ls "$PSITE_DIR"/*_psites.hd5 2>/dev/null | wc -l)
nc=$(( $(wc -l < "$PSITE_DIR/${DS}_collapsed.txt" 2>/dev/null || echo 1) - 1 ))
echo "  psites hd5: $nh    collapsed ORF calls: $nc"
[[ "$nh" -ge 1 ]] || { echo "ABORT: no psites hd5 written" >&2; exit 6; }
du -sh "$PSITE_DIR" | sed 's/^/  /'
echo "=== $DS DONE=$(date -Iseconds)"
