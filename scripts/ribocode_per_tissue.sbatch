#!/usr/bin/env bash
#SBATCH --job-name=rc_tissue
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --array=0-8
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/ribocode_per_tissue/%x_%A_%a.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/ribocode_per_tissue/%x_%A_%a.err
#
# Regenerate the per-tissue RiboCode ORF calls (<Tissue>_collapsed.txt) that were lost with the
# biotype_probe tree on 2026-08-15. 15 files depend on them and no copy survived anywhere.
#
# Mirrors scripts/heldout/ribocode_heldout.sh, whose own header says it mirrors the LOST
# expression_context_human phase30 run_ribocode_per_tissue.sh. Same two steps: metaplots to derive
# per-read-length P-site offsets over the tissue's pooled transcriptome BAMs, then RiboCode for the
# joint ORF call + per-nt *_psites.hd5.
#
# WRITES IN-PROJECT, under data/ribocode_per_tissue/. This is the whole point: the original lived in
# ANOTHER experiment's data dir, which is exactly why it vanished (see docs/PENDING_FOLLOWUPS.md,
# task #93). Nothing this project needs to re-run may live outside riboseq_signal_model/.
#
#   sbatch scripts/ribocode_per_tissue.sbatch
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
RIBOCODE_BIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
SAMTOOLS=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/samtools
# GENCODE v49 primary, matching star_index_grch38_v49 which every Chothani BAM was aligned against.
# RiboCode joins on VERSIONED transcript ids, so a release mismatch silently drops transcripts.
ANNOT=$NEW/data/human_ribocode_annot_primary

TISSUES=(ES Fat Fibroblast HA_EC HCAEC Hepatocytes HUVEC VSMC Brain)
T=${TISSUES[${SLURM_ARRAY_TASK_ID:-0}]}
# Brain was dropped from training on low periodicity, so its BAMs were never staged into
# ribo_bam_e2e_human/; they are still in the raw chothani_bam/ tree.
if [[ "$T" == "Brain" ]]; then BAM_DIR=$NEW/data/chothani_bam/Brain
else                          BAM_DIR=$NEW/data/ribo_bam_e2e_human/$T; fi

OUT=$NEW/data/ribocode_per_tissue/$T
mkdir -p "$OUT" "$NEW/logs/ribocode_per_tissue"
export PYTHONNOUSERSITE=1

BAM_LIST=$OUT/${T}_bam_list.txt; : > "$BAM_LIST"
shopt -s nullglob
for bam in "$BAM_DIR"/*.Aligned.toTranscriptome.out.bam; do
  [[ -f "$bam.bai" ]] || "$SAMTOOLS" index "$bam"
  echo "$bam" >> "$BAM_LIST"
done
n=$(wc -l < "$BAM_LIST")
[[ "$n" -gt 0 ]] || { echo "FATAL: no BAMs under $BAM_DIR" >&2; exit 2; }
echo "[$(date +%T)] $T: $n pooled BAMs from $BAM_DIR"

cd "$OUT"
echo "[$(date +%T)] $T: metaplots"
"$RIBOCODE_BIN/metaplots" -a "$ANNOT" -i "$BAM_LIST" -o "$T"
CONFIG=$OUT/${T}_pre_config.txt
[[ -f "$CONFIG" ]] || { echo "FATAL: metaplots produced no config: $CONFIG" >&2; exit 3; }
echo "[$(date +%T)] $T: pre_config has $(wc -l < "$CONFIG") rows"

echo "[$(date +%T)] $T: RiboCode"
"$RIBOCODE_BIN/RiboCode" -a "$ANNOT" -c "$CONFIG" -l no -g -o "$T"

COLLAPSED=$OUT/${T}_collapsed.txt
[[ -s "$COLLAPSED" ]] || { echo "FATAL: no collapsed calls written" >&2; exit 4; }
echo "[$(date +%T)] $T: DONE -- $(( $(wc -l < "$COLLAPSED") - 1 )) collapsed ORF calls"
