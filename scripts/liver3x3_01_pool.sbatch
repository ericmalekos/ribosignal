#!/usr/bin/env bash
#SBATCH --job-name=l3x3_pool
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/l3x3_pool_%A_%a.out
set -euo pipefail
# STAGE 1 of the mouse-liver 3x3 factorial: pool each dataset from BAMs onto the SHARED UNION
# universe (22,974 transcripts).
#
# ALL THREE DATASETS GO THROUGH THIS ONE PATH. That is the entire point. The step-0 gate showed that
# pooling one dataset differently from the others shifts P-site totals ~1.7% -- the same order as the
# effects the factorial measures -- so Wang is no longer sliced from its legacy full-annotation pack.
# Its Ribo runs were re-fetched and its RNA runs aligned specifically so it could join this path.
#
# WHAT THIS PRODUCES. Each task writes a diagonal pack (Ribo_X + RNA_X) and, more importantly,
# reusable per-sample hd5 under <out>/_work/{psites,coverage}/. Stage 3 assembles all NINE cells from
# those cached hd5, so the three expensive BAM passes happen once and every cell is assembled by an
# identical code path.
#
# PER-DATASET metaplots OFFSETS ARE KEPT, not pinned to a shared read-length set. That is the standard
# RiboCode workflow and the correct choice here: the three datasets have genuinely different footprint
# length distributions, and forcing one config across them would impose a different bias than the one
# it removed. The consequence -- rebuilt Janich differs ~1.7% from its historical pack -- is recorded
# in docs/mouse_liver_3x3_plan.md and the numbers are not mixed with pre-2026-08-10 results.
#
# RNA DEPTH IS NOT MATCHED across datasets (Janich 7 samples, GSE243134 19, Wang 2). Deliberate: the
# user's call was that unmatched depth is closer to real application. A matched-depth arm is tracked
# separately as task #75.
#
#   sbatch --array=0-2 scripts/liver3x3_01_pool.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
# rnaseq_coverage.py imports pysam, which cas12a does NOT have. prepare_from_bams defaults
# --pysam-python to sys.executable, so the coverage stage died on every dataset after the Ribo stage
# had already completed. riboseq carries pysam 0.23.3.
PYSAM=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/python3
ANNOT=$NEW/data/mouse_ribocode_annot
RBIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
REF=$NEW/data/heldout_refs
UNI=$REF/mouse_liver_union3_tx.txt
# RIBO comes from the CANONICAL aligner (EndToEnd + mm1 + filter on query-grouped bams,
# scripts/riboseq_align.sbatch). The old data/heldout_bam Ribo alignments lacked EndToEnd and the
# cross-gene filter, which together account for ~9.5% of records; the final recipe reproduces
# the reference pack to 3.04% (r=0.957) where the off-recipe one managed 11.69%.
# RNA is unchanged: RNA-seq stays Local, so those alignments remain valid and are not rebuilt.
E2E=$NEW/data/ribo_bam_e2e
BAM=$NEW/data/heldout_bam
L=$NEW/data/liver3x3

DS=(janich gse243134 wang)
D=${DS[$SLURM_ARRAY_TASK_ID]}
OUT=$L/pool_${D}_canon
mkdir -p "$OUT"

case "$D" in
  janich)
    mapfile -t RIBO < <(ls $E2E/janich/*.Aligned.toTranscriptome.out.bam | sort)
    # Match *.toTranscriptome.out.bam explicitly, not `*.bam | grep -v bai`. mouse_janich_liver is
    # in fact clean (7 transcriptome BAMs, no genome-coordinate ones, no Ribo runs -- verified), so
    # this selects the identical 7 files and changes no result. It is tightened because the same
    # idiom on data/heldout_bam/mouse_gse120762_* WOULD have pulled in coordinate-space BAMs and the
    # old Ribo alignments; a glob that is correct only by virtue of what happens to sit in the
    # directory is not a guarantee.
    mapfile -t RNA  < <(ls $BAM/mouse_janich_liver/*.toTranscriptome.out.bam | sort) ;;
  gse243134)
    mapfile -t RIBO < <(while read s; do echo "$E2E/gse243134/${s}.Aligned.toTranscriptome.out.bam"; done < $L/gse243134_ribo_runs.txt)
    mapfile -t RNA  < <(while read s; do echo "$BAM/mouse_gse243134_liver/${s}.Aligned.toTranscriptome.out.bam"; done < $L/gse243134_rna_runs.txt) ;;
  wang)
    mapfile -t RIBO < <(ls $E2E/wang/*.Aligned.toTranscriptome.out.bam | sort)
    mapfile -t RNA  < <(ls $BAM/mouse_wang_liver_rna/*.Aligned.toTranscriptome.out.bam | sort) ;;
esac

for f in "${RIBO[@]}" "${RNA[@]}"; do
  [[ -s "$f" ]] || { echo "missing BAM: $f" >&2; exit 1; }
done
echo "dataset=$D ribo=${#RIBO[@]} rna=${#RNA[@]} universe=$(wc -l < $UNI) host=$(hostname) start=$(date -Iseconds)"
printf '  ribo: %s\n' "${RIBO[@]##*/}"
printf '  rna:  %s\n' "${RNA[@]##*/}"

# RESUMABLE. The Ribo stage (metaplots + RiboCode, ~45 min for GSE243134's 21 BAMs) writes its hd5
# before the coverage stage runs, so a coverage failure should not cost the Ribo work twice. If a
# complete set of cached psites hd5 is already present, feed those instead of the BAMs.
CACHED=("$OUT"/_work/psites/*_psites.hd5)
if [[ -e "${CACHED[0]}" && "${#CACHED[@]}" -eq "${#RIBO[@]}" ]]; then
  echo "  reusing ${#CACHED[@]} cached psites hd5 (skipping the Ribo stage)"
  RIBO_ARGS=(--ribo-psites "${CACHED[@]}")
else
  RIBO_ARGS=(--ribo-bam "${RIBO[@]}")
fi

# --no-ncrna-filter IS CORRECT HERE, and is a declared deviation rather than a skip. The canonical
# aligner (scripts/riboseq_align.sbatch) already applied filter_tx_heldout to the RAW query-grouped
# STAR output, which is the only point at which it CAN be applied -- it needs query grouping, and
# these bams are now coordinate-sorted, so re-running it would (a) double-filter and (b) be refused
# by the sorted-input guard. The flag records in provenance that prepare_from_bams did not filter,
# because the filtering happened upstream.
#
# REUSE the RNA coverage from the pre-canonical pool. RNA-seq alignment is unchanged (it stays
# Local; only Ribo needed EndToEnd), and re-deriving it means ~2.5 h per sample on a paired-end
# transcriptome bam -- 1.48 BILLION records for an 82 M-read library. The universe is identical, so
# the hd5 are directly valid. Falls back to --rna-bam if the cache is absent or incomplete.
mapfile -t RCOV < <(ls $L/pool_${D}/_work/coverage/*_coverage.hd5 2>/dev/null | sort)
if [[ ${#RCOV[@]} -eq ${#RNA[@]} ]]; then
  echo "  reusing ${#RCOV[@]} RNA coverage hd5 from pool_${D} (RNA alignment unchanged)"
  RNA_ARGS=(--rna-coverage "${RCOV[@]}")
else
  echo "  coverage cache incomplete (${#RCOV[@]}/${#RNA[@]}); deriving from BAMs"
  RNA_ARGS=(--rna-bam "${RNA[@]}")
fi

"$PY" "$NEW/scripts/prepare/prepare_from_bams.py" \
  "${RIBO_ARGS[@]}" \
  "${RNA_ARGS[@]}" \
  --pysam-python "$PYSAM" \
  --universe-tx "$UNI" \
  --annot "$ANNOT" --tx2biotype "$NEW/data/tx2biotype_mouse.tsv" \
  --ribocode-bin "$RBIN" \
  --no-ncrna-filter \
  --group "l3x3_$D" --species mouse --out "$OUT" --threads "${SLURM_CPUS_PER_TASK:-12}"

# Gate on the cached hd5 stage 3 depends on, not on exit status alone.
# Count with a nullglob array, NOT `$(ls ... | wc -l)`. Under `set -euo pipefail` a failing `ls` on a
# missing directory propagates through the pipe and kills the script with no message -- which is
# exactly what happened once coverage was REUSED and $OUT/_work/coverage/ therefore never existed.
# The packs had already been built correctly; only this check failed. Same failure mode as an
# unguarded `x=$(md5sum missing)` earlier in this project.
shopt -s nullglob
_ps=("$OUT"/_work/psites/*_psites.hd5); NP=${#_ps[@]}
_cv=("$OUT"/_work/coverage/*_coverage.hd5); NC=${#_cv[@]}
shopt -u nullglob
echo "  cached hd5: psites=$NP (want ${#RIBO[@]})  coverage=$NC (want ${#RNA[@]})"
[[ "$NP" -eq "${#RIBO[@]}" ]] || { echo "ABORT: $NP psites hd5 for ${#RIBO[@]} Ribo BAMs" >&2; exit 1; }
[[ "$NC" -eq "${#RNA[@]}" || ${#RCOV[@]} -eq ${#RNA[@]} ]] || { echo "ABORT: $NC coverage hd5 for ${#RNA[@]} RNA BAMs" >&2; exit 1; }
"$PY" -c "
import numpy as np
t=np.load('$OUT/target_counts.npy'); c=np.load('$OUT/coverage.npy')
n=sum(1 for _ in open('$OUT/tx_order.txt'))
print(f'  pack: {n:,} tx, {t.sum():,} P-sites, {c.sum():,} coverage over {len(t):,} nt')
assert n==$(wc -l < $UNI), 'pack tx != union tx'
"
echo "done=$(date -Iseconds) -> $OUT"
