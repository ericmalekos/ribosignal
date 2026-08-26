#!/usr/bin/env bash
#SBATCH --job-name=cov_from_bam
#SBATCH --partition=long
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/cov_bam_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/cov_bam_%j.err
set -euo pipefail
# FALLBACK: per-nt coverage from an ALREADY-ALIGNED transcriptome BAM. Skips STAR entirely.
#
# Why this exists: align_line_multirun.sbatch removes its scratch only AFTER a successful
# `mv "$TBAM" "$BAMDIR/"`, so a wall-clock kill during the coverage step leaves the BAM intact in
# data/_scratch_align/<LINE>_mm<MM>_<jobid>/. Re-running STAR would waste ~9 h; this picks up from
# the BAM. On `long` (14 d) rather than `medium` (12 h) precisely because the coverage step is the
# slow tail that caused the timeout.
#
# Only works for jobs whose scratch was on the GROUP FS. Jobs submitted before that change wrote to
# node-local /data/tmp, which is unreachable from the head node and may be cleared on job death --
# those need a full STAR re-run.
#
# Env: LINE, MM, SCRDIR (the surviving scratch dir).
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
LINE=${LINE:?set LINE}; MM=${MM:?set MM=1 or 10}; SCRDIR=${SCRDIR:?set SCRDIR=<surviving scratch dir>}
P=$NEW/proteogenomics/data/${LINE}_pilot
COVDIR=$P/coverage_mm${MM}_full; BAMDIR=$P/bam_mm${MM}_full
TBAM=$SCRDIR/${LINE}.Aligned.toTranscriptome.out.bam
mkdir -p "$COVDIR" "$BAMDIR"
echo "host=$(hostname) LINE=$LINE MM=$MM start=$(date -Iseconds)"
[ -s "$TBAM" ] || { echo "no BAM at $TBAM" >&2; exit 1; }
echo "  BAM $(du -h "$TBAM" | cut -f1)"
# a BAM truncated by the job's death keeps a VALID HEADER -- quickcheck verifies the EOF block
"$ENV/bin/samtools" quickcheck "$TBAM" || { echo "TRUNCATED $TBAM -- STAR must be re-run" >&2; exit 2; }
[ -s "$SCRDIR/${LINE}.Log.final.out" ] && cp "$SCRDIR/${LINE}.Log.final.out" "$COVDIR/${LINE}.Log.final.out"
grep -E "Number of input reads|Uniquely mapped reads %|% of reads mapped to multiple loci" \
  "$COVDIR/${LINE}.Log.final.out" 2>/dev/null | sed 's/^/  /'
"$ENV/bin/python" "$NEW/scripts/rnaseq_coverage.py" "$TBAM" "$COVDIR/${LINE}_coverage.hd5" "$LINE"
mv "$TBAM" "$BAMDIR/" && rm -rf "$SCRDIR"
echo "done=$(date -Iseconds) -> $COVDIR/${LINE}_coverage.hd5"
