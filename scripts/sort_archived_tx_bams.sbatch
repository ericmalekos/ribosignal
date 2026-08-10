#!/usr/bin/env bash
#SBATCH --job-name=sort_tx_bams
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=8:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/sort_tx_bams_%j.out
set -euo pipefail
# Coordinate-sort + index the archived STAR transcriptome BAMs IN PLACE.
#
# WHY THEY ARE UNSORTED. `--outSAMtype BAM SortedByCoordinate` sorts the GENOME BAM
# (Aligned.sortedByCoord.out.bam). The transcriptome BAM produced by `--quantMode TranscriptomeSAM`
# is a SEPARATE output that STAR always emits in READ order, with no @HD SO: line at all -- its
# intended consumers (RSEM, salmon) want reads grouped, not coordinate-ordered. Anything that needs
# random access (samtools index, RiboCode metaplots via pysam fetch) therefore fails on them with
# "cannot be indexed". The B721 BAMs are sorted only because filter_tx_heldout.py sorts as a side
# effect; BAMs archived straight from STAR are not.
#
# Verified by header scan (2026-08-07): gse243134_liver 0/40 sorted, janich_liver 0/7,
# gse155087_tcell 5/8, gse120762_nt 2/5, gse120762_lps 2/4.
#
# SAFETY. The archive is the deliverable (feedback_keep_alignment_bams), so this never overwrites in
# a way that can lose data:
#   1. sort to NODE-LOCAL scratch (this job is the only reader; group fs is at ~93% of a 15 TB quota)
#   2. verify the sorted copy: header says SO:coordinate, samtools quickcheck passes, and the record
#      count MATCHES the original
#   3. only then move it over the original and index
# A failure at any step leaves the original untouched.
#
#   sbatch scripts/sort_archived_tx_bams.sbatch            # all archived dirs
#   DIRS="mouse_gse243134_liver" sbatch scripts/...        # one dir
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
S=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/samtools
BAMROOT=$NEW/data/heldout_bam
T=${SLURM_CPUS_PER_TASK:-8}
DIRS=${DIRS:-"mouse_gse243134_liver mouse_janich_liver mouse_gse155087_tcell mouse_gse120762_nt mouse_gse120762_lps"}
SCR=/data/tmp/emalekos/sort_tx_${SLURM_JOB_ID}; mkdir -p "$SCR"
trap 'rm -rf "$SCR"' EXIT
echo "host=$(hostname) start=$(date -Iseconds)"

done_n=0; skip_n=0; fail_n=0
for d in $DIRS; do
  [[ -d "$BAMROOT/$d" ]] || { echo "no dir $BAMROOT/$d"; continue; }
  echo "=== $d ==="
  for b in "$BAMROOT/$d"/*.Aligned.toTranscriptome.out.bam; do
    [[ -s "$b" ]] || continue
    nm=$(basename "$b")
    # NOTE: `view -H | head -1 | grep -q` is WRONG under `set -o pipefail`. head/grep exit early,
    # SIGPIPE the upstream samtools, pipefail propagates that, and the test reads FALSE even when the
    # tag is present. The first run of this script reported "failed=64, sorted=0, already=0" entirely
    # from that artifact. Capture the line with sed (which drains stdin) instead.
    if [[ "$("$S" view -H "$b" 2>/dev/null | sed -n 1p)" == *SO:coordinate* ]]; then
      [[ -s "$b.bai" ]] || "$S" index -@ "$T" "$b"
      echo "  SKIP $nm (already coordinate-sorted)"; skip_n=$((skip_n+1)); continue
    fi
    tmp=$SCR/$nm
    n_in=$("$S" view -c -@ "$T" "$b")
    "$S" sort -@ "$T" -m 2G -T "$SCR/st_$nm" -o "$tmp" "$b"
    # gate on CONTENT, not on the sort exiting 0
    "$S" quickcheck "$tmp" 2>/dev/null || { echo "  FAIL $nm (quickcheck)"; fail_n=$((fail_n+1)); rm -f "$tmp"; continue; }
    [[ "$("$S" view -H "$tmp" | sed -n 1p)" == *SO:coordinate* ]] || { echo "  FAIL $nm (no SO tag)"; fail_n=$((fail_n+1)); rm -f "$tmp"; continue; }
    n_out=$("$S" view -c -@ "$T" "$tmp")
    [[ "$n_in" == "$n_out" ]] || { echo "  FAIL $nm (records $n_in -> $n_out)"; fail_n=$((fail_n+1)); rm -f "$tmp"; continue; }
    mv -f "$tmp" "$b"
    "$S" index -@ "$T" "$b"
    echo "  SORTED $nm ($n_in records, $(du -h "$b" | cut -f1))"
    done_n=$((done_n+1))
  done
done
echo "done=$(date -Iseconds)  sorted=$done_n  already=$skip_n  failed=$fail_n"
[[ $fail_n -eq 0 ]] || exit 1
