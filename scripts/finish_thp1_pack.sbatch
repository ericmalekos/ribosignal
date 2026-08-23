#!/usr/bin/env bash
#SBATCH --job-name=thp1_finish
#SBATCH --partition=long
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=1-00:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/thp1_finish_%j.out
set -euo pipefail
# RESUME the THP-1 (GSE208041) final-recipe pack. Job 36786310_1 hit its 11:45 wall on the LAST of four
# RNA coverage files, with 5/5 psites and 3/4 coverage already written -- about 95% done.
#
# WHY NOT JUST RESUBMIT. rebuild_heldout_canon.sbatch has NO resume guards: it would regenerate all
# five psites and all four coverage files from scratch, ~11 h of work already on disk, and would very
# likely hit the wall again. `prepare_from_bams.py` explicitly supports passing pre-made hd5 to skip a
# stage, so this finishes the one missing file and assembles the pack directly.
#
# WHY `long` AND 24 h. The RNA BAMs here are ~1.5 BILLION records each and each coverage pass took
# 2-3 h. medium's 12 h ceiling is what killed the first attempt; the marginal cost of a longer wall
# on `long` is zero if the job finishes early.
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
W=$NEW/data/packed_heldout_human_gse208041_canon/_work
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
RPY=/private/groups/carpenterlab/emalekos/conda_envs/riboseq/bin/python3
RNADIR=$NEW/data/processed/GSE208041_thp1
RUNS=$NEW/data/heldout_refs/human_gse208041_rna_runs.txt

echo "host=$(hostname) start=$(date -Iseconds)"
echo "  psites present:   $(ls "$W"/psites/*_psites.hd5 2>/dev/null | wc -l) (expect 5)"
echo "  coverage present: $(ls "$W"/coverage/*_coverage.hd5 2>/dev/null | wc -l) (expect 4)"

# Generate ONLY the coverage files that are missing. A partial file from the killed job would be
# unusable, so anything not matching the expected name is regenerated rather than trusted.
while read -r s; do
  [[ -n "$s" ]] || continue
  out=$W/coverage/${s}_coverage.hd5
  if [[ -s "$out" ]]; then echo "  skip $s (present)"; continue; fi
  bam=$RNADIR/${s}.toTranscriptome.bam
  [[ -s "$bam" ]] || { echo "ABORT: missing $bam" >&2; exit 1; }
  echo "  generating coverage for $s $(date -Iseconds)"
  "$RPY" "$NEW/scripts/rnaseq_coverage.py" "$bam" "$out" "$s"
done < "$RUNS"

n=$(ls "$W"/coverage/*_coverage.hd5 2>/dev/null | wc -l)
[[ "$n" -eq 4 ]] || { echo "ABORT: $n coverage hd5, expected 4" >&2; exit 1; }

echo "--- assembling pack $(date -Iseconds) ---"
"$PY" "$NEW/scripts/prepare/prepare_pack.py" \
  --ribo-psites "$W/psites" --rna-coverage "$W/coverage" \
  --group human_gse208041_canon --species human \
  --out "$NEW/data/packed_heldout_human_gse208041_canon" \
  --universe-tx "$NEW/data/heldout_refs/human_gse208041_universe_tx.txt" \
  --tx2biotype "$NEW/data/tx2biotype.tsv"

echo "done=$(date -Iseconds)"
