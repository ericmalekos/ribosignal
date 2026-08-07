#!/usr/bin/env bash
# RiboCode on the REAL B721.221 Ribo-seq (Ouspenskaia 2022) -- the orthogonal validation.
#
# WHY THIS MATTERS. Every model-vs-null comparison so far scores the model against enumerated or
# coding-potential nulls; none of it checks the model's calls against actual measured translation in
# the same cells. HBL-1 was supposed to provide that and was abandoned at ~2% unique mapping
# (~98% rRNA). B721 aligned at 10.7-20.4% unique for 327 M unique footprints across 4 HLA alleles,
# so it is the first usable same-sample Ribo-seq on an immunopeptidome line in this project.
#
# The 7 runs are 4 biological libraries: A0101 and B1501 and B4402 were each sequenced twice for
# depth (per GEO), A3303 once. Pooling all 7 into one RiboCode call is therefore correct -- they are
# the same parental line, and metaplots derives a per-BAM P-site offset anyway.
#
# NOTE the alleles differ from the RNA-seq (C0401/C0701-transduced, Sarkizova 2020). Same B721.221
# parental line, different HLA transduction. Fine for translation, but it is NOT same-sample matching
# and is recorded as such in DATA_PROVENANCE.md.
#
#   sbatch proteogenomics/scripts/ribocode_b721.sbatch
#SBATCH --job-name=ribocode_b721
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=110G
#SBATCH --time=12:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/proteogenomics/logs/ribocode_b721_%j.out
set -euo pipefail
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ECH=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human
RC=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
ANNOT=$ECH/data/ribocode_annot            # v49 prepared-transcripts, same annotation as training
RB=$NEW/proteogenomics/data/B721_pilot/riboseq
OUT=$NEW/proteogenomics/data/B721_pilot/ribocode; mkdir -p "$OUT"

SRRS=(SRR10846525 SRR10846526 SRR10846527 SRR10846528 SRR10846529 SRR10846530 SRR10846531)
BAMLIST=$OUT/b721_bams.txt; : > "$BAMLIST"
for s in "${SRRS[@]}"; do
  b=$RB/$s.Aligned.toTranscriptome.out.bam
  [[ -f "$b" && -f "$b.bai" ]] || { echo "missing/unindexed $b" >&2; exit 1; }
  echo "$b" >> "$BAMLIST"
done
echo "[$(date +%T)] host=$(hostname)  $(wc -l < "$BAMLIST") BAMs pooled"

cd "$OUT"
echo "[$(date +%T)] metaplots (P-site offsets per BAM)"
"$RC/metaplots" -a "$ANNOT" -i "$BAMLIST" -o B721 > "$OUT/metaplots.log" 2>&1
CONFIG=$OUT/B721_pre_config.txt
[[ -f "$CONFIG" ]] || { echo "no pre_config produced" >&2; tail -20 "$OUT/metaplots.log"; exit 2; }
echo "[$(date +%T)] pre_config:"; cat "$CONFIG"

echo "[$(date +%T)] RiboCode joint ORF call"
"$RC/RiboCode" -a "$ANNOT" -c "$CONFIG" -l no -g -o B721 > "$OUT/ribocode.log" 2>&1
echo "[$(date +%T)] DONE"; ls -lh "$OUT" | grep -i B721 | head
wc -l "$OUT"/B721*collapsed.txt 2>/dev/null || true
