#!/usr/bin/env bash
#SBATCH --job-name=huvec_test
#SBATCH --partition=medium
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH --time=11:45:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/huvec_test_%j.out
set -euo pipefail
# Is the CURRENT pipeline aligned with the training data, or was the 12.6% my own recipe drift?
#
# The Chothani divergence run reported HUVEC at 11.69% against packed_union_HUVEC. But those BAMs
# were aligned with a recipe I wrote on 2026-08-10 that adds two stages the originals cannot have
# had: bowtie2 contaminant depletion, bowtie2 snoRNA depletion (that index was built the same day),
# and --alignEndsType EndToEnd. So 11.69% conflates "pipeline vs training data" with "my recipe vs
# the original recipe", and only the first is the question being asked.
#
# This re-aligns the same 3 HUVEC runs with the MINIMAL recipe -- the closest analogue to
# align_janich_ribo.sbatch, which reproduced its pack to 1.69%:
#     cutadapt --trim-n -m 20 -M 40   (N-padded reads; verified, no ligated adapter present)
#     STAR --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM   (no EndToEnd, no bowtie2)
# then pools with the same code path and compares to the same reference pack.
#
# Reading the result:
#   ~1-2%   -> the current pipeline IS aligned with the training data; 11.69% was my recipe drift,
#              and the depletion/EndToEnd stages account for the difference.
#   ~12%    -> a real mismatch between the repo's pipeline and the data the models were trained on,
#              independent of anything added this week. That would be the serious case.
#
#   sbatch scripts/huvec_pipeline_alignment_test.sbatch
NEW=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
ENV=/private/groups/carpenterlab/emalekos/conda_envs/riboseq
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
RBIN=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin
IDX=/private/groups/carpenterlab/emalekos/STAR_indexes/star_index_grch38_v49
ANNOT=$NEW/data/human_ribocode_annot_primary
D=$NEW/data/external/chothani_ribo_refetch
PACK=$NEW/data/packed_union_HUVEC
W=$NEW/results/chothani_regeneration/_minimal_HUVEC
BAMDIR=$W/bam
mkdir -p "$BAMDIR"
T=${SLURM_CPUS_PER_TASK:-12}
RUNS=(SRR15513213 SRR15513214 SRR15513215)
echo "host=$(hostname) start=$(date -Iseconds)"

for SRR in "${RUNS[@]}"; do
  OUTBAM=$BAMDIR/${SRR}.Aligned.toTranscriptome.out.bam
  [[ -s "$OUTBAM" ]] && { echo "  skip $SRR (bam exists)"; continue; }
  SCR=/data/tmp/emalekos/huvec_${SLURM_JOB_ID}_${SRR}
  mkdir -p "$SCR"
  "$ENV/bin/cutadapt" --trim-n -m 20 -M 40 -j "$T" -o "$SCR/t.fq.gz" "$D/${SRR}.fastq.gz" \
    > "$BAMDIR/${SRR}.cutadapt.log" 2>&1
  "$ENV/bin/STAR" --genomeDir "$IDX" --readFilesIn "$SCR/t.fq.gz" --readFilesCommand zcat \
    --runThreadN "$T" --outFilterMultimapNmax 1 \
    --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM \
    --outSAMattributes NH HI AS nM --outFileNamePrefix "$SCR/${SRR}."
  grep -E "Uniquely mapped reads %" "$SCR/${SRR}.Log.final.out" | sed "s/^/  $SRR /"
  "$ENV/bin/samtools" sort -@ "$T" -m 3G -o "$OUTBAM" "$SCR/${SRR}.Aligned.toTranscriptome.out.bam"
  "$ENV/bin/samtools" index -@ "$T" "$OUTBAM"
  rm -rf "$SCR"
done

# ---- same Ribo stage and same pooling code path as the divergence run ----
if ! ls "$W"/*_psites.hd5 >/dev/null 2>&1; then
  ls "$BAMDIR"/*.Aligned.toTranscriptome.out.bam > "$W/bam_list.txt"
  ( cd "$W" && "$RBIN/metaplots" -a "$ANNOT" -i "$W/bam_list.txt" -o min_HUVEC )
  CFG=$(ls "$W"/min_HUVEC_pre_config.txt | head -1)
  ( cd "$W" && "$RBIN/RiboCode" -a "$ANNOT" -c "$CFG" -l no -g -o min_HUVEC )
fi

"$PY" - <<PYEOF
import numpy as np, glob, json, pathlib, sys
sys.path.insert(0, "$NEW/scripts")
from prepare.packlib import pool_per_nt
pack = pathlib.Path("$PACK")
order = [l.strip() for l in open(pack / "tx_order.txt")]
lens = dict(zip(order, np.load(pack / "lengths.npy").astype(np.int64)))
hd5 = sorted(glob.glob("$W/*_psites.hd5"))
pooled = pool_per_nt(hd5, "p_sites", set(order))
new = np.concatenate([np.asarray(pooled[t], dtype=np.int64) if t in pooled
                      else np.zeros(lens[t], dtype=np.int64) for t in order])
old = np.load(pack / "target_counts.npy").astype(np.int64)
d = np.abs(new - old)
res = dict(tissue="HUVEC", recipe="minimal (no bowtie2, no EndToEnd)", n_bams=len(hd5),
           old_psites=int(old.sum()), new_psites=int(new.sum()),
           pct_of_signal=float(100 * d.sum() / max(old.sum(), 1)),
           pearson=float(np.corrcoef(new, old)[0, 1]))
pathlib.Path("$W/divergence_minimal.json").write_text(json.dumps(res, indent=2))
print(f"  MINIMAL RECIPE: old {res['old_psites']:,} | new {res['new_psites']:,} | "
      f"{res['pct_of_signal']:.2f}% of signal | per-nt r={res['pearson']:.4f}")
print("  compare: my depleted/EndToEnd recipe gave 11.69% and r=0.9350 on this same tissue")
PYEOF
echo "done=$(date -Iseconds)"
