#!/usr/bin/env bash
#SBATCH --job-name=o2_sym_janich
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A5500:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=08:00:00
#SBATCH --output=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/o2_sym_%j.out
#SBATCH --error=/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/logs/janich/o2_sym_%j.err
set -euo pipefail
# O2 SYMMETRIC ARM: rebuild the Janich pack on the SAME universe as Wang, then dump + call.
#
# Why the rebuild: the original Janich pack was built on a TPM>=1 expressed universe and came out at
# 2,434 transcripts, because mouse liver's transcriptome is dominated by a handful of genes (Alb/Ttr/Apo*)
# which compresses every other isoform's TPM below the threshold. Wang's pack, by contrast, is 221,835
# transcripts (essentially the whole annotation). Comparing a 2,434-tx model arm against a 221,835-tx one
# is not a symmetric check -- the Janich arm would be structurally unable to recover most of Wang's calls.
# --ref-pack puts both arms on identical tx_order/lengths, so the only difference left is the data.
BASE=/private/groups/carpenterlab/emalekos/RNAZoo_meta
NEW=$BASE/RNAZoo/experiments/riboseq_signal_model
ECH=$BASE/RNAZoo/experiments/biotype_probe/expression_context_human
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
RIBOPY=/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python3
SIF=/private/groups/carpenterlab/emalekos/singularity_cache/ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img
RUN=$NEW/results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes
REF=$NEW/data/heldout_refs
WANGPACK=$NEW/data/packed_heldout_mouse_wang_liver
PACK=$NEW/data/packed_heldout_mouse_janich_liver_wanguni
OUT=$NEW/results/o2_liver/janich_wanguni
mkdir -p "$NEW/logs/janich" "$OUT"
echo "O2 symmetric start=$(date -Iseconds) host=$(hostname)"

# 1. Janich pack on the Wang universe (Janich Ribo psites + Janich RNA coverage, Wang tx_order)
"$PY" "$NEW/scripts/prepare/prepare_pack.py" \
  --ribo-psites "$ECH/data/ribocode_mouse_liver/Liver_5samp" \
  --rna-coverage "$NEW/data/heldout_rnaseq_coverage/mouse_janich_liver" \
  --ref-pack "$WANGPACK" \
  --tx2biotype "$NEW/data/tx2biotype_mouse.tsv" \
  --build-orf-track --fasta "$REF/mouse_wang_universe.fa" --kozak none \
  --group janich_liver_wanguni --species mouse --out "$PACK"
echo "pack: $(wc -l < "$PACK/tx_order.txt") tx"

# 2. predict on it (nokozak model <- nokozak track; train/inference must match)
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
  --env RIBO_ORF_TRACK="$PACK/orf_track_v2_nokozak.npy" \
  --env RIBO_ONEHOT_FASTA="$REF/mouse_wang_universe.fa" \
  -B "$BASE" "$SIF" \
  python3 "$NEW/scripts/dump_pred_profiles.py" --run "$RUN" \
    --heldout mouse_janich_liver_wanguni --tx_list "$REF/mouse_wang_universe.txt" \
    --fasta "$REF/mouse_wang_universe.fa" --out "$OUT" --device cuda

"$PY" -c "
import numpy as np,sys
d=np.load('$OUT/pred_profiles.npz',allow_pickle=True)
nt=len(d['tx_ids']); ns=int(str(d['meta'][3])); f=ns/(nt+ns) if (nt+ns) else 0
print(f'dump: {nt} tx, {ns} skipped ({f:.1%})'); sys.exit(1 if f>0.05 else 0)"

# 3. RiboCode drop-in (mouse annot): real + the de-novo arm
export PYTHONNOUSERSITE=1
mkdir -p "$OUT/dropin"
for V in real pred_preddepth pred_obsdepth; do
  "$RIBOPY" "$NEW/scripts/ribocode_dropin.py" --profiles "$OUT/pred_profiles.npz" \
    --annot "$ECH/data/ribocode_annot_mouse" --variant "$V" --out "$OUT/dropin" \
    --min_aa 5 --pval 0.05
done

# 4. symmetric comparison: Janich model vs WANG observed (mirror of the landed Wang arm)
"$PY" "$NEW/scripts/heldout/o2_between_dataset.py" \
  --model "$OUT/dropin/pred_preddepth_collapsed.txt" \
  --model-obsdepth "$OUT/dropin/pred_obsdepth_collapsed.txt" \
  --model-name "Janich_model" --other-name "Wang_obs" --swap \
  --out "$NEW/results/o2_liver/o2_summary_symmetric.txt"
echo "O2 symmetric done=$(date -Iseconds)"
