# RiboSignal

Predicts a **per-nucleotide ribosome P-site profile** for a transcript from its mature mRNA
sequence and matched RNA-seq coverage. No ribosome-profiling experiment is needed at inference.

The predicted profile can be fed to an ORF caller in place of real Ribo-seq, which is what it is
for: calling translated ORFs, including upstream and non-canonical ones, in samples where no
Ribo-seq exists.

**Weights:** <https://huggingface.co/emalek/RiboSignal> (MIT)

| checkpoint | mixer | params | device |
|---|---|--:|---|
| `mamba4_best.pt` | dilated CNN + 4 bidirectional Mamba blocks | 7,521,026 | GPU only (`mamba_ssm`) |
| `attn_best.pt` | dilated CNN + 2 transformer layers | 5,071,106 | CPU or GPU |

---

## Minimal example: predict a profile

### 1. Get the weights

```bash
pip install huggingface_hub torch numpy
python -c "
from huggingface_hub import snapshot_download
snapshot_download('emalek/RiboSignal', local_dir='weights')"
```

### 2. Build the three inputs

A checkpoint alone is not enough. The model reads **ten channels per nucleotide**: 4 one-hot
A/C/G/T, 5 ORF-candidate track, 1 RNA-seq coverage. Each comes from one script in `scripts/`.

```bash
# a. transcript FASTA -> ORF-candidate track (frames, start context, stop)
python scripts/build_orf_track.py \
    --fasta transcripts.fa \
    --out   orf_track_v2_nokozak.npy \
    --kozak none

# b. RNA-seq BAM (transcriptome coords) -> per-nucleotide coverage
python scripts/rnaseq_coverage.py \
    --bam  rnaseq.toTranscriptome.bam \
    --out  coverage.hd5
```

> **`--kozak none` is required.** Both released checkpoints were trained without the Kozak
> heuristic. Feeding them a heuristic-Kozak track silently degrades the prediction.

### 3. Predict

```bash
export RIBO_PACK_DIR=packed/          # tx_order, offsets, lengths, coverage
export RIBO_ORF_TRACK=orf_track_v2_nokozak.npy
export RIBO_ONEHOT_FASTA=transcripts.fa

python scripts/dump_pred_profiles.py \
    --run    weights/ \
    --device cpu                      # use cuda for the mamba4 checkpoint
```

Writes `pred_profiles.npz` with, per transcript, a `pred_flat` profile summing to 1 and a
`pred_total` count.

### 4. Optional: call ORFs from the predicted profile

```bash
python scripts/ribocode_dropin.py \
    --profiles pred_profiles.npz \
    --annot    ribocode_annot/ \
    --variant  pred_preddepth \
    --out      calls/ \
    --min_aa 5 --pval 0.05
```

`--variant pred_preddepth` is the fully de novo call: predicted shape and the model's own count
head, using no observed Ribo-seq. Add `--pred_poisson --pred_scale 0.05` to trade non-canonical
yield for non-canonical precision.

---

## What is here

```
scripts/     the prediction path and the scripts that build its inputs
release/     which checkpoint is which, with configs and held-out metrics
containers/  Dockerfile and Singularity definition
env/         pinned environment specs
tests/       invariant tests
docs/        reference and dataset registries
proteogenomics/   the proteogenomic database-search pipeline
```

## Training

Human tissue Ribo-seq (GEO **GSE182371**), leave-one-tissue-out with Hepatocytes held out,
unique-mapper alignments only. `scripts/train.py` is the training entry point; see
`release/README.md` for the exact configuration of each released checkpoint.

## Citation

Manuscript in preparation. Until then cite this repository and
<https://huggingface.co/emalek/RiboSignal>.

## License

MIT.
