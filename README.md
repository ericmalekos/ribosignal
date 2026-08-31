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

### 2. Prepare the RNA-seq

The model reads **ten channels per nucleotide**: 4 one-hot A/C/G/T, 5 ORF-candidate track, and
1 RNA-seq coverage. The coverage channel comes from a STAR alignment in **transcriptome
coordinates**.

```bash
# a. STAR index over the genome, with the annotation the transcripts come from
STAR --runMode genomeGenerate --genomeDir star_index/ \
     --genomeFastaFiles genome.fa --sjdbGTFfile annotation.gtf \
     --sjdbOverhang $((READ_LEN - 1)) \
     --genomeSAindexNbases $NB          # min(14, log2(genome_len)/2 - 1)

# b. trim adapters. -m 20 only; do NOT pass --maximum-length or --discard-untrimmed,
#    which are correct for ribosome footprints and wrong for RNA-seq
cutadapt -a "$ADAPTER" -m 20 -j "$T" -o t1.fq.gz reads.fq.gz        # single-end
cutadapt -a "$ADAPTER" -A "$ADAPTER" -m 20 -j "$T" \
         -o t1.fq.gz -p t2.fq.gz R1.fq.gz R2.fq.gz                  # paired-end

# c. align to the TRANSCRIPTOME
STAR --genomeDir star_index/ --readFilesIn t1.fq.gz [t2.fq.gz] --readFilesCommand zcat \
     --runThreadN "$T" --outSAMtype None --quantMode TranscriptomeSAM \
     --outFilterMultimapNmax 10 --outSAMattributes NH HI AS nM \
     --outFileNamePrefix out/sample.

# d. verify the BAM, then convert to per-nucleotide coverage
samtools quickcheck out/sample.Aligned.toTranscriptome.out.bam
python scripts/rnaseq_coverage.py \
       out/sample.Aligned.toTranscriptome.out.bam  coverage.hd5  SAMPLE_ID
```

`scripts/xspecies/align_rna_xspecies.sbatch` runs a-d end to end as a SLURM array;
`scripts/xspecies/build_species_refs.sbatch` builds the index and the RiboCode annotation.

Four things that are easy to get wrong, each of which has cost this project real time:

- **`--quantMode TranscriptomeSAM`.** The model works in transcript coordinates. A genome BAM
  will not do.
- **`--sjdbOverhang` and `--genomeSAindexNbases` are per-genome.** The mammalian default of 14 is
  wrong for small genomes (about 10 for a 12 Mb yeast genome), and STAR does **not** error on a
  bad value, it silently builds a poor index.
- **`samtools quickcheck`, not a header grep.** A BAM truncated mid-write keeps a perfectly valid
  header, so `view -H | grep SO:coordinate` accepts it.
- **RNA-seq should be poly(A)-selected.** A ribo-depleted total-RNA arm retains tRNA and 7SL and
  once collapsed a transcript universe from tens of thousands to 3,321. `library_selection = cDNA`
  in SRA metadata does not prove poly(A); check the protocol text.

### 2b. Build the ORF-candidate track

```bash
python scripts/build_orf_track.py \
       --fasta transcripts.fa \
       --out   orf_track_v2_nokozak.npy \
       --mode atg --kozak none
```

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
scripts/xspecies/  reference build and RNA-seq alignment drivers
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

MIT, see `LICENSE`.
