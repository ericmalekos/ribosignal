# RiboSignal

Predicts a **per-nucleotide ribosome P-site profile** for a transcript from its mature mRNA
sequence and matched RNA-seq coverage. No ribosome-profiling experiment is needed at inference.

The predicted profile can be fed to an ORF caller in place of real Ribo-seq, which is what it is
for: calling translated ORFs, including upstream and non-canonical ones, in samples where no
Ribo-seq exists.

**Weights:** <https://huggingface.co/emalek/RiboSignal>

| checkpoint | mixer | params | device |
|---|---|--:|---|
| `mamba4_best.pt` | dilated CNN + 4 bidirectional Mamba blocks | 7,521,026 | GPU only (`mamba_ssm`) |
| `attn_best.pt` | dilated CNN + 2 transformer layers | 5,071,106 | CPU or GPU |

---

## Minimal example: predict a profile

### 1. Get the weights

The published container already has them baked in at `$RIBO_WEIGHTS`:

```bash
docker run --rm -it ghcr.io/ericmalekos/riboseq-model:latest \
       ls $RIBO_WEIGHTS
```

Or fetch them directly:

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

### 4. Optional: tighten precision with the Poisson dial

The fully de novo call over-calls short non-canonical ORFs, because the caller gates on absolute
P-site thresholds while the model's count head is scaled to its training depth. Sampling the
predicted density as Poisson at a reduced effective depth gives it the detection noise a real
experiment has, so weak diffuse ORFs fail the frame test and drop out.

Add two flags to the call in step 5:

```bash
    --pred_poisson --pred_scale 0.05
```

`--pred_scale` (theta) simulates an experiment of that fraction of the predicted depth. Lower is
stricter. At theta = 0.05 on held-out data this costs about six points of annotated-CDS recall and
buys 0.10 to 0.15 of uORF precision and 0.10 to 0.29 of non-canonical precision, roughly halving
the number of non-canonical calls.

Anchor the choice on CDS, the one class the annotation makes trustworthy: pick the most permissive
theta whose CDS precision and recall both stay above 0.90. Note that CDS precision is nearly flat
in theta, so recall is the binding constraint.

### 5. Call ORFs from the predicted profile

```bash
python scripts/ribocode_dropin.py \
    --profiles pred_profiles.npz \
    --annot    ribocode_annot/ \
    --variant  pred_preddepth \
    --out      calls/ \
    --min_aa 30 --pval 0.05
```

`--variant pred_preddepth` is the fully de novo call: predicted shape and the model's own count
head, using no observed Ribo-seq. `--min_aa 30` is the floor used throughout this work; for MHC
immunopeptidomics, where peptides are 8 to 11 aa, the floor is 7 instead.

Two other variants exist for evaluation rather than deployment: `real` runs the caller on an
observed profile, and `pred_obsdepth` uses the predicted shape scaled to an observed depth, which
isolates whether the model places ribosomes correctly.

### 6. Optional: call ORFs with a different caller

The predicted profile is not tied to RiboCode. Four callers have been run on it, and running more
than one separates a model failure from a caller failure. What each needs differs:

| caller | how it takes the prediction | needs a bam? |
|---|---|---|
| RiboCode | native density injection (step 5) | no |
| Ribo-TISH | `predict --inprofile` | no |
| RiboTaper | synthesized genome bam | yes |
| ribotricer | synthesized genome bam | yes |

PRICE was tested and **does not work** on predicted profiles: it fits its cleavage model to the
read mismatch distribution, which a density array does not have. See `docs/price_not_usable.md`
for the evidence and the two fixes that failed.

#### 6a. Ribo-TISH

RiboCode is not the only caller the predicted profile can drive. Ribo-TISH accepts a
per-transcript P-site profile through `--inprofile`, so it needs no bam either, and running both
separates a model failure from a caller failure.

```bash
# the GTF must be RESTRICTED to the transcripts you are scoring: Ribo-TISH only skips the
# bam path when the profile covers every transcript of a gene
awk -v F=tx_ids.txt 'BEGIN{while((getline l < F)>0) k[l]=1}
     !/^#/ { if (match($0,/transcript_id "[^"]+"/)) {
       t=substr($0,RSTART+15,RLENGTH-16); if (t in k) print } }' \
    gencode.annotation.gtf > scored.gtf

python scripts/orfcallers/ribotish_dropin.py \
    --profiles pred_profiles.npz \
    --gtf      scored.gtf \
    --genome   genome.fa \
    --variant  pred_preddepth \
    --out      calls_ribotish/ \
    --longest --minaalen 5 --fpth 0.05 --numproc 8
```

`--variant`, `--pred_scale` and `--pred_poisson` mean exactly what they do in step 5: the density
is built by the same `build_density` the RiboCode drop-in uses, imported rather than
reimplemented, so both callers see identical arrays.

**Pass `--longest`.** Without it Ribo-TISH reports every in-frame downstream ATG as a separate
`Truncated` ORF, which was 86% of calls in testing. `--longest` keeps one ORF per stop codon, the
same rule as RiboCode's `(gene_id, ORF_gstop)` collapsing.

Run it under a python that can import `ribocode_dropin`, since that is where `build_density`
lives. The script shells out to the `ribotish` binary, so the two do not need to share an
environment:

```bash
/path/to/ribocode-env/bin/python scripts/orfcallers/ribotish_dropin.py ... \
    --ribotish /path/to/ribotish-env/bin/ribotish
```

The script checks every transcript's cDNA length against the GTF and drops mismatches rather than
writing a shifted profile; it reports the count, which should be zero if the GTF and the pack came
from the same annotation release.

Output is Ribo-TISH's native table with a `TisType` column (`Annotated`, `5'UTR`, `3'UTR`,
`Internal`, `Novel`, `Truncated`, `Extended`, and `:Known` / `:CDSFrameOverlap` variants). These
do not map one-to-one onto RiboCode's `ORF_type`, and the two callers use different significance
gates, so do not compare their counts directly without first imposing a common threshold.
See `docs/ribotish_dropin_results.md` for per-class precision, recall and F1 on held-out
hepatocytes.

#### 6b. RiboTaper

RiboTaper has no profile hook. It reads bam files at every stage, and two of the four coverage
tracks its ORF finder consumes are built by `coverageBed -abam`, so the predicted profile has to
be turned back into reads. `scripts/orfcallers/synth_ribo_bam.py` does that, and RiboTaper then
runs unmodified.

```bash
# 1. a genome-coordinate Ribo-seq bam whose P-sites reproduce the predicted profile.
#    One read length with one offset, so RiboTaper's P-site recovery is deterministic.
python scripts/orfcallers/synth_ribo_bam.py \
    --profiles pred_profiles.npz \
    --gtf      scored.gtf \
    --fai      genome.fa.fai \
    --variant  pred_preddepth \
    --read_len 29 --offset 12 --selfcheck 200000 \
    --out ribo.unsorted.bam
samtools sort -o ribo.bam ribo.unsorted.bam && samtools index ribo.bam

# 2. one-time annotation build (GTF + samtools-faidx'd genome)
create_annotations_files.bash scored.gtf genome.fa false false annot_dir/

# 3. RiboTaper. n_cores must be > 1; read_lengths and cutoffs must match step 1.
Ribotaper.sh ribo.bam rna.bam annot_dir/ 29 12 12
```

**`--selfcheck` is not optional.** It re-derives the P-site from each emitted read exactly as
RiboTaper's own awk does and asserts it lands on the intended base. It earned its keep: the first
version measured the minus-strand offset from the wrong end of the read and 36,198 of 77,682
reads landed on the wrong base. A silent version of that error shifts every ORF's frame. The
script exits non-zero on any mismatch.

**RiboTaper requires an RNA-seq bam.** RiboCode and Ribo-TISH do not. Use the real matched
RNA-seq; do not synthesise it. `scripts/orfcallers/synth_rna_bam.py` exists only for the case
where no RNA bam can be recovered, and is marked superseded for that reason.

**Cost.** RiboTaper is the slowest of the three and by far the heaviest on disk: about 1 to 2.5
hours per arm on 12 cores for one chromosome, and roughly 9 GB of `P_sites_all` per arm, one line
per P-site. Budget accordingly, or restrict the GTF to a chromosome subset.

Output is `ORFs_max_filt`, with a `category` column (`ORFs_ccds`, `uORF`, `dORF`, `ncORFS`,
`Overl_uORF`, `Overl_dORF`). Its `uORF` and `dORF` are strictly non-overlapping: an upstream ORF
that runs into the CDS is relabelled `Overl_uORF` (`CCDS_orf_finder.R:989`), verified here on real
data as 0 of 192 `uORF` calls extending past the annotated start.

#### 6c. ribotricer

The lightest of the four. It builds its own candidate-ORF index once, then reads a genome bam and
infers P-site offsets itself by cross-correlation.

```bash
ribotricer prepare-orfs --gtf scored.gtf --fasta genome.fa \
    --min_orf_length 90 --start_codons ATG --prefix idx/chr1

ribotricer detect-orfs --bam ribo.bam --ribotricer_index idx/chr1_candidate_orfs.tsv \
    --prefix calls_ribotricer/chr1 --stranded yes
```

The index build takes seconds (36,392 candidate ORFs over 6,864 transcripts) and `detect-orfs`
tens of minutes. Output is `<prefix>_translating_ORFs.tsv` with a `phase_score`, plus metagene
profiles, the inferred P-site offsets and wig tracks. Use the same synthesized bam as RiboTaper.

#### 6d. Projecting a real transcriptome bam onto the genome

RiboTaper and ribotricer need genome coordinates, and this project keeps transcriptome bams. This
converts one to the other, so the **observed** arm can use real reads rather than synthesized
ones, and so genome-based QC tools become reachable.

```bash
python scripts/orfcallers/tx_bam_to_genome.py \
    --in-bam sample.Aligned.toTranscriptome.out.bam \
    --gtf annotation.gtf --fai genome.fa.fai \
    --out genome.bam --name-sort
samtools sort -o genome.sorted.bam genome.bam && samtools index genome.sorted.bam
```

**It deduplicates isoform expansion, which is the whole reason it is not a one-liner.**
`--quantMode TranscriptomeSAM` writes one genomic alignment once per compatible isoform, so a
genomically unique read routinely carries `NH:i:17`. Projecting every record would write that read
seventeen times at one locus. Secondary records are dropped before projection, alignments are
deduplicated per read name on the projected coordinate, and **NH is recomputed** as the number of
distinct genome loci, which is the genomic multimapping the transcriptome NH never measured. On
five Hepatocyte runs: 3.54 billion records in, 105.4 million genome records out, 0 cDNA-length
mismatches, 0 unprojectable.

Input must be name-grouped. STAR's raw output is; anything coordinate-sorted is not, so
`--name-sort` re-sorts first. Verified against the reference: 93.1% of projected reads match the
genome exactly and 98.7% at 90% identity or better, with junction-spanning reads scoring slightly
**better** than ungapped ones (0.9976 vs 0.9947 mean identity), which is the check that the exon
walk and the N insertion are right.

### Comparing the callers

`scripts/orfcallers/compare_three_callers.py` puts the callers on one class vocabulary and one key,
and `export_caller_tables.py` writes the three summary tables from its JSON.

The callers use **different, strand-dependent coordinate conventions**. Keyed on the ORF's genomic
start, the offsets against RiboCode's `ORF_gstart` are Ribo-TISH `+1` on the plus strand and `0`
on the minus, RiboTaper the reverse. Keying on the stop coordinate does not work at all, because
the three disagree over whether the stop codon is inside the ORF. The comparison script re-derives
these offsets on every run rather than hardcoding them, so a change in any caller's output format
shows up as a shifted offset instead of silently as lost overlap.

Do not read cross-caller count differences as biology without checking this. On observed
hepatocyte data, with no model involved, three callers agree on only 914 canonical ORFs out of
1,160 to 1,852, and on 5 dORFs out of 14 to 25.

## Training

Human tissue Ribo-seq (GEO **GSE182371**), leave-one-tissue-out with Hepatocytes held out,
unique-mapper alignments only. `scripts/train.py` is the training entry point. The exact
configuration and held-out metrics for each checkpoint ship with the weights as
`<arch>_config.json` and `<arch>_test_metrics.json`.

## Citation

Manuscript in preparation. Until then cite this repository and
<https://huggingface.co/emalek/RiboSignal>.

## License

Not yet specified.
