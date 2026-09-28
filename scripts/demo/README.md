# End-to-end demo: public RNA-seq in, ORF calls out, scored against real Ribo-seq

`run_demo.sh` executes the whole pipeline on public data and scores the result. It is the
worked example behind the tutorial at <https://ribosignal.readthedocs.io>. The Nextflow
pipeline runs its prediction half on RNA-seq alone; the Ribo-seq steps and the scoring are
here only.

    THREADS=32 CHROM=all bash scripts/demo/run_demo.sh

**Run it on Linux.** It aborts early on macOS by design; see *Known traps* below.

---

## What it does, and why each choice

One donor, both assays, from Chothani et al. (GEO **GSE182371** Ribo-seq + **GSE182372**
RNA-seq), **Hepatocytes_1**:

| run | assay | | role |
|---|---|---|---|
| `SRR15513269` | RNA-seq | paired, 22,062,488 pairs, 76 nt | the model's input |
| `SRR15513208` | Ribo-seq | single, 182,236,906 reads | the ground truth |

Hepatocytes is the tissue **held out of training** for both released checkpoints, so this
is a genuine held-out test rather than a memorised one, and the matched Ribo-seq means the
run ends with a number rather than an assertion that it ran.

| step | what | note |
|---|---|---|
| 1 | GENCODE v49 + GRCh38 primary assembly, `gffread` transcriptome, STAR index | `CHROM=chr22` (default) for a small run, `CHROM=all` for the whole genome (~32 GB RAM at index time) |
| 2 | fetch FASTQ from ENA | no SRA toolkit needed |
| 3 | **measure** both adapters | `scripts/demo/measure_adapter.py`; see below |
| 4 | RNA-seq → `cutadapt` → STAR `--quantMode TranscriptomeSAM` → `rnaseq_coverage.py` | `--outFilterMultimapNmax 10` |
| 5 | Ribo-seq → `cutadapt` → STAR → `ribo_psites.py` | `--outFilterMultimapNmax 1` and `--alignEndsType EndToEnd`, matching the released checkpoints |
| 6 | `build_pack.py --coverage … --psites …` → `build_orf_track.py --mode ext --kozak none` | the track setting both checkpoints were trained with |
| 7-8 | fetch weights, checksum, `dump_pred_profiles.py` for **both** architectures | |
| 9-10 | `ribocode_dropin.py` (`real`, `pred_preddepth`, `+ --pred_poisson --pred_scale 0.05`) and `ribotish_dropin.py` | two callers separate a model failure from a caller failure |
| 11 | `score_demo.py` | profile Pearson and ORF-call precision/recall/F1, predicted vs observed |

Resumable: each stage skips when its output exists. Delete the output to redo a stage.

## Known traps, all of them hit for real

**The adapter is measured, never defaulted.** `SRR15513269` is *already adapter-trimmed*:
Illumina TruSeq appears in 4 reads out of 400,000, and lengths are ragged 74-76. Passing
`-a` would trim what is not there, and `--discard-untrimmed`, correct for footprints and
catastrophic for RNA, would have destroyed the library. Step 3 prints the evidence and
step 4 asserts the verdict has not changed. `SRR15513208` measured as already trimmed too
(35 nt reads, TruSeq in 0.11%); step 3 still measures it on every run and step 5 acts on
the answer.

**RNA and Ribo need opposite cutadapt flags.** `-M`/`--discard-untrimmed` are right for
length-selected footprints and wrong for RNA. The script uses each in exactly one place.

**macOS cannot run this.** bioconda's osx-64 STAR **2.7.11b aborts on any `--quantMode`**
(`could not open input file /geneInfo.tab`: it resolves `genomeDir` to empty), and
**2.7.10b accepts `--readFilesCommand` but silently delivers ZERO reads and exits 0**,
producing a valid empty BAM. The second is the dangerous one, so step 4 asserts STAR's
input-read count is non-zero regardless of platform.

**The coverage guard is annotation-scope-dependent.** `rnaseq_coverage.py` refuses to
write an hd5 covering fewer than 5,000 transcripts, which is right for a whole
transcriptome and fires spuriously on a chromosome subset (chr22 has 11,615 transcripts in
total). The script scales the floor to `NTX/10` rather than removing it; the guard exists
because a corrupt or untrimmed FASTQ yields a small, valid-looking file that then poisons
the pack.

**`--pred_scale` only means something under `pred_preddepth`.** Under `pred_obsdepth` the
depth is the observed total by definition, and the flag is now a hard error rather than
silently ignored.

**Force the fast Mamba path on a GPU box.** Set `RIBO_MAMBA_IMPL=cuda` so a missing
`mamba_ssm` is an error instead of a silent fall back to the pure-PyTorch reference and a
large slowdown.

## Two scripts this demo needed that the repo did not have

- **`scripts/build_pack.py`**: pools `rnaseq_coverage.py` hd5s into the pack the model
  reads. Without it, `rnaseq_coverage.py`'s output was a dead end.
- **`scripts/ribo_psites.py`**: Ribo-seq BAM to the P-site hd5 `build_pack.py --psites`
  consumes. Without it, the observed arm, the training target, and any
  predicted-vs-observed comparison were unreachable from a public BAM.

`ribo_psites.py` calibrates the 5'-end offset per read length against annotated start
codons and refuses to write a file whose frame-0 fraction is below 0.4.

## Verification status, stated plainly

**Verified.** The whole script has run end to end with `CHROM=chr22` in a fresh working
directory: `ribo_psites.py` calibrated on the real Ribo-seq BAM (frame-0 fraction 0.726),
and the scores are the ones on the tutorial's score page. `ribo_psites.py`'s
transcript-space CDS walk: 0 disagreements against the GTF's independent `start_codon`
features over 5,485 GENCODE v49 chr22 transcripts, and 99.818% (5,475/5,485) land on `ATG`
in the extracted sequence. That check also caught a real defect: 272 CDS-bearing
transcripts are 5'-incomplete with the CDS at transcript position 0, so calibrating against
them calibrates against an annotation artifact; calibration is now restricted to annotated
start codons. `score_demo.py`'s tie-averaged ranks are unit-tested.

**Not verified.** `CHROM=all` has not been run, and the verified run used a conda
environment on the host, not the published image.

## Expected cost

`CHROM=chr22`: about 14 GB disk. On 16 CPU threads prediction alone takes 30 min for `attn`
and 3 h 10 min for `mamba4`; on an A100 both finish in about 95 s.
`CHROM=all`: ~150 GB disk, ~32 GB RAM for `genomeGenerate`, several hours; the Ribo-seq
alignment dominates at 182M reads.
