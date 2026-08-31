# Methods: Ribo-seq signal prediction model

> Current vs superseded results: **`docs/STATUS_CURRENT_VS_ARCHIVED.md`**.

Goal: train a transformer to predict per-nucleotide Ribo-seq P-site signal directly
from (a) foundation-model per-token embeddings of the transcript (RiNALMo and/or
Orthrus) and (b) per-nucleotide binned RNAseq coverage, using the Chothani matched
primary-tissue Ribo-seq + deep RNAseq (GENCODE v49 / GRCh38). First evaluation is
within-tissue, held-out-chromosome (gene-disjoint); the real target is a model that
generalizes to new datasets given deep RNAseq.

Data reused from the sibling project `expression_context_human` (referred to as `$ECH`):
`/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human`.
This project lives at
`/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model`
(referred to as `$NEW`). Start tissue: Fibroblast (deepest, 32 Ribo-seq samples).

---

## 1. Coordinate frame and axis alignment (why transcriptome coords)

Everything is built on the transcript-nucleotide axis (spliced transcript, 0-based),
so the three signals align 1:1 position-for-position with no splicing bookkeeping:

- Ribo-seq P-site target: per-nt array of length L per transcript (this document).
- RiNALMo per-token embedding `{tx}_tokens.npy` is shape `(L, 1280)` where `L` equals
  the transcript length in nt exactly (CLS/EOS special tokens stripped by the extractor),
  row i == nucleotide i. Verified in planning against the GENCODE v49 transcript FASTA
  length field on the 444-transcript pilot.
- Future per-nt RNAseq coverage will be produced on the same axis.

Transcript identity is the versioned Ensembl transcript id (e.g. `ENST00000832824.1`).
The RiboCode transcriptome BAM `@SQ` names, the P-site hd5 `transcript_ids`, the
RiboCode ORF `.txt` `transcript_id`, and GENCODE FASTA header field 1 are all this same
versioned id, so joins are direct with no stripping.

Storage: one variable-length array per transcript (HDF5 vlen), never an `(N, L, D)`
dense/memmap block. On this parallel filesystem, slicing the length axis across many
rows of an `(N, L, D)` memmap is pathologically slow (scattered reads); the vlen
one-row-per-transcript layout is read whole-row, which is fast and is exactly the schema
RiboCode already uses for the inputs.

---

## 2. Ribo-seq target: input provenance and cleaning posture (verified)

### 2.1 Source hd5 (RiboCode per-nt P-sites)

RiboCode (`process_bam.py`, RiboCode 1.2.15) already wrote per-nucleotide P-site arrays
in transcriptome coordinates, one HDF5 per Fibroblast Ribo-seq sample:
`$ECH/data/ribocode_per_tissue/Fibroblast/SRR155131{65..96}.Aligned.toTranscriptome.out_psites.hd5`
(32 files, ~3.7 GB each). Fibroblast Ribo-seq samples = SRR15513165 to SRR15513196.

Each hd5 schema (all 32 share an identical 509,650-transcript axis in identical order,
the full GENCODE v49 annotation: protein_coding + lncRNA + pseudogenes + others):

- `transcript_ids` : vlen-str, versioned ENST, shape (509650,).
- `p_sites` : vlen-int32, shape (509650,); `p_sites[i]` is a per-nucleotide array of
  length equal to that transcript's length, 0-based transcript coordinates.
- attr `psites_number` : total distinct P-site reads in that sample (for normalization).

P-site assignment inside RiboCode: for each read on the sense strand (strand flag `yes`),
of a read length listed in that sample's config, P-site position =
`read.reference_start + offset[read_length]`. The per-sample, per-read-length offsets
and strand are in `$ECH/data/ribocode_per_tissue/Fibroblast/Fibroblast_pre_config.txt`
(one config block per sample; 28/29/30-nt reads to offset 12, 26-nt to 9, some 33/34-nt
to 6 or 15; strand `yes` for all 32). These are already baked into the hd5, so no
offset work is needed when pooling.

### 2.2 Cleaning posture of the per-sample Ribo-seq BAMs (Task 1 verification)

There are three per-sample BAM flavors per SRR under `$ECH/data/riboseq_bam/`, in
DIFFERENT states. This distinction matters:

| per-sample file | posture | evidence |
|---|---|---|
| `{SRR}.Aligned.sortedByCoord.out.bam` (genomic) | CLEANED in place: NH==1 kept + rRNA/tRNA/miRNA loci dropped | 0 secondary records, NH all == 1; read count matches sum going into the per-tissue merge exactly |
| `{SRR}.Aligned.toTranscriptome.out.bam` (transcriptome) | CLEANED in place: ncRNA transcripts dropped + cross-gene multimappers dropped; within-gene isoform multimappers RETAINED (each counted) | primary count roughly halved vs raw STAR output; coordinate-sorted + indexed |
| `{SRR}.toTranscriptome.sorted.bam` (transcriptome) | RAW / dirty: never filtered (120M+ records, NH up to 1000+, rRNA/tRNA laden) | DO NOT USE |

Alignment posture: STAR ran with `--outFilterMismatchNoverLmax 0.05`,
`--outFilterMatchNminOverLread 0.7`, `--outFilterMultimapNmax 20`,
`--alignEndsType EndToEnd`, `--quantMode TranscriptomeSAM` (cutadapt `--trim-n
--minimum-length 20`, no 3' adapter for human). Multimappers were ALLOWED at align time;
the NH==1 filter is a separate Phase-30 in-place step applied to the GENOMIC BAMs only.
Per-tissue merges (`$ECH/data/ribo_per_tissue_bam/{TISSUE}.merged.bam`) are a lossless
union of the already-cleaned genomic BAMs.

Provenance of the target inputs: the Fibroblast hd5 are timestamped Jun 28 00:13 and
later; the in-place clean of `Aligned.toTranscriptome.out.bam` finished Jun 27 17:42.
RiboCode therefore read the CLEANED transcriptome BAM (rRNA/tRNA/miRNA transcripts and
cross-gene multimappers already removed). This is confirmed empirically by
`scripts/verify_target_inputs.py` (Step 1): drop-list ncRNA transcripts on the axis pool
to approximately zero P-sites while housekeeping controls (ACTB, GAPDH) pool to large
counts (see results.md).

### 2.3 Multimapper posture for the target (decision A, user-confirmed)

On the transcriptome axis, a footprint that is genomically unique still maps to every
compatible isoform of its gene. The target counts these within-gene isoform multimappers
with posture A: reuse the 32 existing RiboCode hd5 as-is, each footprint counted weight 1
on every compatible isoform (RiboCode default). Rationale: fastest path (data exists,
offsets/strand baked, ncRNA + cross-gene already clean); consistent with the sibling
project's `_ribo_coverage` feature, which is built from these same hd5; and the isoform
magnitude inflation is absorbed by per-transcript normalization at train/eval time, since
the model is scored with scale-invariant per-transcript correlation and 3-nt periodicity.
A primary-alignment-only rebuild (posture B) is pre-registered as a later one-off
sensitivity check; fractional 1/N (posture C) is a possible future extension.

---

## 3. Target build (Step 2)

Script: `scripts/build_fibroblast_psite_target.py` (runner
`scripts/build_fibroblast_psite_target.sbatch`, SLURM partition medium, 6 h, 2 cpu, 32 G,
python `conda_envs/cas12a/bin/python3`, numpy + h5py). Adapts the sibling project's
`$ECH/scripts/phase31_7_build_ribo_coverage.py` (expression_context_human), replacing its scalar
`p_sites[i].sum()` with per-nucleotide element-wise accumulation.

Method:

- Positional accumulation. All 32 hd5 share an identical transcript order; the build
  asserts `np.array_equal(tx_ids, canonical)` per file (fail fast) and accumulates
  `acc[i] += p_sites[i]` by index, avoiding per-transcript string hashing.
- int64 accumulation (grand total is order 1e9, near the int32 ceiling). `p_sites`
  is stored int32 only after asserting the global per-nt max is below 2^31 (it is;
  a single pooled nucleotide never approaches that); the guard falls back to int64
  otherwise.
- Chunked reads (10,000 transcripts per slice) so no single 3.7 GB `p_sites[:]` is
  ever materialized; peak memory is the int64 ragged accumulator (about 8 GB) plus the
  int32 write buffer.
- Build for all 509,650 transcripts (cheap, one pass). Selecting an expressed
  protein_coding + lncRNA universe is a downstream filter over the summary TSV, no rebuild.

Outputs (`$NEW/data/target/`):

- `Fibroblast_psites_pooled.hd5` mirroring the RiboCode schema: `transcript_ids`
  (vlen-str) + `p_sites` (vlen-int32, per-nt pooled, gzip-4). Attrs: `tissue`,
  `n_samples`, `n_transcripts`, `psites_number_total`, `psites_number_per_sample`
  (JSON srr->int), `multimapper_posture`, `source_glob`.
- `Fibroblast_psites_summary.tsv` (one row per transcript): `tx_id, length,
  total_psites, n_nonzero_nt, max_psite, gene_id, gene_name, chrom, transcript_type`.

### 3.1 Biotype / chromosome table

`scripts/build_tx2biotype.py` parses the transcript and exon features of
`/private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations/gencode.v49.annotation.gtf`
into `$NEW/data/tx2biotype.tsv` (`tx_id, gene_id, gene_name, chrom, strand,
transcript_type, gene_type, length`; length = summed exon lengths = mature transcript
length). This covers 507,365 transcripts. About 2,285 transcripts on the 509,650-transcript
hd5 axis have no GTF transcript-feature match and join as `unknown`/`NA` in the summary;
protein_coding and lncRNA membership can be independently cross-checked against the two
GENCODE transcript FASTAs (245,535 pc + 197,211 lncRNA). This table is reused for the
later held-out-chromosome split (gene to chromosome).

---

## 4. RNAseq per-nucleotide coverage (second input feature)

The model's second input is per-nucleotide RNAseq read coverage on the same
transcriptome-nucleotide axis as the P-site target, so it aligns 1:1 with the target and
the FM embeddings. It carries RNA abundance/context, which the FM embedding does not
encode (a Ridge probe recovers length/GC from the embedding at R^2 ~ 0.9 but mean-TPM at
R^2 ~ 0), so the two inputs are orthogonal and genuinely additive.

Samples: the 32 Fibroblast RNAseq paired-end samples SRR15513233 to SRR15513264
(`$ECH/data/rnaseq_input/{SRR}_{1,2}.fastq.gz`; a different SRR block from the Ribo-seq
165-196). Reads ~76 nt; 3' adapter content is negligible (about 0.0015% of reads carry
AGATCGGAAGAGC), so no trimming is applied. Strandedness is ISR (reverse-stranded),
salmon-inferred with compatible_fragment_ratio 1.0 and strand_mapping_bias 0.0.

Alignment (`scripts/align_rnaseq_sample.sh`, SLURM array `align_rnaseq_array.sbatch`,
env conda_envs/riboseq, STAR 2.7.11b):

    STAR --genomeDir genomes/star_index_grch38_v49 \
         --readFilesIn {SRR}_1.fastq.gz {SRR}_2.fastq.gz --readFilesCommand zcat \
         --outSAMtype None --quantMode TranscriptomeSAM \
         --outFilterMultimapNmax 1 --outSAMattributes NH HI AS nM

Multimap posture (updated 2026-07-21): `--outFilterMultimapNmax 1` keeps only reads that map
uniquely at the GENOME level. This was `20` originally; it is now `1` for every RNA-seq coverage
alignment, for three reasons: smaller intermediates, faster STAR, and one consistent multimap
posture shared with the Ribo-seq alignment. This is licensed by the DoHH2 mm1-vs-mm20 proxy
(results.md Task 25): the per-nt coverage the model reads is identical on 85.7% of expressed
transcripts (median per-transcript Pearson 1.0000; only the ~5% paralog/repeat/multimap-heavy
tail differs), so single-mapper coverage does not materially change the model input, and the
deployed model (trained on mm20 coverage) accepts mm1 coverage without a retrain. NOTE this
genomic single-mapper filter is ORTHOGONAL to the isoform-level posture A below: a read that is
unique in the genome is still projected onto all compatible isoforms of its gene by
rnaseq_coverage.py (posture A is unchanged).

The index was built from the same gencode.v49 GTF and has exactly 509,650 transcripts,
so the RNAseq toTranscriptome `@SQ` axis is the same transcript set as the Ribo-seq
target (joined by versioned tx-id; the STAR `@SQ` order may differ from the RiboCode
order, so joins are by id, not position). `--outSAMtype None` suppresses the genomic BAM;
only the transcriptome BAM is produced, in node-local scratch, and deleted after coverage
(the BAMs never touch the group filesystem).

Coverage (`scripts/rnaseq_coverage.py`, pysam):
  - Strand: ISR sense-only. Each transcript is + in its own coordinates, so the sense
    mates are read2-forward and read1-reverse, i.e. a mate is kept when
    `is_read2 != is_reverse`. This drops antisense mates so antisense-overlapping loci do
    not get spurious coverage.
  - Multimapper posture A (matches the target): no NH filter; every alignment record
    (primary + secondary) contributes, so a fragment compatible with N isoforms covers
    all N.
  - Depth: each kept mate adds 1 to every transcript position it spans
    (reference_start .. reference_end, exclusive), per-mate depth like `samtools depth`;
    skips unmapped + supplementary.
  Output per sample: `data/rnaseq_coverage/per_sample/{SRR}_coverage.hd5` (vlen schema:
  transcript_ids + coverage int32 gzip-4; attrs srr, libtype, records_seen,
  mapped_sense_counted).

Counter choice (benchmarked on SRR15513233): the pysam counter above is the fleet tool.
Posture A on this deep RNAseq gives about 1.4 billion alignment records per sample (a
55.4M-pair sample expands to 1,406,801,328 records, mean within-gene isoform multiplicity
~12.7x per mate). deepTools `bamCoverage` was evaluated as a multithreaded alternative
(`scripts/align_rnaseq_sample_bamcov.sh`, `bw_to_coverage_hd5.py`) but is pathologically
slow on this workload: 509,650 transcript "chromosomes" at `--binSize 1` over a 1.4B-record
BAM ran >3.5 h on one sample and was killed, versus 81 min for the single-threaded pysam
counter (STAR ~12-20 min + coverage ~60 min). The pysam counter needs no coordinate sort
and streams the BAM once, so it wins here despite being single-threaded. The bamCoverage
scripts are retained for reference but not used. To limit reserved-core waste during the
single-threaded coverage phase, the array requests 4 cpus per task (STAR still fast enough)
at 12 concurrent.

Pooling (`scripts/pool_rnaseq_coverage.py`): element-wise sum of the 32 per-sample
coverage arrays (positional accumulation with an axis-order assert, int64 accumulate,
int32 store-guard) into `data/rnaseq_coverage/Fibroblast_rnaseq_coverage_pooled.hd5`,
the same schema and axis as the pooled target.

Depth normalization of the coverage input (`--cov_norm`, added Task 10). The pooled coverage is
raw per-nt read depth, so `log1p(cov)` carries the dataset's absolute sequencing depth -- fine
within one dataset but not transferable across datasets of different depth. The `global_mean`
mode (default for new runs) divides coverage by the dataset global mean per-nt depth (7,733 for
Fibroblast; `data/packed/coverage_norm.json`) before `log1p`, a CPM-like normalization that a new
dataset reproduces from its own coverage. It is depth-invariant by construction (a 10x-deeper
library with its own mean gives an identical input). `raw` is retained as the legacy path so
Tasks 6 to 9d (all raw) evaluate unchanged; `eval_extra`/plotting read `cov_norm` from `args.json`
(default `raw`). No model change -- the count head sums the now-normalized channel. The 5-fold
within-Fibroblast check (results.md Task 10) shows this costs a small ~0.01 median Pearson (the raw
absolute depth was a mild within-dataset abundance prior) -- the expected, worth-paying price of
making the input transferable for the cross-tissue / cross-dataset objective.

---

## 4B. Expressed universe, per-token FM embeddings, and splits

Expressed universe (`scripts/define_universe_and_fasta.py`): protein_coding + lncRNA
transcripts with Fibroblast salmon meanTPM >= 1 (`scripts/compute_salmon_mean_tpm.py`,
mean over the 32 RNAseq samples) and mature length <= 10000 nt (the RiNALMo cap) =
36,668 transcripts (34,248 pc + 2,420 lncRNA), written to `data/fibroblast_universe.tsv`
(with chrom, gene, mean_tpm, ribo_psites) and `data/fibroblast_universe.fa` (bare
versioned-ENST headers).

Per-token RiNALMo embeddings (`scripts/extract_rinalmo_universe_array.sbatch`): RiNALMo
giga-v1 (1280-d per token) in the ghcr rnazoo-rinalmo SIF on A5500, one {tx}_tokens.npy
per transcript, shape (L, 1280) float16 with L == transcript length (CLS/EOS stripped),
so each row aligns 1:1 to the same transcript-nt position as the P-site target and the
RNAseq coverage. The FASTA is split into 8 length-balanced chunks (`scripts/split_fasta.py`)
run as a GPU array; outputs under `data/rinalmo_token_emb/chunk_{1..8}/`, indexed to
tx_id -> path by `scripts/build_embedding_index.py` (`data/rinalmo_token_emb/tx_index.tsv`).
36,668/36,668 extracted, 0 dropped, ~234 GB. The Orthrus 4-track per-token counterpart
(section 4C) is extracted the same way for the embedding-model comparison.

Held-out-chromosome split (`scripts/build_chrom_split.py`): transcripts grouped by
chromosome (gene-disjoint automatically) into 5 folds balanced by transcript count
(`data/splits/fibroblast_chrom_kfold.json`); fold 0 (chr1/7/8/22, ~7,473 test tx) is the
canonical first held-out-chromosome test.

Mitochondrial exclusion: chrM protein-coding genes are dropped from every split. Mitochondrial
mRNAs are mitoribosome-translated (a different genetic code), leaderless (`utr5=0`, no uORFs), and
lack the cytoplasmic 3-nt periodicity the ORF track and FM embeddings encode, so they are not a
valid target and the model predicts them poorly (13-gene whole-tx Pearson median 0.124 vs 0.618 for
real pc). `Mt_rRNA` / `Mt_tRNA` were already removed by the ncRNA drop-list; the 13 chrM mRNAs
(MT-ND*, MT-CO1/2/3, MT-CYB, MT-ATP6/8) slipped through as `protein_coding`. Two gates:
`dataset.excluded_tx()` filters `load_split` (train + eval; effective for the current pack and all
LOTO runs), and `define_universe_and_fasta.py` skips `chrom in {"chrM"}` (future universe / FASTA /
embedding builds). Impact on the reported Fibroblast matrix (computed with the 13 present) is
negligible: pc profile median 0.6180 -> 0.6181, uORF/dORF unaffected (leaderless), rankings unchanged.

---

## 4C. Orthrus 4-track per-token embeddings (FM comparison)

To test whether the choice of foundation model matters for the profile, the same universe
is embedded with a second FM, Orthrus (a Mamba-SSM mature-mRNA model, 512-d), and the
model is trained on RiNALMo alone, Orthrus alone, and the two concatenated, on the
identical fold-0 split.

- **Variant: 4-track, sequence-only.** Orthrus ships a 4-track model (A/C/G/T one-hot) and
  a 6-track model (4 nt + per-nt CDS + 5' splice indicator channels built from a GTF). The
  6-track channels would hand the model the annotated CDS boundaries, which is exactly the
  ORF location that shapes the Ribo-seq profile, so 6-track is a different, easier task
  (annotation-informed) rather than a fair sequence-only FM. The comparison uses the
  **4-track** model (`quietflamingo/orthrus-large-4-track`), matching RiNALMo's
  sequence-only input.
- **Per-token, not pooled.** The existing RNAZoo Orthrus extractor calls
  `model.representation(x, lengths, channel_last=True)`, which mean-pools over the length
  axis to a single (512,) vector. The signal model needs the pre-pool per-token tensor:
  `model.forward(x, channel_last=True)` returns (1, L, 512). Orthrus one-hot input has no
  CLS/EOS special tokens (unlike RiNALMo), so row i is nucleotide i with no slicing -- it
  aligns 1:1 to the P-site target, RNAseq coverage, and RiNALMo per-token rows. Verified on
  GPU that `forward(...).mean(axis=0)` reproduces `representation(...)` to ~1e-7 (so
  `forward` is the true pre-pool tensor) and that Orthrus L == RiNALMo L for cross-checked
  transcripts.
- **Method** (`scripts/extract_orthrus4t_pertoken.py`, `scripts/extract_orthrus4t_array.sbatch`):
  Orthrus SIF on A5500, `--max-length 10000` and the same 8 length-balanced FASTA chunks as
  RiNALMo, so the Orthrus tx set matches the RiNALMo set 1:1 and the split is identical
  across backends. One `{tx}_tokens.npy` (L, 512) float16 per transcript under
  `data/orthrus4t_token_emb/chunk_{1..8}/`, indexed by
  `scripts/build_embedding_index.py --emb-dir data/orthrus4t_token_emb`. Skip-if-exists
  resume + per-record OOM guard + length-sorted iteration, as for RiNALMo.
- **Training** (`scripts/train_backends.sbatch`): `train.py --emb_backend {rinalmo,orthrus,concat}`.
  The loader (`dataset.py` `PackedStore` + `BACKENDS`) takes one index (single FM) or two
  (concat, feature-axis stack -> 1792-d), and `usable()` intersects the chosen indices so
  every backend trains and tests on the same transcripts. All three use identical
  hyperparameters, seed, and fold-0 split; only the input embedding differs. Runs land in
  `results/backend_cmp/{rinalmo,orthrus,concat}_f0/`; `scripts/aggregate_backends.py`
  collates them.

---

## 5. Verification method

Step 1 (`scripts/verify_target_inputs.py`, light, head-node safe):

- ncRNA posture: partition the 7,587-transcript drop list
  (`$ECH/data/phase30/ncrna_filter_tx.txt`) into present-on-axis vs absent; fancy-index
  the present rows across all 32 hd5 and confirm pooled P-sites are approximately zero
  (rRNA/tRNA/miRNA reads were dropped upstream). Report residual and its biotype makeup.
- Positive controls: pool all ACTB and GAPDH isoform rows across the 32 hd5; expect very
  large counts (proves the reader sees real signal, so the ncRNA zero is real, not an
  all-zero bug). Also confirm hd5 per-nt array length equals the GENCODE mature length
  from tx2biotype for those transcripts (proves the per-nt axis is the transcript-nt axis).

Step 2 sanity (Tier 2, reductions over the summary TSV, no second heavy pass): grand
total pooled P-sites, number of transcripts with nonzero signal, and the
protein_coding vs lncRNA vs other breakdown of signal-bearing transcripts.

Cross-check (post-build): for spot-checked transcripts, the pooled hd5 per-nt array
length equals the corresponding RiNALMo `{tx}_tokens.npy` row count L.

---

## 5B. Model design rationale (data-driven)

Design choices are grounded in the empirical target/coverage distributions measured by
`scripts/inspect_distributions.py` (report in `logs/inspect_distributions.txt`; sample of
1,500 expressed+translated transcripts, plus whole-universe per-transcript summaries).
Key measurements (Fibroblast):

- Per-nucleotide P-site counts are sparse and extremely overdispersed: 74.2% of nt are
  zero, per-nt mean 6.66, variance 15,674, **var/mean ~ 2,355**. A Poisson head (which
  assumes var == mean) is therefore inappropriate. The signal is also periodic: P-site
  autocorrelation peaks at lags 3/6/9/12 (0.32/0.24/0.23/0.21) versus ~0.02-0.05 off-frame
  -- the 3-nt ribosome footprint signature is strong and must be reproduced, so it is a
  first-class evaluation metric.
- Cross-transcript dynamic range is huge (median 3,463 P-sites/tx, mean 17,741, max
  2.8M). Absolute Ribo-seq magnitude is also a function of library depth, which differs
  arbitrarily across datasets. Both facts mean the transferable, dataset-invariant target
  is the **per-transcript-normalized P-site profile** (shape), not absolute counts.
- RNAseq coverage is dense (5.0% zero, per-nt mean 147, max ~20k) -- a good dense
  conditioning input; fed as `log1p(coverage)` given the 0-to-20k range.
- The two inputs are orthogonal in exactly the intended way: within a transcript,
  coverage-vs-P-site per-nt Spearman is only ~0.17 (coverage does NOT localize where
  ribosomes sit), while across transcripts the totals correlate (log1p Pearson ~0.78,
  coverage-vs-P-site). So the embedding must supply intra-transcript localization
  (CDS vs UTR, reading frame) and coverage supplies transcript-level magnitude. (The
  weaker total-coverage-vs-salmon-TPM r ~ 0.41 is the expected transcript-length confound
  -- total coverage scales with TPM*length; the length factor cancels in the
  coverage-vs-P-site total comparison, hence its cleaner 0.78.)

These lead to a **BPNet-style dual-head design** (Avsec et al. 2021, base-resolution
count profiles) over the per-nt FM embeddings:

- Input per nt: `[ FM emb (1280 for RiNALMo) ; log1p(pooled RNAseq coverage) ]`, projected
  to a working width.
- Body: a stack of residual **dilated 1-D convolutions** (dilations 1,2,4,...,512) rather
  than full self-attention over up to 10,000 positions -- O(L) not O(L^2), and the FM
  embedding already encodes long-range sequence context, so the head's job is to mix in
  coverage and sharpen to a base-resolution, frame-aware profile (a local operation).
  A local-attention transformer variant is kept as an alternative to test.
- Profile head: per-nt logits -> log-softmax over the transcript -> **multinomial NLL**
  against the per-nt counts. Scale-free, captures periodicity and CDS localization, and is
  the transferable target for new datasets.
- Count head: a scalar total (from pooled coverage + embedding) with a
  negative-binomial / log-count regression, capturing magnitude. Dataset-depth-dependent,
  so it is the non-transferable head; profile is what carries across datasets.
- Per-transcript normalization is intrinsic (multinomial profile is per-transcript by
  construction), directly addressing the cross-transcript dynamic range.

Evaluation (held-out chromosome, fold 0): per-transcript Pearson and Spearman between
predicted and observed per-nt profile (scale-invariant, so it isolates shape), plus
recovery of the lag-3/6/9 periodicity peak, reported separately for protein_coding and
lncRNA transcripts.

---

## 5C. Implementation and training protocol

Packed training data (`scripts/pack_target_coverage.py` -> `data/packed/`): for the 36,668
universe transcripts (sorted tx id = canonical order), the per-nt P-site counts (from the
target hd5) and pooled RNAseq coverage (from the pooled coverage hd5) are joined BY id and
concatenated into ragged 1-D int32 arrays (`target_counts.npy`, `coverage.npy`) with an
`offsets.npy` index (row k = `arr[offsets[k]:offsets[k+1]]`), plus `lengths.npy`,
`tx_order.txt`, and `pack_meta.tsv`. This makes the loader numpy-only (so it runs in the
CUDA SIF, which lacks h5py) and reads contiguous per-transcript slices (fast, memmap, no
`(N,L,D)` middle-axis pathology, far cheaper than random hd5 vlen-row reads per epoch). The
packer asserts every universe tx is present in both hd5 and that the target and coverage
per-nt lengths agree per transcript, so structural misalignment fails loudly.

Data loader (`scripts/dataset.py`): per transcript, feats = `[ emb (L, d_emb) float32 ;
log1p(coverage) (L, 1) ]`; embeddings are read from the per-tx `.npy` via `tx_index.tsv`,
target/coverage sliced from the packed arrays. Splits are gene-disjoint by chromosome from
`fibroblast_chrom_kfold.json`: fold 0 = test, fold 1 = val, folds 2/3/4 = train. Batching
is by a token budget (summed length per batch, length-bucketed) so variable-length
transcripts (up to 10,000 nt) pack with minimal padding and even GPU load; padded positions
are masked everywhere (conv inputs zeroed, profile logits set to -1e9 before softmax).

Model (`scripts/model.py`, `RiboSignalModel`): input 1x1 conv (d_emb+1 -> C) -> N residual
dilated conv1d blocks (dilations 1,2,...,512; GroupNorm + GELU) -> two heads. Profile head:
1x1 conv to per-nt logits. Count head: masked-mean body features + log1p(total coverage) ->
MLP -> scalar log-count. Default C=256, 10 blocks (~4.7M params at d_emb=1280).

Loss: `profile_multinomial_nll` = per-transcript multinomial NLL, `-sum_i (c_i/N) log
softmax(logits)_i`, equal-weighted across transcripts (scale-free); plus `count_weight`
(default 0.1) times `count_mse` = MSE of predicted log-count vs `log1p(total P-sites)`.

Training (`scripts/train.py`): AdamW, cosine schedule with warmup, grad-clip 1.0, early
stopping on val profile Pearson (patience 6). Evaluation is per-transcript and
scale-invariant: Pearson and Spearman between the predicted position distribution and the
observed profile `c/N`, plus 3-nt periodicity (autocorr lag 3/6/9/12 minus off-frame) of
the predicted and observed profiles, all reported for `all`, `protein_coding`, and
`lncRNA` separately (transcripts with < 50 total P-sites skipped in metrics). Runs inside
the CUDA SIF (`ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img`, torch 2.2.0+cu118 + numpy)
on one A5500 via `scripts/train_fold0.sbatch`, bind-mounting the group fs at the same path
so the scripts' absolute paths resolve read-write. The whole chain end-to-end (pack ->
loader -> model -> train -> eval -> checkpoint -> early stop) was validated on CPU against
a 1,500-transcript dry-run pack before launch.

Automated pipeline chain (SLURM `afterok` dependencies): coverage array -> pool ->
post-pool (`scripts/validate_coverage.py` then the full pack) -> `train_fold0.sbatch`, so
the remaining steps run unattended once the coverage fleet finishes.

Extra evaluation (`scripts/eval_extra.py`, `eval_extra.sbatch`): loads a run's `best.pt`
and, on its held-out test fold, adds two metrics the profile Pearson does not cover.
(1) Magnitude: per-transcript Pearson between the count head's predicted log-count and the
observed `log1p(total P-sites)` -- this is where the RNAseq coverage input is credited
(coverage feeds the count head). (2) Posture-B proxy: profile Pearson split by whether a
transcript's gene has a single annotated isoform (posture A == posture B, no multimap
ambiguity) vs multiple isoforms (which bear the multimapping) -- a rebuild-free check that
the posture-A inflation is not distorting shape. Writes `<run>/extra_metrics.json`.

---

## 5D. Multi-fold confirmation of the improvement matrix

The Task 8 improvement matrix was fold 0 only, where the lncRNA test set is n=264 (see the
funnel: 189,152 annotated lncRNA -> 2,420 expressed at TPM>=1 & len<=10k -> 583 on the fold-0
chromosomes -> 264 clearing the >=50 P-site floor). The +0.02 to +0.06 lncRNA deltas therefore
rest on a small, noisy median. This step extends every config to all 5 gene-disjoint chromosome
folds. Because each scorable transcript is held out in exactly one fold, pooling the 5 folds'
per-transcript Pearson gives a lncRNA median over ~1,239 transcripts (and ~33k pc), a firmer
basis for the deltas.

- **Cells trained** (`scripts/train_multifold.sbatch`, SLURM array 0-24, throttled to 12
  concurrent): the 5 non-baseline configs (orf, orf_v2, attn, orf_attn, orf_v2_attn) for folds
  1-4 (fold 0 REUSED from the existing `results/improve/*_f0` runs so the documented fold-0
  column stays exact), plus the baseline for folds 0-4. The baseline is retrained fresh in the
  improve lineage because the pre-existing `fold{1..4}_rinalmo` runs predate the `emb_backend`
  refactor (their fold-0 scores 0.604 vs the improve-matrix baseline `backend_cmp/rinalmo_f0`
  at 0.612) and are not apples-to-apples with the improve configs. `baseline_f0` also serves as
  a reproducibility check against `backend_cmp/rinalmo_f0`. All settings match that baseline
  exactly (emb_backend rinalmo, budget 24000 / 16000 for attention, seed 0, identical CNN
  hyperparameters); `orf_v2*` set `RIBO_ORF_TRACK` to the non-AUG-aware track.
- **Aggregation** (`scripts/aggregate_multifold.py`) reports two views per config, for pc AND
  lncRNA: (1) per-fold, the mean +/- std of the 5 per-fold `pearson_median` (standard k-fold CV
  reporting); (2) pooled, the median over the concatenated per-transcript `profile_pearson` from
  all 5 folds. The pooled view is the headline small-n fix. The per-transcript values come from
  `pertx.tsv`, which `eval_extra.py` now writes (same `MIN_SIGNAL=50` filter as training-time
  evaluation, so the pooled set matches the per-fold medians). `pertx.tsv` for the 5 reused
  fold-0 runs is backfilled eval-only by `scripts/backfill_pertx_f0.sbatch`.

## 5E. uORF-aware evaluation (5'UTR / 3'UTR profile Pearson)

The whole-transcript pc profile Pearson is dominated by the main CDS (the bulk of the P-sites),
so it barely reflects prediction of uORFs -- which are the most common class of alternative ORF.
To score uORF prediction directly, the evaluation restricts the profile Pearson to the 5'UTR
window `[0, cds_start)` (uORF region) and, for symmetry, the 3'UTR window (dORF region), among pc
transcripts whose window carries a translated ORF's signal.

- **Windows** come from `scripts/build_tx2cds.py` -> `data/tx2cds.tsv` (per pc transcript:
  utr5_len, cds_len, utr3_len, has_start_codon), derived from the GENCODE v49 GTF. v49 annotates
  undifferentiated `UTR`, so 5' vs 3' is resolved by strand relative to the CDS extent (5'UTR =
  exonic bases 5' of the CDS; verified on ACTB = 87 nt). `tx_len` matches the mature-transcript L
  exactly (0 mismatches vs `tx2biotype`), so the window aligns 1:1 to the target/embeddings.
- **Annotation-free training preserved:** the CDS coordinates are used ONLY to define eval
  windows. The model never sees the CDS boundary; it predicts a whole-transcript profile, and the
  metric asks how well it did inside the 5'UTR. A correct 5'UTR prediction is therefore genuine
  (the model was not told where the 5'UTR ends).
- **Scored set:** pc transcripts with a 5'-complete CDS, window length >= 30 nt, and >= 20 P-sites
  in the window (a translated uORF/dORF). ~74% of 5'-complete 5'UTRs clear the P-site floor, so
  the uORF metric pools ~22,000 transcripts (~4,400/fold) -- about 18x the lncRNA n, so it is far
  better powered. Caveat: some 5'UTR signal is scanning/spillover rather than uORF 80S elongation;
  the metric captures 5'UTR ribosome distribution, which uORFs dominate where present. A stricter
  variant restricted to sequence-defined uORF ATG..stop windows can be added if needed.
- **Implementation:** `eval_extra.py` computes per-transcript 5'UTR and 3'UTR profile Pearson,
  emits `uorf_5utr` / `dorf_3utr` median blocks in `extra_metrics.json`, and adds
  `utr5_pearson` / `utr3_pearson` columns to `pertx.tsv`. `aggregate_multifold.py` pools the 5'UTR
  Pearson across folds per config, so every improvement (ORF track, CTG, attention) is judged on
  uORF prediction, not just whole-transcript shape. Applied to already-trained models by re-eval
  (`reeval_uorf.sbatch` for the 6 core configs, `reeval_ctg_uorf.sbatch` for the CTG configs); no
  retraining.

## 5F. Leave-one-tissue-out (cross-tissue transfer)

The project's real objective (results.md Task 12; design `LOTO_PLAN.md`): train on several tissues,
test on a held-out tissue. Nine Chothani tissues share the fixed 36,668-tx universe -- the FM
embeddings and ORF track are sequence-based / tissue-independent, so they are reused with zero
re-extraction; per-tissue expression enters through the RNA-seq coverage input, not the transcript
set. Each tissue has its own pack (`data/packed_<Tissue>/`: target_counts + coverage + a per-tissue
`coverage_norm.json`); `tx_order`/`offsets`/`lengths`/`orf_track` are identical across tissues, so the
shared embeddings + ORF track align to every pack. Build chain: `build_psite_target.py` (pool the
tissue's RiboCode P-site hd5), `align_rnaseq_sample.sh` -> `pool_rnaseq_coverage_tissue.py` (STAR the
tissue's RNA-seq, per-nt coverage, pool), `pack_target_coverage_tissue.py`.

Split (`dataset.loto_split`): train = the held-in tissues concatenated (one PackedStore each ->
`ConcatDataset`, token-budget sampler over combined lengths), val = a gene-disjoint chromosome fold
across the training tissues (for honest early stopping), test = the held-out tissue. A transcript is
scorable in a tissue with >= 50 pooled P-sites there; each training tissue is capped
(`--max_tx_per_tissue 6000`) so deep tissues do not dominate. `cov_norm=global_mean` (each tissue by
its own mean depth) -- exactly the transferability fix from Task 10. chrM excluded (sec 7.1).
`train_loto.py` reuses the model, loss, and evaluate() from `train.py`; it writes a per-epoch `last.pt`
(model + optimizer + epoch) and supports `--resume` so a mid-run interruption (the first pilot died on
a filesystem-quota write failure) recovers instead of restarting. Eval reuses `eval_extra.py` (LOTO
branch: build the held-out tissue's store + `loto_split` test set); the consolidated table is built
from pertx by `loto_table.py`, and `robustness_concentration_loto.py` checks the pc number is not
spike-inflated. Backends compared via `--emb_backend {rinalmo,orthrus}` (packs are backend-agnostic).

## 5G. Localization eval (does the model concentrate signal on translated ORFs?)

Profile Pearson (whole-transcript or windowed) scores the *shape* of the predicted profile but never
asks the sharpest question: on a held-out transcript, does the model localize predicted Ribo-seq
signal to the specific ORFs that are genuinely translated -- especially non-canonical ones (uORF /
dORF / novel) -- rather than to the many untranslated candidate AUG ORFs on the same transcript?
`scripts/eval_localization.py` answers this against a RiboCode truth set.

- **Truth set:** RiboCode ORF calls, `<Tissue>_collapsed.txt` (from the sibling
  `expression_context_human` project, per tissue). Coordinate conventions verified against the
  universe FASTA: `ORF_tstart` is 1-based (the AUG sits at 0-based `[tstart-1, tstart+2)`, 100% of
  calls), `ORF_tstop` is 1-based inclusive of the last stop-codon base (stop at `[tstop-3, tstop)`,
  100%), and ORF length `tstop-(tstart-1)` is divisible by 3 (100%). RiboCode here was run
  **AUG-only** (13,445/13,445 in-universe calls start at AUG), so this is a translation-localization
  test, **not** a non-AUG-start discovery test -- the non-AUG information in the ext/ctg tracks is an
  input, but there is no non-AUG *truth* to score against (a non-AUG-enabled ORF caller is future
  work). Calls are richly typed: annotated (canonical CDS) vs non-canonical (uORF, Overlap_uORF,
  dORF, Overlap_dORF, novel, internal).
- **Candidates:** every maximal AUG..in-frame-stop ORF on a test transcript (the same construction
  as the orf_track occupancy channels; every in-frame AUG is a candidate start, so nested / internal
  AUGs are included), length >= 30 nt. A candidate is POSITIVE if its 0-based start matches a
  RiboCode start on that transcript (inheriting the RiboCode ORF type), else NEGATIVE. Critically the
  orf_track marks candidate positions + reading frames but NOT which candidates translate, so
  separating positives from negatives must come from the model (FM embedding + Ribo-seq-trained
  profile), not the track the model was handed.
- **Scores** (per candidate, from the predicted per-nt softmax profile `p` over the transcript):
  `density = mean_{i in body} p[i]` (per-nt predicted P-site density) and
  `frame0 = sum_{i in body, (i-start)%3==0} p[i] / sum_{i in body} p[i]` (periodicity signature: ~1/3
  if the model predicts no translation there, ->1 for a clean in-frame ORF), over the body incl. stop
  `[tstart-1, tstop)`.
- **Metric 1 (discrimination):** AUROC (rank-sum, ties averaged) of each score separating translated
  from untranslated candidates, overall and per ORF type, with an ORF-length AUROC floor and an
  observed-profile AUROC ceiling. Because canonical CDS are long, the length floor is high on the ALL
  stratum (annotated ORFs are trivially separable by length); the informative test is the
  **non-canonical stratum**, where length is uninformative, plus a **length-controlled** AUROC
  (`auroc_stratified`: Mann-Whitney within ORF-length quantile bins, so only translated-vs-untranslated
  pairs of similar length are compared) for both the model scores and the observed ceiling.
- **Metric 2 (fidelity, positives only):** median predicted vs observed `frame0` and the
  predicted-vs-observed `frame0` Pearson -- does the model reproduce the 3-nt periodicity it was never
  explicitly trained on (the loss is multinomial NLL on raw per-nt counts, with no periodicity term)?
- **Held-out only:** within-tissue runs use the fold's test transcripts + Fibroblast calls; LOTO runs
  use the held-out tissue's test transcripts + that tissue's own calls (the mature-mRNA sequence, and
  hence every candidate ORF, is tissue-independent, so one universe FASTA serves all tissues). Reuses
  `eval_extra.py`'s model load + batched inference; writes `<run>/localization_metrics.json` +
  `localization_orfs.tsv`; `loc_table.py` renders the cross-run comparison. Runners:
  `eval_localization.sbatch` (within-tissue fold-0 baseline + orf_v2_attn) and
  `eval_localization_loto.sbatch` (held-out-Hepatocytes LOTO models).
- **Per-candidate diagnostic columns:** `localization_orfs.tsv` also carries `tx_id`, `orf_start`,
  `orf_end`, and `cds_rel` -- each candidate's relationship to the annotated CDS (`data/tx2cds.tsv`),
  classified by where its start sits relative to the CDS and whether it shares the CDS frame:
  `cds_canonical` (start == CDS start), `cds_inframe` (in-frame start strictly inside the CDS: an
  internal alt start / 5'-truncation), `dorf_inframe` (in-frame start at or past the CDS stop),
  `ext_inframe` (in-frame start 5' of the CDS start that still overlaps it: an N-terminal extension),
  `cds_offframe` (out-of-frame CDS overlap), `utr5`, `utr3`, `no_cds` (lncRNA). This resolves why the
  untranslated-candidate frame0 distribution is bimodal: the high cluster (predicted frame0 > 0.7) is
  ~94% `cds_inframe` -- internal in-frame AUGs sharing a translated CDS's frame (RiboCode collapses each
  locus to one representative call, so the internal fragments are labelled untranslated yet sit in real
  periodic translation), which the model correctly predicts as periodic.
  `scripts/make_untranslated_decomp_figure.py` visualizes it.
- **Internal in-frame alt-start AUG ORFs are collapsed by default** (`--no_collapse_inframe_cds` to
  disable): an AUG candidate that starts strictly inside the annotated CDS in the CDS reading frame
  (`cds_inframe`) is a 5'-truncation re-initiating the same protein, not a distinct ORF, so it is dropped
  from the metric candidate set. This is the ONLY collapsed class. Everything else is kept as a
  candidate: in-frame dORFs at or past the stop (`dorf_inframe`), out-of-frame CDS overlaps
  (`cds_offframe`), in-frame N-terminal extensions (`ext_inframe`), and all UTR / lncRNA ORFs. The one
  kept exception inside `cds_inframe` is a candidate whose start coincides with an annotated RiboCode
  call (144 in the within-Fibroblast fold: a RiboCode-vs-`tx2cds` CDS-start coordinate offset lands a
  handful of annotated calls at an internal in-frame position, and those are real translation that must
  not be dropped). This removes ~62k internal in-frame candidates from the within-Fibroblast fold and 0
  non-canonical positives -- every non-canonical positive is a uORF / dORF / novel / off-frame ORF, none
  is an internal in-frame alt start. The full candidate list is still dumped to the tsv so the
  decomposition diagnostic survives. **dORFs (in-frame past the stop and off-frame overlapping) and
  out-of-frame overlapping ORFs are genuine distinct ORFs and are kept** -- only internal, in-frame alt
  starts are excluded. **The collapse is also AUG-specific, and this matters for alt-ORF detection.** An
  in-frame *non-AUG* start 5' of the annotated AUG is not a truncation but a genuine N-terminal EXTENSION
  via alternative initiation -- the CUG-initiated extensions of MYC (L-Myc/p67) and VEGFA (L-VEGF) are
  the textbook cases, and near-cognate initiation happens at meaningful rates precisely at such sites. So
  an in-frame CUG extension must be kept as a distinct candidate, never collapsed into the CDS.
  Candidates are AUG-only for now (RiboCode was run AUG-only), so the collapse is safe today; when a
  non-AUG-enabled caller is wired in, the collapse must stay gated on an AUG start codon, and in-frame
  near-cognate extensions become first-class alt-ORF candidates (harder, and exactly the biology this
  project targets).

## 5H. RiboCode drop-in test (do the predicted profiles reproduce the ORF caller's own calls?)

The localization eval (5G) scores the model's periodicity with a custom AUROC. A stronger, more concrete
question: if the model's predicted per-nt profile is fed into the ACTUAL ORF caller (RiboCode) in place
of the experimental Ribo-seq, do the ORF calls track the calls RiboCode makes on the real data? This is
the "predicted profile as a drop-in replacement for the experiment" test.

- **Injection point:** RiboCode's per-transcript P-site density is a plain `dict {transcript_id:
  np.ndarray}` (0-based, length == transcript length), built from the BAM by `process_bam.psites_count`
  and consumed by `detectORF.main` (RiboCode v1.2.15, `conda_envs/ribocode`). Because it is just a
  positional argument, `scripts/ribocode_dropin.py` bypasses BAM reading entirely: it builds that dict
  from a predicted profile and calls `detectORF.main` directly, reusing the prebuilt annotation DB
  (`expression_context_human/data/ribocode_annot`, the same GENCODE v49 universe -- 509,650 tx, versioned
  ENST, lengths verified to match the model universe 1:1) and RiboCode's default parameters (start ATG,
  min_AA 5, pval 0.05, longest_orf, stouffer none, fdr_bh). The identical statistical machinery
  (Wilcoxon f0>f1 and f0>f2, Stouffer-combined, ORF typing, collapsing) runs; only the density array is
  swapped.
- **Profiles:** `scripts/dump_pred_profiles.py` reruns the held-out-Hepatocytes orf_v2_attn LOTO model
  (the exact 5G inference) and dumps, per test transcript, the predicted softmax profile (the SHAPE),
  the real pooled P-site counts (the reference), and the count head's predicted total `expm1(pred_
  logcount)` (the predicted DEPTH), flat-packed into `pred_profiles.npz` (numpy-only; the rinalmo SIF
  lacks h5py).
- **Three density variants, identical caller + params, ONLY the density differs:** `real` = the real
  pooled counts (reference; validated to reproduce the official `Hepatocytes_collapsed.txt`);
  `pred_obsdepth` = predicted shape x real per-tx total, rounded to int (isolates the profile shape at
  correct depth); `pred_preddepth` = predicted shape x count-head predicted total, rounded (fully
  standalone -- no experimental data used at all). **Two magnitude rules matter:** RiboCode hard-gates
  per-tx sum>=5, per-ORF frame0 sum>=5, and >=5 nonzero frame0 codons (`detectORF.py:306`), so the
  density must be scaled to realistic depth; and the pred profiles are rounded to integer counts because
  a smooth density has almost no exact zeros and would pass RiboCode's per-codon coverage/tie logic more
  easily than real integer counts, looking artificially more periodic.
- **Metric:** `scripts/compare_dropin_calls.py` matches calls on the **genomic ORF locus**
  `(gene_id, ORF_gstop)` by default (`--key genomic`), restricted to the model's test-tx genes, filtered
  on raw `pval_combined <= 0.05` and `--min_len` (default **90 nt / 30 codons**).
  **FDR posture (matched ORFs are BH-controlled, matched after BH filtering):** RiboCode already applies
  Benjamini-Hochberg (`pval_adj=fdr_bh`, cutoff 0.05) WITHIN each run at call generation, so every ORF in
  every `_collapsed.txt` (real, both predicted variants, official) is already FDR-controlled -- verified
  `adjusted_pval <= 0.05` for all rows (max 0.0499). Matching is therefore a set-intersection over
  BH-passed survivors; multiple-testing correction lives in the caller (where it belongs), and F1 of
  set overlap is not itself a hypothesis test, so there is no second BH to apply in the compare. The
  compare re-filters on the RAW per-ORF `pval_combined`, NOT `adjusted_pval`, as a deliberate near-no-op:
  it drops nothing on top of RiboCode's gate (adjusted <= 0.05 implies raw <= 0.05), but keeping it raw
  makes it tx-set-independent, so any tightening of `--max_pval` below RiboCode's gate applies an identical
  bar to every set. `adjusted_pval` must NOT be the compare threshold because its BH burden scales with the
  tx-set size (34k held-out vs 509k full annotation), so the same ORF gets a different adjusted value per
  set. Consequently the primary `pred_*_vs_real` comparisons (both run over the identical held-out universe)
  have matched FDR burden that cancels; only the `*_vs_official` comparisons carry a residual universe-size
  asymmetry, hence their harness-validation framing. Genomic keying is essential because RiboCode's per-locus collapse
  picks the longest ORF per genomic stop (`detectORF.py:440-442`, density only breaking exact-start
  ties), so restricting the transcript universe merely reshuffles which isoform carries a call;
  `(transcript_id, ORF_tstart)` matching (`--key transcript`) is available but understates every number
  by penalizing representative-isoform disagreements (it is why the isoform-level real-vs-official
  precision is only 0.61 while the genomic-locus precision is 0.99 -- the same calls, different
  representative isoform). The 90 nt floor drops the sub-90 nt regime where the predicted density
  over-calls (a smooth softmax never zeros the 3'UTR, and depth-scaling turns tiny 3'UTR probabilities
  into callable counts; 37% of the pre-floor dORF false positives sit in 3'UTR windows with zero real
  P-sites). `--min_enrichment` (default **0.5**) additionally drops PREDICTED calls whose mean predicted
  density over the ORF is below 0.5x uniform (mean probability x transcript length) -- the diffuse 3'UTR
  softmax leak that yields the spurious dORFs; the real/official reference calls are never
  enrichment-filtered. The uniform-enrichment score separates genuine uORF from spurious dORF at AUROC
  0.947 (a CDS-relative denominator is worse, 0.78-0.84), so the operating point is set purely by the
  threshold; 0.5 was chosen on the held-out Hepatocytes sweep (`ribocode_dropin.py --floor_mult`
  per-position density floor + `analyze_floor_sweep.py`) and should be re-validated per tissue. It reports
  precision / recall / F1 of the predicted calls against the real calls (overall +
  recall by ORF type), the real-vs-official overlap (harness validation), and on matched calls the
  agreement of RiboCode's own confidence (-log10 pval_combined), frame0 magnitude, and ORF_type. Runners:
  `dump_pred_profiles.sbatch` (GPU), `ribocode_dropin.sbatch` (+`ribocode_dropin_floor.sbatch` for the
  floor sweep), then `compare_dropin_calls.py`.

## 5I. External held-out dataset (Ruiz-Orera 2024 human iPSC-CM)

Motivation: every eval to this point (fold-0 held-out chromosome, LOTO, the drop-in) stays within the
Chothani smORF-atlas data. The remaining validation gap is a fully independent dataset from a different
lab and protocol, to test whether the model and the drop-in operating point transfer across studies.

Dataset: Ruiz-Orera et al. 2024, Nature Cardiovascular Research, "Evolution of translational control and
the emergence of genes and open reading frames in human and non-human primate hearts"
(DOI 10.1038/s44161-024-00544-7; PMC11473369; code github.com/jorruior/RuizOrera_etal_2024). Archive
ENA PRJEB65856 (a 4-species iPSC-CM + primate-LV panel: human, chimp, macaque, gorilla). Standard
monosome Ribo-seq, single-end 1x51 with a clipped 3' adapter, matched paired-end mRNA-seq. The paper
reports 82 to 84 percent of 29-nt RPFs in the primary reading frame (strong 3-nt periodicity; the only
candidate surveyed with periodicity confirmed in-paper rather than inferred from depth).

Cell type and held-out status: the human portion is iPSC-derived cardiomyocytes (line BIHi242, 5 clones)
plus an iPSC baseline; it does NOT contain the adult left-ventricle tissue (human LV in this study is EGA
gated; only chimp and macaque LV are public here). Cardiomyocyte is absent from the 9-tissue training
panel (Brain, ES, Fat, Fibroblast, HA_EC, HCAEC, Hepatocytes, HUVEC, VSMC), so it is a genuinely
out-of-distribution cell type; the iPSC baseline overlaps the training "ES" tissue conceptually and can
serve as a sanity anchor (same role as Day-0 in the hESC differentiation sets).

Accession correction (source-study provenance): the sibling project doc
`expression_context_human/methods.md` (Phase 28) cites `GSE183670` as the Chothani source. That GSE is
actually "Transcriptomic profiling of oncogenic KRAS and TP53 in lung adenocarcinoma cells" (Berger,
PMID 36522653), an unrelated RNA-only cell-line study. The Chothani smORF-atlas Ribo-seq runs actually
used here (`SRR15513xxx`; verified `SRR15513165` = "Fibroblast_21 (RiboSeq)") belong to
PRJNA756018 / SRP333128 ("smORF atlas part 1-RiboSeq"), GEO GSE182371 (Ribo) + GSE182372 (RNA),
SuperSeries GSE182377. The training DATA is correctly Chothani; only the GSE label in the sibling doc is
wrong. Flagged here; the sibling doc is not edited without explicit instruction.

What was pulled (human subset only): cross-species arms were deliberately skipped to avoid non-human
primate transcript-annotation overhead for the FM embeddings; cross-species generalization will instead
use mouse (vM38 annotation + embedding pipeline already exist). Pulled 21 files, 44 GB:
- 5 Ribo-seq (SE, iPSC-CM clones C1 to C5): ERR12549926 to ERR12549930.
- 5 matched RNA-seq (PE, same clones C1 to C5): ERR12549891 to ERR12549895.
- 3 iPSC RNA-seq (PE, pluripotent baseline): ERR12549905 to ERR12549907.
The usable held-out is the 5 clone-matched iPSC-CM Ribo+RNA pairs.

Provenance and paths: `data/external/ruizorera2024/`. Serial single-stream head-node pull
(`download.sh`, `wget -c` one file at a time per the project download rule), each file md5-verified
against the ENA filereport checksum (`logs/verify_manifest_human.tsv`: ok=21, fail=0). `manifest_human.tsv`
is the pulled subset; `manifest_all.tsv` retains the full 4-species archive listing for reference;
`manifest_raw_ena.tsv` is the raw ENA filereport.

Next steps for this held-out: (1) adapter + read-length + periodicity QC on the RPFs (check for residual
3' adapter per the project Ribo-seq rule); (2) STAR align with `--outFilterMultimapNmax 1` and drop
rRNA/tRNA/miRNA/Mt loci (project Ribo-seq rule) against GRCh38 + GENCODE v49, transcriptome-coordinate
P-sites via the same RiboCode config path as the Chothani tissues; (3) run the RiboCode drop-in (section
5H) and per-nt model prediction against these calls as a true external, out-of-distribution held-out.

## 5J. Cross-species footprint held-out (mouse liver, Wang 2021) + held-out processing pipeline

Motivation: a mouse cross-species held-out (human-trained model applied to mouse), chosen over non-human
primate because the GRCm39 / GENCODE vM38 annotation + mouse FM-embedding + RiboCode-annot pipeline already
exist (no new-annotation overhead).

### Dataset: Wang et al. 2021 mouse liver (GSE94982)
Wang et al. 2021 Nucleic Acids Research, "A spatio-temporal landscape of the mouse translatome"
(GEO GSE94982 / SRA PRJNA375080). Genuine CHX-arrested monosome ribosome FOOTPRINTING, verified by read
inspection (short ~17-30 nt inserts followed by the TruSeq adapter `AGATCGGAAGAGCACACGTCT` in ~92% of
reads), paper-reported ~79% in-frame periodicity, with matched mRNA-seq. Adult (P42) LIVER is the closest
analog to the human "Hepatocytes" training tissue, isolating species as the main variable. Pulled the P42
liver arm: 2 RPF (SRR5262890/91) + 2 matched mRNA (SRR5262874/75), 8.1 GB, md5-verified. Path
`data/external/wang2021_mouse/`. Brain + heart RPF+mRNA exist in the same study for a multi-tissue mouse
held-out if wanted (not pulled).

Rejected dataset (recorded so it is not re-tried): GSE157050 / PRJNA659995 ("Ribosome profiling of
metformin in primary murine hepatocytes") was pulled first but is actually POLYSOME PROFILING, not
footprinting (GEO extract protocol: sucrose-gradient polysome fractionation + TruSeq mRNA-seq; 51 nt
full-length reads, no footprint adapter; the GEO library-strategy field is mislabeled "Ribosome
profiling"). No P-site periodicity, cannot feed the per-nt model. It was briefly retained as a
translation-efficiency (polysome/input) substrate for a treatment-response layer, but that layer was
DROPPED (user decision 2026-07-13, "not doing TE") and the 34 GB dataset was DELETED. Lesson: verify
footprint-vs-polysome/TRAP by read structure + GEO protocol before pulling, not by the title.

### Human cross-STUDY held-out (section 5I) runs the same pipeline
The human held-out (section 5I, Ruiz-Orera PRJEB65856 iPSC-CM) is genuine footprint in the SAME format as
the Chothani training reads (35 nt, N-padded, pre-clipped), so it runs through the identical pipeline below.

### Processing pipeline (`scripts/heldout/`, posture A -- matches the Chothani training P-sites)
Replicates expression_context_human phase28 (align) + phase30 (filter) + RiboCode, parameterized per
dataset in `scripts/heldout/heldout_config.sh` (datasets `human_ruizorera`, `mouse_wang_liver`):
1. `align_and_filter.sh` (array `align_and_filter.sbatch`, medium/16cpu/40G, one task per Ribo SRR):
   cutadapt (human `--trim-n --minimum-length 20` on the already-clipped reads; mouse `-a AGATCGGAAGAGC
   --trim-n --minimum-length 20 --maximum-length 40`) -> STAR `--quantMode TranscriptomeSAM
   --outFilterMultimapNmax 20 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNminOverLread 0.7
   --alignEndsType EndToEnd` (v49 index human / vM38 index mouse) -> `filter_tx_heldout.py` drops ncRNA-tx
   records + cross-gene multimappers, keeps within-gene isoform multimappers (posture A), sorts + indexes
   the transcriptome BAM in place. Output `data/heldout_bam/<dataset>/<SRR>.Aligned.toTranscriptome.out.bam`.
2. `ribocode_heldout.sh` (`ribocode_heldout.sbatch`, chained via `--dependency=afterok`): pooled `metaplots`
   (per-read-length P-site offsets over the dataset's BAMs) -> `RiboCode -l no -g` -> per-nt
   `<SRR>...toTranscriptome.out_psites.hd5` (same vlen schema as the Chothani training P-sites: transcript_ids
   + p_sites int32, 0-based, length == transcript length) + `<dataset>_collapsed.txt` ORF calls, in
   `data/heldout_psites/<dataset>/`.
Mouse cross-gene filter needs a vM38 tx->gene map, built once to `data/heldout_refs/mouse_tx_to_gene.tsv`
(278,326 tx from the vM38 GTF). Reference infra (both species, on disk, no building): STAR indexes
`genomes/star_index_{grch38_v49,grcm39_vM38}`; RiboCode annot `ribocode_annot` + `ribocode_annot_mouse`;
ncRNA-tx lists `phase30/ncrna_filter_tx.txt` (human, 7,587) + `phase31/mouse_ncrna_tx.txt` (mouse, 2,579).

Launched 2026-07-13: human `human_ruizorera` align array 35250015 (5 tasks) -> RiboCode 35250016
(dependency). Mouse `mouse_wang_liver` runs the same two-step chain once its FASTQs finish downloading.
FASTQ deletion (standing rule) is deferred until the psites.hd5 are validated on the first run. Next:
cross-study drop-in (section 5H) + per-nt prediction on both held-out sets; RNA-seq per-nt coverage for the
matched mRNA is a separate track.

## 5K. Held-out pack + eval implementation (cross-study + cross-species, 2026-07-14)

The external held-out P-sites (section 5I/5J) + matched RNA-seq coverage are pooled into a
model-input pack and scored by the SAME three evals as the LOTO models (profile Pearson,
localization AUROC, RiboCode drop-in), with a Chothani-trained orf_v2_attn model applied
without retraining.

- **Held-out pack** (`scripts/heldout/build_heldout_pack.py`, runner `build_heldout_pack.sbatch`):
  pools an external dataset's per-SRR RiboCode P-sites + per-SRR RNAseq coverage directly into
  a pack (posture A), joining strictly by versioned tx id (coverage hd5 @SQ order differs from the
  RiboCode order), WITHOUT writing the 3.7 GB pooled-hd5 intermediates (ceph is tight). Two universe
  modes: (a) human Ruiz-Orera reuses the Fibroblast pack's tx_order / offsets / lengths verbatim, so
  the shared `orf_track_v2.npy`, the one-hot FASTA, and the RiNALMo/Orthrus per-token embeddings all
  align 1:1 and ALL THREE backends run; the build asserts every pooled per-nt length equals the
  reference length (v49 annotation match). (b) mouse Wang uses a species-specific universe with fresh
  offsets/lengths, one-hot backend only (no mouse FM embeddings -- the point of the one-hot control).
  Human pack: 36,668 tx, 33,185 scorable (>=50 P-sites), global-mean depth 743 (vs Fibroblast 7,733,
  so ~10x shallower -- exactly why the `global_mean` coverage normalization is needed for transfer).
  Mouse pack: 221,835 tx, 41,096 scorable, global-mean depth 43.9.
- **Mouse universe** (`scripts/heldout/build_mouse_universe.py`): all vM38 protein_coding + lncRNA
  transcripts (from the two vM38 transcript FASTAs) with length <= 10,000 nt present on the Wang
  psite+coverage axis = 221,835 tx (65,936 pc + 155,899 lncRNA). No expression prefilter -- the
  >=50-P-site scorable floor defines the scored set; untranslated tx are harmless. Writes
  `data/heldout_refs/mouse_wang_universe.{txt,fa,tsv}`.
- **Mouse ORF track**: `build_orf_track.py` gained `--pack` / `--fasta` overrides; the mouse
  `orf_track_v2.npy` (ext mode, near-cognate) is built on the mouse pack's tx_order + mouse FASTA
  (2.9 GB, 5 channels, aligned to the mouse pack sum_L 291,028,013).
- **Mouse CDS table**: `build_tx2cds.py` gained `--gtf` / `--out`; `data/mouse_tx2cds.tsv` (66,618
  coding tx) is parsed from the vM38 GTF with the identical logic, so the mouse localization eval's
  in-frame-CDS collapse + cds_rel diagnostic mirror the human eval. `load_cds_bounds()` reads
  `RIBO_TX2CDS` to select it.
- **Held-out eval mode** added to `dump_pred_profiles.py` + `eval_localization.py` via `--heldout
  <name>` (+ `--min_signal`, `--fasta`, `--out`): the pack becomes `packed_heldout_<name>`, the test
  set is every scorable tx (no train/val split), and the localization truth calls come from
  `data/heldout_psites/<name>/<name>_collapsed.txt`. `dataset.py` gained `heldout_pack_dir` /
  `heldout_test_tx` and an `RIBO_ONEHOT_FASTA` env override so the one-hot backend runs on a
  cross-species FASTA. The trained run's `cov_norm=global_mean` uses the held-out pack's OWN mean depth
  (from its `coverage_norm.json`), which is the transferability mechanism. `RIBO_ORF_TRACK` points at
  the shared `orf_track_v2.npy` for human (offsets identical) or the mouse track for mouse.
- **Runners** (`scripts/heldout/`): `eval_heldout_predict.sbatch` (SIF/GPU: dump profiles +
  localization, one array task per backend, env ORFTRACK/FASTA_OVERRIDE/TX2CDS for cross-species) and
  `eval_heldout_dropin.sbatch` (ribocode env: 3 density variants through detectORF + compare against
  the held-out calls, ANNOT=ribocode_annot for human / ribocode_annot_mouse for mouse).
  `assemble_heldout.py` collates the three-metric table (profile Pearson computed directly from
  `pred_profiles.npz`). Launched 2026-07-14: human predict/dropin 35314728/35314729 (3 backends),
  mouse orf_track->predict->dropin 35314742/35314746/35314747 (one-hot).

## 5L. Cross-study landing, architecture ablations, and the bidirectional Mamba mixer (2026-07-15)

- **Cross-study (Ruiz-Orera human) landed** (results.md Task 16). The predicted profile transfers to
  an unseen human study: one-hot protein-coding profile Pearson 0.425 (ties/beats Orthrus/RiNALMo),
  drop-in F1 ~0.93 (one-hot best 0.931), and the fully standalone `pred_preddepth` path (predicted
  shape x count-head total, no target Ribo-seq) matches the `pred_obsdepth` diagnostic (F1 gap ~0:
  +0.002/-0.001/-0.002 for onehot/orthrus/rinalmo). Notably the predicted-profile F1 (0.93) EXCEEDS
  the real shallow held-out profile vs the official deep calls (0.878) -- the model denoises. This
  **resolves the count-head transferability decision**: the absolute-count head transfers fine for
  ORF calling, so the TE / per-million refactor is NOT needed for single-dataset deployment; it only
  matters for coherent multi-dataset TRAINING (`design_count_magnitude_transferability.md` Section 9).
  Mouse cross-SPECIES drop-in (Wang liver, vM38, one-hot) LANDED 2026-07-15: pred_obsdepth F1 0.929 /
  pred_preddepth 0.919 (gap -0.010), real_vs_official 0.981 (harness validated) -- the human-trained
  model recovers 92-93% of RiboCode's own mouse calls from predicted profiles alone, ~= the human
  cross-study 0.931, and the count-head resolution holds across species. This required fixing a
  hardcoded human tx->gene map in `compare_dropin_calls.py` (added `--tx2gene` + a loud 0-match guard;
  built `data/mouse_tx2biotype.tsv` from vM38 via the now-parameterized `build_tx2biotype.py --gtf/--out`).
- **AMENDED 2026-08-08.** Everything in the bullet above was measured on
  `orf_v2_*_onehot_holdout_Hepatocytes`, not on the shipping union models, and the staleness was
  invisible because run directories are named for the DATASET rather than the CHECKPOINT. Re-dumped both
  released models over Ruiz-Orera and three mouse-liver studies and re-scored identically. The
  `pred_obsdepth` (shape) numbers survive: Ruiz-Orera 0.929/0.934, Wang 0.923/0.927 for attn/mamba4,
  against 0.931 and 0.929 before. The **`pred_preddepth` (standalone) numbers do not**: Ruiz-Orera
  0.876/0.880 and Wang 0.867/0.867, against 0.930 and 0.919. The obsdepth-to-preddepth gap is therefore
  NOT ~0; it is -0.024 to -0.062. The loss is pure precision (Ruiz-Orera 0.942 -> 0.832) with recall
  unchanged or higher, i.e. the count head over-calls on the 84,472-tx union universe. The count-head
  transferability *decision* stands -- the head still transfers across study and species without a TE
  refactor -- but the Ribo-seq-free Poisson calibration is load-bearing for the standalone arm rather
  than an optional refinement. Two consequences for method text elsewhere: quote standalone numbers only
  with the calibration arm stated, and read per-dataset numbers from the generated tutorial pages
  (`tutorial/make_heldout_{human,mouse}.py`), which print the source checkpoint next to every row.
- **Input-modality ablation** (results.md Task 17A), Hepatocytes hold-out, converged held-out test-set
  medians. `--input_mode` in {both, emb (sequence-only), cov (RNA-seq-only)}. Confirms the dual-head
  design premise measured cross-tissue for the first time: sequence-only recovers 91% of the
  both-modality pc profile Pearson (0.579 of 0.639) and reproduces (sharpens) the 3-nt periodicity
  (+0.471 vs +0.391), while RNA-seq-only collapses to 0.122 with periodicity destroyed (-0.011). So
  the profile SHAPE (and all periodicity) is sequence-borne; RNA-seq supplies MAGNITUDE (via the count
  head) plus a +0.060 whole-transcript-envelope lift on top of sequence.
- **Repeated on the DEPLOYED union recipe** (results.md Task 61, 2026-08-08), because the numbers above
  are from the pre-union recipe and the claim they support is load-bearing for the proteogenomics
  argument. `scripts/train_union_inputablation.sbatch` clones `train_loto_union.sbatch` and changes only
  `--input_mode`; the `both` arm is the deployed run itself, so the contrast carries no extra seed noise.
  The dissociation reproduces and SHARPENS: sequence-only now recovers **98.5%** of the pc profile
  Pearson (0.6601 of 0.6699, was 91%), periodicity again sharpens without RNA-seq (0.3704 vs 0.3359),
  RNA-seq-only still collapses to 0.1226 with periodicity destroyed (-0.0155), and the count head loses
  **-0.193** Pearson (0.8990 -> 0.7063) when RNA-seq is removed. Consequence for how the claim is
  worded: the model's cell-type specificity runs through the COUNT HEAD (which ORFs clear the caller's
  depth threshold), not through the profile shape, which is nearly cell-type-invariant. Quote -0.193,
  never a shape number, when supporting a cell-type-specificity claim.
- **Capacity ablation** (results.md Task 17B). Neither deeper (4 attn layers: pc Pearson 0.603) nor
  wider (ch 384: 0.633) beats the baseline (0.639); deeper HURTS whole-tx Pearson while sharpening
  periodicity (+0.430). The 2-attn/256-ch body is at the capacity sweet spot for whole-tx shape. (This
  revised an earlier mid-flight val-Pearson read that had wrongly reported attn=4 as a small lift.)
  DROP-IN TEST of attn=4 (jobs 35357991-993): the +0.039 periodicity DOES convert to calling gains on the
  non-canonical strata -- non-canonical recall +0.064 (0.614->0.678), uORF recall +0.120 (0.689->0.809) --
  but at a precision cost (-0.042), so aggregate drop-in F1 is slightly lower (0.904 vs 0.922). So whole-tx
  Pearson proxies aggregate F1 well (both favor 2-layer) but MASKS the recall/precision tradeoff: 4-layer is
  the better uORF / non-canonical DISCOVERY model, 2-layer the better-calibrated general caller.
- **Bidirectional Mamba mixer** (results.md Task 17C; `manuscript/related_work.md`). `model.py` gained
  `BiMambaBlock` (forward Mamba + reversed-sequence Mamba, summed, pre-norm residual) and `MambaBody`;
  `RiboSignalModel(mixer="transformer"|"mamba", mamba_d_state/d_conv/expand)` branches the `n_attn_layers>0`
  block, leaving the transformer path's state-dict keys unchanged (backward compatible). `train_loto.py`
  exposes `--mixer mamba` (+ `--mamba_d_state/d_conv/expand`); all four model-reconstruction sites
  (`eval_extra.py`, `dump_pred_profiles.py`, `eval_localization.py`, `plot_example_uorf.py`) read `mixer`
  from `args.json`. Adapted from the seq2ribo polisher (Kaynar & Kingsford 2026) but made BIDIRECTIONAL:
  that model is causal only because it refines an sTASEP simulation that already carries downstream
  context, which this model lacks. Trained in the Orthrus SIF (mamba_ssm 1.2.0.post1). 5,768,962 params.
  CONVERGED (job 35349301, early-stopped ~epoch 9): held-out test pc Pearson 0.599 vs the transformer
  anchor's 0.639 (-0.040), but higher periodicity (0.458 vs 0.391) -- the same Pearson-for-periodicity
  tradeoff as deeper attention (17B). Competitive, not better, on this hold-out; O(L) vs O(L^2) is the
  standing edge, so keep it as an option for very long transcripts, not the default. (The 1-epoch smoke's
  periodicity 0.623 was an undertraining artifact; it settled to 0.458.)
- **Full 9-tissue LOTO** (results.md Task 17D). Every Chothani tissue held out once (previously only
  Hepatocytes), one-hot `orf_v2_attn`. 9/9 DONE, mean pc Pearson 0.505 (0.539 excl Brain); range
  Hepatocytes 0.638 / Fibroblast 0.629 down to Brain 0.234. KEY FINDING: the cross-tissue spread is a
  held-out-TARGET data-quality ceiling, not a generalization gradient -- the model's predicted periodicity
  is ~constant (0.30-0.45) across tissues, but pc Pearson broadly tracks the OBSERVED periodicity of the
  held-out data (period_OBS 0.044-0.260, strong at the low end). Brain is the extreme (period_OBS 0.044,
  near noise -> pc capped at 0.234), same shallow-data ceiling as Task 16.
- **New scripts**: `scripts/train_loto_ablation.sbatch` (parameterized modality/capacity ablation via
  `--export`), `scripts/train_loto_mamba.sbatch` (mamba variant), `scripts/heldout/dump_finalize_gpu.sbatch`
  (self-contained un-sharded dump + drop-in, replacing the timeout-prone sharded chain),
  `scripts/heldout/mouse_annot_and_dropin.sbatch` (vM38 prepare_transcripts + drop-in).

## 5M. Depth crossover: predict vs measure (2026-07-16, results.md Task 18)

Concrete, quantitative form of the 5L / Task 17D "target-quality ceiling" claim and the Task 16 denoising
finding. Binomially thin the REAL Hepatocytes P-site profile to a fraction f of depth (`ribocode_dropin.py
--subsample f --seed s`, each footprint kept i.i.d. w.p. f; array runner `subsample_depth.sbatch`), run
RiboCode `detectORF` on the thinned counts, score vs the deep OFFICIAL calls (genomic key). The MEASUREMENT
curve F1(f) is compared to the depth-INDEPENDENT prediction line `pred_preddepth` (predicted shape x
count-head predicted depth; RNA-seq + sequence, no Ribo-seq). `subsample_depth_curve.py` aggregates +
interpolates the crossover; `plot_depth_crossover.py` renders `figures/depth_crossover/depth_crossover.png`
(+ `FIGURE_DATA_INPUTS.md`). RESULT: standalone prediction F1 0.818; measurement crosses it at ~2.9e7
test-set P-sites (f~0.047) -- below ~29M P-sites on these 33,918 tx, predicting from RNA-seq beats measuring
Ribo-seq. At 6.3e6 P-sites (normal-depth library) measuring is 0.733 vs predicting 0.818. Non-canonical
recall crosses at ~4e7. So a Pearson/F1 measured against shallow Ribo-seq is a MEASUREMENT ceiling, not a
model ceiling. (`--subsample`/`--seed` added to `ribocode_dropin.py`, backward compatible: 1.0 = off.)

## 5M2. Non-AUG (ATG + CTG) RiboCode drop-in (2026-07-16, results.md Task 19)

Extends the Task 14 drop-in (section 5H) to non-AUG initiation: can the model's predicted profile help
call CTG-initiated ORFs, or only canonical ATG? The predicted / observed per-nt profiles already dumped
for the drop-in (`dump_pred_profiles.py` -> `pred_profiles.npz`) are reused with NO re-inference;
`ribocode_dropin.py --alt_start_codons CTG` re-runs RiboCode `detectORF` with
ALTERNATIVE_START_CODON_LIST=["CTG"] on each density variant (real / predicted-at-observed-depth /
predicted-at-predicted-depth), writing to a separate `dropin_ctg/` dir. RiboCode treats alternative
starts as FALLBACK-only: a CTG opens an ORF only where its in-frame stop has no ATG, so CTG calls are a
DISJOINT set on top of the unchanged ATG calls (ATG F1 is therefore an invariant sanity check).

Comparison (`compare_dropin_ctg.py`) is start-codon-stratified on the genomic-locus key (gene_id,
ORF_gstop -- invariant to isoform collapse), with the standard filters (min_len 90, raw pval_combined
<= 0.05, predicted enrichment >= 0.5). It reports strict precision/recall/F1 SEPARATELY for ATG- and
CTG-initiated ORFs, plus a locus-level recall/precision for CTG (where the collapse representative may
record a different codon). Start codons are looked up from the RiboCode annotation FASTA.

Generalization across backend / dataset / species / architecture (`ribocode_dropin_ctg.sbatch` per run,
consolidated by `aggregate_ctg.py` -> `results/ctg_across_runs.json`): the same drop-in is run on the
Hepatocytes deployment model (onehot / rinalmo / orthrus), the independent human Ruiz-Orera set, the
cross-species mouse Wang set, and the attn4 / mamba mixer variants, so CTG-detection quality can be read
against the fixed ATG-F1 line on every axis. `kozak_context_alt_orfs.py` separately checks whether the
Kozak context of the real CTG calls differs from ATG calls (it is a uORF-vs-CDS effect, not a codon
effect -- feeds the Task 20 pre-registration).

## 5N. Kozak start-context ablation (remove / empirical / learned) (2026-07-16, results.md Task 20)

The ORF-track `ext` (v2) start-propensity channel (channel 3, `build_orf_track.py`) multiplies a
per-codon start weight by a hand-picked Kozak factor `kz = 0.5*[purine at -3] + 0.5*[G at +4]`. The
`{0.5,0.5}` weights at exactly -3/+4 are a heuristic. This experiment asks three questions and isolates
the start channel by changing ONLY channel 3 (occupancy channels 0-2 and the stop channel 4 are ATG-based
and byte-identical across all arms): (1) how do metrics hold if the heuristic is removed, (2) does an
empirically-fit Kozak beat it, (3) can the model learn the context weights itself. All arms are one-hot,
trained on Fibroblast (`data/packed/`) and tested on the held-out Hepatocytes tissue
(`data/packed_Hepatocytes/`, offsets md5-identical so a new orf_track drops in via `RIBO_ORF_TRACK`);
`orf_v2_attn` config (`--n_attn_layers 2 --emb_backend onehot --use_orf_track`, budget 16000, 40 epochs,
patience 8). The anchor V0 is trained fresh under this single-tissue protocol (the deployment orf_v2_attn
is train-on-8, not comparable).

Four arms (channel 3 only differs):
- **V0 heuristic** -- `START_W[codon] x (0.5 + 0.5*kz_heuristic)`; the existing `orf_track_v2.npy`
  (reproduced byte-identical by the extended builder as a check).
- **V1 no-Kozak** -- `START_W[codon]` alone (Kozak factor dropped); `build_orf_track.py --kozak none`.
- **V2 empirical PWM** -- `START_W[codon] x m_pwm`. `build_kozak_pwm.py` fits a position x nucleotide
  log2-odds PWM from 32,332 annotated CDS ATG starts (`tx2cds.tsv` has_start_codon=1, start at
  `utr5_len`, flanks from `fibroblast_universe.fa`) at context positions {-6,-5,-4,-3,-2,-1,+4} (the ATG
  triplet is excluded -- constant, captured by START_W) versus the universe background nt frequency
  (pseudocount 1). The multiplier `m_pwm = clip(0.5 + 0.5*(score - p5)/(p95 - p5), 0.5, 1.0)` uses the
  summed log-odds and the p5/p95 of the annotated-start score distribution, matching the heuristic's
  [0.5,1] range so V2-vs-V0 is a clean swap. NO test-label leakage: the PWM is fit from annotation
  SEQUENCE, not Hepatocytes Ribo-seq (same posture as the literature heuristic). The fitted PWM recovers
  the canonical Kozak (GCCRCCAUGG): -3 purine 85.8% (its strongest position), +4 G only 51.5%.
  `data/kozak_pwm.json`; tracks `data/packed/orf_track_v2_{nokozak,pwm}.npy`.
- **V3 learned** -- input is the V1 (no-Kozak) track; inside the model (`model.py
  --learn_start_context`, one-hot only) a `Conv1d(4 -> 1, kernel_size=10)` reads the one-hot slice with
  asymmetric pad (6 left, 3 right) so output[i] depends on one-hot[i-6 .. i+3] = Kozak positions
  [-6..-1, +1..+4]; `gate = sigmoid(conv)`; the start channel is recomputed `w x (0.5 + 0.5*gate)` before
  `in_proj`. The 4x10 kernel is a LEARNED position x nucleotide Kozak matrix (41 params incl. bias),
  extracted post-training and plotted against the empirical PWM and the heuristic
  (`plot_kozak_weights.py`).

Eval per arm (foreground what the start channel touches; bulk CDS metrics are a flat-sanity control):
1. **CTG non-AUG drop-in** (the sharpest test): `eval_kozak_arms.sbatch` dumps the per-nt predicted
   profiles, then `ribocode_dropin_ctg.sbatch` re-runs RiboCode `detectORF` with
   ALTERNATIVE_START_CODON_LIST=["CTG"] on real / predicted-at-observed-depth / predicted-at-predicted-
   depth density; `compare_dropin_ctg.py` reports start-codon-stratified precision/recall/F1 (genomic
   locus key, min_len 90, pval<=0.05, enrichment>=0.5).
2. **uORF / non-canonical localization** AUROC (length-controlled pred_frame0) + positives-only fidelity
   Pearson (`eval_localization.py`).
3. **Sanity** (expected flat): annotated-ORF AUROC, ATG drop-in F1, whole-pc + 5'UTR profile Pearson
   (`extra_metrics.json`, written at train time).
Orchestration `run_kozak_eval_chain.sh` (per-arm eval -> dropin -> compare SLURM dependency chain);
consolidation `aggregate_kozak.py` -> `results/kozak/kozak_summary.md`.

Pre-registered prediction (`KOZAK_PLAN.md`, from `kozak_context_alt_orfs.py`): RiboCode's real CTG calls
are ~92% uORFs, and both ATG-uORF and CTG-uORF starts sit in a weak, G-rich 5'UTR context far from the
strong C-rich canonical Kozak (the real dichotomy is uORF-vs-CDS, not CTG-vs-ATG). So the annotated-CDS
PWM (V2) is a miscalibrated prior for the non-canonical starts this project targets (it pushes ~15-25% of
uORFs to the 0.5 floor), and V2 may NOT beat V1 on uORF/CTG detection.

Outcome + action (2026-07-16): confirmed. In-distribution the four arms tie (best Fibroblast-val Pearson
0.6418-0.6458); on the cross-tissue Hepatocytes test V1 (no-Kozak) is directionally best on ~7 metrics
(CTG drop-in F1 0.369 vs 0.340, uORF frame0 Pearson 0.503 vs 0.461); V2 (empirical) is a wash; V3 (learned)
rediscovers the Kozak PWM at r=0.807 yet its explicit gate is the most conservative detector (worst CTG F1
0.320). The explicit Kozak factor is redundant with what the sequence backbone learns. ACTED ON:
`build_orf_track.py --kozak` default flipped `heuristic` -> `none` (function + argparse + docstring); the
existing deployment `orf_v2_attn` keeps its heuristic track `orf_track_v2.npy` (untouched, no mismatch), so
only new `ext` builds are no-Kozak. Deployment retrain deferred pending a 3-seed confirmation of the V1
edge. Full result: results.md Task 20.

## 5O. Posture-B multimap sensitivity check (2026-07-16, results.md Task 21)

The pre-registered heavy multimap check (plan section 2.3). Posture A (the target in use) counts each
Ribo-seq footprint weight 1 on every within-gene isoform it is compatible with (RiboCode default); the
resulting ~18x magnitude inflation is absorbed by per-transcript normalization, but a barely-expressed
isoform can inherit an expressed sibling's footprints in shared exons and so get a "phantom" profile it
did not earn. Posture B commits each footprint to ONE representative isoform.

Key equivalence that makes this cheap (no target re-derivation): for a gene's HIGHEST-EXPRESSED isoform,
posture-A per-nt counts already EQUAL its posture-B counts. A footprint on that isoform lands in an exon
it contains, so committing it to that isoform (the assignment target) does not change the
representative's per-nt counts; only low-expressed siblings lose inherited signal, and those fall below
min_signal and drop out of every split anyway. So posture B reduces to "restrict train+eval to the
highest-expressed isoform per gene," using the existing posture-A counts, embeddings, ORF track, and
coverage unchanged. `make_representative_tx.py` picks the argmax-mean-TPM isoform per gene (tie: longest,
then tx_id) from the matched Fibroblast salmon quant (`fibroblast_salmon_mean_tpm.tsv`) among the 36,668
packed tx -> `data/fibroblast_representative_tx.txt` (12,485 genes = 34% of packed tx; 11,241 pc + 1,244
lncRNA). `dataset.representative_tx()` reads it via env `RIBO_REPRESENTATIVE_TX` and every split
(`load_split`, `loto_split`, `heldout_test_tx`) keeps only those tx, mirroring the chrM exclusion.

Design (`train_postureB.sbatch`): a clean 2x2, {baseline, orf_v2_attn} x {A = full posture-A universe,
B = representative-only}, one-hot, fold-0 held-out-chromosome, identical code + seed 0, so the ONLY
variable is the transcript universe. Both orf configs use the no-Kozak track (the new default). Retraining
BOTH postures fresh (not reusing an old posture-A run) makes the comparison airtight. `eval_extra.py`
scores each on its own test set; the pc/lncRNA/uORF/dORF pooled medians are compared A vs B. Expected to
CONFIRM (Task 6's single-isoform-gene proxy already showed multi-isoform transcripts were not better,
i.e. no phantom inflation); a match closes the multimap question.

## 5P. mm1 (unique-mapper) RNA coverage + no-Kozak retrain + CLI data-prep hardening (2026-07-23)

Three coupled changes to the deployment model + the data-prep interface.

**mm1 RNA coverage (unique mappers).** All RNA-seq is now aligned with STAR
`--outFilterMultimapNmax 1` (unique genomic mappers), matching the Ribo-seq mapping posture
for consistency and smaller intermediates. The mm20-vs-mm1 proxy check (DoHH2, Task 25) showed
85.7% of transcripts identical, so this is posture-invariant. The 55 Chothani RNA SRRs were
re-aligned (BAMs archived to `data/rnaseq_bam_mm1/`, not deleted -- long-term storage), and each
tissue pack's `coverage.npy` was rebuilt from the mm1 per-sample coverage. The mm20 coverage is
preserved per pack as `coverage_mm20.npy`.

**No-Kozak ORF track.** All retrains now use the no-Kozak track
(`data/packed/orf_track_v2_nokozak.npy`, `build_orf_track.py --kozak none`), never the heuristic.
The Task-20 ablation showed the hand-picked Kozak factor is redundant with what the sequence
backbone learns and mildly harmful cross-tissue. The old deployed model stays heuristic-pinned
(train/inference consistency); its GSE120762 prediction therefore also stays heuristic.

**Brain dropped.** The retrain holds out Hepatocytes (same eval fold as the deployed model) and
trains on 7 tissues, dropping Brain (period_obs 0.044, noise). Final model:
`results/loto/orf_v2_attn_onehot_noBrain_nokozak_mm1_holdout_Hepatocytes` (all other
hyperparameters identical to the deployed model; see `train_loto_noBrain.sbatch`).

**CLI data-prep hardening (`scripts/prepare/`).** The pack-build chain was Chothani-specific and
fragile (absolute project root baked into every script, hardcoded tissue lists, `assert len==32`,
`SRR15513*` globs, name-keyed `DATASETS` dict / `heldout_config.sh` case blocks). A new path-driven
CLI toolkit replaces it, reusing the validated building blocks:
- `paths.py` -- project root from `$RIBOSEQ_SIGNAL_MODEL_ROOT` or auto-detect (no baked path).
- `packlib.py` -- shared, path-free pool + 7-file-pack writer (lifted from `build_heldout_pack.py`).
- `prepare_pack.py` -- build one pack from `*_psites.hd5` + `*_coverage.hd5`; `--universe-tx` OR
  `--ref-pack`; `--reuse-target` (swap coverage against a fixed target/universe, the mm1 path);
  optional `--build-orf-track --kozak none`.
- `prepare_from_bams.py` -- BYO-BAM entry: STAR toTranscriptome BAM paths -> filter_tx_heldout +
  RiboCode metaplots/RiboCode (auto P-site offsets) + rnaseq_coverage -> `prepare_pack`. All
  external tools resolved from `$RIBOCODE_BIN`/`$RIBO_PYSAM_PYTHON`/`$RIBO_SAMTOOLS`, never hardcoded.
- `prepare_packs.py` -- samplesheet batch (`group,assay,path`; BAM or hd5); one pack per group.
- `README.md` documents the interface. The 6 legacy `pool_*`/`pack_*`/`build_*psite_target*` scripts
  carry deprecation banners pointing here.

Validation: `prepare_pack.py` pooling HUVEC's Ribo target from the raw RiboCode hd5 reproduced the
existing pack's `target_counts.npy` **byte-for-byte**, with mm1 coverage differing as expected. The
production rebuild (`prepare/rebuild_coverage_mm1.sbatch`, `--reuse-target`, per-tissue) asserts the
target stays byte-identical before swapping coverage, and backs up mm20. The trainer/eval consume
these packs unchanged via `$RIBO_PACK_DIR`/`$RIBO_ORF_TRACK`/`$RIBO_ONEHOT_FASTA`.

## 5Q. ORF-call evaluation + count-head calibration methodology (2026-07-24, results.md "De novo ORF calling")

The standing way to evaluate any Ribo-seq-signal model: not just profile Pearson / localization AUROC, but
how the downstream **RiboCode ORF calls** change. Full plan + status: `docs/count_head_calibration_checks.md`.

- **Drop-in** (`scripts/ribocode_dropin.py`): runs RiboCode's `detectORF` on a SUBSTITUTED per-nt density --
  the model's predicted profile instead of experimental Ribo-seq -- with the prebuilt annotation + default
  params (`--min_aa 5 --pval 0.05`). Variants: `real` (observed P-sites; the reference), `pred_obsdepth`
  (predicted shape x REAL per-tx total; isolates shape), `pred_preddepth` (predicted shape x count-head total;
  standalone, no Ribo). Emits RiboCode-native `<variant>_collapsed.txt` with `ORF_type`.
- **Calibration knobs added:** `--pred_scale theta` (scale the predicted effective depth) and `--pred_poisson`
  (draw `Poisson(rate x theta)` instead of rounding = inject the shot noise real Ribo-seq has). The over-calling
  is a smoothness problem, so Poisson is the fix (Check 4).
- **Two-arm protocol (STANDING):** every ORF-calling test reports BOTH `pred_preddepth --pred_scale 1` (standard)
  AND `--pred_poisson --pred_scale theta` at the CDS-anchored operating point.
- **CDS-anchored dial** (`dropin_sweep_*.sbatch` + `cds_anchored_curve.py`): sweep theta; per theta compute CDS
  recall/precision (CDS = the trusted anchor) + non-canonical (uORF/novel/dORF) count + precision vs observed.
  Operating point = theta giving a target CDS recall (~90%). Lever comparison: `compare_levers.py` (depth-scale
  vs Poisson vs p-value at matched CDS recall). Seed stability: `poisson_seed_stability.py` (3 seeds).
- **Over-call validation:** a model's "extra" calls (predicted, absent from that sample's observed) are
  validated against a DEEPER independent observed reference (`gse120762_overcall_check.py`); real calls reproduce
  ~87% across conditions, spurious extras ~13%. Precision must NEVER be judged against a shallow same-dataset
  reference (it penalizes real recovery) -- always a deep independent one.
- **O2 reproducibility ceiling** (`o2_between_dataset.py`): coord-key-matched (genomic coords, annotation-version
  robust) `F1(datasetA_obs, datasetB_obs)` between two independent same-tissue Ribo-seq datasets = the de novo
  bar; then `F1(model, datasetB_obs)`. model/ceiling ~1.0 => as good as a second real experiment.
- **Periodicity QC:** RiboCode metaplots `f0_percent` (frame-0 enrichment) per read length from the sample's
  `*_pre_config.txt`; >70% on the dominant footprint lengths = good.

## 5R. O3 -- FP-filter + anti-smoothing model for novel-ORF over-calling (2026-07-24, results.md "O3")
Two tracks addressing RiboCode's over-calling of short non-canonical ORFs on the smooth predicted density.
Survey verdict (`docs/o3_orf_caller_miniplan.md`): every periodicity-based caller inherits the bug (a smooth
density is trivially periodic wherever the model puts in-frame P-sites), so keep RiboCode for candidate
generation and add a SHAPE-based FP-filter; the discriminating signal is read distribution, not periodicity.

- **Feature extraction (`scripts/heldout/orf_features.py`):** per candidate ORF, slice the model's predicted
  per-nt density over the ORF (start..stop, in-frame at position 0), reshape to codons, and compute shape
  features -- `frac_f0` (in-frame fraction), `f1max`, `pme` (per-codon in-frame percent-max-entropy = codon
  uniformity), `gini`, `cv`, `max_to_mean`, `frac_zero`, `ramp5` (5' start ramp), `drop3` (3'-of-stop drop),
  `len_codons`. Computed on the SMOOTH predicted density (not the Poisson sample). Label = 1 if the ORF's
  genomic coord key is in the deep observed calls. **Phase-1 finding:** PME alone AUC ~0.50 for novel (the
  model smooths REAL ORFs too), so the shape signal is largely absent from the output; best multivariate
  5-fold CV = logreg = RF = 0.669.
- **FP-filter (`scripts/heldout/fp_filter.py`):** `train` fits a per-ORF-class (+ pooled non-canonical)
  balanced `LogisticRegression` on the feature table with StandardScaler, 5-fold stratified CV for honest
  precision-at-recall, and stores per-target-recall thresholds; `apply` re-featurizes a `*_collapsed.txt`
  against `pred_profiles.npz`, scores each candidate, and drops those below the class threshold (CDS/annotated
  pass through untouched by default) -- this IS the post-RiboCode `--fp-filter` step. `pr_at_recalls` validated
  against `sklearn.precision_recall_curve`. In-sample lift +0.10 novel / +0.12 uORF precision at fixed CDS
  recall; cross-dataset (Wang/Janich) is the honest test, pending GPU.
- **Anti-smoothing loss (`model.py::profile_entropy_gap`, `train_loto.py --peakiness_weight`):** adds
  `w * mean_tx( relu( H(softmax(logits)) - H(counts/N) ) )` to the profile multinomial NLL. The NLL alone is
  minimized by matching the empirical profile but tolerates a smooth hedge; this term penalizes only
  OVER-smoothing (relu), referenced to the deep pooled target's entropy (real periodicity, not shot noise).
  `train_loto_noBrain_peaky.sbatch` trains an exact A/B variant (`..._peaky0p2_...`) alongside the current
  no-Kozak mm1 model. Default weight 0 keeps the deployed behaviour unchanged.

## 5S. `pgx` -- RiboCode-called proteogenomics search pipeline (2026-08-01)

Replaces the `pred_frame0`-threshold database build (`enumerate_score_orfs.py` +
`build_a549_dbs.py`) with the validated ORF-calling stack. Package:
`proteogenomics/scripts/pgx/`; plan and results in `docs/proteogenomics_pipeline_miniplan.md`;
CLI reference in `proteogenomics/scripts/pgx/README.md`.

### Why the f0 threshold was retired

`f0 = pm[0::3].sum() / pm.sum()` is a ratio over the ORF interval, so prepending a zero-signal
upstream region changes neither numerator nor denominator: `f0` is INVARIANT to N-terminal
extension, not merely diluted. Measured on BMDM over 12,719 stop-codon groups containing an
annotated CDS start plus upstream in-frame starts, **96.2%** of 51,646 upstream non-AUG ORFs fell
within 0.05 of the true CDS's `f0`; **median |delta f0| = 0.0026**. The old build also applied no
Poisson calibration, no periodicity test and no significance filter: the entire gate was
`f0 >= thresh` (`build_a549_dbs.py:103`), with `pred_density` computed and never read. Enabling
non-AUG enumeration multiplied candidates ~10x (625,698 -> 5,956,246), which `f0` passed at the
same ~7.4% rate, inflating the model DB 13,835 -> 142,248 without adding evidence.

### Pipeline

1. **Universe** (`pgx.universe`). Mean-TPM merge of the replicates' salmon quants, then the
   validated `build_line_universe` rule: TPM >= 1, `{protein_coding, lncRNA}`, chrom != chrM,
   mature length <= 10,000. One universe feeds the model arm and BOTH nulls. Verified to
   reproduce the existing BMDM universe byte-identically (21,393 tx).
2. **Calibration** (`pgx.calibrate`, implements Check 6). theta scales the predicted effective
   depth (`ribocode_dropin.py --pred_scale`) with Poisson sampling (`--pred_poisson`), the lever
   Check 4 identified as winning. Anchored on **relative** CDS recall, `hits / hits at the
   saturating theta`, NOT absolute recall of the annotation: only ~53% of expressed CDS-bearing
   transcripts are detectably translated (BMDM plateau 10,236 of 19,190), so an absolute 0.90
   target is unreachable by construction and the dial always saturates. theta* is the smallest
   theta clearing the relative floor, since tightening improves novel precision monotonically
   (Check 3). On BMDM this returns theta = 0.05, independently reproducing the operating point
   Check 5 validated as transferable.
3. **Two-arm calling** (`pgx.call_orfs`), per the standing two-arm rule: standard
   (theta = 1, deterministic) and Poisson (theta*). The Poisson arm at theta* is REUSED from the
   calibration sweep rather than recomputed.
4. **Extensions** (`pgx.extensions`). See below.
5. **Databases** (`pgx.build_dbs`): `db_gencode`, `db_model_<arm>`, `db_null_atg`, `db_null_nc`,
   sharing one GENCODE base plus inline `REV_` decoys.
6. **Search** (`pgx.search` + `search.sbatch`): MSFragger, `--enzyme tryptic|nonspecific`,
   keyed by `sha256(db) + sha256(params)` so nothing is searched twice.
7. **Report** (`pgx.report`): global 1% FDR for canonical, class-specific 1% FDR for novel.

### N-terminal extensions require their own code path

RiboCode cannot supply them. `orf_finder.orf_find` sets `alt_flag = 0` as soon as any in-frame ATG
precedes the stop, so **a near-cognate extension of an annotated CDS is unreachable at any
parameter setting**. For ATG extensions it does emit a call (`classfy_orf` labels same-stop ORFs
`annotated`, and `start_check` with `only_longest_orf=True` takes the most 5' start), but it tests
the WHOLE ORF and never asks whether that start is used.

`pgx.extensions` applies RiboCode's own statistic (`extract_frame` + `test_frame`: Wilcoxon
f0>f1 and f0>f2, Stouffer-combined, then BH) to the **extension region alone**, from the candidate
start to the annotated CDS start. Density is the arm's rate construction at theta* with
deterministic rounding rather than a Poisson draw, so it keeps real zeros (no smoothness artefact
through RiboCode's tie logic) while remaining exactly reproducible and independent of the caller's
RNG stream. Only CDSs the arm called translated are extended.

BMDM: 30,243 candidate upstream starts -> 5,396 too short, 23,963 with fewer than 5 nonzero
frame-0 codons in the extension region -> 884 testable -> **564 passing q <= 0.05** (CTG 152,
GTG 95, AAG 72, AGG 67, TTG 46, ACG 37, ATC 34, ATT 27, ATG 18, ATA 16). Of the 884 that reached
the test, the retired whole-ORF `f0 >= 0.5` criterion would have admitted **884 of 884**, and by
inheritance essentially all 30,243: a **54x** reduction on evidence rather than threshold. RiboCode
itself reported 236 `nterm_ext` among its `annotated` calls but only 18 ATG extensions survive the
region-specific test, so its most-5'-start label is unsupported roughly 92% of the time.

### Significance gate

RiboCode filters per-ORF on `--pval` BEFORE BH adjustment, so on the surviving set BH moves p by at
most 4e-4 (measured, 10,786 BMDM calls) and `--orf-qvalue 0.05` is a no-op. **`--pval` is the real
significance lever**; the q-filter only tightens further.

### Peptide-to-GENCODE assignment

MSFragger delimits its protein list with `;`, not `,`. Splitting on `,` collapses the list to one
token, and the classification then survives only because MSFragger happens to emit the canonical
protein first. Measured: 3,489 of 191,358 rank-1 PSMs in the BMDM model search map to BOTH a novel
ORF and a GENCODE protein, all ordered canonical-first, so the pre-existing scripts were correct by
luck. `pgx.report` splits on both characters, so the no-double-count rule holds regardless of order.

### MSFragger index sizing

MSFragger 4.2 sizes its own fragment-index slices from the JVM heap and REJECTS `num_slices`
("Unknown parameters"). The 2.36M-target near-cognate null built 672,835,480 fragments as 6.27 GB
in a single slice, so the only lever for a large database is `-Xmx` (`XMX` in `search.sbatch`).

## 6. Roadmap (subsequent steps)

Done: per-token RiNALMo embeddings for the 36,668-transcript expressed universe
(section 4B); held-out-chromosome gene-disjoint split (section 4B); RNAseq per-nt coverage
fleet + pooling (section 4); packed training data + dual-head model + trainer/eval
(section 5C); the fold-0 held-out-chromosome RiNALMo baseline (Task 6 in results.md:
protein_coding median profile Pearson 0.60, periodicity reproduced).

Next (updated 2026-07-11; the original items 1-3 are now DONE -- 5-fold CV, the RiNALMo vs
Orthrus comparison, and the input ablations are in results.md Tasks 6-9, and the improvement
matrix (ORF track / attention / CTG), the uORF/dORF-aware eval, and the concentration robustness
check are Tasks 9 to 9d):

1. **Non-canonical-ORF truth-set localization eval (PARKED 2026-07-11, highest priority).** The
   shape metrics (whole-tx / uORF / dORF profile Pearson) score how well predicted signal matches
   the OBSERVED Ribo-seq; they do not test the project's actual goal -- CALLING a translated ORF
   (start + reading frame), especially non-AUG uORF / dORF / lncRNA ORFs, against INDEPENDENT
   truth. This reframes eval from "predict the profile" to "localize/detect the ORF."
   - Truth set: Ribo-seq ORFs Consortium (Mudge et al., Nat. Biotech. 2022; ~7,200 human
     translated ORFs with GENCODE coordinates) as primary; SmProt and MS-validated peptides (from
     the sibling `expression_context_human` project) as orthogonal cross-checks. Map onto the
     Fibroblast universe transcripts. Negatives: sequence-derived candidate ORFs (ATG/CTG..stop)
     with no translation evidence.
   - Formulations: (a) start-site localization -- rank candidate starts (ATG + near-cognates) by
     the model's predicted local signal; is the true start top-k. (b) ORF detection -- score each
     candidate ORF by predicted in-ORF / in-frame signal; report AUROC / AUPRC / precision@k of
     translated-vs-not. (c) frame localization -- predicted periodicity phased to the true ORF
     frame. Stratify by ORF class (uORF/dORF/lncRNA) and start codon (AUG vs non-AUG): the key
     question is whether the non-AUG track actually helps LOCALIZE non-AUG ORFs, not just correlate
     profiles (the profile metric already says orf_v2_attn wins uORF shape; this tests calling).
   - Circularity guard (essential): a translated ORF has observed signal, so "predicted signal
     inside it" is partly circular. The eval must show the MODEL's prediction beats the raw
     observed coverage AND a sequence-only baseline (CPAT), especially on ORFs that are hard to see
     in raw Ribo-seq -- that is the real use case: predict translation from SEQUENCE where clean
     Ribo-seq is absent (a new lncRNA, a designed mRNA).
   - First build: map the Consortium set to the universe, score Consortium ORFs vs decoy ORFs with
     each config's predicted signal, report start-site top-k + detection AUROC by class and
     AUG/non-AUG, with raw-coverage and CPAT baselines alongside.
2. **Clean CTG isolation.** Task 9b's v3 track conflated occupancy-extension with start-channel
   shrink; run ATG-only occupancy + a graded start channel restricted to {ATG,CTG} to separate the
   two effects (does the CTG start CHANNEL help, distinct from CTG occupancy hurting).
3. **posture-B (primary-only) target rebuild (PARKED 2026-07-11).** The pre-registered heavy
   multimap sensitivity check (see 2.3). On the transcriptome axis a genomically-unique footprint
   still maps to every isoform of its gene that contains the exon, so per-transcript assignment is
   ambiguous within a gene. Posture A (current) counts each footprint weight 1 on every compatible
   isoform (RiboCode default); the resulting ~18x magnitude inflation is absorbed by scale-invariant
   per-transcript normalization, so it cannot change a single transcript's profile Pearson.
   - What posture B would catch: not magnitude (normalization handles that) but PHANTOM signal --
     a barely-expressed isoform inherits its expressed sibling's footprints in shared exons and so
     gets a plausible profile it did not earn. It is fundamentally a data-quality / isoform-
     resolution check, distinct from the shape and localization questions.
   - Posture B assignment options: commit each footprint to ONE representative isoform -- MANE /
     canonical, or the highest-expressed isoform per gene from the matched RNA-seq quant (salmon /
     RSEM), or STAR's primary-alignment flag. Posture C (fractional 1/N across the N compatible
     isoforms, EM-style) is the softer alternative.
   - Why heavy: posture B changes the LABELS, not just the eval, so it is not a re-eval -- it needs
     the target re-derived (re-assign footprints -> new per-nt hd5), re-packed
     (`target_counts.npy`/`offsets.npy`), all models RE-TRAINED, then re-scored. That is why it is a
     one-off rather than routine.
   - Already largely passed by proxy: Task 6 split the test set into single-isoform genes (where
     posture A == posture B, no ambiguity) vs multi-isoform genes, and the multi-isoform transcripts
     (the ones affected by posture A) were NOT worse -- so the rebuild is expected to CONFIRM, not
     change, the conclusions. It is the definitive version of a check the proxy already mostly cleared.
   - First build: pick the assignment (highest-expressed-isoform from the existing Fibroblast RNA-seq
     quant is the most defensible), rebuild target + pack, re-train the baseline + orf_v2_attn on
     the posture-B target, and compare pc/lncRNA/uORF/dORF pooled medians against the posture-A
     numbers; a match closes the multimap question.
4. **8 matched tissues: leave-one-tissue-out**, then cross-dataset generalization given deep
   RNAseq (the real objective).

---

## 7. Conventions

Cluster `prism` (SLURM): heavy work via sbatch (partitions short 1 h / medium 12 h /
long 14 d / gpu); head node `mustard` for light inspection and serial downloads only.
Outputs on the group filesystem; scratch under `/data/tmp/emalekos`. Writing style:
ASCII punctuation only (no en/em dashes), specific descriptors instead of first-person
possessive, no AI attribution anywhere, no commit or push without explicit instruction.

### 7.1 Project rule: exclude mitochondrial genes (chrM)

Standing rule for this project. Mitochondrial genes are NEVER part of the transcript
universe, target, training set, or evaluation -- for any tissue, dataset, or rebuild.
Mitochondrial mRNAs are mitoribosome-translated (a different genetic code), leaderless
(no 5'UTR / uORFs), and lack the cytoplasmic 3-nt periodicity the ORF track and FM
embeddings encode, so they are not a valid target and are predicted poorly (see Task 11).
`Mt_rRNA` / `Mt_tRNA` are handled by the ncRNA drop-list; this rule additionally removes
the 13 chrM protein-coding mRNAs that otherwise pass as `protein_coding`. Enforced by two
gates that every new build must keep: `EXCLUDE_CHROMS = {"chrM"}` in
`scripts/dataset.py` (`excluded_tx()` filters `load_split` -- train + eval) and in
`scripts/define_universe_and_fasta.py` (universe / FASTA / embedding builds). Any new
universe, split, or eval path added later must route through these, or apply the same
`chrom in EXCLUDE_CHROMS` filter.

### 7.2 Project rule: the RNA arm must be poly(A)-selected (2026-07-31)

Standing rule for this project. Any RNA-seq used to build an expressed-transcript universe
or to supply the model's RNA-coverage input must be **poly(A)-selected**, verified before
the download is started.

**SRA `library_selection` is not acceptable evidence.** Ribo-Zero total-RNA libraries are
routinely filed as `library_selection = cDNA`. Accept only `library_selection =
PolyA`/`Oligo-dT`, GEO `!Sample_molecule_ch1 = polyA RNA`, or an explicit oligo-dT /
poly(A)-enrichment sentence in `!Sample_extract_protocol_ch1`.

Rationale, measured on GSE243134 mouse liver (the third liver arm). Ribo-Zero removes rRNA
but **not tRNA, 7SL/SRP or snRNA**, and total RNA additionally retains intronic pre-mRNA
that cannot map to a transcriptome index:

| | GSE243134 totalRNA | Wang `P42_Liver_mRNA` poly(A) |
|---|--:|--:|
| salmon mapping rate | 9.75-23.64% (median **17.73%**) | 91-93% |
| transcripts detected (TPM > 0) | 62,206 | 52,843 |
| transcripts at TPM >= 1 | **3,534** | 14,331 |
| top-10 share of total TPM | **93.6%** | 32.9% |
| resulting universe | **3,321 tx (aborted)** | 13,041 tx |

Detection was not the problem; 62,206 transcripts received reads. The failure is that ~94%
of all TPM is absorbed by ten structural transcripts (`n-TKctt14`, a tRNA-derived locus,
alone accounts for 82.8% of the library; `Gm59647` is the 7SL/SRP duplicate that also
dominated Janich), which pushes every real transcript below TPM >= 1. Both are typed
`lncRNA` in GENCODE, so the protein_coding+lncRNA biotype filter cannot exclude them.
Alignment and trimming were correct; only the library chemistry was wrong.

Two diagnostics worth reusing:
- A low **transcriptome** mapping rate beside a normal STAR **genome** unique-mapping rate
  (18% vs 28-39% on the same reads) is the signature of this failure, not a processing bug.
- A universe-size guard (`[ "$NTX" -gt 6000 ]`) in the pack build is what caught it before a
  bad universe propagated into the eval. Keep such a guard on every new arm.

Consequence: GSE243134 contributes its 21 WT **RPF** libraries only (kept: 40/40 aligned,
median 24.5M unique reads). Its 19 totalRNA libraries are replaced by an external
poly(A)-selected liver set, screened on one run per candidate against the same three
statistics before committing to a full download
(`scripts/heldout/screen_liver_polya.sh`). That external RNA is **not donor-matched** to the
Ribo-seq, unlike the Wang and Janich arms; since liver is strongly circadian and the RPF is
ZT5/ZT12, this is stated as a caveat on the model-prediction sub-arm. The arm's primary
contribution, an independent set of *observed* RiboCode calls, needs the universe but not
donor-matched per-nt RNA.

Near-miss caught by this rule: **GSE73554 (Atger 2015)**, circadian-matched mouse liver,
filed as `library_selection = cDNA`, but GEO states "TruSeq Stranded Total RNA ... Ribo-Zero
Gold depletion set". It would have reproduced this failure exactly.

### 7.3 Project rule: FDR-filter every column of an MS DB-comparison table (2026-08-01)

Standing rule for this project. In any mass-spec search-database comparison, **every reported
column must be FDR-filtered by the same procedure. Never place a filtered column beside an
unfiltered one.**

Raw rank-1 counts measure how many spectra a database ABSORBED, not how many it identified, and a
larger database absorbs more by chance. A table that filters the novel column but not the canonical
and net columns will therefore reward database size and can inverst its own conclusion.

Which FDR where:

- **canonical and net**: GLOBAL 1% FDR (all targets vs all decoys). This is what the search
  reports and it keeps the target-decoy competition intact.
- **novel**: CLASS-SPECIFIC 1% FDR (novel targets vs novel decoys only). A global FDR inflates
  non-canonical discovery ~10-13x.
- **always** also report discovery rate per 100k database sequences; raw counts are not comparable
  across databases of different size.

Do NOT apply class-specific FDR to the canonical side. A large novel space cannibalises canonical
DECOYS faster than final-recipe targets (macrophage BMDM: `other_d` -29.5% vs `canon_t` -13.2%), which
deflates the estimated canonical FDR, relaxes its threshold, and makes a pure-junk database appear to
GAIN canonical identifications (+3,317 peptides over a canonical-only search, which is impossible).
Both the unfiltered and the naive class-specific versions are biased toward the larger database, by
different mechanisms.

Measured consequences on the 12 mouse macrophage populations, which is why this rule exists:

| claim from the unfiltered table | what survives FDR |
|---|---|
| null displaces 1,480,032 canonical PSMs (2.8x the model DB) | canonical PSMs at 1% FDR are 335,548 / 335,804 / 335,392 / 336,284 for canonical-only / model@0.5 / model@0.74 / null -- a **0.27% spread** (-0.6% to +1.0% across all 12) |
| net PSM favours the null (+1.0%) | the null's lead is 1,592,864 unfiltered novel matches yielding 374 FDR-surviving peptides (4,259 raw per real peptide, vs 399 for the CDS-anchored arm) |
| model DB protects previously-confident canonical IDs | WITHDRAWN -- the displaced PSMs were sub-threshold and were never confident |
| model 440 vs null 374 novel peptides | UNCHANGED -- same decoy class, procedure and denominator across arms |

Sanity check to run on any new table: **canonical counts should be near-flat across databases.** A
large apparent canonical churn means the FDR filter was skipped somewhere.

Reference implementations: `proteogenomics/scripts/net_psm_fdr.py` (corrected computation),
`proteogenomics/scripts/macro_churn_aggregate.py` (now global-FDR canonical/net + class-specific
novel). Superseded tables preserved under
`proteogenomics/data/macrophage_tissue/superseded_unfiltered/` with a README naming the wrong claims.

### 7.4 Project rule: the mandatory proteogenomics table columns (2026-08-01)

Standing rule. Every proteogenomics DB-comparison table reports, **per population** (never collapsed
to a TOTAL row alone), exactly these columns:

| column | definition | FDR |
|---|---|---|
| `novel` | unique novel peptides, MODEL DB | 1% **class-specific** (vs `REV_nuORF|` decoys only) |
| `null` | unique novel peptides, NULL DB | 1% class-specific |
| `DB seqs` | model DB novel target sequences | -- |
| `null seqs` | null DB novel target sequences | -- |
| `canon base` | canonical peptides, canonical-ONLY DB | 1% **global** |
| `dPC model` | change in canonical vs baseline, model DB | 1% global |
| `dPC null` | change in canonical vs baseline, null DB | 1% global |
| `ncStart` | model peptides absent from the null amino-acid space | 1% class-specific + substring scan |

**Why both dPC arms.** Adding sequences to a search perturbs the target-decoy competition, so
canonical identifications can be lost purely because the database GREW, independent of whether the
added ORFs are real. Reporting the model arm alone cannot separate a database-size cost from a
database-quality cost. Showing both makes it separable, and it is important to know when canonical
identifications are lost even if the cause is only database size. FDR control is the central
methodological hazard in proteogenomics; this table is where it is made visible rather than assumed.
The two DB-size columns are mandatory for the same reason -- discovery counts are uninterpretable
without their denominator.

**FDR assignment** follows 7.3: global for canonical/dPC (keeps the target-decoy competition intact),
class-specific for novel (a global FDR inflates non-canonical discovery ~10-13x). Raw rank-1 counts
are never reported.

**`ncStart`** is the model-only-discovery column: peptides that pass FDR in the model search and occur
NOWHERE in the null DB's amino-acid space, found by substring scan against every null target sequence
(a peptide is a protein fragment, so set membership is insufficient). The expected mechanism is a
non-canonical start codon extending an ORF N-terminally past the first in-frame AUG, yielding tryptic
peptides that exist in no AUG-only database. This is the cleanest demonstration of model value,
because no naive enumeration can produce these at any threshold.

**Current status: `ncStart` is structurally 0.** `enumerate_score_orfs.candidate_orfs()` matches
`is_atg = (c0 == A) & (c1 == T) & (c2 == G)` only -- verified empirically, 400,000/400,000 BMDM
candidate ORFs begin ATG -- and `db_model` is a strict SUBSET of `db_null` (both are built from a
single `candidates.faa`), so no model peptide can be absent from the null space. Making the column
meaningful requires extending enumeration to non-AUG starts. The model's ORF track already SCORES
them (CUG 0.5, GUG/ACG 0.35, UUG 0.3, AUA 0.25, AUU/AUC 0.2) and Task 19 studied CTG initiation, but
the enumerator emits AUG only. The column is COMPUTED rather than hard-coded so it activates
automatically once that lands.

Reference impl: `proteogenomics/scripts/proteomics_table.py --tag --reuse-tag --out`.

### 7.5 Project rule: the ORF-length floor is set BY ASSAY (2026-08-15)

**The rule.** Tryptic whole-cell-lysate proteomics uses a **30 amino acid** minimum ORF length. MHC /
HLA immunopeptidomics uses **7**. Nothing else sets it -- not the dataset, not convenience, not what
a previous run happened to use.

**Why the two differ.** 30 aa is not a proteomics threshold at all: it is `build_loader`'s ORF-call
filter (>= 90 nt, and `ORF_length` excludes the stop codon, so 90/3 = exactly 30 residues), applied
identically to model calls and reference calls. A proteogenomics database built below it therefore
searches ORFs the ORF-call track excludes **by rule**, and the two halves of the project stop being
comparable. HLA-I is the one genuine exception: its peptides are 8-11 residues and are not tryptic
products, so a 30-aa floor would exclude the biology being measured.

**Classify from the search parameters, never from a name.** `fragger.params` with `trypsin` and
`termini = 2` is tryptic. In this project: the 12 mouse macrophage populations and A549 are tryptic;
HBL-1, SU-DHL-4, DoHH2, THP-1 and B721.221 are MHC.

**Enforcement.** `pgx/build_dbs.py --assay {tryptic,mhc}` SETS `min_aa` (30 / 7) and aborts if an
explicit `--min-aa` contradicts it. `pgx/filter_db_min_aa.py` length-filters an existing database and
asserts target/decoy balance, which is exact: `min_aa` is a pure length filter, so filtering a built
FASTA is equivalent to rebuilding at the higher floor.

**Why it matters quantitatively.** At 7 aa the macrophage AUG null carried 235,686 sequences, ~65% of
them sub-30-aa. Two P11 headline ratios were inflated as a result -- database size 156x -> **78x**,
discovery density 111x -> **43x**. In the human immunopeptidome data the floor is load-bearing in the
other direction: 27-33% of model discoveries come only from sub-30-aa ORFs, versus 0% for CPAT/CPC2,
so applying 30 aa there would delete the model's advantage along with the biology.

## 5S. RNA-quality factorial: separating the two paths RNA-seq takes into the model (Task 68, 2026-08-08)

**Motivation.** Every mouse-liver arm in the project varies RNA-seq and Ribo-seq together, so "better
RNA" has never been separated from "different experiment". RNA-seq enters the model along two distinct
paths, and they can be varied independently:

1. **UNIVERSE** -- salmon TPM >= 1 decides which transcripts exist at all.
2. **COVERAGE** -- the per-nt RNA depth track the model consumes as an input channel.

**Design.** Ribo-seq is held FIXED across all three arms: the same 5 Janich P-site hd5, read out of arm
a's `provenance.json` rather than re-listed, so the arms cannot drift apart through a typo.

| arm | universe | coverage | pack |
|---|---|---|---|
| a | Janich-decontaminated (14,887 tx) | Janich | `data/packed_heldout_mouse_janich_liver_decon` (pre-existing) |
| b | Janich-decontaminated (14,887 tx) | PRJEB34766 | `data/packed_heldout_rnaq_b_uniJanich_cov34766` |
| c | PRJEB34766 (16,474 tx) | PRJEB34766 | `data/packed_heldout_rnaq_c_uni34766_cov34766` |

- **a vs b = the pure COVERAGE effect.** Arm b is built with `prepare_pack.py --ref-pack <arm a>`, so
  `tx_order` and `lengths` are copied verbatim and arm a's `orf_track_v2_nokozak.npy` stays valid
  (symlinked in). The job ASSERTS this rather than assuming it: tx_order equal, lengths equal,
  `target_counts.npy` byte-identical (proving the Ribo arm really did stay fixed), and `coverage.npy`
  NOT equal (proving the RNA swap really happened). Measured coverage mean 152.54 -> 842.90.
- **b vs c = the pure UNIVERSE effect**, and it is NOT a like-for-like F1 comparison: arm c's reference
  call set comes from RiboCode over a larger transcript space, so its F1 has a different denominator.
  The scorer prints the two contrasts in separate blocks and says so, because reading a b-vs-c F1 delta
  as an effect size is the obvious way to misuse this table.

**Comparator choice.** Both existing mouse-liver RNA arms are SINGLE-end (Janich 50-51 nt, GSE243134
75-80 nt). PRJEB34766 is PAIRED, `library_selection=PolyA` per ENA, 7 CreNeg control runs / 223 M
pairs, C57BL/6 male adult liver taken from the ENA sample `strain` field rather than inferred from a
filename -- matching the Ribo-seq strain. Verified on arrival: 90.6-92.4% salmon mapping, 81.1% STAR
unique, top single transcript 4.3-5.0% of TPM (Janich before decontamination was 52-70%). It is a
circadian study, so pooling the 7 controls averages over a rhythmic transcriptome; that is acceptable
for a coverage track and is stated rather than hidden. CrePos/KO/WT runs are excluded because Reverba
deletion perturbs the circadian liver transcriptome.

**PRJEB86747 was rejected and deleted**, not merely unused. It was the depth extreme (481 M pairs,
NovaSeq X) but carried two confounds fatal to an experiment whose entire purpose is isolating RNA
quality: strain C3H/HeNRj against C57BL/6 Ribo-seq, and a mutation-accumulation design (Heredity 2025,
doi 10.1038/s41437-025-00819-0, E-MTAB-14914) whose subject is between-line expression VARIANCE, so
pooling replicates works against what the data was collected to show. "Better RNA" would have been
inseparable from "different mouse".

**Compute note.** mamba4 cannot run on CPU -- `mamba_ssm` dispatches to `causal_conv1d_cuda` and raises
`Expected x.is_cuda() to be true` at the first block, with no CPU fallback in the installed build. The
attn arms run either way, so when the gpu partition is saturated the attn dumps go to `medium` with
`--export=ALL,DEVICE=cpu` and only the mamba4 dumps wait for a GPU. Same script both ways, so the two
paths cannot drift; the script now aborts immediately on `mamba4 + DEVICE=cpu` rather than 15 s in.

Reference impl: `scripts/rna_quality_factorial_packs.sbatch` (stage 1, packs),
`scripts/rna_quality_factorial_dump.sbatch` (stage 2, dumps + 3 drop-in variants),
`scripts/score_rna_quality_factorial.py` (stage 3, the two contrasts).

## 5T. Canonical alignment recipe, its enforcement, and the redo it forced (2026-08-13)

### The recipe

`scripts/riboseq_align.sbatch` is the single Ribo-seq aligner for both species, samplesheet-driven
(`dataset, run, fastq, adapter, species`). Species selects only the STAR index, the ncRNA transcript
list, and the tx->gene map. Fixed for every run:

- `--alignEndsType EndToEnd`. Soft-clipping shifts the inferred P-site, which is the exact quantity
  the model predicts. Measured on HUVEC against its reference pack with alignment as the only
  variable: EndToEnd r = 0.935, soft-clipped r = 0.583.
- `--outFilterMultimapNmax 1`. RiboCode does not drop multi-mappers itself.
- `--quantMode TranscriptomeSAM`, then coordinate-sort, index, and `filter_tx_heldout`
  (ncRNA-transcript drop + cross-gene read drop).

Four other STAR filters were swept on HUVEC across 13 arms
(`--outFilterMismatchNoverLmax`, `--outFilterMatchNminOverLread`, `--outFilterScoreMinOverLread`,
`--seedSearchStartLmax`). In-frame fraction was 55.95-56.05% in every arm, a 0.10-point spread
across a 3x change in allowed mismatch rate. They stay at STAR defaults: tightening costs up to 5.3%
of reads, buys no frame purity, and would widen divergence from the training packs.

### Two guards, both of which fired

- **`filter_tx_heldout.py` refuses coordinate-sorted input.** Cross-gene detection needs
  query-grouped records; fed a coordinate-sorted BAM it silently reports ~0.30% dropped instead of
  the true 9-22%. That false negative had already produced a wrong conclusion ("the filter is a
  no-op") and, downstream of it, a wrong claim that the training data could not be reproduced to
  better than 13%.
- **`prepare_from_bams.py` requires the ncRNA filter arguments** whenever `--ribo-bam` is passed,
  unless `--no-ncrna-filter` is given explicitly. Deviation by declaration, never by omission.

### Reproducibility, resolved

Under the final recipe, `packed_union_HUVEC` reproduces at **3.04% P-site divergence, r = 0.957**
(reference Janich: 1.69%). Earlier reconstructions reached 11.69% and 12.82%. The gap was recipe
drift across three compounding errors, not a property of the data: the pipeline and the training
data were aligned all along.

### The redo

Everything built off-recipe was rebuilt rather than annotated: 40 Ribo alignments -> 3 pools ->
9 packs -> 18 dumps -> rescore, plus the three interpretability analyses. Outputs went to `*_canon`
paths alongside the archived originals (`results/_archive_offrecipe_2026_08_13/`) so each finding
could be compared old-vs-new instead of silently replaced, and expectations were written down before
the comparison. Final-recipe pooling removes 3.5-4.9% of P-sites and 3.8-4.3% of observed ORF calls.
Every qualitative conclusion held; see results.md.

### New scripts this sprint

| script | purpose |
|---|---|
| `scripts/riboseq_align.sbatch` | THE canonical aligner; replaces 5 one-off `align_*` scripts |
| `scripts/score_channel_ablation.py` | scores the ORF-track channel ablation; previously an inline shell heredoc |
| `scripts/prepare/extract_universe_fasta.py` | subset a universe FASTA from GENCODE; generalizes the Fibroblast-hardcoded `define_universe_and_fasta.py` |
| `scripts/rebuild_heldout_canon.sbatch` | rebuild a held-out pack on canonical Ribo alignments, keeping its universe and ORF track |
| `docs/PIPELINE_POLICY.md` | the three standing rules, with the actual failure as the worked example |

### A near-miss worth recording

The leukocyte rebuild initially globbed `data/heldout_bam/<arm>/*.bam` for RNA inputs. Those
directories are COMBINED: they hold the arm's RNA alignments, the old off-recipe Ribo alignments,
and (for the GSE120762 arms) genome-coordinate BAMs beside transcriptome ones. The glob would have
fed Ribo runs in as RNA coverage and double-counted samples present as both BAM types, producing a
plausible pack rather than a crash. Caught by reading the job's own echoed input list within a
minute of submission, and cancelled.

The fix is explicit per-arm RNA run lists (`data/heldout_refs/<arm>_rna_runs.txt`) plus a guard that
aborts if any Ribo run appears among the RNA inputs. Those lists were then checked against
`scripts/heldout/heldout_config.sh`, the authoritative record of the original splits, and match its
`RNA_SRRS` exactly (lps 41/42, nt 38/39/40, tcell 20/21/22) -- confirming the ORIGINAL packs were
built on the correct RNA runs, so their only deviation is the alignment recipe. The general lesson
matches `feedback_blast_radius_scan_by_artifact`: the scan was run over the artifact (every RNA BAM
directory, every script using that glob idiom) rather than the one arm where it was noticed, which
is what showed `mouse_janich_liver` to be clean and bounded the problem to two directories.

### The aligner grew two samplesheet columns, and why that is the point

`scripts/riboseq_align.sbatch` gained optional `umi` and `cutadapt_extra` columns (5-column
samplesheets still parse unchanged). With them, all four datasets this project handles are expressible
as samplesheet ROWS rather than as branches inside a script, which is the operational content of
policy rule 1. `docs/PIPELINE_POLICY.md` had already declared UMI a legitimate per-dataset knob; the
script simply had not implemented it, and that gap was what blocked the human rebuild.

UMI ordering is forced, not chosen: `filter_tx_heldout` needs query-grouped input, `umi_tools dedup`
needs coordinate-sorted + indexed input, and the filter leaves the latter behind. STAR -> filter ->
dedup is the only order that satisfies both, and running dedup first now trips the sorted-input guard
rather than silently under-filtering.

### CAR-T: deduplication that ran on a file the pack never read

Found while scoping the human rebuild. `process_cart_gse304796.sbatch` runs `umi_tools dedup`, but on
the GENOME BAM. It then copies the raw STAR transcriptome BAM to `${RUN}.toTranscriptome.bam`, and the
pack is built from exactly that file (confirmed via `_work/psites/SRR34901743.toTranscriptome_psites.hd5`).
The deduplication was real and the log looked correct; it just applied to a different BAM.

Duplicate rate, from the pipeline's own `read_counts.txt`:

| run | aligned | deduped | PCR duplicates |
|---|--:|--:|--:|
| SRR34901743 | 33,114,136 | 23,888,001 | 27.9% |
| SRR34901744 | 43,372,505 | 30,637,962 | 29.4% |
| SRR34901745 | 45,122,621 | 31,747,439 | 29.6% |

So the CAR-T pack carries ~29% PCR duplicates plus the ~9.7% the filter would have removed.
Duplicates concentrate at specific positions rather than spreading evenly, so profile SHAPE is
affected and not merely depth. Every CAR-T-derived number is provisional pending the rebuild,
including the human codon-occupancy value of 0.548, which pools THP-1 with CAR-T; the mouse codon
number, which is what the transfer argument rests on, is untouched. Flagged in
`docs/PIPELINE_POLICY.md` and in the tutorial's codon page rather than left implicit.

This is the same failure shape as the coordinate-sorted filter earlier in the sprint: a step that ran,
logged success, and had no effect on the artifact downstream actually consumed. Both are invisible to
exit status and visible only by asking which file the next stage opens.

### Status of the human rebuild

Ribo arms only (RNA-seq is unchanged and not re-aligned): 5 GSE208041 + 3 GSE39561 + 3 CAR-T = 11
runs, samplesheet `data/human_ribo_canon_samplesheet.tsv` with its deviations recorded in the sibling
README. Two deviations from the previous human processing are declared there: no bowtie2 snoRNA
depletion (the final recipe has none and mouse never did; measured 0.01-0.87% of reads, a null),
and CAR-T dedup now acting on the BAM the pack reads. The first means a rebuilt GSE208041 pack differs
from the current one by two changes rather than one -- the filter dominates by two orders of
magnitude, but the comparison is not single-variable and should not be described as one.

## 5U. What the final-recipe rebuild actually changed (2026-08-13/14)

The rebuild was undertaken for reproducibility, but it produced a substantive result about the
reference data that any non-canonical number depends on.

### The model did not change; the reference did

The deployed checkpoint is frozen and its inputs are sequence, the ORF track, and RNA-seq coverage.
RNA-seq was never re-aligned, so none of those moved. Verified rather than assumed, on the matched
3x3 cell:

| file | off-recipe md5 | canonical md5 | |
|---|---|---|---|
| `pred_preddepth_collapsed.txt` | `27d5084b` | `27d5084b` | **byte-identical**, 19,601 calls both |
| `pred_obsdepth_collapsed.txt` | `556409b1` | `15c6ab49` | differs -- scaled by OBSERVED depth |
| `real_collapsed.txt` | `488e747a` | `88c4e3cf` | differs, 15,771 -> 14,877 |

So the standalone arm's F1 change (0.8762 -> 0.8661) is entirely the reference moving.

### The reference lost non-canonical calls, disproportionately

| class | off-recipe | canonical | change |
|---|--:|--:|--:|
| annotated | 10,535 | 10,385 | -1.4% |
| uORF | 2,587 | 2,263 | **-12.5%** |
| novel | 1,228 | 954 | **-22.3%** |
| dORF | 414 | 362 | -12.6% |
| Overlap_uORF | 591 | 534 | -9.6% |
| internal | 355 | 324 | -8.7% |
| TOTAL | 15,771 | 14,877 | -5.7% |

Canonical alignment strips 22.3% of novel and 12.5% of uORF reference calls while touching annotated
CDS by 1.4%. That asymmetry is the signature of soft-clipped and unfiltered reads manufacturing
spurious non-canonical ORFs, which is what `EndToEnd` and `filter_tx_heldout` exist to remove: deep,
unambiguous CDS is unaffected, and marginal low-count calls are where artifacts land.

**Therefore the off-recipe non-canonical F1 was INFLATED** -- the model was credited for matching
alignment artifacts. Canonical non-canonical numbers are the honest ones and must never be compared
against previously reported values, which were measured against a reference containing ~22% more
novel calls. It also means the recipe's value exceeds what the P-site totals suggest: a 4.9% drop in
P-sites removed 22% of novel calls.

### Reproducibility is heterogeneous across tissues

Extending the single-tissue test to all eight (74 runs re-aligned, per-tissue re-pool, same metric;
the generalised script reproduces the HUVEC number exactly at 3.04% / r 0.9567): median **3.86%**,
range 3.04-11.66%, per-nt r 0.432-0.980. **HUVEC is the best of the eight**, so the single-tissue
number understated the typical case and understated the worst by nearly 4x. Divergence and
correlation are partly decoupled (Fat moved most signal but coherently; HA_EC less but incoherently),
and there is NO depth relationship. HA_EC (r = 0.432) remains an unexplained outlier with four
hypotheses tested and rejected; no causal story about it belongs in the manuscript.

Table `results/chothani_regeneration/divergence_by_tissue.tsv`, figure `figures/S_reproducibility/`.

## 5V. Merged 28-library reference and the over-call validation (2026-08-14, results.md 2026-08-14)

**Purpose.** Every observed reference in this project is one dataset deep, so none of them can
distinguish a model false positive from a real ORF that dataset was too shallow to call. This builds
a deeper arbiter and, critically, the control that makes its verdict interpretable.

### Building the merged reference

`scripts/merged_liver_ribocode.sbatch` (job 36780045, medium partition, 12 cpus, 180 GB, 2h01m).

- **Inputs:** 28 aligned with the final recipe mouse-liver Ribo-seq transcriptome BAMs -- Janich 5,
  GSE243134 21, Wang 2. GSE243134's members come from `data/liver3x3/gse243134_ribo_runs.txt`, the
  same run list `liver3x3_01_pool.sbatch` uses, NOT a directory glob: that submission ships more
  libraries than the 21 WT RPF ones, and a merge over a different set would not be comparable to the
  per-dataset references it exists to be compared against. The script asserts exactly 28 and writes
  `bam_list.txt`.
- **Calling:** RiboCode `metaplots` derives a P-site offset per LIBRARY (28 config rows, one per BAM,
  each with its own selected read lengths), then `RiboCode -a <annot> -c <config> -l no -g` tests
  periodicity on the POOLED signal.
- **This is a joint call, not a union of three call sets.** An ORF with weak but consistent signal in
  all three datasets can clear significance here while failing in each dataset alone. Unioning three
  separately-called sets cannot do that, and would inherit each set's threshold.
- **Result:** 26,388 raw ORF calls, vs 14,877 (GSE243134, +77.4%), 14,322 (Janich, +84.2%),
  13,060 (Wang, +102.1%).
- Resumable: the psites hd5 are written before anything reads them.

### Scoring the over-call (`scripts/score_merged_liver_overcall.py`)

Filtering and keying are `compare_dropin_calls.build_loader` -- the merged reference is restricted to
the SAME 22,974-transcript space as every other drop-in number here and matched on
`(gene_id, ORF_gstop)`. The model arm is `pred_preddepth` (standalone) with GSE243134 RNA, matching
`figures/P9_upset_mouse3x3`.

Four strata, and the design point is that the middle two are useless without each other:

| stratum | definition | what it establishes |
|---|---|---|
| **anchor** | model calls that ARE supported by some experiment | that the merged reference is not the limitation (99.3%) |
| **control** | observed calls made by exactly ONE experiment, missed by the other two | how often "unsupported by the other arms" just means "too deep for them" (74.5%) |
| **question** | model calls supported by NO experiment | the number of interest (7.7% / 8.6%) |
| **shared** | called by BOTH architectures, no experiment | whether cross-architecture agreement is evidence (8.7% -- it is not) |

**Why the control is mandatory.** The merged reference is deeper than any single dataset, so it
recovers more of everything; a bare recovery rate for the model's extras would be confounded with
depth. The control is the closest available matched comparison: equally unsupported by the other 3x3
arms, but known-real because an experiment saw them.

**Why the gene-coverage split is mandatory.** A model-unique call the merged reference misses is only
FP evidence if the merged reference had signal at that locus. Splitting on whether merged called ANY
ORF on the same gene separates "the specific ORF was rejected on its merits" (881 calls, 16.9%
corroborated) from "the gene shows no translation at all" (1,061 calls). The 0% recovery in the
second stratum is definitional -- recovery requires a merged call on the gene -- so its `n` is the
measurement and its rate must never be quoted.

Per-class rates are emitted separately because lncRNA and uORF fail in opposite directions, and
because `internal` is the one class where the control itself is weak (33-54%), making a low model
number there genuinely ambiguous.

**Outputs:** `results/merged_liver_ribocode/overcall/{overcall_recovery.tsv,
overcall_recovery_by_class.tsv, overcall_meta.json}`. `figures/P9_upset_mouse3x3` reads these at draw
time so the plotted FP rate cannot drift from the scorer, and falls back to the old "candidate FPs"
label if the scorer has not run rather than plotting a stale number.

**Scope.** This is a fourth, deeper validation arm. It does NOT replace the three independent
observed arms the 3x3 factorial needs, and the same construction has not yet been done for any human
dataset.

---

## 2026-08-16 -- RNA census columns, attn exemplars, portable architecture spec

### RNA facts in the dataset census (`scripts/build_dataset_census.py`)

`rna_facts(pack_rel)` resolves `(n_rna_libraries, rna_coverage_total, source)` in three tiers,
because no single field is trustworthy across pack generations:

1. `provenance.json` -> `rna_coverage_inputs` (the actual file list). Authoritative.
2. `coverage_norm.json` -> `n_rna_samples`, **only when > 0**. It is 0 in every `_canon` pack,
   which reused its union twin's coverage without carrying the count.
3. The union twin, accepted only when `coverage_total` matches exactly, which proves the same
   pooled data.

Unresolved rows are written **empty, never 0** -- an empty count is a missing measurement, a zero
is a claim. `RNA_PACK_OVERRIDE` points B721.221 at `data/packed_heldout_human_b721` because its
registry row has `pack=None` (depth comes from an external reference) while RNA is its whole
input. `RNA_FROM_SRR_MAP` fills Brain from `data/loto_rnaseq_srr_tissue.tsv`.

`rna_coverage_total` was validated as a raw total by comparing against `coverage.npy.sum()` on the
int32 arrays, and `global_mean_coverage` confirmed to be exactly `total / sum_L`.

### ATTN profile exemplars (`select_profile_exemplars.py` / `plot_profile_exemplars.py`)

Re-ran the selector with `--model attn --top 0`. **`--top 0` is required** whenever the output
feeds a class-resolved figure: the default `--top 200` keeps the highest-scoring 200 rows per
dataset across all classes, and annotated CDS crowd out every scarce class. Under the default,
human Hepatocytes showed no `novel` / lncRNA row at all; with `--top 0` there are 413 lncRNA
transcripts clearing the same gates (>= 200 P-sites, ORF >= 90 nt). The gates never excluded them.

Per-nt exemplars regenerated at `--per-class 3` over all 7 classes, at both `--zoom-mode start`
and `--zoom-mode center`. The new file is a strict superset of the previous 4-class one: same
class winners, and 19,799 shared per-nt rows verified value-identical.

### Portable architecture spec (`scripts/export_arch_spec.py`, new)

Emits `figures/arch_attn/arch_spec_attn.json` and `figures/arch_attn/orf_track_spec.json` so a
build that cannot read the cluster filesystem can trace every printed structural number to a file
rather than to a prose table. Everything is **read from artifacts, never transcribed**: parameter
count summed from `best.pt`'s state_dict (5,071,106), channel and kernel shapes from the tensor
shapes, hyperparameters from the run's `args.json`, `START_W` imported from `build_orf_track.py`,
and the display-window ORF track computed by calling `orf_track()` rather than copied.

Two things the spec records that the prose table did not:

- **The 4,093 nt receptive field is the CONVOLUTIONAL one only.** `1 + 2*(k-1)*sum(dilations)` =
  `1 + 4*1023`. The mixer is full self-attention masked for padding alone
  (`src_key_padding_mask=~mask`, no causal and no local window), so the attn model's effective
  context is the **entire transcript**. "Receptive field ~4 kb" understates it.
- **The RNA coverage channel is `log1p(coverage / global_mean_coverage)`**, not raw depth
  (`cov_norm=global_mean`). It is always the last of the 10 input channels; the layout is
  one-hot(4) | orf_track(5) | coverage(1).

`START_W` has **10** graded start codons, not the 6 the figure doc recorded: ATG 1.0, CTG 0.5,
GTG 0.35, ACG 0.35, TTG 0.3, ATA 0.25, ATT 0.2, ATC 0.2, AAG 0.15, AGG 0.15. The deployed model
uses `--kozak none`, so the start channel is the raw codon weight with no context multiplier.

## Dataset citations: accession to source publication (2026-08-18)

Until now every dataset in this project was tracked by accession only, with **no source publication
recorded anywhere**. That gap was found while filling two `[author] et al. [year]` placeholders on
the poster and is closed here. Full detail, including drop-in reference text, is in
`docs/POSTER_MSG_REFS_4_AND_5.md`.

| dataset | accession | citation |
|---|---|---|
| Human Ribo-seq compendium (training, 9 tissues) | GSE182371 / GSE182372 (RNA) | Chothani et al. 2022 |
| Human iPSC-cardiomyocyte | ruizorera2024 | Ruiz-Orera et al. 2024 |
| Mouse liver (Wang) | GSE94982 | Wang H, ... Xie Z. 2021. *Nucleic Acids Res* 49(11):6165-6180. doi:10.1093/nar/gkab482. PMID 34107020 |
| Mouse liver (Janich) | -- | Janich P, ... Gatfield D. 2015. *Genome Res* 25(12):1848-1859. doi:10.1101/gr.195404.115. PMID 26486724 |
| Mouse liver (third) | GSE243134 | Nadimpalli HP, ... Gatfield D. 2024. *Genome Biol* 25(1):128. doi:10.1186/s13059-024-03270-2. PMID 38773499 |
| Human THP-1 | GSE208041 | Ansari SA, ... Uhlenhaut NH. 2022. *Comput Struct Biotechnol J* 20:5622-5638. doi:10.1016/j.csbj.2022.09.042. PMID 36284713 |
| Human CAR-T | GSE304796 | Shi Y, ... Sadelain M. 2026. *Cell* 189(2):401-417.e21. doi:10.1016/j.cell.2025.11.005 |
| Mouse CD4+ T cell | GSE155087 | Matheson LS, ... Turner M. 2022. *Sci Rep* 12(1):19657. doi:10.1038/s41598-022-24132-6. PMID 36385275 |
| Mouse BMDM | GSE120762 | Jackson R, ... Flavell RA. 2018. *Nature* 564(7736):434-438. doi:10.1038/s41586-018-0794-7. PMID 30542152 |
| Mouse tissue macrophage proteome (MS) | iProX IPX0001245001 = PRIDE PXD021583 | Qie J, ... Ding C. 2022. *Nat Commun* 13:7389. doi:10.1038/s41467-022-35095-7 |

GSE39561 (THP-1) is dropped and needs no citation.

**Two of these are inferred, not declared.** GEO lists a Citation PMID for GSE208041, GSE155087,
GSE120762 and GSE94982, so those four links are asserted by the submitters. **GSE243134 and
GSE304796 have no GEO citation**; both were matched on a verbatim-identical series title plus author
overlap (Nadimpalli/Gatfield; Shi/Sadelain). Strong matches, but flag them as inferred if either is
used to support a claim.

**Two caveats worth carrying:**

- **The CAR-T paper is 2026, not 2025.** GEO submission was Aug 2025; *Cell* publication was Jan 2026.
- **Two of the three mouse livers share a senior author.** Janich 2015 and GSE243134 2024 are both
  Gatfield lab, nine years apart; only Wang is independent of the other two. "Cross-study" is
  accurate, "three independent experiments" is not. This matters where the mouse-liver comparison is
  framed as independent replication: the arm the model wins on recall (Wang) is both the shallowest
  and the only lab-independent one.

## Per-tissue Chothani TPM recovery (2026-08-20)

`scripts/salmon_chothani_alnmode.sbatch` + `scripts/aggregate_chothani_tpm.py`.

**Inputs.** 55 STAR transcriptome BAMs in `data/rnaseq_bam_mm1/` (within-gene isoform multimappers
retained -- required, salmon's EM has nothing to apportion without them), against
`data/human_ribocode_annot_primary/transcripts_sequence.fa`, which matches the BAMs' `@SQ` space
exactly at 509,650 sequences. SRR-to-tissue map at
`data/tpm/chothani_alnmode/srr_tissue_have.tsv`.

**Command.** `salmon quant -t <txfa> -l A -a <bam> -o <out> -p 8 --noErrorModel`, one array task per
run, 8 concurrent. Per-run `quant.sf` pooled to an unweighted per-tissue mean.

**Method difference from the original, which must be carried with any use of these numbers.**
CORRECTED 2026-08-21: an earlier version of this section blamed the decoys alone. That understated
it. There are THREE differences, and the read-population one is larger:

1. **STAR ran with `--outFilterMultimapNmax 1`, so genomically multi-mapping reads were DISCARDED
   before salmon existed. Measured loss: 18.6% of reads** (143.6M of 176.4M over six sampled runs,
   range 77-86% kept). Salmon-from-FASTQ would have seen and apportioned those.
2. Reads STAR failed to map are absent; from FASTQ, salmon gets its own attempt at them.
3. No decoys.

NOTE the NH tags in these BAMs run to 17, which looks like multimappers were kept. They were not:
in a TRANSCRIPTOME BAM, `--quantMode TranscriptomeSAM` expands ONE genomic alignment across all
compatible isoforms, so NH counts transcripts of one gene, not genomic loci. Isoform multimapping
survives (which is what the within-gene EM needs); genomic multimapping does not.

These are therefore **STAR-unique-read TPMs**, not a reproduction of the original. The lost
`salmon_quant_chothani_decoy` was **decoy-aware selective alignment from FASTQ**. This is
**alignment-based mode against the transcriptome only, no decoys.** Decoys exist to stop intronic
and repeat-derived reads being assigned to transcripts; without them those reads land on real
transcripts instead. Validated on Fibroblast (the only tissue whose decoy-aware table survived),
same 32 libraries both ways: **Spearman 0.755 on expressed transcripts, median ratio 0.84x, no
length dependence (-0.093)**. Use for **isoform ranking within a gene**; do NOT quote as absolute
TPM and do NOT compare against numbers produced by the original pipeline. **Treat paralog-rich loci
with extra suspicion** -- GTF2I sits in the chr7 Williams-syndrome segmental duplication alongside
GTF2IRD1/GTF2IRD2, exactly where `Nmax 1` drops the most reads.

SUPERSEDED 2026-08-21 -- confirmed to mis-assign isoforms at paralog-rich loci; see results.md
"chothani_alnmode is SUPERSEDED". Gene-level totals usable, isoform-level NOT.
interim recovery, retained for comparison.

**Coverage.** 7 of 8 tissues. Brain has no surviving RNA BAM.

**Why not just re-download.** Measured 2026-08-20: ENA 0.30 MB/s, AWS S3 0.6-4.3 MB/s, and four
concurrent streams gave 2.4 MB/s **aggregate** -- worse than one, so the limit is our shared egress
and no download strategy avoids it. 67 GB would have taken 30-50 hours. **Check
`data/rnaseq_bam_mm1/` for a surviving transcriptome BAM before ever planning an RNA re-download.**

## mm10 RNA coverage: two retrain arms, A and B (2026-08-23)

Follows the multimap audit in results.md (2026-08-23). That audit established two separate facts,
and this experiment separates the two fixes they imply rather than bundling them.

**Fact 1, the model INPUT is corrupted at multimapped loci.** The RNA-seq coverage channel was
built from mm1 alignments. Where a transcript overlaps a repeat or a paralog, mm1 discards the
reads and the channel reads near-zero even though the RNA is abundant. For GTF2I
`ENST00000901263.1` a 1,617 nt window (transcript positions 1682:3299) reads 8.9% of that
transcript's CDS mean at mm1 and 138.3% at mm10. The earlier "the model memorised a blank CDS"
reading of this was wrong: the model read a corrupted input correctly.

**Fact 2, the training LABEL is corrupted at the same loci.** The Ribo target is also mm1, and
must stay mm1 -- RiboCode does not drop multimappers and they destroy the 3-nt periodicity P-site
calling depends on. So the hole is in both channels at the same coordinates by construction. The
same GTF2I window holds 1 observed P-site out of 13,258 on that transcript.

### The two arms

Both use the union universe (84,472 tx), one-hot backend, `--kozak none`, noBrain, holdout
Hepatocytes, mamba4 (`--mixer mamba --n_attn_layers 4`), seed 0 -- identical to the standing
baseline `results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes`
(best val Pearson 0.6190 @ epoch 23) in every respect except what each arm names.

| arm | RNA coverage | Ribo target | training set | isolates |
|---|---|---|---|---|
| baseline | mm1 | mm1 | full | -- |
| **A** | **mm10** | mm1 | full | the input fix |
| **B** | **mm10** | mm1 | **7,640 tx dropped** | the label fix, on top of A |

A vs baseline changes one thing (the coverage channel). B vs A changes one other thing (the
training set). The test split is identical in all three, so all three are directly comparable.

### Arm A: mm10 coverage

57 Chothani RNA libraries re-aligned with `--outFilterMultimapNmax 10` (`--outSAMtype None
--quantMode TranscriptomeSAM`), then `rnaseq_coverage.py` per library
(`scripts/rebuild_coverage_mm10.sbatch`, job 36859009, 57/57 COMPLETED, 213.6 GB of hd5). Packs
built by `scripts/rebuild_packs_mm10cov.sbatch` via `prepare_pack.py --reuse-target`, which reuses
the source pack's universe and `target_counts.npy` verbatim and writes only new coverage. Each job
asserts with `cmp` that `target_counts.npy` is byte-identical to its `packed_union` source and
fails hard otherwise, so a moved label cannot be mistaken for an input effect.

mm10 rather than mm25: at GTF2I, mm10, mm25, and best-alignment-only give identical numbers
(565.91/nt in the gap, 138.3% of CDS), so mm10 already recovers everything mm25 would and keeps
the count inflation lower. Mean RNA multimap rate across the 57 libraries is 3.8%, so this is a
small global change concentrated at specific loci -- aggregate metrics are expected to move little
and locus-level behaviour is where the effect should appear.

### Arm B: drop multimap-corrupted transcripts from training

Exclusion set: `frac_discarded > 0.5` in `data/mm25_diagnostic/training_multimap_cost.tsv.gz`
(`bamlib.multimap_cost`, mm1 vs mm25 total records per transcript over the 74 Chothani Ribo
libraries), restricted to the union universe. 7,640 tx of 81,611 (9.36%), listed in
`data/mm25_diagnostic/train_exclude_mmcost_gt0.5.txt`. These are transcripts whose mm1 Ribo label
is missing more than half its reads.

Implemented as `dataset.train_excluded_tx()`, env `RIBO_TRAIN_EXCLUDE_TX`, applied in `loto_split`
and `load_split`. It is deliberately **not** the existing `representative_tx()` hook, which is an
include-list applied to every split:

- **Dropped from train and val.** Val too, because checkpoint selection is on val Pearson; leaving
  them in val would make arm B select checkpoints against a different criterion and introduce a
  second difference between the arms. Val metrics are therefore NOT comparable across arms.
- **Kept in test.** They are exactly what the resulting model must be scored on, and keeping them
  holds the test split identical across all three arms. Measured effect: train 420,348 -> 387,959
  (tissue, tx) pairs (-7.71%), val 45,102 -> 41,295 (-8.44%), test 70,883 -> 70,883 (unchanged).
- `heldout_test_tx()` is untouched -- it is test-only by design.

### Prediction registered before the runs

Arm A supplies true coverage across windows where the label still says approximately zero, so at
multimap sites it trains the association "abundant RNA, no ribosome". Over the other ~90% of the
universe the corrected input is strictly better information. The net is genuinely uncertain and A
may regress. Arm B removes the transcripts that teach that association and is expected to improve.
A regression in A alongside a gain in B is a clean result, not a failure.

### Success criteria, registered before the runs

1. The standing ORF-call evaluation: RiboCode drop-in on both arms (deterministic theta=1 and
   `--pred_poisson` at CDS-anchored theta), per-condition composition and F1, expression-restricted
   and stating `n_ref`.
2. The locus check: GTF2I `ENST00000901263.1`, predicted density in positions 1682:3299 relative to
   its CDS mean. A model that still paints that window empty has not been reached by the input fix,
   whatever the aggregate numbers say.

### Provenance note on the reused target

`packed_union*/provenance.json` for several tissues is a reconstructed record with
`inputs_known: false` -- those packs predate provenance recording and their source hd5 are partly
deleted. `--reuse-target` copies those bytes forward, so all three arms share one fixed inherited
Ribo target. The `cmp` assert proves the target did not change between arms; it does not prove the
target is reproducible from source. That is a known and accepted limitation of this comparison.

## Cross-species expansion: datasets, references, annotation (2026-08-24/25)

Context: the per-nt model has only been trained and held out on human and mouse. The expansion adds
gorilla, chimp, macaque, zebrafish, C. elegans, D. melanogaster and S. cerevisiae. Plan:
`~/.claude/plans/abstract-hatching-peach.md`. Scope, endpoint and matrix were fixed with the user:
all seven species, endpoint = built packs, full 5-way alignment matrix, and the 12 undifferentiated
iPSC RNA runs of PRJEB65856 dropped because they have no Ribo partner.

Sections below, in the order the work runs:
**A** scope, **B** dataset selection, **C** download and integrity, **D** adapter determination,
**E** annotation normalization (the bulk), **F** reference construction.
Operational companion, including what is running and how to resume:
`docs/HANDOFF_XSPECIES_2026_08_24.md`.

### A. Scope

144 runs, 305.3 GiB, seven species, all matched Ribo + RNA. Endpoint is built packs and the full
5-way alignment matrix (RNA mm1 + mm10 + salmon, Ribo mm1 + mm25).

| species | Ribo | RNA | GiB | project |
|---|---|---|---|---|
| human (Ribo backfill only) | 5 | 0 | 14.5 | PRJEB65856 |
| gorilla | 3 | 3 | 35.4 | PRJEB65856 |
| chimp | 8 | 8 | 66.8 | PRJEB65856 |
| macaque | 7 | 7 | 61.0 | PRJEB65856 |
| zebrafish | 8 | 24 | 69.0 | PRJNA200706 + PRJNA288987 |
| C. elegans | 11 | 10 | 26.4 | PRJNA230441 + PRJNA230374 |
| yeast | 16 | 16 | 19.2 | PRJNA726548 |
| fly | 9 | 9 | 12.9 | PRJNA390134 |

The **human arm is a backfill**: its five Ribo runs existed only as BAMs and P-sites, with no FASTQ
locally and none in the warm archive, which would have left human as the one iPSC-CM species with no
mm25 diagnostic and no ability to re-trim.

### B. Dataset selection

Authoritative list: `data/external/xspecies/download_manifest.tsv`, written by
`scripts/xspecies/select_runs.py` from the full ENA run tables in `data/external/xspecies/PRJ*_runs.tsv`
(fetched by `build_run_tables.sh`, 18 columns including `fastq_md5`, which the old 4-column
`PRJEB65856_runs.tsv` lacked entirely, so nothing could be verified from it).

Every exclusion is expressed as a **rule in the selection script**, not a typed accession list, so the
reasoning survives and an ENA metadata correction is picked up by re-running. `--check` asserts the
expected per-species counts and exits nonzero on drift.

The script deliberately does not trust `library_strategy` to separate Ribo from RNA: that vocabulary
term postdates most of these submissions, and zero worm / fly / yeast / zebrafish runs carry it.
Assay is decided per project by whatever is actually reliable there, and the basis is recorded in an
`assay_basis` column of the manifest.

#### The PRJEB65856 design, recovered

`sample_title` is empty for every run in this study, which is why the earlier handoffs had no design
table. `sample_alias` and `sample_description` are populated and give it in full, as
`<species>_<tissue>_<individual>_{mR,Ri}`:

| species | tissue | Ribo | RNA | matched pairs | Ribo len | RNA len |
|---|---|---|---|---|---|---|
| gorilla | iPSC-CM | 3 | 3 | 3 | 51 nt | 75 nt |
| chimp | left ventricle | 5 | 5 | 5 | 51 nt | 75 nt |
| chimp | iPSC-CM | 3 | 3 | 3 | 51 nt | 75 nt |
| macaque | left ventricle | 4 | 4 | 4 | 51 nt | 75 nt |
| macaque | iPSC-CM | 3 | 3 | 3 | 35 nt | 100 nt |
| human | iPSC-CM | 5 | 5 | 5 | 35 nt | 100 nt |
| gorilla / chimp / macaque / human | iPSC | 0 | 3 each | 0 | - | - |

23 matched Ribo/RNA pairs. The 12 undifferentiated-iPSC RNA runs are RNA-only with no Ribo partner
and were dropped (user decision), saving ~43 GiB.

Two consequences that shape everything downstream. **`library_selection = Oligo-dT` on every RNA run**,
so the poly(A) requirement is verified rather than assumed, with no GEO protocol reading needed. And
the study contains **two trimming regimes split by read length inside one project**, which is why the
adapter policy is per run (section D). The clean four-species comparison is iPSC-CM (human 5,
gorilla 3, chimp 3, macaque 3); left ventricle is a bonus two-species set (chimp 5, macaque 4) with
no gorilla and no human arm.

#### The four organisms `SPECIES_EXPANSION_PLAN.md` section D left unresolved

- **zebrafish** -- GSE70549's own GEO summary states verbatim that it is "sample matched to a
  previous ribosome profiling dataset GSE46512". Eight developmental stages, matched 1:1 by stage
  label, 24 RNA (PE 2x76, 3 reps x 8 stages) against 8 RPF. The Ribo arm is labelled `RNA-Seq`; it is
  identifiable by **`library_selection = "size fractionation"`**, which is a better within-study
  discriminator than `library_strategy` and is worth adding to the metadata-trap notes.
- **fly** -- PRJNA390134 / GSE99920, three genotypes x three replicates, both arms in one project.
  The three TRAP runs are excluded: affinity-tagged ribosome pulldown is a different assay and cannot
  serve as a per-nt Ribo target.
- **yeast** -- PRJNA726548 / GSE173654, eight strains x two replicates, `RNA_<strain>-<rep>` /
  `Ribo_<strain>-<rep>` pairing exactly.
- **C. elegans** -- the plan's PRJNA1183817 was **rejected**. Its RNA arm is Total RNA-seq, the exact
  failure mode that crushed the GSE243134 liver universe to 3,321 tx, and its RNA libraries are 0.11
  to 0.30 GB against 0.5 to 4.6 GB for its Ribo arm, so the RNA is shallow as well as wrong.
  Replaced by **GSE52861 + GSE52905** (Hendriks and Gaidatzis), whose RNA series title is literally
  "RNA-seq for polyA enriched mRNAs" and whose protocol names the TruSeq stranded mRNA kit; Ribo is
  gel-excised 28-30 nt monosome RPFs. Eight shared hourly timepoints, 22 to 36 hr.
  **Caveat to carry: the arms are replicate timecourses A (RNA) and B (Ribo), matched by stage rather
  than by lysate**, so no downstream analysis may assume paired samples.

Rejected while searching, recorded so they are not re-tried: **GSE269587/589** (worm hibernation,
attractive same-lysate design but the RNA kit is ScriptSeq v2 on rRNA-depleted total RNA),
**GSE80157/80133** (let-7/LIN41, same problem), **GSE300687** (disome-seq, not monosome footprinting).

### C. Download and integrity

`scripts/xspecies/fetch_from_manifest.sh`, serial on the head node, one connection, `wget -nv -c`,
md5-verified against the manifest, idempotent (a file already matching its md5 is skipped, so an
interruption resumes). Derived from `fetch_heldout_rna.sh` with three changes, each fixing a real
failure mode in that script: no per-accession ENA round trip (it calls the API inside
`set -euo pipefail`, so one transient curl failure aborts the chain and loses every remaining
accession, and it parses only the first response line so a study accession silently fetches one run);
the redundant merged `<ACC>.fastq.gz` is skipped when `_1`/`_2` are present; and an empty `fastq_ftp`
is a loud error rather than a silent skip.

**A CRLF bug caught by the dry run, worth recording.** Python's `csv` writer defaults to
`lineterminator="\r\n"`, which makes the last column of every manifest row end in a carriage return.
The last column is `fastq_md5`, so **every checksum comparison would have failed on every file**.
Fixed at the writer, and the fetcher now refuses a manifest with CRLF endings rather than fail
100% of its verifications.

**ENA's FTP endpoint truncates.** Two files exhausted all three FTP retries, and the failure mode was
truncation rather than corruption -- each attempt returned a different, shorter length:

| run | got | expected | |
|---|---|---|---|
| ERR12549929 | 1,816,666,112 | 3,227,905,296 | 56% |
| SRR5667268 | 2,016,301,056 | 2,148,789,443 | 94% |

The same files' HTTPS endpoint returns a correct `Content-Length`, so FTP is the flaky leg.
`scripts/xspecies/verify_and_repair.sh` md5-checks everything on disk and refetches failures over
**HTTPS**, deleting first (never resuming onto truncated bytes, which `wget -c` would do and could
never recover from). Both files were repaired on the first HTTPS attempt. It skips whichever label
the fetch chain is currently inside, detected from the driver log, to keep one writer per path.

Run it once at the end regardless of whether failures were seen: an unverified FASTQ becomes a
library with silently missing reads, not an error.

### D. Adapter determination is per run, and measured

`scripts/xspecies/probe_adapters.sh` -> `data/external/xspecies/adapter_probe.tsv`, one row per run
with the measured adapter percentages, read-length profile, N-padding fraction, verdict and the exact
cutadapt argument. It probes the **complete** file via `zcat | head` (a streamed byte-range slice once
reported 0% adapter where the real file was 93%, because a mid-file gzip offset does not decompress to
record boundaries).

**It tests a panel, not just TruSeq, and that is not a formality.** The standing cluster rule names
`AGATCGGAAGAGC`, and every Ribo-seq dataset this project had processed used it. The worm arm does not:

| worm GSE52905 | hit rate |
|---|---|
| `TGGAATTCTCGGGTGCCAAGG` (Illumina small-RNA) | **57 to 84%** |
| `AGATCGGAAGAGC` (TruSeq) | **0.0%** |

Trimming those with the standing default would have removed nothing, left ~20 nt of adapter on every
51 nt read, and then `-m 20 -M 40` would have **discarded** them for being too long -- a near-empty
BAM with no error anywhere, the same shape as the Janich library that once mapped 0.03%.

Two details the probe encodes, both learned by getting them wrong first:

- It selects the **longest** panel member above threshold, not the highest-scoring. A 13 nt prefix
  always out-scores the full 21 nt adapter because it also matches reads where the adapter runs off
  the read end (measured: 95.9% vs 83.2%, same adapter). cutadapt matches partial 3' adapters itself,
  so the full sequence is strictly the better argument.
- It decides `--trim-n` on **N-padding alone**. The pre-clipped human arm is 99.7% N-padded but its
  read lengths run 35-51, so an initial fixed-length test missed it.

Measured verdicts, all seven species:

| label | verdict | evidence |
|---|---|---|
| human iPSC-CM Ribo | `--trim-n` | 98.6-99.7% N-padded, 0% adapter |
| gorilla Ribo | TruSeq | 72.7-79.5% |
| macaque Ribo | **MIXED** | 4 runs TruSeq 75.4-76.3% (51 nt LV), 3 runs `--trim-n` 99.3-99.5% N-pad (35 nt CM) |
| chimp Ribo | TruSeq | pending download completion |
| worm Ribo | small-RNA | 57.3-83.5% |
| fly Ribo | already trimmed | 0% adapter, modal 30 nt (footprint-sized), variable length |
| yeast Ribo | already trimmed | 0% adapter, modal 24-25 nt |
| all RNA arms | none needed | 0% adapter at these insert lengths |

`primate_rm_ribo` is **the only label with mixed verdicts**, confirming empirically what the read
lengths had only implied: a per-dataset adapter setting would have been wrong for 3 or 4 of its 7
runs. The human result independently agrees with `heldout_config.sh`'s separately-recorded
`ADAPTER_MODE=trimn` for `human_ruizorera`.

### E. Annotation normalization

#### The problem: `gene_type` is a GENCODE-only spelling

Every non-GENCODE reference in scope uses `gene_biotype` / `transcript_biotype`. Measured:

| GTF | `gene_type` hits | `gene_biotype` hits |
|---|---|---|
| `genomes/xspecies_refs/yeast/genomic.gtf` (RefSeq) | **0** | 6,477 |
| `genomes/xspecies_refs/gorilla/genomic.gtf` (RefSeq) | **0** | 31,925 |
| `Saccharomyces_cerevisiae.R64-1-1.116.gtf.gz` (Ensembl) | **0** | 41,879 |
| `Danio_rerio.GRCz11.116.gtf.gz` (Ensembl) | **0** | 1,161,865 |

`scripts/build_contaminant_index.sbatch` matches `$9 ~ ("gene_type \"(rRNA|Mt_rRNA|...)\"")`, so on
every new reference it selects **zero** lines. Its `DEPLETE` BED has a `[[ -s ]]` guard; its
`POSTFILT` BED has no non-empty check at all and would be written as a 0-byte file. The standing
rule that rRNA/tRNA/miRNA/Mt loci never reach an ORF caller would silently not execute, on all seven
species, with no error anywhere.

#### Ensembl was evaluated as an alternative and rejected

Checked against Ensembl release 116 (current). Ensembl shares RefSeq's attribute key, so it does not
solve the problem; its vocabulary is an older GENCODE variant (`lincRNA`, `antisense`,
`processed_transcript`) that still needs mapping; and for the three primates its assemblies are far
behind the RefSeq T2T builds already on disk:

| species | Ensembl 116 | RefSeq on disk |
|---|---|---|
| gorilla | gorGor4, 2014-12 | T2T NHGRI_mGorGor1-v2.1, 2025 |
| chimp | Pan_tro_3.0, 2016-05 | T2T NHGRI_mPanTro3-v2.1, 2026 |
| macaque | Mmul_10, 2019-02 | T2T-MMU8v2.0, 2025 |
| zebrafish | GRCz11 | GRCz11 (same assembly) |
| C. elegans | WBcel235 | WBcel235 / WS298 (same assembly) |
| yeast | R64-1-1, GTF frozen at release 63 | R64 / SGD R64-5-1, 2026 |

Decision: keep RefSeq for the six fetched species and FlyBase r6.67 for the fly, and normalize.

#### `scripts/normalize_annotation.py`

One general script, four conventions, GENCODE as the target vocabulary.

| convention | key | notes |
|---|---|---|
| gencode | `gene_type` / `transcript_type` | the target |
| refseq | `gene_biotype` / `transcript_biotype` | no `Mt_*` classes at all; `mRNA` not `protein_coding` at transcript level; generic `transcript` defers to the gene; spelling is not self-consistent (`lnc_RNA` in gorilla and zebrafish, `lncRNA` in chimp, macaque, worm) |
| ensembl | `gene_biotype` / `transcript_biotype` | older GENCODE vocabulary; does carry `Mt_*` |
| flybase | none | biotype is the feature-type column, and `mRNA` appears where a GTF normally writes `transcript` |

Source detection is automatic (probe for the attribute key, then RefSeq vs Ensembl by transcript-ID
prefix). **Any source term with no explicit mapping is a hard error, never a silent pass-through**,
which is what makes the script safe to point at an eighth species. Terms with no GENCODE equivalent
(`piRNA`, `V_gene_segment`, `transposable_element`) are marked PASSTHROUGH as a recorded decision,
not a fallback, and are in no drop set.

Outputs per species into `data/annot/<species>/`: `ncrna_tx.txt`, `tx_to_gene.tsv`,
`tx2biotype.tsv`, `ncrna_deplete.bed`, `ncrna_postfilter.bed`, `normalization_report.json`, and
optionally `<species>.normalized.gtf` (`--write-gtf`, needed by STAR, RiboCode and gffread; for
FlyBase it also re-emits each transcript feature as `transcript` so those tools find it).

#### Two drop sets, and why they differ deliberately

`core` = `{rRNA, Mt_rRNA, rRNA_pseudogene, Mt_tRNA, miRNA}`. This is the set that **exactly**
reproduces the hand-built ground truth already on disk, which had no generator script:

| derived from | reproduces | rows | result |
|---|---|---|---|
| `gencode.v49.annotation.gtf` | `data/ncrna_tx_human_v49.txt` | 2,447 | exact match, 0 missing, 0 extra |
| `gencode.vM38.annotation.gtf` | `data/ncrna_tx_mouse_vM38.txt` | 2,579 | exact match, 0 missing, 0 extra |
| `gencode.v49.annotation.gtf` | `data/tx_to_gene_human_v49.tsv` | 507,365 | exact match |
| `gencode.vM38.annotation.gtf` | `data/tx_to_gene_mouse_vM38.tsv` | 278,326 | exact match |

Run it with `--selftest`. This gate comes before any new species: a normalizer that cannot
regenerate the two lists already on disk has no business being pointed at one where no ground truth
exists.

`extended` = core plus `{tRNA, tRNA_pseudogene, SRP_RNA, RNase_P_RNA, RNase_MRP_RNA}`, and is the
default for every non-GENCODE source. The difference is intentional, not a discrepancy: GENCODE's
main GTF contains **no cytoplasmic tRNA** (they ship separately in `gencode.*.tRNAs.gtf.gz` and are
handled by the bowtie2 contaminant index) and types 7SL, RNase-P and RNase-MRP as plain `misc_RNA`.
RefSeq and FlyBase annotate all of them inline, so if they are not dropped here they are not dropped
anywhere. The magnitude is not marginal: **zebrafish alone has 8,839 tRNA transcripts**, 68% of its
entire drop list. 7SL is specifically the RNA that collapsed the GSE243134 liver universe to 3,321
transcripts when it survived a biotype filter as lncRNA.

#### Result on all seven references, zero unmapped terms

| species | source | mito contig | tx | ncRNA dropped | composition |
|---|---|---|---|---|---|
| yeast | refseq | NC_001224.1 | 6,477 | 316 | tRNA 275, Mt_tRNA 24, rRNA 12, Mt_rRNA 2, RNase_P 1, SRP 1, RNase_MRP 1 |
| C. elegans | refseq | NC_001328.1 | 56,731 | 1,586 | miRNA 720, tRNA 612, tRNA_pseudogene 209, Mt_tRNA 22, rRNA 20, Mt_rRNA 2, rRNA_pseudogene 1 |
| fly | flybase | mitochondrion_genome | 35,736 | 1,174 | miRNA 747, tRNA 290, rRNA 113, Mt_tRNA 22, Mt_rRNA 2 |
| zebrafish | refseq | NC_002333.2 | 82,520 | 12,980 | tRNA 8,839, rRNA 3,049, miRNA 1,066, Mt_tRNA 22, RNase_P 2, Mt_rRNA 2 |
| gorilla | refseq | NC_011120.1 | 99,500 | 1,475 | miRNA 740, tRNA 472, rRNA 239, Mt_tRNA 22, Mt_rRNA 2 |
| chimp | refseq | NC_001643.1 | 135,577 | 2,065 | miRNA 1,309, tRNA 500, rRNA 232, Mt_tRNA 22, Mt_rRNA 2 |
| macaque | refseq | **none** | 140,820 | 2,387 | miRNA 1,494, rRNA 464, tRNA 429 |

Every species shows the expected 22 Mt_tRNA + 2 Mt_rRNA (the standard metazoan mitochondrial gene
complement) except macaque.

#### The macaque assembly has no mitochondrion, and the guard is what caught it

`GCF_049350105.2` T2T-MMU8v2.0 ships **22 contigs: chr1-20, X, Y, and no MT sequence**. There is no
mitochondrial contig to detect, nothing to re-type, and correspondingly zero Mt_tRNA / Mt_rRNA in
its annotation. Its 464 rRNA and 429 tRNA entries are all nuclear, so the drop set is still correct.

The consequence is for alignment rather than annotation, and it is not cosmetic. With no chrM in the
index, mitochondrial footprints (a real and often large fraction of a Ribo-seq library) have nowhere
correct to map, and the likely destination is a NUMT, which a T2T assembly resolves in unusual
completeness. That is silent signal corruption on one of the four species in the iPSC-CM comparison,
and it also makes macaque the only species whose mito reads are handled differently from the rest.

Recorded as an open decision rather than patched, because adding `NC_005943.1` to the FASTA and GTF
modifies a reference shared across projects. Options: (a) accept and document, quantifying the
unmapped fraction after alignment; (b) supplement the reference with the RefSeq macaque
mitochondrion, which restores consistency with the other six species.

#### Four annotation defects the normalizer had to absorb (2026-08-24)

Each was found by a gate rather than by reading the GTF, and each would have produced a
plausible-looking but wrong result. They are recorded because the next species will have its own.

**1. `transcript_id ""` on every RefSeq gene line breaks gffread.** RefSeq writes a literal empty
`transcript_id` on `gene` features (6,477 of them in yeast alone); GENCODE omits the attribute.
gffread reads the empty string as the record's ID and aborts with
`Error: no valid ID found for GFF record`. This happens on the **raw** RefSeq GTF, before
normalization touches it, so without stripping the attribute there is no transcriptome FASTA for any
RefSeq species, and therefore no salmon index. The normalizer now removes it.

**2. FlyBase transcripts were extracted twice.** FlyBase writes `mRNA` where a GTF normally writes
`transcript`, and RiboCode and gffread both look for `transcript`. The first implementation emitted
the feature a second time under the new name, which made gffread extract every fly transcript twice:
**71,107 sequences for 35,736 transcripts, 35,371 duplicate IDs**. The salmon index built from it
reported `num_refs=58133`, a number that looks entirely reasonable in isolation. Fixed by RENAMING
the feature instead of duplicating it; the biotype the feature name carried is preserved in the
`transcript_type` attribute the normalizer adds. Caught by the registry's transcript-count column,
not by any error.

**3. Nine FlyBase transcripts have strand `.`, and RiboCode dies on the first one.** All nine are
`mod(mdg4)`, the trans-spliced Drosophila gene whose mature transcripts are assembled from both
strands, so FlyBase legitimately declines to give them a strand. `prepare_transcripts` raises
`ValueError: strand is neither "+" nor "-"` and the annotation build fails outright. They are now
dropped with an explicit count and ID list in the report, since a single-strand transcriptome model
cannot represent them anyway.

**4. RefSeq C. elegans reuses transcript IDs across loci.** `unassigned_transcript_572` and
`unassigned_transcript_574` each appear at two distinct loci. `tx_to_gene` is keyed by transcript id,
so the second occurrence silently OVERWROTE the first and the map came out with 56,731 rows, while
gffread emitted both and the FASTA had 56,733. The build's join check compared **sets**, so it
subtracted to zero and passed. Both copies of a reused id are now dropped (which locus is genuine is
not knowable from the GTF, and a silently-wrong pick is worse than losing two placeholder
transcripts), and the check is being tightened to compare counts, not just membership.
Only C. elegans is affected: gorilla, chimp, macaque, zebrafish and yeast have zero duplicated ids.

Consequence for the work: the worm and fly STAR and salmon indexes built before fixes 2 and 4 were
discarded and rebuilt. Nothing had been aligned against them yet.

#### The macaque reference was supplemented with its mitochondrion (2026-08-25)

`GCF_049350105.2` T2T-MMU8v2.0 ships 22 contigs (chr1-20, X, Y) and **no mitochondrial sequence**.
Every other species in the expansion carries the standard metazoan 22 Mt_tRNA + 2 Mt_rRNA; macaque
had zero, which is how the normalizer's mito guard surfaced it. The drop set was still correct (its
464 rRNA and 429 tRNA are all nuclear), but the alignment consequence is real: with no chrM in the
index, mitochondrial footprints have nowhere correct to map and the likely destination is a NUMT,
which a T2T assembly resolves in unusual completeness. That is silent signal corruption on one of
the four species in the iPSC-CM comparison, and it would have made macaque the only species whose
mito reads were handled differently from the rest.

Resolved on user decision by appending **NC_005943.1**, the RefSeq reference mitochondrial genome
for *Macaca mulatta* (16,564 bp, 37 genes = 13 protein-coding + 22 tRNA + 2 rRNA, TaxId 9544,
updated 2023-04-03). That is the same gene complement the chimp mitochondrion `NC_001643.1` already
contributes, so macaque is now structurally at parity with the other primates.

Done by `scripts/xspecies/add_mito_contig.py --species macaque --accession NC_005943.1`, which is
general rather than macaque-specific, is idempotent (it refuses to append a contig already present),
and records the change in `PROVENANCE.json` with before/after md5s of both files and the exact
command to reproduce it.

**Two caveats recorded in the provenance**, because macaque's reference is no longer a single
unmodified NCBI download: NC_005943.1 is a separate RefSeq record rather than part of the T2T
assembly, and it derives from a different individual than the assembly (isolate MMU2019108-1), so a
few footprints will carry mismatches against it. At Ribo-seq read lengths and STAR's default
mismatch tolerance that is immaterial next to having no mitochondrion at all.

##### The incomplete stop codon, again, and how the GTF handles it

**6 of the 13 macaque mito CDS have no stop codon in the genome.** NCBI annotates this explicitly
with `transl_except=(pos:N,aa:TERM)` and the note "TAA stop codon is completed by the addition of
3' A residues to the mRNA"; their CDS lengths are 1 mod 3. This is the same phenomenon that killed
the fly build on mt:ND2 / mt:CoII / mt:ND4 / mt:ND5.

The converter therefore emits a `stop_codon` only for the 7 CDS that actually have one (taking the
final 3 bases and excluding them from the CDS feature, which is the RefSeq GTF convention), and for
the other 6 emits the CDS whole with no stop_codon. Writing a stop codon for those would place it
outside the transcript's exon, which is precisely the condition that makes RiboCode's
`prepare_transcripts` die with "Can't transform the genomic interval".

The script asserts every emitted codon feature lies inside its transcript's exons before it writes
anything, and the output was validated against RiboCode in isolation first: `prepare_transcripts`
on the mito-only FASTA and GTF accepts all 37 transcripts, including the 6 without a stop codon.
Emitted: 37 gene, 37 transcript, 37 exon, 13 CDS, 13 start_codon, 7 stop_codon.

The macaque STAR index, salmon index and RiboCode annotation were rebuilt, since the originals
predate the mitochondrion.

### F. Reference construction

`scripts/xspecies/build_species_refs.sbatch`, one SLURM array task per species, five steps each with
an idempotent skip guard so a re-run resumes: normalize -> transcriptome FASTA -> STAR index ->
decoy-aware salmon index -> RiboCode annotation.

Array index: `0 yeast, 1 celegans, 2 fly, 3 zebrafish, 4 macaque, 5 chimp, 6 gorilla`.
Parameters come from `scripts/xspecies/species.tsv`, the single table also read by
`riboseq_align.sbatch`, so an annotation version cannot drift between the index that was built and
the index that gets used. `riboseq_align.sbatch` now resolves any species with a row in that table
automatically; human and mouse keep their historical hardcoded paths, because those packs were built
against those exact files and re-pointing them is a re-validation, not a refactor.

The reference tree was renamed `genomes/primates_t2t` -> **`genomes/xspecies_refs`**, since it holds
worm, fish, fly and yeast as well as the primates, and fly was added as symlinks to the canonical
FlyBase copy. **Never point `archive_to_warm.sh` at `xspecies_refs/fly`**: it resolves with
`readlink -f` and would archive and then delete the FlyBase originals.

#### Two index parameters that must be per species

`--sjdbOverhang` is max RNA read length minus 1, taken from the download manifest.
`--genomeSAindexNbases` is `min(14, log2(genome_bp)/2 - 1)`; **14 is the mammalian default and is
wrong for the small genomes**. STAR does not error on a bad value, it silently builds a worse index,
so the registry reads both back out of `genomeParameters.txt` rather than trusting that they were
passed.

| species | genome bp | sjdbOverhang | genomeSAindexNbases |
|---|---|---|---|
| yeast | 12,157,105 | 132 | **10** |
| C. elegans | 100,286,401 | 49 | **12** |
| fly | 143,726,002 | 72 | **12** |
| zebrafish | 1,679,203,469 | 75 | 14 |
| macaque | 3,115,120,525 | 99 | 14 |
| chimp | 3,177,756,316 | 74 | 14 |
| gorilla | 3,545,850,636 | 74 | 14 |

Two further STAR settings are non-default and required at this scale: `--limitGenomeGenerateRAM
90000000000` (the 31 GB default is too low for a 3.5 Gb genome) and `--limitSjdbInsertNsj 4000000`
(the 1e6 default is exceeded by RefSeq mammalian annotation).

#### The transcriptome FASTA comes from gffread, not NCBI's `_rna.fna`

NCBI's `rna_from_genomic` headers are `lcl|NC_001133.9_mrna_NM_001180043.1_1`, not bare accessions,
so they do not join to the GTF or to tx2gene without rewriting -- a silent failure for salmon, whose
quantification would simply not aggregate. Deriving with `gffread -w` from the **normalized** GTF
guarantees the transcript IDs are byte-identical to the ones STAR, RiboCode and tx2gene use, and the
build asserts that every transcriptome ID joins to `tx_to_gene.tsv` before continuing.
`gffread` is not in the `riboseq` env; it lives in `conda_envs/{crisprware,sqanti3}`.

The salmon gentrome and `decoys.txt` are **kept beside the index** this time. They were not kept for
the human and mouse indexes, and there is no `cmd_info.json` either, so those builds cannot be
reproduced from what is on disk.

#### Verification

`scripts/xspecies/build_reference_registry.py` writes `docs/SPECIES_REFERENCE_REGISTRY.md` and
`data/species_reference_registry.tsv`, and with `--check` exits nonzero if any species fails a gate,
so it can block the alignment phase. It checks: every artifact present and non-empty; the STAR index
recorded the intended per-species parameters; the ncRNA drop list is above a floor (an empty one is
the silent failure the normalizer exists to prevent); the salmon index has decoys; RiboCode kept most
of the transcripts; and every transcriptome ID joins to tx2gene.

Compute note: `--partition=gpu` with no `--gres` was used throughout. The `medium` partition was
repeatedly backed up on fair-share priority while the gpu nodes had ~900 to 1100 idle CPUs and over a
terabyte of free memory.

#### Script inventory for this expansion

| script | does | gate it enforces |
|---|---|---|
| `scripts/xspecies/build_run_tables.sh` | pull full ENA run tables, 18 fields incl. `fastq_md5` | warns on any run with no md5 |
| `scripts/xspecies/select_runs.py` | apply the selection rules, write the download manifest | `--check` asserts expected per-species counts |
| `scripts/xspecies/preflight_survey.sh` | full-depth local + archive survey before any fetch | nonzero if anything in scope already exists |
| `scripts/xspecies/fetch_from_manifest.sh` | serial md5-verified ENA fetch | refuses a CRLF manifest; loud on empty `fastq_ftp` |
| `scripts/xspecies/verify_and_repair.sh` | md5-check everything, refetch failures over HTTPS | nonzero if anything is still bad |
| `scripts/xspecies/probe_adapters.sh` | per-run adapter panel + length profile | flags labels with mixed verdicts |
| `scripts/normalize_annotation.py` | four GTF conventions -> GENCODE vocabulary + drop lists | `--selftest` exact-match vs human/mouse ground truth; hard error on any unmapped biotype term |
| `scripts/xspecies/add_mito_contig.py` | append a mitochondrion from a separate RefSeq record | asserts every codon feature is inside its exons before writing |
| `scripts/xspecies/build_species_refs.sbatch` | the five-step per-species reference build | per-step asserts; RiboCode transcript-ratio floor |
| `scripts/xspecies/build_reference_registry.py` | write and verify the reference registry | `--check` blocks alignment on any failure |
| `scripts/xspecies/species.tsv` | the single per-species parameter table | read by both the build and the aligner, so versions cannot drift |

The recurring pattern is worth stating once: **every one of these fails loudly rather than producing
a plausible empty or truncated result.** The four annotation defects, the CRLF md5 bug, the FTP
truncations and the worm adapter were each caught by a gate, not by reading output, and every one of
them would otherwise have produced a number that looked entirely reasonable.

#### Ensembl re-examined on annotation CONTENT (2026-08-25)

The original Ensembl assessment (above) rejected it on two axes only: it uses `gene_biotype` exactly
like RefSeq, so it does not solve the attribute-key problem, and its primate assemblies are
gorGor4 (2014) / Pan_tro_3.0 (2016) / Mmul_10 (2019) against the RefSeq T2T builds. Both points
stand, but the write-up implied a broader evaluation than had actually been done: **annotation
content was never compared.** That gap was worth closing, because the content differences are
exactly where the biotype census found its oddities.

Ensembl release 116 GTFs were fetched for the four species where Ensembl and RefSeq share an
assembly, normalized through the SAME classifier, and counted the same way. (Note the file naming:
for non-vertebrates the main Ensembl FTP tree carries the Ensembl **Genomes** release number, so the
current fly/worm/yeast GTFs are `.63.`, not `.116.`.)

| species | source | total tx | coding | lncRNA | small ncRNA | dropped ncRNA | pseudogene | other |
|---|---|---|---|---|---|---|---|---|
| fly | FlyBase r6.67 | 35,723 | 30,823 | **0** | 3,361 | 1,174 | 365 | 0 |
| fly | Ensembl BDGP6.54 | 41,590 | 30,792 | **0** | 3,362 | 1,174 | 366 | 5,896 |
| C. elegans | RefSeq WS298 | 56,729 | 31,219 | **306** | 23,618 | 1,586 | 0 | 0 |
| C. elegans | Ensembl WBcel235 | 60,000 | 31,865 | **306** | 24,278 | 1,596 | 1,955 | 0 |
| zebrafish | RefSeq GRCz11 | 82,520 | 59,013 | 8,815 | 1,649 | **12,980** | 32 | 31 |
| zebrafish | Ensembl GRCz11 | 59,876 | 49,447 | 6,636 | 901 | **2,425** | 317 | 150 |
| yeast | RefSeq R64-5-1 | 6,477 | 6,039 | 8 | 114 | 316 | 0 | 0 |
| yeast | Ensembl R64-1-1 | 7,127 | 6,600 | 0 | 101 | 323 | 12 | 91 |

**The fly and worm lncRNA gap is NOT a RefSeq artifact and Ensembl does not fix it.** Ensembl fly
reports **0 lncRNA** and Ensembl worm reports **306**, identical to the other registry in both cases.
Both inherit the catch-all `ncRNA` biotype from FlyBase and WormBase respectively, which is where the
distinction is genuinely not drawn. No choice of registry recovers it, which is the stronger form of
the claim made in the census write-up.

**For zebrafish, RefSeq is substantially the richer annotation**, contrary to what might be assumed:
82,520 transcripts against 59,876, more lncRNA (8,815 vs 6,636), and critically **8,839 inline tRNA
genes that Ensembl does not annotate at all** (its whole drop set is 2,425 against RefSeq's 12,980).
For this pipeline the tRNA annotation is an asset, not noise: an annotated tRNA transcript can be
dropped from the transcriptome, whereas an unannotated tRNA locus still absorbs reads in the genome
and has nowhere to be filtered.

**Two places Ensembl is genuinely better, recorded honestly:**
- **Pseudogene typing.** Worm 1,955 vs 0 transcripts, zebrafish 317 vs 32, yeast 12 vs 0. RefSeq
  under-annotates pseudogenes in these species. This is not cosmetic for a Ribo-seq pipeline: an
  annotated pseudogene absorbs its own multimapping reads instead of donating them to the parent
  gene, which is precisely the cross-gene paralog ambiguity `filter_tx_heldout.py` exists to remove.
- **Transposable elements.** Ensembl types 5,896 TE transcripts in fly and 91 in yeast; neither
  RefSeq nor FlyBase's GTF carries them as transcripts.

**Decision unchanged: RefSeq for the six, FlyBase for the fly.** Per species: fly gains only TEs from
Ensembl and FlyBase r6.67 is newer than BDGP6.54; worm is near-identical apart from pseudogene
typing; zebrafish is clearly better under RefSeq; yeast's RefSeq annotation (SGD R64-5-1, 2026) is
fifteen years newer than Ensembl's R64-1-1. The pseudogene deficit is the one real cost, and it is
recorded here rather than discovered later.

This comparison also exercised the normalizer's unmapped-term gate, which refused two Ensembl
biotypes it had never seen (`pre_miRNA` x262, `miRNA_primary_transcript` x10, both miRNA precursors)
rather than passing them through. They are now mapped to `miRNA`, matching how GENCODE types
pre-miRNA and how RefSeq's equivalent `primary_transcript` is already handled. The selftest still
reproduces the human and mouse ground truth exactly.

#### Ensembl for the primates: content, not just assembly age

The primates were the case where switching registries would have cost the T2T assemblies, so it was
worth checking whether Ensembl's annotation was rich enough to justify that. It is not: Ensembl is
worse on essentially every axis, so there is no trade-off to weigh.

| species | source | tx | genes | tx/gene | coding | lncRNA | small nc | dropped nc | pseudogene |
|---|---|---|---|---|---|---|---|---|---|
| human | GENCODE v49 | 507,365 | 78,691 | 6.4 | 294,386 | 190,272 | 5,116 | 2,447 | 14,704 |
| gorilla | RefSeq T2T 2025 | 99,500 | 34,172 | 2.9 | 84,606 | **10,298** | 2,862 | 1,475 | 62 |
| gorilla | Ensembl gorGor4 2014 | 53,705 | 30,084 | 1.8 | 44,987 | **703** | 4,717 | 2,569 | 522 |
| chimp | RefSeq T2T 2026 | 135,577 | 34,784 | 3.9 | 119,088 | **11,353** | 2,876 | 2,065 | 111 |
| chimp | Ensembl Pan_tro_3.0 2016 | 61,457 | 33,729 | 1.8 | 49,717 | **3,091** | 5,120 | 2,812 | 485 |
| macaque | RefSeq T2T 2025 | 140,857 | 38,839 | 3.6 | 120,554 | **13,744** | 2,681 | 2,411 | 1,383 |
| macaque | Ensembl Mmul_10 2019 | 64,228 | 35,432 | 1.8 | 48,600 | **6,556** | 7,340 | 795 | 767 |

Ensembl carries **roughly half the transcripts** (45 to 54%) and far less long non-coding annotation,
most starkly for gorilla: **703 lncRNA transcripts against RefSeq's 10,298**, a factor of 15. Its
only advantage is pseudogene typing in gorilla and chimp (522 and 485 against 62 and 111), and even
that reverses for macaque (767 against 1,383). So for the primates the older assembly comes with the
poorer annotation, and RefSeq T2T wins outright.

**The isoform-depth confound is real and no registry fixes it.** Transcripts per gene: human 6.4,
but 2.9 to 3.9 under RefSeq and a flat **1.8 under Ensembl** for all three primates. Great ape
genomes are near-identical to human, so this is annotation effort, not biology, and it is worse
under Ensembl than under RefSeq. Since `--quantMode TranscriptomeSAM` expands one genomic alignment
across every compatible isoform, the same library yields systematically different per-transcript
record counts in human than in gorilla for annotation reasons alone.

That has a direct consequence for the cross-species evaluation: **per-transcript coverage is not
comparable between human and the non-human primates as-is.** It has to be handled analytically, by
aggregating to the gene level, or by restricting to 1:1 orthologs with comparable isoform counts,
rather than by choosing a different reference. Recorded here so the comparison is not made naively
later.

The unmapped-term gate fired once more on this data, refusing Ensembl macaque's `Y_RNA` (x666, the
Ro60 ribonucleoprotein component) rather than passing it through; it is now mapped to `misc_RNA`,
matching GENCODE, which places it in the postfilter BED rather than the drop set.

## Running the Mamba mixer on CPU (2026-08-25)

mamba_ssm's fast path is CUDA-only -- `causal_conv1d_fwd` asserts `x.is_cuda` -- which had been
recorded as "no CPU fallback". That is wrong. The package ships pure-PyTorch reference twins for
every CUDA kernel it uses, and swapping them in makes the whole model CPU-runnable:

    selective_scan_fn -> selective_scan_ref
    mamba_inner_fn    -> mamba_inner_ref
    causal_conv1d_fn  -> causal_conv1d_ref

**The non-obvious part.** Each name was resolved at import time in SEVERAL modules, so patching one
namespace is not enough. Patching only `causal_conv1d_interface` still dies, because
`mamba_inner_ref` calls the `causal_conv1d_fn` bound inside `selective_scan_interface`. All three
namespaces (`selective_scan_interface`, `causal_conv1d_interface`, `modules.mamba_simple`) must be
patched before the model is constructed.

Implemented as `_enable_mamba_cpu()` in `scripts/dump_pred_profiles.py`, gated on `RIBO_MAMBA_CPU=1`
(off by default -- on a GPU the fast path is much quicker). Verified to produce finite output on the
real `RiboSignalModel` (mamba4, 10 conv blocks, d_state 16) at L=1000 and L=3000.

**Cost, and a correction.** The microbenchmark is ~7.9 s per 3,000 nt transcript single-threaded.
That was then extrapolated to "a 41,096-tx dump at 16 cores lands in the same few-hour range as the
attn CPU dump" -- and that extrapolation was WRONG. The real run produced no dump at all in 12 h and
was killed by the partition wall, i.e. roughly an order of magnitude worse than predicted. Either
the dump does not parallelise the way the estimate assumed, or per-transcript cost scales worse at
real transcript lengths than at the 3,000 nt test point.

**So: CPU mamba is for smoke tests and small transcript sets.** It is NOT a general substitute for a
GPU. For a full dump, shard hard across many jobs or wait for a GPU (~10-20 min there). The attn
mixer, by contrast, genuinely does complete a 41,096-tx CPU dump in ~7 h.

**Why it matters operationally:** the GPU partition on this cluster sits at 48/48 allocated for long
stretches (8 A100 + 40 A5500). When that happens, a CPU run started immediately finishes before a
GPU job that is still queued. Used on 2026-08-25 to run both Wang mm10 cross-species arms (attn and
mamba4) concurrently on CPU rather than waiting. Training still wants a GPU -- the reference path is
far too slow for a backward pass.

### G. Drosophila is DROPPED from the modelling arm: no usable 3-nt periodicity (2026-08-26)

Fly aligned cleanly and then failed at P-site calling. RiboCode reported "No obviously periodicity
are detected" for all 9 BAMs and exited with "Error, can not determine the P-site locations". The
alignment was not the problem: 472.7M input reads, 98.6M uniquely mapped, 31-38 nt, landing on CDS.

**The cause is the nuclease.** GSE99920's GEO protocol reads: *"240 ul of lysate was incubated with
anti-Flag antibody coated magnetic beads and 10000 units of RNase T1 ... to perform digestion of
exposed mRNA and ribosome purification simultaneously."* RNase T1 cleaves specifically after G, so
footprint ends are ragged and sequence-dependent, and the sub-codon precision that periodicity
depends on is destroyed. RNase I, the standard, cuts without base preference and gives sharp ends.

This is the same class of error as the GSE157050 rejection recorded in section 5J (polysome
profiling mislabelled as ribosome profiling). The screening at dataset-selection time checked title,
matched-RNA design, replicate count and poly(A) status, but **not the nuclease**. Add it: for any
Ribo-seq candidate, read the digestion step out of the GEO extract protocol before downloading.

#### The measurement, and a methodological correction worth keeping

The first frame test pooled all read lengths under a single +12 P-site offset and reported
38.0 / 35.1 / 27.0 for GSE99920. That number is right but the method was wrong: the correct offset
varies with read length, so pooling 30-38 nt reads under one offset BLURS a real signal and can make
a periodic library look flat. **Always compute the frame histogram per read length, and validate the
method on an arm already known to be periodic.** Redone that way, with yeast as the positive control:

| dataset | dominant read length | n | f0 / f1 / f2 | max |
|---|---|---|---|---|
| yeast GSE173654 (RiboCode PASSED) | 28 nt | 581,729 | **93.9** / 2.6 / 3.6 | **93.9%** |
| fly GSE99920 (rejected, RNase T1) | 32 nt | 244,481 | 15.2 / 56.2 / 28.7 | 56.2% |
| fly GSE147619 (candidate, RNase I) | 32 nt | 131,623 | 21.5 / 19.5 / 59.0 | 59.0% |

Yeast concentrates 582k reads into ONE length class at 93.9% in frame. That is what sub-codon
resolution looks like. Both fly datasets sit near 56-59% and, more tellingly, are **incoherent
across read lengths** (GSE147619 gives 59.0 / 42.7 / 44.4 / 44.7 over its top four classes with no
consistent phase), which is what a ragged-ended library looks like.

#### The replacement candidate was tested and also rejected

A literature search for an untreated fly study with matched RNA-seq and RNase I footprinting
returned **GSE147619** (PRJNA615808, Drosophila embryo smORF translation) as the best candidate:
wild-type Oregon-R with no treatment or genotype arm, matched cytoplasmic RNA-seq from the same
lysate, 6 Ribo + 6 RNA over three developmental windows x 2 replicates, and RNase I confirmed in the
GEO protocol (*"Footprinting was performed overnight at 4C with RNaseI (Invitrogen)"*).

It was tested empirically on a 2M-read subset before committing to the 43.1 GB download, and
rejected on three counts: periodicity is weak and phase-incoherent (table above); only **12.07%**
of reads mapped uniquely with **87.31%** to too many loci; and the digestion follows polysome
purification with monosomes deliberately discarded (*"We purified mRNAs in polysomes, away from
monosomes (80S)"*), so the footprints are not the standard total-lysate monosome fraction.

**GSE49197** (Dunn and Weissman, stop-codon readthrough) is disqualified without testing: it used
**MNase**, and the paper states plainly that MNase's 3' A/T bias *"prevents us from achieving the
sort of sub-codon resolution seen in ribosome profiling datasets generated with RNase I."*

#### Decision

Fly is dropped from the modelling arm. Its references, alignments, BAMs, coverage and QC are RETAINED
and documented so the arm reads as processed-but-unusable rather than as a gap, and so the negative
result is not re-derived. The expansion proceeds with **eight species**. The four-species iPSC-CM
comparison, which is the primary generalisation test, is unaffected.

Caveat on the search: the literature workflow completed only 58 of 93 agents before hitting a
session limit, and the synthesis step failed. GSE60384, GSE168879/GSE168878 and GSE233555 were named
as candidates but never fully verified. If fly is wanted later, those three are where to resume, and
the nuclease is the first thing to check for each.

---

## Post-hoc rescoring of the HLA-I immunopeptidome arms (MS2Rescore, 2026-08-27)

### Why

The pgx tables report novel peptides from MSFragger alone: rank-1 PSM per spectrum, thresholded by
MSFragger's `hyperscore` at 1% class-specific FDR against `REV_nuORF|` decoys. Hyperscore is a
fragment-count statistic and is weakest exactly where these searches live, in nonspecific digestion
over a large novel-ORF space. Rescoring adds orthogonal evidence: predicted fragment-ion intensities
and predicted retention time.

The purpose here is **confirmation of an ordering, not a larger headline count**. The pgx claim is a
database-composition contrast (a model-derived database beats an every-AUG null at 25 to 150 times
fewer sequences). Rescoring changes the scoring function, applied identically to every arm. If the
model-over-null ordering survives a scoring function that shares nothing with hyperscore, the claim
is robust; if it flips on a line, that line was fragile.

### Configuration

`proteogenomics/msfragger/ms2rescore_pgx_hla.json`, MS2Rescore 3.2.1:

| component | setting | reason |
|---|---|---|
| feature generator | `basic` | charge, mass error, peptide length |
| feature generator | `ms2pip`, model `Immuno-HCD`, tol 0.02 | fragment-intensity model trained on non-tryptic HLA peptides |
| feature generator | `deeplc`, calibration 0.15 | retention-time prediction; absent from the earlier run in this project, and among the strongest discriminators for nonspecific searches |
| rescoring engine | `mokapot` 0.10 | its output is used only for the feature table and as a sanity check; the reported counts come from the class-aware pass below |
| decoy pattern | `^REV_` | matches the pgx database convention |

### Work units

`pgx/rescore_manifest.py` enumerates (line, arm, fraction). The `gencode` and `null_atg` databases do
not depend on which signal model produced the ORF calls and are byte-identical between the attn and
mamba4 pgx directories (verified by fasta md5), so each is rescored once rather than twice. That
gives 6 unique arms per line: `gencode`, `null_atg`, and `{attn,mamba4}_model_{standard,poisson}`.

126 units: 12 each for HBL-1, DoHH2 and SU-DHL-4 (2 fractions), 90 for THP-1 (15 fractions).
`pgx/rescore.sbatch` runs one unit per array task, indexed by manifest data-line number so a rerun of
an individual failure is exact. Output is gzipped (the feature table is 107 columns and compresses
about 5x; the group filesystem is at 93% of its 15 TB quota).

### Counting: three corrections, each forced by a measured failure

`pgx/rescore_report.py`.

1. **Class-aware.** Global mokapot optimises the canonical-dominated target/decoy separation and
   interleaves novel targets with novel decoys, which collapses the novel-class FDR. Measured
   earlier in this project on A549: 10/13 novel peptides (model/null) on raw hyperscore became
   **1/6** under global rescoring. The fix, already written in `compare_rescored_classaware.py`, is
   to retrain mokapot on the novel class only (`novel_t` vs `novel_d`) using the same MS2Rescore
   feature columns. The novel class is too sparse to bootstrap at mokapot's default 1% `train_fdr`,
   so the trainer walks a ladder of 0.05, 0.10, 0.25 and reports which one trained.

2. **Seed-ensembled.** At two-fraction depth a single mokapot run is bimodal, not noisy. Across 8
   seeds (`results/mokapot_stochastic.json`) SU-DHL-4 returned 2, 2, 18, 2, 18, 23, 23, 2 novel
   peptides (SD 10.1) and DoHH2's null returned 7 to 29 (SD 9.0). That spread is the entire effect
   size, so a single roll is not reportable. The fix is to average the **posterior error
   probability** over N seeds. Averaging the PEP rather than the raw discriminant is deliberate:
   mokapot's discriminant is not calibrated to a common scale across runs, while the PEP is a 0-1
   quantity with the same meaning in every run. A PSM that mokapot dropped in every seed is scored
   as PEP 1.0, never given a free pass.

3. **Deterministic final cut.** The 1% threshold is a target/decoy cut on the ensembled score using
   the same `pgx.report.cut` the hyperscore table uses, not mokapot's own q-values. This is what
   makes the two columns comparable: the score changes and nothing else does.

### What a rescored row counts, and how it differs from a hyperscore row

MS2Rescore defaults to `max_psm_rank_input=10` and `max_psm_rank_output=1`: it rescores MSFragger's
top-10 candidates per spectrum and writes back only the winner. The pgx hyperscore pipeline only
ever looked at MSFragger's rank 1. Rescoring can therefore promote a rank-2-to-10 novel peptide over
a rank-1 canonical one, so a changed count reflects **both** a changed threshold and a changed
spectrum-to-peptide assignment. Separating those two contributions would require a rerun at
`max_psm_rank_output` > 1, roughly 4x the output volume (the probe run discarded 30,658 non-winning
PSMs against about 10,000 kept), and is not done here. It is a diagnostic, not the claim being
tested: the hyperscore table already reports hyperscore selection under a hyperscore threshold, and
the rescored table reports rescored selection under a rescored threshold, which is the comparison
that decides whether the ordering is robust.

### Known limitation

These are shallow datasets: 1,415 to 5,527 baseline canonical PSMs per line, and 2 fractions for
three of the four lines. The standing project finding is that Percolator-class rescoring helps at 30
mzML and hurts at 10. THP-1, with 15 fractions, is the only line comfortably above that line, and it
is also the line the earlier rescoring work never covered. Where the class is too sparse to train at
any `train_fdr`, the report says `untrainable` rather than substituting a number.

### H. The RefSeq C. elegans GTF mislabels its CDS rows, and what it cost (2026-08-27)

NCBI's RefSeq annotation for C. elegans (GCF_000002985.6, WormBase WS298) writes the
isoform's own `standard_name` into the feature-type column where `CDS` belongs:

```
NC_003284.9  RefSeq  F31B12.1k  10838615  10838650  .  -  0  gene_id "CELE_F31B12.1"; ...
                     ^^^^^^^^^ should be CDS
```

204,530 of 204,542 CDS rows are affected, spread over 30,548 distinct names, leaving
**12** rows in the whole file that actually say `CDS`. None of the other six RefSeq
annotations used here (yeast, zebrafish, gorilla, chimp, macaque, human) contains a
single non-standard feature row, so this is one file's defect rather than a RefSeq
convention.

**Why it survived every existing check.** `exon`, `transcript`, `gene`, `start_codon`
and `stop_codon` rows are all untouched, so the file builds a correct STAR index, a
correct transcriptome FASTA, a correct ncRNA drop list and a correct universe. The
transcript counts all reconcile. Only a CDS-consuming step notices, and RiboCode does
not fail on it: with no annotated CDS to compare against, it types **every** called ORF
`novel`. That is exactly what it did, and the tell was in the output all along:

| arm | annotated | novel |
|---|---|---|
| yeast | 5,344 | 4 |
| zebrafish | 10,947 | 273 |
| human | 23,144 | 6,976 |
| **C. elegans** | **0** | **17,745** |

A worm arm reporting 17,745 novel ORFs and zero annotated ones is not a biological
result. It was caught by comparing the ORF-type composition across species, not by any
assertion in the pipeline.

**The repair, in `normalize_annotation.py`.** `repair_feature()` rewrites a row to `CDS`
only when all three of these hold, and all three held for 204,530/204,530 rows:

1. the feature type is not a known GTF feature,
2. the frame column is a valid `0`/`1`/`2` (only CDS carries a frame), and
3. the feature string equals the row's own `standard_name` attribute.

Anything else non-standard is a **hard error**, matching the treatment of an unmapped
biotype term: the normalizer is meant to be pointed at an eighth species safely, and a
silent pass-through is what caused this.

**A second silent failure found while fixing the first.** `--genome` was accepted
without an existence check, so an unmatched shell glob was passed through as the literal
pattern `*genomic.fna`, `detect_mito()` found no header to read, and the only symptom was
24 worm tRNAs keeping `tRNA` instead of `Mt_tRNA`. A missing `--genome` now aborts.

**Scope of the rebuild.** After re-normalizing, `ncrna_tx.txt`, `tx_to_gene.tsv` and
`tx2biotype.tsv` are byte-identical to the pre-fix versions and CDS rows go 12 to
204,536. The STAR index, salmon index, transcriptome FASTA, universe and ncRNA filtering
are therefore all unaffected and were not rebuilt. Only `data/ribocode_annot_celegans`
and the worm P-site / ORF-call outputs were regenerated. The superseded outputs are kept
alongside as `*.broken_cds.2026-08-27` rather than deleted, per the intermediate-data
policy.

The GENCODE selftest (exact set equality against `data/ncrna_tx_human_v49.txt` = 2,447
and `ncrna_tx_mouse_vM38.txt` = 2,579) still passes after both changes.

### I. Genome-coordinate Ribo-seq BAMs and bigwig tracks (2026-08-27)

`scripts/xspecies/ribo_genome_bigwig.sbatch`, driven by
`data/external/xspecies/samplesheets/ribo_ALL_genome.tsv` (all 67 Ribo runs across the
seven active species plus the dropped fly arm).

**Why this is a re-alignment and not a conversion.** Every Ribo BAM this project had was
`Aligned.toTranscriptome.out.bam`. `scripts/riboseq_align.sbatch` runs STAR with
`--outSAMtype BAM Unsorted --quantMode TranscriptomeSAM` and copies out only the
transcriptome BAM; the genomic one is written into node-local `$SCR` and removed by that
job's own `trap ... EXIT`. `bamCoverage` cannot be pointed at what was on disk, because
those BAMs' references are TRANSCRIPTS (201,583 `@SQ` lines for human, 6,477 for yeast),
which no genome browser can display against an assembly. The FASTQs were all still local,
so the fix was to re-align rather than to restore anything.

**The posture is deliberately identical** to the transcriptome run, or the tracks would not
describe the same reads as the ORF calls: `--outFilterMultimapNmax 1`,
`--alignEndsType EndToEnd`, and the same per-run cutadapt arguments read from the same
samplesheet column. The single intentional difference is the ncRNA drop. In transcript
space `filter_tx_heldout.py` drops reads by TRANSCRIPT ID; in genome space the equivalent
is dropping reads that fall in rRNA / tRNA / miRNA / Mt loci, using
`data/annot/<species>/ncrna_postfilter.bed` via `samtools view -L <bed> -U <keep.bam>`
(`-L` selects reads overlapping the intervals, `-U` writes the ones that do not, which is
the set to keep). The job refuses to run on an interval file with fewer than 10 lines,
since an empty one would filter nothing and report success.

**bigwig settings**: `--binSize 1`, no normalization, so the track is raw per-base read
depth and is directly comparable to the P-site profiles the model is scored against. The
consequence to remember is that tracks are NOT comparable ACROSS runs without accounting
for library size.

**Validation, on the yeast smoke test before the other 66 were launched.** The genome and
transcriptome runs of `SRR26680426` agree exactly: 25,378,479 input reads and 1.69%
uniquely mapped in both. That equality is the evidence the posture really is the same.

Two numbers that look wrong and are not:

- **1.69% unique mapping.** STAR reports 96.16% of reads going to "too many loci". That is
  `--outFilterMultimapNmax 1` discarding rRNA-repeat multimappers in an undepleted yeast
  library, it is identical in the pre-existing transcriptome BAMs, and that arm still
  yields 4,973 CDS calls at precision 1.000 with clean periodicity.
- **The genomic BAM keeps MORE records than the transcriptome one** (300,660 vs 248,777),
  inverting the usual isoform-expansion expectation. The two ncRNA filters are not
  equivalent: `filter_tx_heldout.py` can only retain reads that align to an annotated
  transcript at all, whereas the genomic BED filter also keeps intergenic and intronic
  reads. Worth knowing before comparing a browser track against a pack profile.

**Contig naming.** These are RefSeq assemblies, so contigs are `NC_000001.11` style rather
than UCSC `chr1`. Load `genomes/xspecies_refs/<species>/genome.fna` as the IGV reference;
a stock hg38 will not match.

**deeptools** was installed into `conda_envs/riboseq` (3.5.6). The env was snapshotted first
to `logs/env_snapshots/riboseq_before_deeptools_2026-08-27.txt` and the solve dry-run
checked: 38 pure additions, no downgrade of numpy, pysam, python, STAR, salmon, samtools
or cutadapt, all of which were re-verified after the install.

Outputs: `data/xspecies_ribo_bigwig_genome/<arm>/<RUN>.bw` and
`data/xspecies_ribo_bam_mm1_genome/<arm>/<RUN>.genome.bam` (+ `.bai`, STAR and cutadapt
logs). The genomic BAMs are retained rather than deleted in-job, per the
intermediate-data policy; they are the only copy and re-deriving them costs 67 STAR runs.

---

## HLA-I immunopeptidome rescoring and grading -- see the dedicated methods document

The full manuscript-grade methods for the 2026-08-27/28 HLA-I analysis (database construction,
MSFragger parameters, MS2Rescore configuration, mokapot training and FDR, the spectrum-grading
rubric, software versions and the script map) are in
**`docs/METHODS_hla_rescoring_and_grading.md`**, with the index of all results, figures and data at
`docs/HLA_IMMUNOPEPTIDOME_README.md`.

Kept there rather than inline because it is manuscript copy: every parameter is the value actually
used, and the document is meant to be lifted into a methods section largely as written.

### J. The mm25 multimap diagnostic arm (2026-08-28)

The plan's fifth matrix cell: Ribo at `--outFilterMultimapNmax 25`, to measure per species
what the standard mm1 posture discards.

**Alignment.** `scripts/riboseq_align.sbatch` hardcoded `--outFilterMultimapNmax 1`. It now
takes `MM=${MM:-1}`, so the default is unchanged for every existing caller and the mm25 arm
is `MM=25` with `OUTROOT=data/xspecies_ribo_bam_mm25`. All 67 runs, same adapters, same
`filter_tx_heldout.py` posture-A filter (isoform multimapping kept, cross-gene dropped).

**ORF calling.** `scripts/xspecies/ribocode_xspecies.sbatch` hardcoded its output to
`data/xspecies_psites/$DS`, which holds the mm1 P-sites, the packs, and every call set the
model was scored against. Running a second posture through it unchanged would have silently
overwritten them. It now takes `PSITEROOT`, defaulting to the old path; mm25 writes to
`data/xspecies_psites_mm25/`.

**Comparison.** `scripts/xspecies/compare_mm1_mm25_calls.py`, keyed `(gene_id, ORF_gstop)`.

**The confound that had to be measured, not assumed.** A raw call-count ratio is
uninterpretable on its own, because RiboCode silently drops any run whose 3-nt periodicity
fails its metaplots gate: that run contributes no uncommented row to `<ds>_pre_config.txt`
and simply is not in the pooled call. Chimp passed 8/8 runs at mm1 and 3/8 at mm25, so its
0.884x looks like "mm25 calls fewer ORFs" when the real effect is five lost libraries.
`runs_used()` counts contributing runs per posture and marks any arm where they differ
`NO_depth_differs`. Six of seven arms are like-for-like; chimp is not.

**Scheduling, and a mistake worth not repeating.** These jobs sat PENDING for hours with
`Reason=Nodes_required_for_job_are_DOWN,_DRAINED_or_reserved...` while the partition had
2,106 idle CPUs, 1.5-2 TB free memory per node, no drained nodes and no reservations. The
message is misleading: the real cause is fairshare (priority 1369 after ~1,100 jobs that
day) plus oversized walltime, since backfill only starts a job in a gap it can finish
inside. Walltimes had been set from worst-case guesses rather than from the mm1 runtimes
already on disk. Cutting them started a job within seconds:

| arm | mm1 actual | requested | cut to |
|---|---|---|---|
| zebrafish | 0:11 | 10:00 | 1:30 |
| gorilla | 0:40 | 11:45 | 3:00 |
| chimp | 1:00 | 11:45 | 4:00 |
| human | 1:51 | 11:45 | 6:00 |
| macaque | 1:53 | 11:45 | 6:00 |

Chimp began running on phoenix-10 immediately after its cut. Switching partition is NOT the
lever (short/medium/long draw on the same nodes, and this account has no `high_priority`
association); multi-partition submission is also unavailable because the association pins a
single partition. **Size walltime from a measured prior run, not from fear.**

**Cost notes.** Disk inflation is negligible after the posture-A filter (yeast final mm25
BAM 1.09x mm1), contrary to the 4.1x suggested by `janich_ribo_mm25`, which came from a
different, unfiltered pipeline. Time is the real cost: `filter_tx_heldout.py` is
single-threaded over every record and human hit 2.02 BILLION records, taking 87 minutes for
the filter alone. The first human attempt died at a 2h limit; the successful re-run took
2:04:48.

**A trap re-confirmed.** The timed-out human job left a 19 GB BAM that PASSED
`samtools quickcheck` but was the unfiltered intermediate, because `filter_tx_heldout.py`
overwrites in place and the job died during the sort that follows. quickcheck validates
structure, not correctness or completeness. Its mtime, not its validity, is what exposed it.
The file was deleted and the run redone.

**Fly.** RiboCode fails on the fly mm25 BAMs exactly as it does at mm1: no read length in
any of the 9 runs passes the periodicity test, so no P-site offsets can be assigned. This is
the documented RNase T1 problem from section G, not a pipeline fault, and no attempt was
made to lower `frame0_percent` to force calls out of a library with no periodicity.

### K. The RNA mm1 posture completes the 5-way matrix (2026-08-29)

The last unbuilt cell of the alignment matrix from the expansion plan: RNA at
`--outFilterMultimapNmax 1`, for comparability with the project's historical posture.
72 runs across 7 arms (fly excluded with the rest of that dropped arm), via
`scripts/xspecies/align_rna_xspecies.sbatch MM=1`, which was already parameterized and
writes to `data/xspecies_rna_bam_mm${MM}` with no code change needed.

All 72 written, all 72 pass `samtools quickcheck`, per-arm counts match the samplesheets
exactly (worm 10, yeast 16, zebrafish 23, chimp 8, macaque 7, human 5, gorilla 3).

**The matrix is now complete:**

| assay | posture | purpose | state |
|---|---|---|---|
| RNA | STAR mm1 | historical-posture comparability | archived |
| RNA | STAR mm10 | corrected input posture; the packs are built from this | archived |
| RNA | salmon decoy | TPM, universe construction | local (`data/xspecies_tpm/`) |
| Ribo | STAR mm1 | the only posture used as a model target | local |
| Ribo | STAR mm25 | multimap diagnostic (section J) | archived |

**Walltime, again.** The first worm RNA mm1 run finished in 6 minutes 14 seconds against a
6-hour request, and the arms sat `(Priority)` behind 425 cluster-pending jobs until the
limits were cut to 1 h (yeast, worm), 2 h (zebrafish) and 3 h (primates), at which point
every arm moved to `Reason=None`. This is the same error recorded in section J and it was
repeated here: **size walltime from a measured prior run.** The measurement existed both
times.

### L. Archive state after the expansion (2026-08-29)

Archived to `/private/warm-archive/carpenterlab/RNAZoo_meta/...` via
`scripts/archive_to_warm.sh` (copy, checksum-verify, index, stub, then delete the source
only after `rsync --checksum --dry-run` reports zero differing paths):

| path | size |
|---|---|
| `data/external/xspecies/fastq` | 308 GB |
| `data/xspecies_rna_bam_mm1` | 257 GB |
| `data/xspecies_rna_bam_mm10/*` (8 arms) | 252 GB |
| `data/xspecies_ribo_bam_mm25` | 60 GB |
| `data/xspecies_ribo_bam_mm1_genome` | 17 GB |

Kept local, because they are what everything downstream reads: the mm1 Ribo BAMs, all
P-site and coverage hd5, the 7 packs, `data/xspecies_tpm/`, the genome bigwigs (1.1 GB),
and every call set and results table.

**An ordering constraint worth stating.** The FASTQ archive had to run AFTER the RNA mm1
alignment, because that alignment reads the FASTQs. The first archive chain was launched
with the FASTQ as its first item and was stopped mid-copy by its recorded PID; because
`archive_to_warm.sh` deletes only after a verified copy, all 186 files were still intact
and no stub had been written. Check `squeue` and the pending work list before archiving an
input, not just before archiving an output.

---

# Poisson stringency sweep on the cross-species arms (2026-08-30)

**Why.** The cross-species ORF-call tables reported RECALL by class and an overall F1. Neither
answers the question a user of the predicted call set has: if the model hands me a uORF call, is
it right? That is precision, it is the number the stringency dial moves, and on the non-canonical
classes it is far from recall (human uORF 0.632 precision against 0.712 recall; zebrafish 0.084
against 0.645).

**The sweep.** `scripts/xspecies/sweep_poisson_xspecies.sbatch`, job 37342702, a 140-task array:
7 species x 2 architectures x theta in {0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1.0, 1.5, 2.5},
seed 0 throughout. Each task runs `ribocode_dropin.py --variant pred_preddepth --pred_poisson
--pred_scale THETA`, which draws the per-nucleotide density as Poisson(theta x predicted_total x
p) instead of rounding it. Outputs under `results/xspecies_poisson_sweep/<sp>_<arch>/theta_<t>/`.

**Scoring.** `scripts/xspecies/analyze_poisson_sweep.py` scores each sweep point against that
species' observed call set through `compare_dropin_calls.build_loader` and the CDS / uORF / ncORF
class definitions from `orf_classes.py`, so membership rules and genomic `(gene_id, ORF_gstop)`
keying are identical to the published tables and a sweep row at theta = 1 is directly comparable
to the deterministic row. Writes `results/xspecies_poisson_sweep_curve.tsv` and
`results/xspecies_poisson_operating_points.json`.

**The operating point took two constraints, and the first rule was wrong.** The rule implemented
first was "the most permissive theta whose CDS precision clears 0.90". The sweep itself showed
that is useless: **CDS precision is essentially flat in theta**, above 0.96 across the entire
swept range in six of seven species, so the rule is non-binding and degenerates to selecting the
top of the grid, which is the point with the *worst* non-canonical precision. What actually
degrades under tightening is CDS RECALL. The rule is now the **strictest theta** (best
non-canonical precision, which is what the dial exists to buy) keeping **CDS precision >= 0.90
AND CDS recall >= 0.90**. A species satisfying neither anywhere is reported as having no
operating point rather than being given one; zebrafish is that case.

**Sbatch hygiene.** The array checks the exit code of `ribocode_dropin.py` AND the existence and
non-emptiness of its output, because `set -uo pipefail` without `-e` previously let a driver
report COMPLETED for tasks whose real work had died.
