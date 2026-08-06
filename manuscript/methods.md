# Methods

Manuscript methods for the per-nucleotide Ribo-seq signal model. This document is the
clean, submittable synthesis; the working development log with per-task provenance and
internal decisions lives in `../methods.md` and `../results.md`. It is maintained
alongside the work: every new configuration, dataset, or evaluation is added here with a
dated entry in the Changelog (Section 13).

Writing conventions: ASCII punctuation only; specific descriptors rather than
first-person possessive; software versions and accessions stated where they affect
reproducibility.

---

## 1. Overview and task definition

The model predicts, for a mature mRNA sequence, the per-nucleotide distribution of
ribosome P-sites that a Ribo-seq experiment would measure. Formally, given a transcript
of length L, the model outputs a length-L probability vector over positions (the profile,
which encodes ORF location and 3-nucleotide reading-frame periodicity) and a scalar total
(the expected footprint depth). The profile is the transferable quantity: it is
per-transcript normalized and therefore invariant to sequencing depth, so a model trained
on one set of tissues can be applied to a new dataset or species given only that dataset's
own RNA-seq. The intended use is to predict translation, and in particular the location of
non-canonical open reading frames (upstream, downstream, and lncRNA ORFs), from sequence
where clean experimental Ribo-seq is unavailable.

Three signals are aligned position-for-position on the spliced-transcript nucleotide axis
(0-based), so no splicing bookkeeping is needed at train or inference time:

1. a per-position sequence representation (a foundation-model embedding, a one-hot
   encoding, or a sequence-derived ORF-candidate track; Section 5),
2. per-position RNA-seq read coverage from a matched RNA-seq library (Section 4),
3. the per-position Ribo-seq P-site profile, which is the prediction target (Section 3).

Transcript identity throughout is the versioned Ensembl transcript identifier (for
example `ENST00000832824.1`); the alignment `@SQ` names, the P-site HDF5 keys, the ORF
caller output, and the GENCODE FASTA headers are all this same identifier, so every join
is direct.

---

## 2. Datasets

### 2.1 Training data (Chothani human primary-tissue atlas)

Ribo-seq and matched RNA-seq for nine human primary cell types were taken from the
Chothani et al. smORF atlas (BioProject PRJNA756018 / SRA study SRP333128; GEO GSE182371
Ribo-seq and GSE182372 RNA-seq; SuperSeries GSE182377). The nine cell types are Brain,
embryonic stem (ES), Fat, Fibroblast, HA_EC (human aortic endothelial), HCAEC (human
coronary artery endothelial), Hepatocytes, HUVEC (umbilical vein endothelial), and VSMC
(vascular smooth muscle). Fibroblast is the deepest (32 Ribo-seq samples,
SRR15513165 to SRR15513196, plus 32 matched RNA-seq samples SRR15513233 to SRR15513264)
and was used for all single-tissue development; the remaining eight tissues enter through
the leave-one-tissue-out protocol (Section 9.4).

### 2.2 External held-out human data (Ruiz-Orera cardiomyocyte)

Cross-study generalization was tested on Ruiz-Orera et al. 2024 (Nature Cardiovascular
Research, DOI 10.1038/s44161-024-00544-7; ENA PRJEB65856), specifically the human
induced-pluripotent-stem-cell-derived cardiomyocyte (iPSC-CM, line BIHi242) arm: five
clone-matched pairs of monosome Ribo-seq (single-end, ERR12549926 to ERR12549930) and
mRNA-seq (paired-end, ERR12549891 to ERR12549895). Cardiomyocyte is absent from the
nine-tissue training panel, making this a genuinely out-of-distribution cell type from an
independent laboratory and protocol. The paper reports 82 to 84 percent of 29-nucleotide
footprints in the primary reading frame.

### 2.3 External held-out mouse data (Wang liver)

Cross-species generalization was tested on Wang et al. 2021 (Nucleic Acids Research; GEO
GSE94982 / SRA PRJNA375080), the adult (postnatal day 42) liver arm: two CHX-arrested
monosome footprint samples (SRR5262890, SRR5262891) and two matched mRNA-seq samples
(SRR5262874, SRR5262875). Footprint identity was confirmed by read structure (short
17 to 30 nucleotide inserts followed by the TruSeq adapter `AGATCGGAAGAGCACACGTCT` in
about 92 percent of reads); the paper reports about 79 percent in-frame periodicity. Adult
liver is the closest analog to the human Hepatocytes training tissue, isolating species as
the primary variable.

### 2.4 Genome annotations

Human data used GENCODE v49 on GRCh38; mouse data used GENCODE vM38 on GRCm39. All
transcriptome axes are the full primary-assembly annotation (509,650 human transcripts;
278,326 mouse transcripts), so protein-coding, lncRNA, and pseudogene loci are all present
and joins between the P-site target, the RNA-seq coverage, and the sequence representations
are exact.

---

## 3. Ribo-seq P-site profile (prediction target)

### 3.1 Read processing and alignment

Footprint reads were adapter-trimmed with cutadapt (human reads arrive 3'-clipped and were
passed `--trim-n --minimum-length 20`; mouse reads carry the TruSeq adapter and were passed
`-a AGATCGGAAGAGC --trim-n --minimum-length 20 --maximum-length 40`). Trimmed reads were
aligned with STAR 2.7.11b to the matched genome and GENCODE annotation, producing
transcriptome-coordinate alignments (`--quantMode TranscriptomeSAM`) under
`--outFilterMismatchNoverLmax 0.05 --outFilterMatchNminOverLread 0.7
--alignEndsType EndToEnd`. Multimapping was permitted at alignment (Section 3.3 explains
the retained posture).

### 3.2 Non-coding and cross-gene filtering

Before P-site assignment, transcriptome alignments to rRNA, tRNA, miRNA, snoRNA, and
mitochondrial rRNA/tRNA transcripts were dropped using a per-annotation drop list (7,587
human transcripts; 2,579 mouse transcripts), and alignments multimapping across distinct
genes were removed. Within-gene isoform multimappers (a genomically unique footprint that
maps to several isoforms of the same gene) were retained (Section 3.3).

### 3.3 Multimapper posture

The target uses posture A: each footprint is counted with weight one on every compatible
within-gene isoform (the default of the P-site caller). This is consistent with the matched
RNA-seq coverage input, which uses the same posture, and the resulting per-transcript
magnitude inflation is absorbed by the per-transcript normalization intrinsic to the
profile target and to every scale-invariant evaluation metric (Section 9). A
primary-alignment-only rebuild (posture B) is pre-registered as a sensitivity check; because
posture A affects magnitude rather than profile shape, and multi-isoform transcripts were
not predicted worse than single-isoform transcripts in a held-out split, posture B is
expected to confirm rather than change the conclusions.

### 3.4 P-site assignment and pooling

Per-nucleotide P-site profiles were computed with RiboCode 1.2.15 (`process_bam.py`). For
each sense-strand footprint of a read length present in that sample's metaplot, the P-site
position is `reference_start + offset(read_length)`, with per-sample, per-read-length
offsets estimated by RiboCode's metaplot step (typically 28/29/30-nucleotide reads to
offset 12). This yields, per sample, one per-nucleotide P-site array per transcript
(length equal to the transcript length, 0-based). Per-transcript arrays are pooled across
a tissue's samples by element-wise summation (positional accumulation with an axis-order
assertion, int64 accumulation, int32 storage), giving one pooled per-nucleotide target per
tissue on the shared transcript axis. Storage is one variable-length array per transcript
(HDF5 vlen), never a dense (N, L, D) block, because slicing the length axis of such a block
is pathologically slow on the parallel filesystem.

### 3.5 Mitochondrial exclusion

Mitochondrial protein-coding mRNAs (the 13 chrM genes) are excluded from every universe,
target, training set, and evaluation. They are mitoribosome-translated under a different
genetic code, are leaderless (no 5'UTR or uORFs), and lack the cytoplasmic 3-nucleotide
periodicity the model targets; they are predicted poorly (whole-transcript Pearson median
about 0.12 versus about 0.62 for cytoplasmic protein-coding mRNAs) and are not a valid
target. Exclusion is enforced at both the universe-construction and split-loading stages.

---

## 4. RNA-seq coverage (input feature)

### 4.1 Rationale and alignment

The second input is per-nucleotide RNA-seq read coverage on the same transcriptome axis. It
carries transcript-level RNA abundance and context, which the sequence representation does
not encode (a ridge probe recovers transcript length and GC content from a foundation-model
embedding at R^2 about 0.9 but mean expression at R^2 about 0); the two inputs are therefore
orthogonal and additive. Matched RNA-seq was aligned with STAR 2.7.11b to the same index used
for the footprints (`--quantMode TranscriptomeSAM --outFilterMultimapNmax 20 --outSAMtype
None`, so only the transcriptome BAM is produced, in node-local scratch, and deleted after
coverage extraction).

### 4.2 Library-strand auto-detection

Coverage counting is strand-aware and the library strand is detected per sample rather than
assumed, because external datasets differ in library orientation (the Chothani training
libraries are reverse-stranded ISR; the Ruiz-Orera samples are a mixture of forward- and
reverse-labeled; the Wang mouse samples are effectively unstranded). For each sample, up to
two million mapped records are inspected and the sense fraction is computed (a mate is sense
to its transcript when `is_read2 != is_reverse`); the library is called ISR when the sense
fraction is at least 0.7, ISF when it is at most 0.3, and unstranded otherwise. Counting
then keeps sense mates for ISR, antisense mates for ISF, and all mates for unstranded. This
per-sample detection is essential: a fixed-strand assumption drops the majority of reads on
a mislabeled library (the initial fixed-ISR implementation kept only about 1 percent of
reads on the forward-stranded external samples).

### 4.3 Per-nucleotide coverage and pooling

Each retained mate adds one to every transcript position it spans (per-mate depth, as in
`samtools depth`); supplementary and unmapped records are skipped. Multimapper posture
matches the target (posture A: every alignment record contributes). Per-sample coverage is
pooled across a dataset's samples by element-wise summation into one per-nucleotide coverage
array per transcript, on the same axis as the target.

### 4.4 Depth normalization

Raw pooled coverage carries the dataset's absolute sequencing depth, which is not
transferable across datasets. For transferable evaluation, coverage is divided by the
dataset's global mean per-nucleotide depth before the `log1p` transform (the `global_mean`
mode; a CPM-like normalization each dataset reproduces from its own coverage), so a
deeper library gives the same normalized input. Within a single dataset this costs a small
(about 0.01) median profile Pearson relative to raw depth, the expected price of
transferability. The model input at each position is the sequence representation
concatenated with this single `log1p(normalized coverage)` channel.

---

## 5. Sequence representations and auxiliary track (input configurations)

Four interchangeable per-position sequence representations were compared; all are frozen
(no gradient flows into a foundation model) and all align 1:1 to the target and coverage.

### 5.1 RiNALMo embedding (1280-dimensional)

RiNALMo giga-v1 per-token embeddings, one array of shape (L, 1280) per transcript with the
CLS/EOS special tokens stripped so row i is nucleotide i. Extracted on GPU in the RNAZoo
RiNALMo container; transcripts capped at 10,000 nucleotides (the model's context limit).

### 5.2 Orthrus embedding (512-dimensional)

Orthrus large 4-track per-token embeddings (a Mamba state-space mature-mRNA model), one
array of shape (L, 512) per transcript. The 4-track (A/C/G/T) variant was used rather than
the 6-track variant, because the 6-track CDS and splice channels would hand the model the
annotated ORF location and make the task annotation-informed rather than sequence-only. The
pre-pool per-token tensor is used (`forward(...)`, not the mean-pooled `representation(...)`),
verified to reproduce the pooled vector on averaging; Orthrus one-hot input has no special
tokens so row i is nucleotide i directly.

### 5.3 One-hot encoding (4-dimensional)

A raw one-hot encoding (A/C/G/T, U folded to T; unknown bases all-zero), read directly from
the universe FASTA with no learned model. This is the control that measures how much lift a
foundation model provides over the raw sequence the model already receives.

### 5.4 ORF-candidate track

An optional sequence-derived track marks candidate ORF starts and their reading frames:
every maximal start..in-frame-stop ORF on the transcript contributes occupancy and
frame-phase channels. Two versions were evaluated: v1 marks AUG-start candidates only; v2
additionally marks near-cognate (non-AUG) starts. The track marks candidate positions and
frames but not which candidates translate, so distinguishing translated from untranslated
ORFs must still come from the model. The track is used only in the ORF-track configurations
(Section 10).

---

## 6. Expressed universe and data splits

The transcript universe for training and within-study evaluation is the set of
protein-coding and lncRNA transcripts with matched-tissue mean expression at least 1 TPM
(salmon, averaged over the tissue's RNA-seq samples) and mature length at most 10,000
nucleotides: 36,668 human transcripts (34,248 protein-coding, 2,420 lncRNA) for Fibroblast.
The nine Chothani tissues share this fixed universe, because the sequence representations and
ORF track are sequence-based and tissue-independent; per-tissue biology enters only through
the RNA-seq coverage input and the P-site target.

Splits are gene-disjoint by chromosome. Transcripts are grouped by chromosome (which is
automatically gene-disjoint) into five count-balanced folds; for within-study evaluation,
fold 0 is the test set, fold 1 is validation, and folds 2 to 4 are training. This guarantees
that no gene appears in more than one split, so held-out performance is not inflated by
paralog or isoform leakage. Transcripts with fewer than 50 pooled P-sites are excluded from
metrics (too little signal to score a shape). Mitochondrial genes are excluded from every
split (Section 3.5).

---

## 7. Model architecture

The model (`RiboSignalModel`, PyTorch) is a base-resolution dual-head convolutional network
in the BPNet family (Avsec et al. 2021), chosen from the empirical target statistics: the
per-nucleotide P-site counts are sparse (74 percent zero) and extremely overdispersed
(variance-to-mean ratio about 2,350, ruling out a Poisson head) and strongly 3-nucleotide
periodic (P-site autocorrelation peaks at lags 3/6/9/12).

- Input: at each position, the sequence representation (Section 5) concatenated with the
  `log1p` normalized RNA-seq coverage channel (Section 4.4), projected by a 1x1 convolution
  to a working width of 256.
- Body: ten residual dilated 1-D convolution blocks (kernel size 3; dilations
  1, 2, 4, ..., 512; GroupNorm(8) and GELU; dropout 0.1), giving an approximately
  4-kilobase receptive field at O(L) cost. The foundation-model embedding already encodes
  long-range sequence context, so the body's task is to mix in coverage and sharpen to a
  frame-aware base-resolution profile, which is a local operation.
- Optional attention: a two-layer pre-norm Transformer encoder (`nn.TransformerEncoder`,
  8 heads, model width 256, feed-forward width 512, GELU, `batch_first`), inserted after the
  convolution body in the attention configurations, with padded positions masked.
- Optional state-space mixer (alternative to attention): a two-layer bidirectional Mamba body
  (`MambaBody` of `BiMambaBlock`, each block a forward Mamba plus a second Mamba over the
  reversed sequence, summed with a pre-norm residual; model width 256, `d_state` 16, `d_conv` 4,
  expand 2). It supplies global context at O(L) rather than the Transformer's O(L^2), and is made
  bidirectional deliberately: the seq2ribo polisher it adapts (Kaynar and Kingsford 2026) is
  causal because it refines a mechanistic sTASEP simulation that already carries downstream
  context, whereas this model has no simulator and must read both directions. Selected by
  `--mixer mamba`; the Transformer path is the default and its state-dict keys are unchanged.
- Profile head: a 1x1 convolution to per-nucleotide logits, masked at padded positions and
  converted to a per-transcript distribution by log-softmax over the transcript length.
- Count head: the masked-mean of the body features concatenated with `log1p(total coverage)`,
  passed through a small MLP to a scalar log-count.

The model is about 5 million parameters (4.7 million at the 1280-dimensional RiNALMo input;
5.8 million for the bidirectional-Mamba variant), small enough to run on CPU when GPUs are
saturated.

---

## 8. Training

Training data are packed into ragged 1-D int32 arrays (`target_counts`, `coverage`) with an
`offsets` index (row k is `arr[offsets[k]:offsets[k+1]]`), joined by transcript id and
asserted length-consistent, so the loader is numpy-only (it runs inside the CUDA container,
which lacks h5py) and reads contiguous per-transcript slices. Batches are formed by a token
budget (summed length per batch, length-bucketed) so variable-length transcripts pack with
minimal padding; padded positions are masked throughout (convolution inputs zeroed, profile
logits set to a large negative before softmax).

The loss is the sum of (i) a per-transcript multinomial negative log-likelihood on the
profile, `-sum_i (c_i / N) log softmax(logits)_i`, equal-weighted across transcripts and
therefore scale-free, and (ii) `count_weight` (default 0.1) times the mean squared error of
the predicted log-count against `log1p(total P-sites)`. Optimization is AdamW (learning rate
3e-4, weight decay 1e-2) with a cosine schedule and warmup, gradient clipping at 1.0, and
early stopping on validation profile Pearson (patience 6). Periodicity is never in the loss;
recovering the 3-nucleotide phase is therefore an emergent property that the evaluation tests
directly.

---

## 9. Evaluation

Three complementary metrics probe the profile from progressively sharper angles, plus two
transfer protocols (cross-tissue and cross-study/species).

### 9.1 Profile Pearson (shape)

Per-transcript Pearson correlation between the predicted position distribution and the
observed profile `c/N`, reported separately for protein-coding whole-transcript, lncRNA
whole-transcript, and the 5'UTR (uORF) and 3'UTR (dORF) windows. The windows are defined
from GENCODE CDS coordinates and used only to restrict the metric; the model never sees the
CDS boundary, so a correct 5'UTR prediction is genuine. This scores the shape of the
predicted profile.

### 9.2 Localization AUROC (does the model concentrate signal on translated ORFs?)

Against a RiboCode ORF-call truth set, every maximal AUG..in-frame-stop candidate ORF on a
held-out transcript (length at least 30 nucleotides) is labeled positive if its start matches
a RiboCode call and negative otherwise; internal in-frame alternative-start AUGs that merely
re-initiate an already-translated CDS are collapsed out (the only collapsed class). Each
candidate is scored from the predicted profile by its in-frame fraction
(`sum of predicted probability at in-frame positions / total over the ORF body`), and the
AUROC of separating translated from untranslated candidates is reported. The informative
stratum is non-canonical ORFs (uORF/dORF/novel), where ORF length is uninformative; a
length-controlled AUROC (Mann-Whitney within length-quantile bins) and an observed-profile
AUROC ceiling accompany it. Because the truth set was called AUG-only, this is a
translation-localization test, not a non-AUG-discovery test.

### 9.3 RiboCode drop-in concordance (predicted profile as an experiment replacement)

The sharpest test feeds the predicted per-nucleotide profile into the actual ORF caller in
place of the experimental Ribo-seq. RiboCode consumes per-transcript P-site density as a
plain dictionary of arrays, so the predicted profile (scaled to realistic depth and rounded
to integer counts) is passed directly to `detectORF.main` with the prebuilt annotation
database and default parameters (start AUG, minimum 5 codons, Benjamini-Hochberg
FDR-adjusted p < 0.05); only the density array is swapped. Each call set (real, both
predicted variants, and the external reference) is thus an independently FDR-controlled list:
RiboCode applies the BH adjustment within each run at generation, so every retained ORF has
adjusted p < 0.05. Predicted calls are matched to real-data calls by set intersection on the
genomic ORF locus (gene id and genomic stop) over these BH-passed survivors, additionally
filtered at 90 nucleotides and at a 0.5x-uniform in-ORF enrichment floor (which removes the
diffuse 3'UTR softmax leak), and precision, recall, and F1 are reported overall and by ORF
type, including novel-ORF recall. The comparison deliberately does not re-threshold on the
adjusted p-value: the BH burden scales with the number of transcripts in a run and the sets
are generated over different-size universes (the held-out test set for the predicted and real
variants, the full annotation for the external reference), so the adjusted value is not
comparable across sets; the retained raw per-ORF p is set-independent and in practice removes
nothing on top of RiboCode's own FDR gate. Because the two predicted variants and the real
variant are run over the identical held-out universe, their FDR burden is matched and cancels;
only the comparisons against the full-annotation reference carry a residual universe-size
asymmetry, which is why those are reported as harness validation rather than the headline. Two
density variants isolate the contributions: predicted shape at the real per-transcript depth,
and the fully standalone predicted shape at the count head's predicted depth (no experimental
data at all).

### 9.4 Cross-tissue transfer (leave-one-tissue-out)

Each tissue is held out in turn: the model trains on the other eight tissues concatenated
(token-budget sampler over the combined transcripts, each training tissue capped at 6,000
scorable transcripts so deep tissues do not dominate), validates on a gene-disjoint
chromosome fold across the training tissues, and is tested on the held-out tissue's own
P-sites and RNA-seq. The sequence representation and ORF track are reused unchanged (they are
tissue-independent); only the coverage and target differ per tissue, and coverage uses the
per-dataset `global_mean` normalization. Hepatocytes is the primary held-out tissue.

### 9.5 Cross-study and cross-species transfer

The external held-out datasets (Sections 2.2 and 2.3) are processed through the identical
footprint pipeline (Section 3) and RNA-seq pipeline (Section 4), yielding a pooled P-site
target, pooled RNA-seq coverage, and a RiboCode call set on the external data. A
Chothani-trained model is then applied without retraining: the sequence representation and
ORF track come from the external dataset's own transcripts, the coverage is normalized by the
external dataset's own global mean depth, and the three evaluations (Sections 9.1 to 9.3) are
run against the external call set. For the mouse cross-species test the one-hot configuration
is used directly, because it requires no species-specific foundation-model embeddings; the
human-trained model is applied to mouse sequence and mouse RNA-seq. This is the strictest
test of the profile's transferability. (Status: pipeline and inputs complete for both
external datasets; the held-out prediction and drop-in runs are in progress; see Changelog.)

---

## 10. Model configurations compared

Configurations are the cross product of an architecture setting and a sequence backend.
Architecture settings, in order of increasing capacity:

| Setting        | ORF track | Attention | Notes |
|----------------|-----------|-----------|-------|
| `baseline`     | none      | none      | dilated CNN only |
| `orf`          | v1 (AUG)  | none      | AUG candidate track |
| `orf_v2`       | v2 (+non-AUG) | none  | near-cognate candidate track |
| `attn`         | none      | 2-layer   | Transformer encoder, no track |
| `orf_attn`     | v1        | 2-layer   | |
| `orf_v2_attn`  | v2        | 2-layer   | full configuration |

Sequence backends: `rinalmo` (1280-d), `orthrus` (Orthrus 4-track, 512-d), `onehot` (4-d),
and `concat` (RiNALMo and Orthrus stacked, 1792-d). All backends train on the identical
split, seed, and hyperparameters; only the input representation differs, and the loader
intersects the available transcripts so every backend is scored on the same set. The two
headline configurations are `baseline` (minimal) and `orf_v2_attn` (full).

Beyond this cross product, three single-variable ablations off the full `orf_v2_attn_onehot`
configuration probe the design axes directly, each on the same Hepatocytes hold-out and the
same profile-Pearson evaluation (Section 9.1): (i) input modality -- both inputs versus
sequence-only (`--input_mode emb`) versus RNA-seq-coverage-only (`--input_mode cov`), which
isolates that sequence carries the profile shape and periodicity while RNA-seq carries
magnitude; (ii) body capacity -- two versus four Transformer layers, and width 256 versus 384;
and (iii) the sequence mixer -- the two Transformer layers replaced by the bidirectional-Mamba
body (Section 7). Results in results.md Task 17.

---

## 11. Result: foundation-model embeddings versus one-hot

Across all three evaluations, on the held-out Hepatocytes cross-tissue task, the one-hot
backend matches or exceeds both foundation-model backends; the foundation-model embeddings
provide no lift for per-nucleotide P-site prediction. This is the central methodological
finding: with per-nucleotide supervision across eight tissues, the convolutional body learns
the task-relevant sequence features (codon and reading-frame identity) directly from raw
sequence, and a frozen foundation-model embedding is a lossy compression optimized for the
foundation model's own objective rather than for this one. The one-hot advantage is largest
in the full `orf_v2_attn` configuration, consistent with the ORF track and attention
supplying exactly what the embedding otherwise contributed.

Held-out Hepatocytes, `orf_v2_attn` unless noted (medians; higher is better):

| Backend | Profile Pearson (pc) | Profile Pearson (lncRNA) | uORF | dORF | Localization non-canonical AUROC (% of ceiling) | Drop-in F1 | Novel-ORF recall |
|---------|------|------|------|------|------|------|------|
| one-hot | 0.640 | 0.447 | 0.559 | 0.245 | 0.832 (94.6%) | 0.923 | 0.862 |
| Orthrus | 0.621 | 0.428 | 0.530 | 0.239 | 0.822 (93.6%) | 0.916 | 0.862 |
| RiNALMo | 0.611 | 0.438 | 0.534 | 0.234 | 0.827 (94.1%) | 0.915 | 0.875 |

Baseline (no ORF track, no attention) profile Pearson (pc): one-hot 0.587, Orthrus 0.585,
RiNALMo 0.565. Baseline localization non-canonical AUROC: one-hot 0.815, Orthrus 0.809,
RiNALMo 0.774. The one-hot lead is present at baseline and widens with the full configuration.
The localization AUROC reaches 93 to 95 percent of the observed-profile ceiling (the AUROC of
the real Ribo-seq profile itself), so the remaining gap to a perfect score is mostly
irreducible label noise rather than model error. The drop-in test recovers about 92 percent of
RiboCode's own calls at about 92 percent precision, including about 86 percent of novel ORFs,
from the predicted profile alone.

---

## 12. Compute environment and software

Work ran on the prism SLURM cluster (heavy jobs via sbatch; the head node for light
inspection and serial downloads only), with outputs on the group filesystem. Training and
inference run inside a Singularity/Apptainer container (`rnazoo-rinalmo`, torch 2.2.0 with
CUDA 11.8 plus numpy) on a single GPU (A5500 or A100); the model is small enough to fall back
to CPU when GPUs are saturated. Key external tools: STAR 2.7.11b (alignment), cutadapt
(adapter trimming), RiboCode 1.2.15 (P-site assignment and ORF calling), salmon (expression
quantification for the universe), RiNALMo giga-v1 and Orthrus large 4-track (foundation-model
embeddings). External datasets were downloaded serially single-stream (`wget -c`) and
md5-verified against the ENA filereport checksums.

---

## 13. Changelog

Maintained going forward: one dated entry per new configuration, dataset, or evaluation.

- 2026-07-08 -- Fibroblast per-nucleotide P-site target built (posture A, chrM excluded);
  cleaning posture verified (non-coding loci pool to zero, housekeeping controls large).
- 2026-07-09 -- RNA-seq per-nucleotide coverage fleet and pooling; expressed universe
  (36,668 transcripts); gene-disjoint chromosome five-fold split; ORF-candidate tracks v1/v2.
- 2026-07-10 -- Dual-head model and trainer; fold-0 RiNALMo baseline (protein-coding profile
  Pearson about 0.60, periodicity reproduced); coverage depth-normalization (`global_mean`)
  added for transferability.
- 2026-07-11 -- Improvement matrix (ORF track, attention) across all five folds;
  uORF/dORF-windowed profile Pearson; localization AUROC eval; Orthrus 4-track backend and the
  RiNALMo-vs-Orthrus-vs-concat comparison.
- 2026-07-12 -- Leave-one-tissue-out protocol (held-out Hepatocytes); RiboCode drop-in
  concordance eval.
- 2026-07-13 -- One-hot backend added (measures foundation-model lift); external held-out
  data pulled and processed: Ruiz-Orera human iPSC-CM (5 Ribo + 5 RNA pairs) and Wang mouse
  liver (2 Ribo + 2 RNA); library-strand auto-detection added to the coverage counter (fixes
  mislabeled external libraries).
- 2026-07-14 -- One-hot leave-one-tissue-out complete across all three evaluations: one-hot
  matches or exceeds both foundation-model backends (Section 11). Manuscript methods document
  created. Cross-study (Ruiz-Orera) and cross-species (Wang) held-out prediction started.
- 2026-07-15 -- Cross-study (Ruiz-Orera human) held-out transfer landed: the profile transfers
  (one-hot protein-coding Pearson 0.425, ties or beats the foundation-model backends) and the
  RiboCode drop-in reaches F1 about 0.93 from the fully standalone predicted profile, with the
  predicted-depth path matching the observed-depth path (F1 gap about 0 across all three
  backends). This resolves the count-head transferability question: no target refactor is needed
  for single-dataset deployment (`design_count_magnitude_transferability.md` Section 9).
  Input-modality and body-capacity ablations added (Sections 7, 10; results.md Task 17):
  sequence-only recovers 91 percent of the profile Pearson and all of the periodicity, while
  RNA-seq-only collapses to 0.12 with periodicity destroyed, and deeper or wider bodies do not
  beat the baseline. Bidirectional-Mamba sequence mixer implemented (Section 7) and a full
  nine-tissue leave-one-tissue-out launched; both in progress.
