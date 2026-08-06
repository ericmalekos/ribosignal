# Related work: sequence-to-ribosome-profile models

This note positions the per-nucleotide Ribo-seq P-site model developed here (referred to
below as the profile model) against the two most directly comparable published methods,
Translatomer (He et al. 2024) and seq2ribo (Kaynar and Kingsford 2026). All three map RNA
sequence to a ribosome-density signal, but they sit at different points on the
data-availability spectrum and differ in coordinate system, resolution, and the way they
decompose shape from magnitude. The short version: Translatomer is the closest analog to
the profile model (both are multimodal, taking RNA-seq coverage plus sequence and
predicting a Ribo-seq profile), while seq2ribo is the sequence-only end of the spectrum,
intended for de novo design where no expression data exist.

## Three-tool comparison

| Dimension | Translatomer (He 2024) | seq2ribo (Kaynar 2026) | Profile model (this work) |
|---|---|---|---|
| Input | RNA-seq coverage + one-hot genomic DNA (6 channels) | mRNA sequence only | RNA-seq coverage + sequence (FM or one-hot) + ORF-candidate track |
| Coordinate frame | genomic, 65,536 bp window (introns included) | transcript CDS (codon) | mature spliced transcript |
| Resolution | about 128 bp per token (512 tokens over 65 kb) | codon | per-nucleotide |
| Method | 12-layer transformer + 1D-conv tokenizer | mechanistic sTASEP simulator + Mamba polisher | dilated CNN + 2-layer transformer OR bidirectional Mamba mixer |
| Sequence encoding | one-hot DNA (no FM tried) | codon identity + ViennaRNA 2D-structure features | FM embeddings (RiNALMo / Orthrus) vs one-hot |
| Output / loss | log(x+1) Ribo-seq coverage, MSE (shape and magnitude coupled) | A-site counts, Poisson NLL, plus a TE head | multinomial profile (shape) + count head (magnitude), decoupled |
| Structure awareness | none explicit | explicit (base-pairing, backbone-angle, positional bucket) | implicit via FM embeddings |
| Sequence contribution | +5.3% over RNA-seq-only on the envelope metric (0.731 to 0.784) | the only input | sequence drives SHAPE (seq-only 0.579 of 0.639 pc Pearson, keeps periodicity); RNA-seq drives MAGNITUDE (RNA-seq-only 0.122, periodicity destroyed); one-hot ties FM (encoding) |
| Transfer demonstrated | cross-tissue 0.72 to 0.80; mouse (human-trained) 0.53 to 0.65 | cross-cell-line polishers retain predictive power | cross-study human: drop-in F1 ~0.93 (one-hot best), localization 0.752 = 94% of ceiling, count head transfers (obs ~ pred depth); cross-species mouse (Wang liver): drop-in F1 0.929, count head transfers across species too (obs 0.929 / pred 0.919) |
| Primary purpose | variant-effect interpretation / in silico mutagenesis | de novo mRNA design | Ribo-seq denoising + non-canonical ORF discovery (RiboCode drop-in) |
| Headline metric | Pearson 0.784 (genomic, 65 kb window, log-coverage) | Shape r 0.186 / Tx-level r 0.920 / TE r 0.732 | profile Pearson + ORF-call F1 + localization AUROC |

## Translatomer (closest analog)

Translatomer is a multimodal transformer that predicts a cell-type-specific ribosome
profiling track from two inputs: a cell-type-matched RNA-seq coverage track and the
genomic DNA sequence. It defines a 65,536 bp window centered on each protein-coding gene
(GENCODE v43, hg38/mm10), one-hot encodes the genome sequence (5-dimensional, including N),
resizes the RNA-seq bigwig to the window, concatenates the two into a 6-channel input, runs
a 1D convolution that tokenizes the window into 512 tokens, applies a 12-layer
self-attention tower, and decodes back to a per-base (binSize 1) ribosome profiling
coverage track. It is trained on paired Ribo-seq and RNA-seq across roughly 33 tissues and
cell lines under an MSE loss on log(x+1) coverage, evaluated by Pearson correlation on a
held-out chromosome.

Reported results: the full model reaches Pearson 0.784 on the held-out chromosome; an
RNA-seq-only ablation (sequence masked to zero) reaches 0.731, so sequence contributes
+5.3% (per-dataset median 0.696 to 0.734), with a 20% MSE reduction (0.060 to 0.048). It
generalizes to unseen human tissues (prostate 0.72, epithelial 0.79, hTERT-RPE-1 0.79,
rhabdomyosarcoma 0.80) and, trained on human, to mouse (brain 0.53, heart 0.59, liver
0.65). Its downstream use is in silico mutagenesis for variant-effect interpretation
(Kozak -3/+4 perturbations, ClinVar and gnomAD disease variants, evolutionary constraint
via phyloP and minor-allele frequency).

Translatomer is the direct predecessor in method class: multimodal RNA-seq plus sequence to
a Ribo-seq profile. The profile model differs on three axes that matter (below).

## seq2ribo (sequence-only design tool)

seq2ribo predicts ribosome A-site locations from mRNA sequence alone, with no RNA-seq and
no genomic context. It is a two-stage hybrid: a structure-aware TASEP simulator (sTASEP)
models ribosome traffic on a 1D codon lattice with fitted wait times that combine codon
identity, base-pairing count, local backbone-angle change, and positional bucket (all
structural features derived from ViennaRNA 2D structure); a Mamba-based "polisher" then
refines the simulated profile. Task-specific heads read translation efficiency (TE) and
protein expression off the polished profile. It is fit per cell line (iPSC, HEK293, LCL,
RPE-1) on GWIPS-viz aggregate A-site counts.

Reported results: transcript-level Pearson up to 0.920, within-transcript Shape r up to
0.186 (with all baselines at or below zero on Shape r), elementwise MAE 30 to 38% below a
sequence-only Translatomer baseline, TE Pearson up to 0.732 (with UTR), protein expression
up to 0.903. sTASEP reduces codon-level MAE 90 to 96% over classical TASEP.

seq2ribo occupies the opposite end of the data-availability spectrum from the profile model:
it needs no experimental covariate, at the cost of a much lower positional ceiling (Shape r
0.186 is the best achievable positional correlation from sequence alone).

## What the profile model shares and where it differs

The profile model shares Translatomer's data niche (RNA-seq available, Ribo-seq wanted) and
its multimodal RNA-seq-plus-sequence design. It differs on three axes:

1. Transcriptomic, per-nucleotide vs genomic, 128-bp-token. Translatomer's 512-token
   bottleneck over a 65 kb genomic window gives an effective internal resolution near 128 bp
   per token, which cannot represent 3-nucleotide codon periodicity. That is adequate for a
   smooth density envelope and for gene-level variant effects, but not for frame-level ORF
   calling. The profile model keeps single-nucleotide resolution on the mature transcript
   (the per-nt transformer attention is the O(L^2) cost that dominates inference), preserving
   the periodicity that RiboCode reads to call ORFs and assign frame. Only the
   per-nt/transcriptomic design supports the ORF-discovery deliverable.

2. Decoupled shape and count vs coupled MSE. Translatomer regresses log(x+1) coverage under
   a single MSE, entangling where ribosomes sit with how many there are. The profile model
   splits a scale-free multinomial profile (shape) from a separate count head (magnitude).
   That split is what enables the "RNA-seq but no Ribo-seq" transfer story and the count-head
   transferability analysis: the shape can be read from any RNA-seq input while the count head
   predicts depth at a fixed training reference (Chothani-scale). seq2ribo independently
   arrives at the same shape/magnitude separation by fitting its simulator on
   abundance-normalized error and reading magnitude off a separate TE head, which corroborates
   the design.

3. Foundation-model vs one-hot sequence ENCODING (distinct from the modality question).
   Translatomer used only one-hot DNA; the profile model tested RiNALMo and Orthrus embeddings
   against one-hot and found one-hot ties the FM. That is an ENCODING axis: how the sequence
   channel is represented. It must not be conflated with the MODALITY axis: how much the sequence
   channel contributes at all. The two models look opposite on the modality axis, but only because
   their metrics measure different things. Translatomer's single envelope Pearson is
   abundance-weighted, so RNA-seq coverage dominates it and sequence adds only +5.3%. The profile
   model's input-modality ablation on the SHAPE metric (results.md 17A) shows the reverse:
   sequence-only recovers 91% of the profile-shape Pearson (0.579 of 0.639 protein-coding) and
   reproduces the 3-nt periodicity (indeed sharpens it), while RNA-seq-only collapses to 0.122 with
   periodicity destroyed. The resolution is
   that RNA-seq drives magnitude (the count head reads total coverage) and sequence drives shape
   (the multinomial profile head); Translatomer's coupled envelope conflates the two. So the
   FM-no-lift result is about encoding within the sequence channel; it is NOT evidence that
   sequence as a modality contributes little -- for the within-CDS shape it contributes almost
   everything.

## Two corroboration points worth carrying into the manuscript

- FM-no-lift (the ENCODING axis) is externally consistent. Translatomer reached its performance
  on one-hot DNA alone, never needing a foundation model, and the profile model finds one-hot ties
  RiNALMo/Orthrus. Both say the sequence channel does not need a learned embedding. This is
  separate from the MODALITY question (how much sequence matters), which the profile model's own
  ablation answers directly: RNA-seq drives magnitude, sequence drives shape (results.md 17A). Do
  not read Translatomer's abundance-weighted "+5.3% from sequence" as saying sequence is a minor
  input in general -- for the within-CDS shape it is the dominant input.

- Shape/magnitude decoupling is independently endorsed. Both seq2ribo (abundance-normalized
  fit + TE head) and the profile model (multinomial profile + count head) separate positional
  shape from total magnitude, and both do so specifically to make the magnitude readout a
  normalized quantity rather than an absolute, depth-dependent count. This supports moving the
  count head toward a translation-efficiency or per-million target for multi-dataset coherence.

## Caveat on cross-paper metric comparability

The three headline numbers are not directly comparable and should not be tabulated as a
ranking. Translatomer's 0.784 is Pearson on log-coverage over a 65 kb genomic window, inflated
by large, easy-to-predict zero stretches in introns and UTRs. seq2ribo's Shape r 0.186 is
codon-level within-transcript correlation with no RNA-seq input. The profile model's profile
Pearson is within-CDS per-nucleotide. Each metric ranks methods within its own paper and
coordinate system, not across them. Comparisons across the three should be made on design
(input modality, resolution, coordinate frame, output decomposition, purpose), not on the
raw correlation values.

## References

- He J, Xiong L, Shi S, Li C, Chen K, Fang Q, Nan J, Ding K, Li J, Mao Y, Boix CA, Hu X,
  Kellis M, Xiong X. Deep learning modeling of ribosome profiling reveals regulatory
  underpinnings of translatome and interprets disease variants. bioRxiv 2024.02.26.582217
  (2024); peer-reviewed version cited by downstream work as He et al. 2024. License:
  CC-BY-NC-ND 4.0 (preprint).
- Kaynar G, Kingsford C. seq2ribo: structure-aware integration of machine learning and
  simulation to predict ribosome location profiles from RNA sequences. Bioinformatics 42,
  btag296 (2026), ISMB/ECCB proceedings supplement. DOI 10.1093/bioinformatics/btag296.
  Code: github.com/Kingsford-Group/seq2ribo (GPU-only; CMU academic/non-commercial license).
- Supporting methods cited by the above: RiboNN (Zheng et al. 2025); Mamba structured state
  space (Gu and Dao 2024); ViennaRNA (Lorenz et al. 2011); RiboCode ORF detection.
