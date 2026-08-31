# Methods: cross-species evaluation of a per-nucleotide Ribo-seq signal model

Eric Malekos. Compiled 2026-08-29.

This document describes the data, references, processing and evaluation for a zero-shot
test of a per-nucleotide ribosome-profiling signal model on seven genomes. The model was
trained on human and mouse only; no parameter was updated for any species described here.

All primary artifacts, scripts and result tables referenced below are on the group
filesystem under
`RNAZoo/experiments/riboseq_signal_model/`. Bulk inputs are in the warm archive and are
indexed in `data/ARCHIVE_INDEX.tsv`.

---

## 1. Scope and design

Seven species were evaluated: human, chimpanzee, gorilla, rhesus macaque, zebrafish,
*Caenorhabditis elegans* and *Saccharomyces cerevisiae*. An eighth, *Drosophila
melanogaster*, was ingested and aligned but dropped before modelling (section 8).

Each species contributes one matched pair of arms: ribosome profiling (Ribo-seq) and
poly(A)-selected RNA-seq over the same biological material. 144 sequencing runs were
newly acquired for this study (67 Ribo, 77 RNA; 305 GiB of compressed FASTQ). The five
human RNA runs were already on disk from an earlier fetch and are therefore absent from
this study's download manifest, although they are used throughout; this is why the human
RNA bar in Figure 1a is zero.

### 1.1 Dataset selection

Selection was implemented as code (`scripts/xspecies/select_runs.py`) rather than as a
list of accessions, so that the criteria are auditable and re-runnable. Requirements:

1. Ribo-seq and RNA-seq from the same study or from explicitly sample-matched studies.
2. RNA-seq must be **poly(A)-selected**, not rRNA-depleted total RNA. This is not a
   stylistic preference: a total-RNA arm previously collapsed a liver transcript universe
   to 3,321 transcripts in this project, and the failure is silent.
3. Ribo-seq must be monosome footprinting, not disome or TRAP.
4. No drug or stress treatment beyond the study's own untreated condition.

Two metadata traps were encountered and are recorded because they generalise.
`library_strategy` is unreliable for identifying ribosome profiling: the zebrafish RPF
project declares `RNA-Seq`, and the discriminating field is
`library_selection = "size fractionation"`. Separately, `sample_title` is empty for every
run of PRJEB65856, and the full experimental design is recoverable only from
`sample_alias` and `sample_description`, which encode
`<species>_<tissue>_<individual>_{mR,Ri}`.

The *C. elegans* candidate originally proposed (PRJNA1183817) was rejected on measurement:
its RNA arm is total RNA-seq, and its RNA libraries are 0.11 to 0.30 GB against 0.5 to
4.6 GB for its Ribo arm. It was replaced by GSE52861 + GSE52905 (Hendriks *et al.*,
Grosshans laboratory), whose RNA arm is poly(A)-selected (TruSeq stranded mRNA) at 1.0 to
2.2 GB per library and whose Ribo arm is gel-excised 28 to 30 nt monosome footprints. A
caveat is recorded in the dataset registry: the two arms are replicate timecourses A and
B, matched by developmental stage rather than by lysate.

### 1.2 Acquisition and integrity

FASTQs were fetched serially from ENA over a single connection
(`scripts/xspecies/fetch_from_manifest.sh`) and verified against ENA's published
`fastq_md5`. Multi-connection downloading was avoided: it silently corrupts large files,
producing a file of exactly the correct size that nonetheless fails `gzip -t`.

During acquisition ENA's FTP layer was unavailable and its HTTPS layer returned ~800 byte
HTML directory listings in place of four zebrafish files.
`scripts/xspecies/verify_and_repair.sh` detects and re-fetches such responses, alternating
HTTPS and FTP per attempt. Four runs could not be retrieved from ENA at all and were taken
from NCBI SRA (`scripts/xspecies/fetch_from_sra.sh`).

An SRA-derived FASTQ cannot reproduce ENA's md5, because `fasterq-dump` regenerates the
file and gzip framing, read order and header formatting all differ. Rather than accept
"md5 differs, probably fine", these files were validated on two independent criteria:
read count equal to ENA's published `read_count` for that run, and, where a mate had
already been ENA-md5-verified, identity of read IDs in order over the first 200,000 reads.
Provenance is recorded per file in `data/external/xspecies/sra_provenance.tsv`, and the
repair script skips them rather than deleting them as corrupt.

### 1.3 Adapter determination

Adapters were determined **per run by measurement**, not assumed per dataset
(`scripts/xspecies/probe_adapters.sh`), by counting matches to a panel of candidate
adapters in the first 400,000 reads of the complete file. A prior failure in this project
motivated the "complete file" requirement: a byte-range slice once reported 0% adapter
where the real figure was 93%.

One override was required and is recorded with its evidence in
`data/external/xspecies/adapter_overrides.tsv`. The zebrafish Ribo arm has fixed 44 nt
reads with no TruSeq adapter; under the default `-M 40` cutadapt window this discards
essentially every read (75 of 200,000 retained). The true adapter is the Ingolia variant
`CTGTAGGCACC`, which retains 96.5 to 98.0%. The probe script now refuses to emit a verdict
for a fixed-length Ribo library longer than 40 nt when no panel member matches, returning
`NEEDS_ADAPTER` instead of a silently wrong default.

---

## 2. Annotation normalisation

### 2.1 The problem

The existing pipeline is written against GENCODE, which annotates biotype in the
`gene_type` / `transcript_type` attributes. None of the new species has a GENCODE
annotation. Four distinct conventions had to be reconciled:

| source | species here | biotype location | notes |
|---|---|---|---|
| GENCODE | human v49, mouse vM38 | `gene_type` / `transcript_type` | the target vocabulary |
| RefSeq | 6 species | `gene_biotype` / `transcript_biotype` | no `Mt_rRNA`/`Mt_tRNA`; mito genes are plain `rRNA`/`tRNA` distinguished only by contig |
| Ensembl | none used | `gene_biotype` / `transcript_biotype` | older GENCODE variant (`lincRNA`, `antisense`) |
| FlyBase | fly r6.67 | **none** | biotype is the feature-type column; `mRNA` appears where a GTF writes `transcript` |

The cost of not handling this is silent rather than loud. The project's existing ncRNA
interval builder matches `gene_type "(rRNA|Mt_rRNA|...)"`, which matches **zero** lines on
every new reference, and its post-filter BED had no non-empty guard. The standing rule
that rRNA, tRNA, miRNA and mitochondrial loci are dropped before any ORF caller sees the
BAM would simply not have happened, on seven species, with no error raised.

### 2.2 Ensembl was evaluated and rejected on measurement

Before building a normaliser, Ensembl release 116 was checked as a possibly
GENCODE-aligned alternative. It is not. `gene_type` occurs **0** times in
`Saccharomyces_cerevisiae.R64-1-1.116.gtf.gz` (41,879 `gene_biotype`) and **0** times in
`Danio_rerio.GRCz11.116.gtf.gz` (1,161,865 `gene_biotype`). Ensembl uses the same
attribute keys as RefSeq, so switching buys nothing on the axis that actually breaks the
tooling. For the three non-human primates Ensembl 116 is also far behind: gorGor4 (2014),
Pan_tro_3.0 (2016) and Mmul_10 (2019) against T2T assemblies from 2025 to 2026. RefSeq was
therefore retained, and the normaliser built.

### 2.3 Implementation

`scripts/normalize_annotation.py` takes any GTF and emits a GENCODE-conformant view plus
every derived artifact the pipeline consumes. It auto-detects the source convention,
identifies the mitochondrial contig from the genome FASTA and re-types mitochondrial rRNA
and tRNA (which RefSeq does not distinguish), maps the source vocabulary onto GENCODE v49
through an explicit checked-in table, and writes the normalised GTF, the ncRNA drop list,
`tx_to_gene.tsv`, `tx2biotype.tsv`, two genomic interval BEDs and a JSON report.

Two design decisions matter for reuse. **An unmapped biotype term is a hard error**, not a
pass-through, which is what makes the script safe to point at an eighth species. And a
`--genome` path that does not exist is a hard error: an unmatched shell glob was silently
accepted as a literal filename, `detect_mito()` found no header to read, and the only
symptom was 24 worm tRNAs retaining `tRNA` instead of `Mt_tRNA`.

### 2.4 Validation

The normaliser must reproduce the two hand-built GENCODE ground-truth files already on
disk, by exact set equality, before being pointed at any species where no ground truth
exists (`--selftest`):

| input | target | rows | result |
|---|---|---|---|
| `gencode.v49.annotation.gtf` | `data/ncrna_tx_human_v49.txt` | 2,447 | exact match |
| `gencode.vM38.annotation.gtf` | `data/ncrna_tx_mouse_vM38.txt` | 2,579 | exact match |
| `gencode.v49.annotation.gtf` | `data/tx_to_gene_human_v49.tsv` | 507,365 | exact match |
| `gencode.vM38.annotation.gtf` | `data/tx_to_gene_mouse_vM38.tsv` | 278,326 | exact match |

The selftest passes after every change described in this document.

### 2.5 A defect in the RefSeq C. elegans annotation

NCBI's RefSeq GTF for *C. elegans* (GCF_000002985.6, WormBase WS298) writes the isoform's
own `standard_name` into the feature-type column where `CDS` belongs:

```
NC_003284.9  RefSeq  F31B12.1k  10838615  10838650  .  -  0  gene_id "CELE_F31B12.1"; ...
                     ^^^^^^^^^ should be CDS
```

204,530 of 204,542 CDS rows are affected, across 30,548 distinct names, leaving **12** rows
in the file that say `CDS`. None of the other six RefSeq annotations used here contains a
single non-standard feature row, so this is one file's defect rather than a RefSeq
convention.

It is invisible to every check that does not consume CDS. Exon, transcript, gene,
start_codon and stop_codon rows are untouched, so the file yields a correct STAR index, a
correct transcriptome FASTA, a correct ncRNA drop list and a correct transcript universe,
and all transcript counts reconcile. RiboCode does not fail on it either: with no
annotated CDS to compare against, it types **every** called ORF `novel`. The worm arm
reported 17,745 novel ORFs and zero annotated ones, against 5,344 annotated of 5,358 for
yeast. The defect was found by comparing ORF-type composition across species, not by any
assertion in the pipeline.

`repair_feature()` rewrites a row to `CDS` only when all three of the following hold, and
all three held for 204,530 of 204,530 rows: the feature type is not a known GTF feature;
the frame column is a valid 0, 1 or 2 (only CDS carries a frame); and the feature string
equals the row's own `standard_name` attribute. Any other non-standard feature is a hard
error.

The repair is narrowly scoped, and this was verified rather than assumed. After
re-normalising, `ncrna_tx.txt`, `tx_to_gene.tsv` and `tx2biotype.tsv` are byte-identical to
the pre-repair versions and CDS rows go from 12 to 204,536. All 11 P-site HDF5 files are
byte-identical before and after, and the ORF call set is the same 17,745 genomic keys with
zero gained and zero lost: detection was always correct, only typing was broken. The STAR
index, salmon index, transcript universe, pack and per-nucleotide profile evaluation are
therefore unaffected and were not rebuilt.

---

## 3. Reference construction

One reference set was built per species (`scripts/xspecies/build_species_refs.sbatch`),
parameterised from a single table (`scripts/xspecies/species.tsv`) read by both the build
script and the aligner so that annotation versions cannot drift between them.

| species | assembly | accession | annotation source |
|---|---|---|---|
| human | GRCh38.p14 | GCF_000001405.40 | RefSeq GCF_000001405.40-RS_2025_08 |
| chimp | NHGRI_mPanTro3-v2.1 | GCF_028858775.2 | RefSeq |
| gorilla | NHGRI_mGorGor1-v2.1 | GCF_029281585.2 | RefSeq |
| macaque | T2T-MMU8v2.0 | GCF_049350105.2 | RefSeq |
| zebrafish | GRCz11 | GCF_000002035.6 | RefSeq |
| *C. elegans* | WBcel235 | GCF_000002985.6 | WormBase WS298 |
| yeast | R64 | GCF_000146045.2 | SGD R64-5-1 |
| fly (dropped) | dmel r6.67 | n/a | FlyBase r6.67 |

Build steps: normalise the annotation; extract the transcriptome with `gffread -w`; build a
STAR index; build a decoy-aware salmon index over an explicitly retained gentrome; and
build a RiboCode annotation with `prepare_transcripts`. Two STAR parameters are set per
species rather than copied: `--sjdbOverhang` is the maximum RNA read length minus one, and
`--genomeSAindexNbases` is `min(14, log2(genome length)/2 - 1)`. The mammalian default of
14 is wrong for the small genomes (approximately 10 for yeast at 12 Mb), and STAR does not
error on a bad value, it simply builds a poor index.

Each build asserts that every output is non-empty and above a per-species floor recorded in
`docs/SPECIES_REFERENCE_REGISTRY.md`, that every source biotype term was mapped, and that
`prepare_transcripts` retained at least 50% of the GTF's transcripts. That last check
exists because FlyBase writes `mRNA` where a GTF writes `transcript`, which the normaliser
renames; an earlier version emitted the feature twice instead, producing 71,107 sequences
for 35,736 transcripts.

---

## 4. Alignment and quantification

A complete five-posture matrix was built for every species.

| assay | tool | posture | purpose |
|---|---|---|---|
| RNA | STAR | `--outFilterMultimapNmax 1` | historical-posture comparability |
| RNA | STAR | `--outFilterMultimapNmax 10` | corrected input posture; packs are built from this |
| RNA | salmon | decoy-aware, `-l A` | TPM, universe construction |
| Ribo | STAR | `--outFilterMultimapNmax 1` | the only posture used as a model target |
| Ribo | STAR | `--outFilterMultimapNmax 25` | multimap diagnostic (section 7) |

Ribo-seq alignment follows the project-wide rule: unique alignments only, because the
common ORF callers do not drop multi-mappers themselves and retained multi-mappers inflate
P-site counts and corrupt the periodicity test. Alignment is end-to-end
(`--alignEndsType EndToEnd`), since footprint ends carry the signal, and
`--quantMode TranscriptomeSAM` is used because downstream P-site calling is in transcript
coordinates. rRNA, tRNA, miRNA and mitochondrial loci are removed before any caller sees
the BAM: in transcript space by transcript ID (`filter_tx_heldout.py`), in genome space by
interval (`samtools view -L`).

RNA-seq alignment omits `--maximum-length` and `--discard-untrimmed`; the latter is correct
for footprints and wrong for RNA, where it once retained 16.4% of an arm against 99.0%
without it.

Measured outcomes per arm (mean over runs):

| arm | assay | n | reads kept after trim | uniquely mapped | mean read length |
|---|---|---|---|---|---|
| `ruizorera_hsCM_ribo` | Ribo | 5 | 98.1% | 22.2% | 27.8 |
| `primate_pt_ribo` | Ribo | 8 | 95.6% | 10.6% | 26.4 |
| `primate_gg_ribo` | Ribo | 3 | 94.8% | 18.5% | 27.0 |
| `primate_rm_ribo` | Ribo | 7 | 96.8% | 14.1% | 26.7 |
| `zf_gse46512_ribo` | Ribo | 8 | 97.0% | 8.0% | 31.1 |
| `worm_gse52905_ribo` | Ribo | 11 | 94.0% | 79.4% | 28.5 |
| `yeast_gse173654_ribo` | Ribo | 16 | 79.4% | 3.3% | 25.9 |
| `ruizorera_hsCM_rna` | RNA | 5 | 99.8% | 83.9% | 200.0 |
| `worm_gse52861_rna` | RNA | 10 | 100.0% | 93.4% | 49.0 |
| `yeast_gse173654_rna` | RNA | 16 | 100.0% | 79.1% | 131.3 |
| `zf_gse70549_rna` | RNA | 23 | 100.0% | 71.5% | 151.0 |

The low Ribo unique-mapping rates are expected and are not a processing fault. In yeast,
96.2% of reads are assigned to "too many loci" under `mm1`; these are rRNA-repeat
multimappers in an undepleted library. That arm nonetheless yields clean periodicity and
4,973 CDS calls at precision 1.000.

Every BAM is validated with `samtools quickcheck`, which verifies the EOF block, rather
than by inspecting the header: a BAM truncated mid-write retains a perfectly valid header,
and this cost 92 minutes of downstream compute in an earlier phase of the project.

---

## 5. Universe, P-sites and packs

The transcript universe per species is protein-coding plus lncRNA transcripts of mature
length at most 10,000 nt that are present on **both** the P-site and RNA-coverage
transcriptome axes (`scripts/xspecies/build_universe_xspecies.py`). The axis intersection
is taken across every run, not just the first: `tx_order` in the pack is this list and the
ORF track, one-hot sequence and model input all index against it positionally, so a
transcript present on one axis but not the other would shift everything downstream of it
with no shape check to catch it.

No expression prefilter is applied, deliberately. The scoring floor at evaluation time
defines the scored set, and an untranslated transcript in the universe is harmless because
it never becomes a test transcript; filtering on TPM here would bake a second, different
expression criterion into the axis itself.

| species | universe | protein-coding | lncRNA |
|---|---|---|---|
| human | 183,028 | 150,635 | 32,393 |
| macaque | 131,981 | 118,457 | 13,524 |
| chimp | 125,558 | 114,692 | 10,866 |
| gorilla | 91,112 | 81,014 | 10,098 |
| zebrafish | 65,673 | 56,887 | 8,786 |
| *C. elegans* | 31,332 | 31,026 | 306 |
| yeast | 6,044 | 6,036 | 8 |

P-sites are called by pooled RiboCode `metaplots` (per-read-length offsets, aggregated over
an arm's BAMs) followed by `RiboCode -l no -g`. The job aborts if no read length passes the
periodicity check, since for a Ribo library that means the trimming or the annotation is
wrong and the failure should surface there rather than downstream.

Packs use `universe_tx` mode with `ref_pack: None`, the one-hot sequence backend. This
requires no foundation-model per-token embeddings, which do not exist for these species and
which would otherwise have been the binding constraint on the whole study. The ORF-candidate
track is computed from the mature mRNA sequence alone, with no GTF or CDS input, so it is
defined identically for annotated CDS, uORFs, dORFs and lncRNA ORFs. This is also why the
*C. elegans* CDS defect of section 2.5 does not propagate into the pack.

---

## 6. Evaluation

### 6.1 Per-nucleotide profile agreement

For each test transcript of at least 100 nt with at least 50 pooled P-sites, the predicted
profile is compared with the observed one (`scripts/xspecies/eval_xspecies_profiles.py`),
and the per-transcript statistics are summarised as a median because the distribution is
long-tailed and short or low-count transcripts give unstable correlations.

`pred_flat` is a per-transcript probability distribution summing to 1; the magnitude lives
in a separate count head, `pred_total`. Correlating `pred_flat` against counts directly
compares a probability with a count and is dominated by the scale difference. The profile
is therefore scaled to the observed depth before the log transform, and the transcript
total is taken from `pred_total`. An earlier version of this script omitted that step and
reported a total correlation of approximately zero for every species, because summing a
normalised distribution yields the constant 1.

**`profile_r` is not the project's `pearson_median`.** The pipeline's own training and
evaluation code computes `pearson(p, c/N)` on the raw profile, and that is the
`pearson_median` tabulated elsewhere in this project. The statistic reported here applies
log1p first, on the grounds that a raw P-site profile spans orders of magnitude and a few
tall peaks otherwise dominate. The two are not interchangeable: per transcript they rank
at Spearman 0.65 to 0.79 and share only 34 to 42% of their top decile. Figures in this
report should not be compared against an existing `pearson_median` value. Which of the two
is the better selection criterion is an open question and is not settled here.

A null is included: `shuffled_r` recomputes the profile correlation after permuting each
observed profile's positions. This distinguishes positional information from magnitude
information, since a model that has learned only a transcript's overall level can still
score a respectable profile correlation.

### 6.2 A tie-handling artifact in the project's Spearman implementation

The project's shared `spearman` (in `scripts/train.py`) ranks via a helper that uses
`argsort(kind="stable")` and assigns **ordinal** ranks 0 to n-1. It does not average tied
ranks. On dense data the two agree; on a sparse P-site profile they do not.

These profiles are 74 to 95% zeros. Under ordinal ranking every tied zero receives a
distinct rank in **positional** order, so the observed rank vector becomes a proxy for
position along the transcript rather than for signal. The model correctly predicts low
density in the 3'UTR; in a species with long 3'UTRs those positions are simultaneously
observed-zero (hence high ordinal rank, being last) and predicted-low, which manufactures
an anti-correlation from nothing. Measured over the full evaluation set:

| species | zero fraction | ordinal rank (as implemented) | average-rank ties (standard) | Pearson on log1p |
|---|---|---|---|---|
| *C. elegans* | 73.6% | 0.152 | 0.312 | 0.327 |
| yeast | 78.1% | 0.156 | 0.285 | 0.269 |
| macaque | 85.5% | **-0.047** | 0.420 | 0.537 |
| chimpanzee | 87.6% | **-0.070** | 0.389 | 0.521 |
| human | 88.8% | **-0.195** | 0.379 | 0.532 |
| gorilla | 89.0% | **-0.140** | 0.352 | 0.487 |
| zebrafish | 95.3% | **-0.405** | 0.142 | 0.137 |

The ordinal statistic falls as the zero fraction rises (the two most zero-sparse species are
the only positive ones, and the most sparse is the most negative), though the relationship is
not perfectly monotonic: human at 88.8% zeros scores -0.195 against gorilla's -0.140 at 89.0%.
With correct tie handling every species is positive and the rank correlation tracks the Pearson
value. The negative `profile_rho`
figures reported by this script before 2026-08-29 were therefore a metric artifact and not
evidence of anti-correlated profiles.

The evaluation script now emits both statistics, labelled. **`train._rank` was corrected on
2026-08-30** to average tied ranks; it now matches `scipy.stats.rankdata` exactly (max absolute
difference 2.2e-16 over 3,000 real transcript profiles), and the corrected median rank correlation
on those profiles is +0.3519 against -0.2345 as previously computed, with the sign disagreeing on
74.9% of transcripts. The function feeds `spearman` only -- it touches neither the loss nor
checkpoint selection, so **no trained model is affected and nothing required retraining**.

Values already written to disk were not recomputed: per-transcript Spearman is stored nowhere, so
correcting a value requires re-evaluating that run's model, and no document, table, or figure reads
the column -- this section is its only citation. Instead the stale key was renamed to
`spearman_median_INVALID_ordinal_ties` in all 255 affected JSONs (6,364 keys), each carrying a note
pointing here, so the wrong number cannot be read as a valid rank correlation. Backup at
`results/spearman_flag_backup_2026-08-30.tgz`; the rename is reproducible via
`scripts/flag_stale_spearman.py`. Any run evaluated after that date gets a correct value.

### 6.3 ORF calling from predicted profiles

To compare ORF calls rather than profiles, RiboCode is run on the model's predicted profile
under the same statistics as on the observed one (`scripts/ribocode_dropin.py`, driven by
`scripts/xspecies/dropin_xspecies.sbatch`). Three variants per arm:

- `real`: the observed profile, the reference call set;
- `pred_obsdepth`: predicted shape scaled to observed depth, isolating whether the model
  places ribosomes correctly;
- `pred_preddepth`: predicted shape and the model's own count head, a fully de novo call
  using no observed data. All results reported here use this variant.

Call sets are compared through the project's shared loader
(`compare_dropin_calls.build_loader`): restriction to the model's test transcripts,
`pval_combined <= 0.05`, ORF length at least 90 nt, and, on the predicted side only, a mean
predicted density over the ORF of at least 0.5x uniform. The enrichment gate removes the
diffuse low-magnitude leak into 3'UTRs that otherwise produces spurious dORFs. Using the
shared loader is not optional: hand-loading the same files previously produced a fourfold
overcount of model-specific calls and inverted a conclusion.

Calls are keyed on `(gene_id, ORF_gstop)`, the genomic ORF locus, **not** `ORF_ID`. One
genomic ORF appears once per transcript it sits on, and RiboCode's collapse step chooses a
representative isoform, so restricting the transcript universe reshuffles which isoform
carries a call; genomic keying neutralises that. Each species uses its own
`tx2biotype.tsv`; the loader raises if no test transcript matches the map, which is the
check that catches a wrong-assembly table.

### 6.4 ORF classes

`ORF_type` alone does not give the biologically interesting classes, because RiboCode's
`novel` label covers both ORFs on non-coding transcripts and novel ORFs on coding
transcripts, which are different claims. Classes are therefore defined on the transcript
biotype as well (`scripts/xspecies/orf_classes.py`):

- **CDS**: `ORF_type == annotated` on a protein-coding transcript;
- **uORF**: `uORF` or `Overlap_uORF` on a protein-coding transcript;
- **ncORF**: any called ORF on a non-coding transcript (lncRNA, antisense_RNA,
  ncRNA_pseudogene);
- **other**: dORF, internal, and novel-on-coding, reported for completeness so the classes
  reconcile against the totals.

Per-class call lists with a `recovered_by_model` column are written to
`results/xspecies_orf_classes/`.

---

## 7. The multimap diagnostic

All 67 Ribo runs were re-aligned at `--outFilterMultimapNmax 25` and ORF-called
independently, to measure what the standard unique-only posture discards. Comparison is by
`scripts/xspecies/compare_mm1_mm25_calls.py`.

Two methodological points govern the interpretation. First, NH tags in a transcriptome BAM
measure **isoform multiplicity**, not genomic multimapping: `--quantMode TranscriptomeSAM`
projects one genomic alignment onto every compatible isoform, so a genomically unique read
routinely carries `NH:i:17`. Posture comparison must therefore use total records per
transcript, and an NH split previously reported "99.8% of reads discarded" against STAR's
own 30%.

Second, and decisively for the primates, RiboCode silently drops any run whose periodicity
fails its metaplots gate: that run contributes no uncommented row to `<ds>_pre_config.txt`
and is absent from the pooled call. Chimpanzee passed 8 of 8 runs under mm1 and only 3 of 8
under mm25, so its two call sets rest on different amounts of data.
`runs_used()` counts contributing runs per posture and flags any arm where they differ, so
a depth artifact cannot be read as an ORF-calling result.

---

## 8. Drosophila was dropped, on measurement

Two independent *D. melanogaster* datasets were tested and both fail the requirement for
sub-codon resolution. The frame distribution of P-sites at each library's dominant read
length:

| dataset | dominant length | n | frame 0 / 1 / 2 |
|---|---|---|---|
| yeast GSE173654 (RNase I, positive control) | 28 nt | 581,729 | **93.9** / 2.6 / 3.6 |
| fly GSE99920 (RNase T1) | 32 nt | 244,481 | 15.2 / 56.2 / 28.7 |
| fly GSE147619 (RNase I) | 32 nt | 131,623 | 21.5 / 19.5 / 59.0 |

Yeast concentrates 582k reads into one length class at 93.9% in frame, which is what
sub-codon resolution looks like. Both fly libraries sit near 56 to 59% on their maximal
frame, and GSE147619 is incoherent across read lengths (59.0 / 42.7 / 44.4 / 44.7 over its
top four classes, with no consistent phase), which is the signature of ragged-ended
footprints. RiboCode assigns no P-site offset for any read length in either dataset and
exits. No attempt was made to lower the periodicity threshold to force calls from a library
with no periodicity. The fly alignments are retained and archived; only the modelling arm
is dropped.

---

## 9. Reproducibility and data availability

Binaries are invoked by absolute path throughout; no step depends on an activated
environment. Scripts, samplesheets and registries are version-controlled. Two registries
record provenance at the level needed to re-derive any number:
`docs/SPECIES_REFERENCE_REGISTRY.md` (one row per species: assembly, annotation source,
biotype convention, resolved paths, assertion floors) and
`docs/XSPECIES_DATASET_REGISTRY.md` (one row per run: accession, study, layout, measured
adapter, alignment metrics, and for each of four postures whether the data is local,
archived or never built).

932 GB of bulk data (FASTQ, and BAMs for four postures) was moved to the warm archive.
`scripts/archive_to_warm.sh` copies, verifies by checksum, indexes, writes a stub in place,
and deletes the source only after `rsync --checksum` reports zero differing paths. Every
archived path is listed in `data/ARCHIVE_INDEX.tsv` with the command to restore it. Retained
locally are the artifacts everything downstream reads: unique-posture Ribo BAMs, all P-site
and coverage HDF5 files, the seven packs, TPM tables, genome-coordinate bigwig tracks and
all call sets and result tables.

### Result tables underlying this report

| file | content |
|---|---|
| `results/xspecies_profile_eval.tsv` | per-nucleotide profile agreement, both architectures |
| `results/xspecies_orf_calls.tsv` | ORF-call precision/recall/F1 and per-type recall |
| `results/xspecies_orf_classes_{attn,mamba4}.tsv` | CDS / uORF / ncORF breakdown |
| `results/xspecies_mm1_vs_mm25_calls.tsv` | multimap posture comparison |
| `results/xspecies_model_vs_mm25_reference_attn.tsv` | model scored against both truth sets |
| `results/xspecies_align_qc.tsv` | per-run trimming and alignment metrics |
| `results/xspecies_biotype_census.tsv` | per-species annotation composition |

---

## Figures

Generated by `figures/X_xspecies_pub/make_xspecies_figures.py`; every value is read from a
result table on disk, and per-panel sources are recorded in
`figures/X_xspecies_pub/X_xspecies_pub_values.json`. PDF and PNG at 400 dpi.

**Figure 1. Study design.** (a) Sequencing runs per species, Ribo-seq and RNA-seq. The five
human RNA runs predate this study's download manifest and so appear as zero. (b) Transcript
universe per species, protein-coding and lncRNA. (c) Frame distribution of P-sites at each
library's dominant read length for the two *Drosophila* datasets tested and the yeast
positive control; the dotted line marks the 33% expected with no periodicity.

**Figure 2. ORF-call agreement between the model and the observed data.** (a) Precision
against recall per species, fully de novo variant. (b) Recall by ORF class; CDS is
recovered essentially completely in every species, and the uORF and ncORF classes carry the
discriminating information. (c) uORF recall ordered by value, annotated with class size and
pooled frame-0 P-site depth; the panel deliberately asserts no trend, because the ordering
follows neither phylogeny nor depth.

**Figure 3. Per-nucleotide profile agreement.** (a) Median per-transcript Pearson
correlation for the per-nucleotide profile and, separately, for the transcript total.
(b) The same profile correlation against a position-shuffled null. (c) The Spearman
tie-handling artifact: ordinal ranks as implemented (open triangles) against average-rank
ties (filled circles), plotted against the zero fraction of the observed profile.

**Figure 4. The multimap diagnostic.** (a) Ratio of ORF calls at up to 25 alignments per
read to calls at unique-only, per species. Chimpanzee is marked with a cross because five
of its eight runs failed the caller's periodicity gate under the permissive posture, so its
two call sets rest on different data. (b) Change in the model's F1 when scored against the
permissive call set instead of the unique-only one; the model's predictions are identical
in both rows and only the reference changes.
