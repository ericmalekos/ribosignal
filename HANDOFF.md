> **HISTORICAL DOCUMENT -- do not use as a guide to the current project.**
> Written 2026-07-08, before any code existed, to kick the project off. It describes the work as not
> yet started and proposes an FM-embedding input that the evidence later rejected (Task 15: one-hot
> matches FM embeddings, so the released models carry no FM dependency).
> For current state read `README.md`, then `results.md`. Kept only for provenance.

# Ribo-seq signal prediction model -- session handoff

Written 2026-07-08 by a prior Claude session (the `expression_context_human`
translation/MS project). This doc hands off everything needed to start a NEW
project: **train a transformer to predict per-nucleotide Ribo-seq signal directly
from (FM embedding of the transcript) + (single-nucleotide-binned RNAseq coverage)**,
using Chothani matched primary-tissue Ribo-seq + RNAseq (GENCODE v49 / GRCh38).

All paths below are verified to exist as of 2026-07-08.

---

## 0. The goal (as stated by the user)

- Input per transcript: **FM embedding (orthrus and/or rinalmo)** + **per-nt binned
  RNAseq coverage**. Output: **per-nt (or per-codon) Ribo-seq P-site signal**.
- Data: **Chothani** (matched Ribo-seq + deep RNAseq, high-quality primary tissue).
- First/easy check: **within-tissue, held-out-chromosome** (gene-disjoint) test.
- Real target: a **general model** that works on entirely new datasets given deep
  RNAseq. So the evaluation must eventually be cross-tissue and ideally cross-dataset,
  not just held-out-chromosome within one tissue.

---

## 1. Project root & where data lives

All existing data is under the sibling project:
`/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/biotype_probe/expression_context_human/`
(referred to below as **`$ECH`**). Put NEW code/outputs in this new dir:
`/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model/`.

Read the sibling project's `results.md` / `methods.md` for how the Ribo-seq and
RNAseq were produced (Phases 28, 28b, 28c, 30a-30j).

---

## 2. Ribo-seq (the prediction TARGET)  -- READY

- **Per-sample BAMs**: `$ECH/data/riboseq_bam/` (237 BAMs, 581 GB). Two coordinate
  systems per sample:
  - `*.toTranscriptome*.bam` -- **transcriptome coordinates** (use THIS for
    per-transcript per-nt signal; no splicing to deal with).
  - `*.Aligned.sortedByCoord.out.bam` -- genomic coordinates.
- **Per-tissue pooled BAMs** (cleaned + merged, used for RiboCode): 
  `$ECH/data/ribo_per_tissue_bam/{TISSUE}.merged.bam` (9 tissues, 14 GB).
- **Cleaning posture** (Phase 30c/d): NH>1 multimappers dropped + rRNA/tRNA/miRNA
  loci dropped. **VERIFY** whether `data/riboseq_bam/` holds the cleaned or the raw
  per-sample BAMs before use (the per-tissue merges are definitely cleaned; per-sample
  cleaning was Phase 30c/d -- check `methods.md` Phase 30 and `samtools view -c` /
  NH-tag distribution on one file). If raw, re-filter or use the merges.
- **P-site offsets** (critical for reads -> P-site position): RiboCode wrote
  per-read-length offsets per sample in `$ECH/data/ribocode_per_tissue/{TISSUE}/{TISSUE}_pre_config.txt`.
  Example (Brain, SRR15513148): read lengths `25,27,28,29,30` -> P-site offsets
  `12,11,12,12,12`, strand `yes`. Use these (or re-derive with riboWaltz/plastid) to
  convert footprint 5' ends -> P-site nt. Periodicity is clean (Chothani is
  high-quality); RiboTISH/RiboCode QC PDFs are in the per-tissue dirs.
- **RiboCode ORF calls** (if you want ORF-level labels too):
  `$ECH/data/ribocode_per_tissue/{TISSUE}/{TISSUE}.txt` (use `.txt`, NOT
  `_collapsed.txt` -- collapse drops ~80% of unique lncRNA tx).
- Ribo-seq input fastqs + SRR->tissue map: `$ECH/data/riboseq_input/manifest.tsv`
  + `manifest_phase28b.tsv`.

## 3. RNAseq (the INPUT feature) -- FASTQs READY, per-nt coverage NOT built yet

- **fastqs**: `$ECH/data/rnaseq_input/{SRR}_1.fastq.gz` + `_2.fastq.gz`
  (**paired-end**, 57 samples). Manifest w/ tissue: `$ECH/data/rnaseq_input/manifest.tsv`
  (cols: `tissue, srr, fastq_url, size_mb`).
- **What exists**: only **salmon decoy-aware quant** (transcript-level TPM) at
  `$ECH/data/salmon_quant_chothani_decoy/{SRR}/quant.sf`. That is NOT per-nt.
- **BUILD STEP REQUIRED**: for "single-nucleotide-binned RNAseq" you must **align the
  RNAseq fastqs** to get per-nt coverage. Two options:
  - STAR to genome (`$STAR_INDEX` below) then project to transcriptome, OR
  - STAR to transcriptome / align directly in transcript coords (simplest to align
    target + input on the same axis).
  - Check **strandedness** first (Chothani RNAseq is likely reverse-stranded Illumina;
    confirm with `infer_experiment.py` or salmon's inferred libtype).
- Do alignment on **SLURM** (multi-sample, heavy); fastqs already downloaded.

## 4. FM embeddings -- SEQUENCE-AVERAGED exist; PER-TOKEN must be extracted

- **Sequence-averaged (one vector/transcript)** -- NOT usable for per-nt prediction:
  `/private/groups/carpenterlab/emalekos/RNAZoo_meta/embeddings/human_{pc,lncRNA}/{orthrus,orthrus4t,orthrus6t,rinalmo}/{fm}_out/`.
- **Per-token / per-position embeddings** (what you need to predict per-nt signal):
  only a **pilot** exists -- `$ECH/data/phase40_token_pilot/rinalmo_out/{tx}_tokens.npy`
  (444 tx). **BUILD STEP REQUIRED**: extract per-token embeddings for all training
  transcripts. Template extractor: `$ECH/scripts/phase40_extract_token_embeddings.sbatch`.
  - **rinalmo**: per-nucleotide transformer, 1280-d/token. Full-length extraction has
    an L cap (~10,000 nt) and needs an OOM guard: peak GPU reserved ~= 2.7 + 5.7e-9*L^2 GB
    (see prior memory `project_rinalmo_lncrna_embeddings`). Use resume + per-tx files.
  - **orthrus**: Mamba over nucleotides (4-track = ACGT, 6-track adds CDS + 5' splice),
    512-d. Per-token = per-nt of the spliced transcript. GPU-only.
  - Both give per-transcript-nt vectors that align 1:1 to the transcript sequence ->
    align to per-nt Ribo-seq on transcriptome coords. That shared axis is the point.

## 5. References / indexes (group fs, prebuilt, reuse -- do NOT rebuild)

- Genome FASTA: `/private/groups/carpenterlab/emalekos/genomes/GRCh38.primary_assembly.genome.fa`
- STAR index (GRCh38 + GENCODE v49): `/private/groups/carpenterlab/emalekos/genomes/star_index_grch38_v49/`
- Salmon decoy index (v49 + GRCh38 decoys): `/private/groups/carpenterlab/emalekos/genomes/salmon_index_decoy_v49/`
- GTF: `/private/groups/carpenterlab/emalekos/RNAZoo_meta/annotations/gencode.v49.annotation.gtf`
- Transcript FASTAs: `.../annotations/gencode.v49.{pc,lncRNA}_transcripts.fa`
- ncRNA drop BED (rRNA/tRNA/miRNA): `$ECH/data/phase30/ncrna_filter.bed`
- For **held-out-chromosome splits**: derive tx->chromosome and gene->chromosome from
  the GTF. Gene-disjoint by chromosome is the clean split; the sibling project used
  gene-grouped splits (`$ECH/data/splits/human/`) you can reuse the logic from.

## 6. Matched tissues (Ribo-seq AND RNAseq both present)

8 matched tissues: **Brain, ES, Fat, Fibroblast, HCAEC, Hepatocytes, HUVEC, VSMC**.
- **HA_EC has Ribo-seq but NO RNAseq** -- exclude from matched-tissue training.
- **Fibroblast is by far the deepest** (~32 RNAseq PE samples, ~many Ribo) -- best
  choice for the first within-tissue, held-out-chromosome check.
- Match Ribo<->RNA by tissue via the two manifests (both tissue-labeled).

---

## 7. Key learnings from the prior work (relevant to THIS task)

1. **FM embeddings do NOT encode expression.** A Ridge probe recovers transcript
   length / GC from the embedding (R^2 ~ 0.9) but mean-TPM at R^2 ~ 0. This is the
   scientific justification for the two-input design: the **FM carries sequence
   grammar, the RNAseq carries abundance/context** -- they are orthogonal, so adding
   RNAseq to the FM is genuinely additive, not redundant.
2. **Ribo-seq STAR posture (if you ever re-align footprints):** `--outFilterMultimapNmax 1`
   (single-mappers only) + drop rRNA/tRNA/miRNA loci, and ALWAYS check raw fastqs for
   3' adapter `AGATCGGAAGAGC` (cutadapt `-a AGATCGGAAGAGC --minimum-length 20
   --maximum-length 40`). The existing BAMs already follow this; RiboCode/RiboTISH
   do NOT drop multimappers themselves.
3. **Coordinate frame:** do everything in **transcriptome (spliced) coordinates** --
   Ribo-seq `toTranscriptome` BAM, orthrus/rinalmo per-token embeddings, and RNAseq
   coverage all live on the transcript-nt axis. Avoids splicing bookkeeping.
4. **Generalization is the hard part, and Ribo-seq transfer is known to be hard.**
   In the sibling project, cross-species ORF-level Ribo-seq transfer collapsed to
   ~0.55-0.62 AUROC on the honest expressed universe (vs 0.87 within-tissue). That
   was classification not signal-regression, so not directly comparable, but treat
   cross-dataset generalization as the real test and expect a gap vs held-out-chromosome.
5. **Storage / compute rules (cluster `prism`, SLURM):**
   - Big outputs on the **group fs** (`/private/groups/carpenterlab/emalekos/...`),
     never under `$HOME` (tight quota, confusing "no space" failures).
   - Scratch: `/data/tmp/emalekos/...` (node-local), never `/tmp`.
   - Heavy work (alignment, embedding extraction, training) -> **SLURM** (`sbatch`,
     partitions: short 1h / medium 12h / long 14d / gpu). Head node `mustard` only for
     light tests + serial downloads. GPUs via `gpu` partition.
   - Downloads: serial on head node (compute nodes throttle to ~0.4 MB/s + trip ENA
     rate limits). fastqs are already local, so N/A here.
   - conda/micromamba envs live at `/private/groups/carpenterlab/emalekos/conda_envs/`.
     A working stdlib+sklearn+matplotlib+xgboost python: `.../conda_envs/cas12a/bin/python`.
     Proteomics/comet tooling: `.../conda_envs/proteomics/` (not needed here).
6. **Per-token storage gotcha:** if you stack per-token embeddings as an `(N, L, D)`
   memmap, slicing the middle (length) axis across all rows is pathologically slow on
   the parallel FS (scattered reads). Read full rows and slice in memory, or store
   one file per transcript (the pilot does the latter).
7. **Singularity/Apptainer (if containerizing):** on phoenix nodes pass
   `--no-home --cleanenv --env PYTHONNOUSERSITE=1` (host `~/.local` site-packages leak
   into SIFs and shadow bundled packages); underlay mode blocks `--writable-tmpfs`
   (use `--env` + bind-mounts).

---

## 8. Suggested first steps for the new session

1. Read this doc + sibling `methods.md`/`results.md` (Phases 28, 28c, 30).
2. Confirm per-sample Ribo-seq BAM cleaning posture (Section 2 VERIFY).
3. Pick **Fibroblast** (deepest). Build the target: transcriptome-coord P-site
   coverage per transcript-nt (apply RiboCode per-read-length offsets), pooled over
   Fibroblast Ribo-seq samples.
4. Build the RNAseq input: align Fibroblast RNAseq fastqs (SLURM), produce per-nt
   transcriptome coverage (check strandedness first). Bin to single-nt.
5. Extract per-token FM embeddings (rinalmo and/or orthrus) for the Fibroblast-
   expressed transcript set (SLURM GPU; OOM guard; one file per tx).
6. Define held-out-chromosome (gene-disjoint) split from the GTF.
7. Train a small transformer: inputs = per-nt [FM emb ; RNAseq coverage], target =
   per-nt P-site signal; loss e.g. Poisson / negative-binomial NLL or correlation-
   based; evaluate per-transcript Pearson/Spearman + 3-nt periodicity recovery on the
   held-out chromosome.
8. Only then scale to all 8 matched tissues + leave-one-tissue-out, then cross-dataset.

## 9. Maintain project docs (global rule)

Create `methods.md`, `results.md`, `logs/` in this dir as you go (per the user's
standing rule for any multi-step project). Figures get a sibling `FIGURE_DATA_INPUTS.md`.
Writing style: no en/em dashes (ASCII only); no "our"; never credit Claude/AI as author;
never commit/push without explicit user instruction.
