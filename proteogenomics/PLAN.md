# Proteogenomics: Ribo-seq model as a translation prior for novel-protein discovery

Kicked off 2026-07-17. Application of riboseq_signal_model. Supersedes the CPTAC plan
(PLAN_cptac_superseded.bak) -- CPTAC RNA-seq is dbGaP-controlled (weeks of approval, multi-TB BAMs),
so we moved to fully OPEN, same-sample, FASTQ-available datasets.

## Hypothesis + design

The model predicts WHICH ORFs are translated (from sequence, with RNA-seq expression context). Those
model-selected ORFs build a custom MS search database; searching matched proteomics should find novel
peptides at a HIGHER-THAN-NULL rate. This tests the sibling project's finding (project_expr_context_
human_ms_fdr): the model's real edge is Ribo-seq TRANSLATION, not MS detectability.

- **Test DB**: canonical proteome + model-selected translated novel ORFs.
- **Null DB**: canonical + 3-frame translation of all >1 TPM expressed transcripts (no model selection).
- **Readout**: novel-peptide IDs at CLASS-SPECIFIC FDR, model-DB rate vs null-DB rate, with PC-churn +
  net-PSM columns (feedback_db_tradeoff_table_columns). Differentiator vs moPepGen (Nat Biotech 2025,
  which built DBs from variants/splicing on these same CCLE proteomes): the ORF selector here is a
  Ribo-seq-trained TRANSLATION model, benchmarked against the 3-frame-expressed null.

## What the model needs from RNA-seq (CORRECTED 2026-07-17 -- earlier draft wrongly called coverage optional)

The model has TWO distinct RNA-seq needs, and they are different data products:
- **Per-nt RNA-seq coverage = a REQUIRED model INPUT** (the coverage channel the model was trained on,
  input_mode "both"). It comes from **STAR alignment** (`--quantMode TranscriptomeSAM` -> per-nt
  transcriptome coverage via `rnaseq_coverage.py`), NOT salmon. This is NOT optional here: it is what
  makes the model's translation prediction A549-SPECIFIC. Run sequence-only and the calls are
  tissue-AGNOSTIC (identical for every cell line), defeating the per-sample search-DB idea. (The Task 8/15
  "sequence-only 91%" number is from a SEPARATELY-trained emb-only model, not the deployment model run
  with zeroed coverage -- do not use it to justify skipping the alignment.) Faithful recipe =
  `align_rnaseq_sample.sh` (STAR v49 `--outFilterMultimapNmax 20`, posture A, auto-strand coverage).
- **Salmon decoy-aware isoform TPM = a SEPARATE product** for the expressed-transcript filter (>1 TPM) +
  the 3-frame null. Runs on FASTQ; != coverage.
Both come from ONE FASTQ revert of the open BAM (`process_rnaseq_a549.sbatch` does revert -> salmon +
STAR->coverage together). Corrected index paths: STAR `STAR_indexes/star_index_grch38_v49`, salmon
`Salmon_indexes/salmon_index_decoy_v49` (CLAUDE.md's `genomes/` paths are stale/missing).

## Datasets (fully open, same-sample)

- **CCLE = workhorse** (Nusinow 2020). 375 cancer cell lines, deep TMT10/SPS-MS3 whole proteome, RAW OPEN
  at MassIVE **MSV000085836**. CCLE RNA-seq OPEN (DepMap TPM + open GDC/AWS BAMs -> FASTQ). Same-sample by
  line. moPepGen-benchmarked. Weakness: tryptic whole-proteome (lower cryptic-ORF yield); TMT search mods.
- **Ouspenskaia 2022 open subset = novel-ORF showcase** (Nat Biotech). MHC-I immunopeptidome (cryptic
  ORFs >20x enriched vs whole proteome) + matched Ribo-seq + RNA-seq. OPEN MS on MassIVE
  (MSV000084787/084442/084172/080527) + Ribo-seq/RNA-seq on GEO **GSE143263** + their **nuORFdb**. Use
  the OPEN cell-line/GEO subset ONLY (patient tier = dbGaP phs001998/001451/001519 -- avoid). Bonus:
  matched Ribo-seq validates the model's translation calls directly.
- **Chothani 2022 (GSE182377)** = the model's own open Ribo-seq+RNA-seq training data -> validation
  reference / smORF truth set, not a matched-MS discovery set.
- Scale-up later: ProCan (PXD030304, 949 lines DIA -- harder custom-DB re-search).

## Pilot (start here): ONE CCLE cell line, end-to-end smoke test

1. Pick a well-characterized line in BOTH proteome + RNA-seq (recon agent recommending; e.g. A549/MCF7/
   K562/HCT116).
2. RNA-seq: DepMap TPM (instant) for the expressed set; open FASTQ/BAM -> salmon decoy-aware isoform TPM
   (+ STAR coverage if wanted).
3. Model -> translated-ORF calls (deployment model + all-alt-ORF caller) -> class-labeled novel ORFs.
4. DBs: canonical / model-novel / 3-frame-null + decoys (reuse build_msfragger_db.py,
   phase76b_build_novel_dbs.py).
5. Download that line's ONE TMT plex from MassIVE -> convert to mzML -> search (MSFragger 4.2 + Comet) ->
   DELETE raw.
6. MS2Rescore (install) -> Percolator; class-specific FDR (phase60_class_specific_fdr.py); model-DB vs
   null-DB novel-peptide rate. Validate the full chain, then scale across lines + add Ouspenskaia.

## Tool locations (verified 2026-07-17)

MSFragger 4.2 `conda_envs/msfragger/share/msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar` (+ own JVM);
Comet `conda_envs/proteomics/bin/comet`; reuse sibling DB/FDR/Percolator (build_msfragger_db.py,
phase76b_build_novel_dbs.py, phase60_class_specific_fdr.py, phase60e/61e_percolator_driver.sh); salmon
decoy-aware human GENCODE v49 index (shared prism). MS2Rescore = fresh install.

## Quota

Group ceph 15 TB SHARED, ~5.2 TB free. Now MUCH easier than CPTAC: a single CCLE TMT plex is tens of GB,
one cell line's RNA-seq is a few GB. Still stream + delete raw per plex; watch getfattr ceph.dir.rbytes.

## Open / next

- Recon agent: CCLE MassIVE plex structure + sizes, DepMap expression files, open RNA-seq raw route for
  one line, cell-line->plex map, pilot line recommendation.
- Then: pull the pilot line's RNA-seq (DepMap TPM first) + run the model + build DBs while the plex
  downloads.

## Standing rules

ASCII only; no "our"; no AI-authorship; no commit/push without ask; scratch to /data/tmp/emalekos; group
ceph 15 TB SHARED (stream + delete raw); salmon decoy-aware; maintain methods.md/results.md/logs.
