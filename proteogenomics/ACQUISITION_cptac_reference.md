# CPTAC3 acquisition -- live-verified 2026-07-17 (PDC + GDC APIs)

Ground truth pulled from the live PDC (`pdc.cancer.gov/graphql`) and GDC (`api.gdc.cancer.gov`) APIs.
Program CPTAC3; GDC project CPTAC-3; dbGaP phs001287.

## Proteomics (PDC) -- fully OPEN (CC-BY, no dbGaP), small, signed CloudFront URLs, no auth

| Cohort | global-proteome study | .raw files | raw total | mzML.gz alt | experiment | cases/aliquots |
|--------|-----------------------|-----------:|----------:|------------:|------------|----------------|
| **GBM** (pilot #1) | **PDC000204** | 264 | **214.5 GB** | 49.7 GB | TMT11 | 111 / 111 |
| CCRCC | **PDC000127** | 575 | 440.2 GB | 118.9 GB | TMT10 | 124 / 208 |
| LUAD | **PDC000153** | 625 | 466.2 GB | 97.0 GB | TMT10 | 115 / 217 |

- LUAD is PDC000153 (NOT the stale PDC000150 in third-party tables -- verified empty on the live API).
- These are the base Discovery global proteomes (NOT phospho/acetyl, NOT the CCRCC DIA PDC000200).
- Format: Thermo `.raw` (`data_category="Raw Mass Spectra"`) or smaller open `.mzML.gz`
  (`data_category="Processed Mass Spectra"`). All 3 raw proteomes together ~1.1 TB -> fit the 5.2 TB
  headroom; proteomics is NOT the bottleneck.
- API: `getPaginatedUIFile(study_name, data_category, offset, limit){total uiFiles{file_name file_size
  md5sum}}` for names/sizes/md5; `filesPerStudy(pdc_study_id, data_category, acceptDUA:true){signedUrl{url}}`
  for wget-able URLs. Portal "Export File Manifest" URLs expire after 168 h -> generate per-cohort just
  before streaming. (Gotcha: `filesPerStudy` scalar fields return null on the current endpoint -- use
  `getPaginatedUIFile` for metadata.)

## RNA-seq (GDC, CPTAC-3) -- the ACCESS FORK

- **OPEN (no dbGaP):** only `Gene Expression Quantification` = STAR-Counts gene-level TSV, ~1-2 GB per
  cohort (GBM 238 files/1.0 GB). Gene-level only -- NO isoform TPM, NO salmon input.
- **CONTROLLED (dbGaP phs001287):** aligned BAMs (no FASTQ in GDC). Genome BAM ~9 GB/sample =
  GBM 2.04 TB / LUAD 4.17 TB / CCRCC 4.34 TB per cohort. Needed for salmon decoy-aware isoform TPM +
  the model's per-nt RNA-seq coverage input.
- API: `POST api.gdc.cancer.gov/files` filter `cases.project.project_id=CPTAC-3` + `experimental_strategy=
  RNA-Seq` + `cases.primary_site` (LUAD also needs `cases.disease_type="Adenomas and Adenocarcinomas"` --
  "bronchus and lung" mixes LUAD+LSCC). Open counts need no token; BAMs need GDC token + dbGaP approval.

## Matching

Join on the CPTAC case id: PDC `case_submitter_id` == GDC `cases.submitter_id` (both `C3L-`/`C3N-`).
PDC `biospecimenPerStudy(pdc_study_id)` gives case<->aliquot<->sample_type; matched cohort = the
case-id intersection (RNA-seq covers more cases than proteome, e.g. 220 brain vs 111 GBM proteome).

## Pilot order: GBM -> CCRCC -> LUAD (smallest raw proteome first).

## The fork (decision needed -- see PLAN.md)

Open-only path (STAR gene counts + sequence-based model, no dbGaP, instant, gene-level) vs controlled
path (dbGaP BAMs -> salmon decoy-aware isoform TPM + per-nt coverage input, multi-TB, weeks of approval).
The user's own null framing ("3-frame translation of expressed GENES >1 TPM") is gene-level -> open path
satisfies it. dbGaP has weeks of lead time, so if the isoform/coverage full study is wanted, start
approval in parallel with the open-only pilot.
