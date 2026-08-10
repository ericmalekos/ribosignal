# DATA PROVENANCE LOG -- proteogenomics/

PROJECT RULE (2026-07-17, user-set): every external dataset touched by this project is logged HERE before
or at download, with full provenance. No data is used without an entry. Be careful and explicit about:
- **raw vs processed**: raw FASTQ (all reads) vs someone's aligned/filtered BAM (reads dropped) vs a
  quant matrix. Reverting a BAM to FASTQ is NOT raw FASTQ -- note it if ever done.
- **access tier**: open vs dbGaP/EGA-controlled.
- **genome/annotation build** (hg19/GRCh37 vs GRCh38/v49) and library strandedness.
- **matching**: same-sample vs same-cell-line-different-project (cross-project) vs cohort-level.
- **checksum**: record md5/size and whether it was verified against a published checksum.
- **how used** + **license**.

Status key: ACTIVE = in use; REJECTED = evaluated + not used (say why); ABANDONED = downloaded then dropped.

| id | what | source (accession / URL) | build/format | access | pulled | size / checksum | matching | how used | status |
|----|------|--------------------------|--------------|--------|--------|-----------------|----------|----------|--------|
| ccle-a549-prot | CCLE A549 whole-proteome TMT10, plex 12 (12 fractions, ch 127n) | MassIVE **MSV000085836** `ccms_peak/raw/d00615..d00626_Prot_12_*.mzML` (Nusinow 2020, Cell) | mzML (CCMS-converted) | OPEN (CC0) | 2026-07-17 | ~10.7 GB; size-checked, md5 NOT yet verified | -- | A549 MS search target (the proteomics being predicted) | ACTIVE |
| depmap-24q4-expr | DepMap 24Q4 expression: gene TPM + Model.csv + OmicsProfiles.csv | figshare+ DOI 10.25452/figshare.plus.27993248 (files 51065489/51065297/51065723) | CSV, log2(TPM+1) | OPEN | 2026-07-17 | 484 MB gene TPM | A549 = ACH-000681 | cross-check of A549 expression (gene-level) | ACTIVE (reference) |
| ccle-a549-rnabam | CCLE A549 RNA-seq **aligned BAM** | `s3://depmap-omics-ccle/data/rna/bam/G20499.A549.2.bam` | BAM, **hg19/GRCh37**, coord-sorted | OPEN | 2026-07-17 | 14.0 GB (byte-verified) | A549 same line | -- | **REJECTED**: aligned+filtered BAM, not raw FASTQ (reverting drops unaligned reads). Deleting to reclaim quota. |
| a549-rnaseq-raw | deep A549 raw RNA-seq FASTQ (cross-project), ENCODE **ENCSR000CON** = untreated A549 polyA+ RNA-seq (Gingeras/CSHL 2011) | ENCODE files ENCFF000EJJ/EJV (rep1) + ENCFF000EJW/EKB (rep2); URLs `encodeproject.org/files/<ENCFF>/@@download/<ENCFF>.fastq.gz` | FASTQ (raw), **PE 2x76, stranded-REVERSE** (salmon ISR / STAR -s2), polyA+ | OPEN | 2026-07-17 | 29 GB (4 files); **all 4 md5 VERIFIED** vs ENCODE (EJJ 76e7.., EJV 7cf5.., EJW ce7d.., EKB 8083..) | A549 line, DIFFERENT project than the CCLE proteome (real-world unmatched) | UC-A: model RNA input (STAR coverage + salmon TPM) to predict ccle-a549-prot | **ACTIVE**; FASTQ DELETED post-STAR (kept: coverage hd5 + salmon quant.sf; re-fetchable via the 4 ENCFF accessions + md5s above) |
| hbl1-imp | HBL-1 (DLBCL) MHC-I immunopeptidome, W6/32 pulldown, 2 reps (Ruiz-Cuevas 2021) | PRIDE **PXD020620** `HBL1_DMSO_200M_050219_{1,2}.raw` | Thermo **.raw** (NO mzML), Q Exactive HF, **non-tryptic label-free** MHC-I | OPEN | 2026-07-18 | 0.89 GB; **SHA-1** verify (1: 89354c8d.., 2: 6eb85483..) | same-sample w/ hbl1-rna | UC-B: immunopeptidome search target (>20x cryptic-ORF enriched); HLA params (no-enzyme, 8-14mer, no TMT, 0.02 Da high-res -- confirmed MS2 FTMS Q Exactive HF, 15.5K+14.2K MS2). **PROCESSING: .raw -> indexed mzML via ThermoRawFileParser 2.0 (proteomics env), default peak-pick=centroid -- msfragger conda build ships the Thermo reader as .exe (needs mono/.NET, absent), so reads .raw as 0-scan; mzML is the search input.** | ACTIVE |
| hbl1-rna | HBL-1 RNA-seq raw FASTQ, run SRR12285172 (paired, Oligo-dT, NextSeq500, 57.5M pairs) | SRA **PRJNA647736** / ENA SRR12285172_{1,2}.fastq.gz | FASTQ (raw), paired | OPEN | 2026-07-18 | 4.87 GB; **MD5** verify (1 d2870784.., 2 7ad1f3f4..) | same-sample w/ hbl1-imp | UC-B: model RNA input (STAR coverage + salmon TPM) for HBL-1; HBL1 universe = 34,775 expressed tx | ACTIVE |
| hbl1-priceorf | authors' own PRICE Ribo-seq ORF DB for HBL-1 (their immunopeptidome search DB) | PRIDE PXD020620 `HBL1_Price_orfs.fasta` | protein FASTA (168,852 ORFs, PRICE ids `<n>_<strand>_NP`) | OPEN | 2026-07-18 | 43 MB | HBL-1 | UC-B: SECONDARY comparison point (their processed PRICE output). db_price = 136,238 unique novel ORFs after dropping 26,580 == canonical; tagged `nuORF|price|`. NOTE: PRICE uses a different ORF caller (non-AUG starts, own ORF taxonomy/boundaries) than the model's enumeration, so model-vs-PRICE conflates SELECTION signal with caller discordance -- superseded as the primary baseline by hbl1-ribo (RiboCode) | ACTIVE (secondary) |
| hbl1-ribo | HBL-1 **DMSO elongating Ribo-seq**, 2 reps (D1+D2), raw FASTQ (Ruiz-Cuevas 2021) -- matched to the DMSO immunopeptidome | SRA PRJNA647736 / ENA **SRR12285192** (HBL1_D1_ribo, 416.9M reads) + **SRR12285191** (HBL1_D2_ribo, 386.2M reads) | FASTQ (raw), SINGLE-end, **3' TruSeq adapter attached** (AGATCGGAAGAGC, confirmed in read 1) | OPEN | 2026-07-18 | 17.1 GB (9.09+8.00); **MD5** verify (192 d2db106a.., 191 a2dd443b..) | same-sample w/ hbl1-imp + hbl1-rna | UC-B PRIMARY baseline: run RiboCode's OWN ORF caller on this (same lineage as the model's training P-sites) to remove PRICE's caller discordance -- isolates measured-vs-predicted translation. Recipe = training's: cutadapt -a AGATCGGAAGAGC -m20 -M40 --trim-n; STAR v49 multimap-20 TranscriptomeSAM (training params); filter_tx_bam_generic.py clean (ncRNA + cross-gene MM drop, keep within-gene = posture A); reuse v49 ribocode_annot; metaplots+RiboCode. H1/H2 (SRR12285174/173) = harringtonine TIS runs, NOT used (elongating only) | **ABANDONED 2026-07-19**: ~2% unique mapping (~98% rRNA -- poor depletion + rRNA absent from GRCh38 primary). Paper-matched processing was correct (u2, offset 12, f0 ~80%) and PRICE ran (19,801 ORFs, 17,322 non-canonical), but a 2%-unique library is too shallow/noisy for a trustworthy re-derived experimental baseline, and RiboCode full-run was anomalous (649 vs 4,755 on the subset). Controlled Ribo-seq baseline dropped. UC-B relies on model-vs-null (no Ribo-seq, 1.73x) + authors' PRICE as a caveated external reference. |
| ruizcuevas-more | other lines (DoHH2, SU-DHL-4) + HBL-1 reps -- scale-up after the pilot | PXD020620 / PRJNA647736 (SRR122851{80,71,75} etc.) | -- | OPEN | pending | DoHH2 31.5 GB, SU-DHL-4 47.2 GB RNA | -- | UC-B scale-up | DEFERRED |
| macro-tissue-prot | IPX0001245001 = 12 tissue-macrophage-population whole-proteome DDA (Microglia, Kupffer, Lung resident+recruited, Peritoneal, Large+Small intestinal, Spleen resident+recruited, Liver recruited, BMDM, RAW264.7), 3 reps x 6 fractions | iProX **IPX0001245000/IPX0001245001** (= PRIDE **PXD021583**), Nat Commun 2022 **s41467-022-35095-7** (art. 7389); `download.iprox.org/IPX0001245000/IPX0001245001/*.raw` | Thermo **.raw**, 216 files, fractionated DDA (SCX) | OPEN | 2026-07-24 | ~173 GB; **size-verified vs Content-Length** (iProX publishes NO md5), local md5 recorded in `downloaded_manifest.tsv`; downstream integrity re-checked at ThermoRawFileParser mzML conversion | same 12 populations as macro-tissue-rna = matched proteome+transcriptome, SAME study/samples | UC (test data): re-search .raw vs a model-informed DB built from macro-tissue-rna. **Skipped the 210 .msf** (authors' Proteome Discoverer results -- replaced by this pipeline, not needed) | DOWNLOADING (`data/macrophage_tissue/`) |
| macro-tissue-rna | PRJNA482293 = RNA-seq of the SAME 12 macrophage populations, 3 reps (36 runs) | SRA **PRJNA482293** / ENA fastq (one `SRR*.fastq.gz` per run) | FASTQ (raw), **SINGLE-END 1x150** (ENA/SRA layout is labelled PAIRED but the data is single-end: `fastq-dump --split-spot` gives 1 read/spot and `--split-3` one file -- R2 never deposited; **align single-end**), Mus musculus -> GRCm39/vM38 | OPEN | 2026-07-24 | 175 GB (36 fastqs); **md5 VERIFIED vs ENA fastq_md5** | same 12 populations as macro-tissue-prot (matched, same samples) | UC (test data): RNA input (STAR coverage + salmon TPM, single-end) -> model ORF calls -> proteogenomics search DB for macro-tissue-prot | DOWNLOADING (`data/macrophage_tissue/`) |

## Use cases

- **UC-A (cross-project, current):** deep A549 raw RNA-seq FASTQ (a549-rnaseq-raw, e.g. ENCODE) -> model
  translation calls -> DB -> search the CCLE A549 proteome (ccle-a549-prot). RNA and proteome from
  DIFFERENT projects, same cell line = the realistic unmatched-omics scenario.
- **UC-B (matched, in parallel):** a single study with open raw RNA FASTQ + deep matched proteome from the
  SAME samples. WINNER (ENA-verified) = **Ruiz-Cuevas 2021** (PRJNA647736 RNA+Ribo FASTQ + PXD020620 deep
  proteome + immunopeptidome, 3 DLBCL lines, same-sample, matched Ribo-seq). Runner-up = Ouspenskaia 2022
  open subset (B721.221 RNA+Ribo+mono-allelic immunopeptidome; MS split across MassIVE MSV000084787 etc.).
  Rejects (fail open-raw-FASTQ): Chong 2020 (EGA), CCLE (BAM-only), CPTAC/GTEx/TCGA (dbGaP).

## MHC-I expansion (Fig 2 pattern) -- added 2026-07-20
| id | what | source | build/format | access | size / checksum | matching | how used | status |
|----|------|--------|--------------|--------|-----------------|----------|----------|--------|
| dohh2-imp | DoHH2 (DLBCL) MHC-I immunopeptidome, 2 reps | PRIDE PXD020620 `DOHH2_400M_050219_{1,2}.raw` | Thermo .raw, Q Exactive HF, non-tryptic label-free | OPEN | 0.61 GB; size-checked | same-line w/ dohh2-rna | MHC-I expansion dataset (DLBCL, same study as HBL-1) | DOWNLOADING |
| dohh2-rna | DoHH2 DMSO RNA-seq (D1), paired | ENA PRJNA647736 SRR12285182_{1,2} | FASTQ paired | OPEN | 7.4 GB; MD5 (99c006.., 76ba26..) | same-line w/ dohh2-imp | model RNA input (expressed universe + coverage) | DOWNLOADING |
| sudhl4-imp | SU-DHL-4 (DLBCL) MHC-I immunopeptidome, 2 reps | PRIDE PXD020620 `SUDHL4_400M_050219_{1,2}.raw` | Thermo .raw, Q Exactive HF | OPEN | 0.68 GB; size-checked | same-line w/ sudhl4-rna | MHC-I expansion dataset | DOWNLOADING |
| sudhl4-rna | SU-DHL-4 DMSO RNA-seq (D1), paired | ENA PRJNA647736 SRR12285190_{1,2} | FASTQ paired | OPEN | 11.4 GB; MD5 (25bfc8.., 066b5d..) | same-line w/ sudhl4-imp | model RNA input | DOWNLOADING |
| b721-imp | B721.221 mono-allelic MHC-I immunopeptidome (Ouspenskaia 2022) | MassIVE MSV000080527 (mono-allelic) / MSV000084172 | .raw/mzML | OPEN | -- | same-line w/ b721-rna | MHC-I expansion (cross-study lymphoblastoid) | ACCESS-CONFIRMED, not yet pulled |
| b721-rna-ribo | B721.221 RNA-seq + Ribo-seq (Ouspenskaia 2022) | GEO GSE143263 (+ GSE131267 RNA) | FASTQ | OPEN | -- | same-line | model RNA input + orthogonal Ribo-seq validation + ribotish metagene QC | ACCESS-CONFIRMED, not yet pulled |

## GSE120762 mouse macrophage +/- LPS (2026-07-22) + Chothani mm1 re-download

| id | what | source | type | access | date | size | matching | use | status |
|---|---|---|---|---|---|---|---|---|---|
| gse120762-mouse | BMDM +/- LPS: 6 RNA-seq (3 WT/NT + 3 LPS) + 4 ribosome-profiling (CHX; 2 NT + 2 LPS) | GEO GSE120762 = SRA SRP163147 / PRJNA494404 (Jackson 2018, "Translation of Non-Canonical ORFs Controls Mucosal Immunity") | raw FASTQ | OPEN | 2026-07-22 | ~25 GB (10 samples; 6 RiboTag skipped) | mouse GRCm39 vM38 | NEW held-out eval: predict Ribo-seq per condition (NT vs LPS) from RNA coverage + sequence, check condition-specific translation | DOWNLOADING (7 MB/s, md5-verified) |
| chothani-rna-mm1 | Chothani training RNA-seq re-download for the mm1 retrain (7 tissues + Hepatocytes eval; Brain DROPPED) | ENA PRJNA756023 (55 of 57 runs) | raw FASTQ | OPEN | 2026-07-22 | ~150 GB (55 PE samples) | human GRCh38 v49 | retrain orf_v2_attn on mm1 coverage, 7 tissues (drop Hepatocytes-holdout + Brain-noise) | STAGED (launches after gse120762; md5-verified) |

## GSE155087 mouse CD4+ T cells (ZFP36/ZFP36L1 study) -- WT held-out (2026-07-25)

Study: GEO **GSE155087** (SRR12318326	PRJNA648398	SRP273454), "Multiomics analysis couples mRNA turnover and translational control of
glutamine metabolism to the differentiation of the activated CD4+ T cell" (Turner lab, Nat Commun 2022).
Mouse C57BL/6 activated CD4+ T cells. **WT/control samples ONLY** were taken (genotype
`Zfp36fl/fl Zfp36l1fl/fl control`, i.e. no CD4-cre); the `CD4-cre dKO` samples, the iCLIP (RIP-Seq)
libraries, and the 4sU-seq metabolic-labelling RNA were all EXCLUDED.

| id | what | source | type | access | date | size | matching | use | status |
|---|---|---|---|---|---|---|---|---|---|
| gse155087-tcell-ribo | Activated CD4+ T-cell **Ribo-seq, WT/control ONLY**, 5 reps (Ribo Ctrl 1-5) | GSE155087 / SRX8819596-8819600 -> **SRR12318326-330** (GSM4694968-972) | raw FASTQ, **SINGLE-end** | OPEN | 2026-07-25 | ~7.1 GB (5 runs; 46.3+101.7+28.3+35.1+32.2 M reads) | mouse GRCm39/vM38 | NEW cross-species + cross-cell-type held-out eval: predict Ribo-seq from RNA coverage+sequence -> RiboCode drop-in + ORF-call metrics (mirrors Wang/GSE120762 mouse held-outs). Adapter-check before STAR (Ribo-seq rule) | DOWNLOADING (md5-verified vs ENA, NPAR=3) |
| gse155087-tcell-rna | matched CD4+ T-cell **mRNA-seq, WT/control ONLY**, 3 reps (mRNA Ctrl 1-3) | GSE155087 / SRX8819590-8819592 -> **SRR12318320-322** (GSM4694962-964) | raw FASTQ, **SINGLE-end** | OPEN | 2026-07-25 | ~11.2 GB (3 runs; 37.2+47.7+44.3 M reads) | mouse GRCm39/vM38; matched to gse155087-tcell-ribo (same study, same WT CD4+ T cells) | model RNA input (STAR coverage + salmon TPM, single-end) -> T-cell universe + pack. **Chose standard mRNA-seq over the 4sU-seq** RNA (4sU is a metabolic-labelling turnover assay, not steady-state coverage) | DOWNLOADING |

## Macrophage search: mouse canonical proteome + mzML conversion QC (2026-07-30)

| id | what | source | type | access | date | size | matching | use | status |
|---|---|---|---|---|---|---|---|---|---|
| gencode-vM38-proteome | GENCODE vM38 protein-coding translations (canonical proteome) | GENCODE EBI FTP `release_M38/gencode.vM38.pc_translations.fa.gz` | protein FASTA, 66,668 seqs | OPEN | 2026-07-30 | 9.7 MB; **gzip -t verified** | mouse GRCm39/vM38 (matches macro-tissue universes) | canonical base + REV_ decoys for db_{canonical,model,null} per macrophage population (`build_a549_dbs.py --canonical`) | PRESENT (`annotations/gencode_proteins/`) |

**mzML conversion QC (macro-tissue-prot -- the integrity re-check anticipated at line 26):** 214/216
`.raw` converted to indexed mzML (ThermoRawFileParser 2.0; Orbitrap-Fusion-Lumos, ITMS/HCD centroid MS2).
**2 files fail conversion:** `Peritoneal_macrophage_Rep3_F6.raw` and `Spleen_recruited_macrophages_Rep1_F1.raw`.
Both are **byte-perfect** (on-disk md5 == `downloaded_manifest.tsv` md5, correct size), so NOT a download
problem; TRFP's vendor RawFileReader throws `Cannot get scan event for 1` at the first scan -- a vendor-library
parse bug on these 2 files (re-download cannot help). **EXCLUDED** from the MSFragger search via the >1 MB
size gate (Peritoneal + SpleenRecruited each retain 17 of 18 fractions). ProteoWizard `msconvert` is the only
recovery path if 216/216 is ever required.

## Third mouse-liver arm: GSE243134 Ribo-seq KEPT, its totalRNA REJECTED (2026-07-31)

Study: GEO **GSE243134** / PRJNA1016416, "Diurnal control of iron responsive element containing mRNAs
through IRP1/IRP2" (2024). Selected as the THIRD independent mouse-liver arm so Wang/Janich
disagreements can be adjudicated by majority rather than coin-flip. Chosen over GSE73554 (Atger)
because Atger shares PubMed 28475894 with Janich's GSE67305 (same Lausanne circadian lab lineage,
so not technically independent).

Scope: WILD-TYPE only (`IRP1_+/+`, `IRP2_+/+`; ZT5 and ZT12; Rep1-3) = 12 biological samples, each
with a matched RPF and totalRNA library. Knockout arms ignored. `_340_`/`_347_` run pairs are the
SAME library on two flowcells and are pooled per sample, not treated as replicates.

| id | what | source | type | access | date | size | matching | use | status |
|---|---|---|---|---|---|---|---|---|---|
| gse243134-liver-rpf | Mouse liver **Ribo-seq (RPF)**, WT only, 21 runs | GSE243134 / PRJNA1016416 | raw FASTQ, SINGLE-end 75 nt | OPEN | 2026-07-31 | 28.8 GB (with RNA); 40/40 aligned | mouse GRCm39/vM38 | **THIRD liver arm** for the three-way Wang/Janich/GSE243134 ORF-call Venn | **KEPT.** 40/40 BAMs, median 24.5M unique reads, 273M total unique RPF. 93.1% adapter, median insert 29 nt -> `cutadapt -a AGATCGGAAGAGC --minimum-length 20 --maximum-length 40`; STAR `--outFilterMultimapNmax 1` |
| gse243134-liver-totalrna | matched liver **totalRNA (Ribo-Zero)**, WT only, 19 runs | GSE243134 / PRJNA1016416 | raw FASTQ, SINGLE-end 80 nt | OPEN | 2026-07-31 | (in the 28.8 GB above) | mouse GRCm39/vM38 | intended as model RNA-coverage input + universe | **REJECTED -- library chemistry.** salmon 9.75-23.64% (median **17.73%**) vs Wang 91-93%; **93.6% of TPM in 10 structural transcripts** (`n-TKctt14` tRNA-derived alone = 82.8%, `Gm59647` = 7SL/SRP dup, both typed `lncRNA` so the biotype filter cannot drop them); TPM>=1 = 3,534; universe **3,321 tx** vs Wang 13,041 / Janich 10,431. The `NTX -gt 6000` guard in `build_gse243134_pack.sbatch` aborted before it propagated. Replaced by an external poly(A) arm |

**Root cause and the rule it produces.** Ribo-Zero removes rRNA but NOT tRNA, 7SL/SRP or snRNA, and
total RNA additionally retains intronic pre-mRNA that cannot map to a transcriptome index (hence 18%
salmon against 28-39% STAR *genome* unique-mapping on the same reads). Poly(A) selection removes this
class by construction. Alignment and trimming were correct; only the chemistry was wrong.

**RULE: the RNA arm must be poly(A)-selected, verified BEFORE download. SRA `library_selection` is not
acceptable evidence** -- Ribo-Zero total-RNA libraries are routinely filed as `library_selection = cDNA`.
Accept only `library_selection = PolyA`/`Oligo-dT`, GEO `!Sample_molecule_ch1 = polyA RNA`, or an
explicit oligo-dT / poly(A)-enrichment sentence in `!Sample_extract_protocol_ch1`.
Near-miss caught by this rule: **GSE73554 (Atger 2015)** is circadian-matched mouse liver filed as
`library_selection = cDNA`, but GEO states *"TruSeq Stranded Total RNA ... Ribo-Zero Gold depletion
set"* and `!Sample_molecule_ch1 = total RNA`. It would have reproduced this failure exactly.

Excluded run: `SRR26055531` (RP125_347_..._RPF), failed run with 45,633 reads vs its `_340_` partner's
16M; the partner covers that sample.

## MHC-I expansion round 2 -- B721.221 + THP-1, access verified 2026-08-05

| id | what | source | build/format | access | size | matching | how used | status |
|----|------|--------|--------------|--------|------|----------|----------|--------|
| b721-rna | B721.221 RNA-seq, 8 reps (4x HLA-C*04:01, 4x C*07:01 transduced) | ENA PRJNA543098 (GEO GSE131267, Sarkizova 2020 NBT, PMID 31844290) | FASTQ, paired | OPEN, ENA-verified | 14 GB on disk, 16/16 files md5-verified | same PARENTAL line as b721-imp/ribo, different HLA transduction | model RNA input (expressed universe + coverage) | **COMPLETE 2026-08-05** |
| b721-ribo | B721.221 Ribo-seq, 7 runs / 4 HLA alleles (A0101 x2, A3303, B1501 x2, B4402 x2 -- GEO notes A0101/B1501/B4402 were sequenced twice for depth) | ENA PRJNA599422 (GEO GSE143263, Ouspenskaia 2022) | FASTQ, single-end | OPEN, ENA-verified | 39.9 GB, 2.20e9 reads (B721 runs only; project total 62.3 GB incl. A375/HCT116/melanocyte, not pulled) | same parental line | ORTHOGONAL Ribo-seq validation + ribotish metagene QC (D5) | PULLING (md5-verified, serial head-node) |
| b721-imp | B721.221 mono-allelic MHC-I immunopeptidome | MassIVE MSV000084172 (Sarkizova) / MSV000080527 (Abelin 2017) / MSV000084787 (Ouspenskaia) | .raw/mzML | OPEN, all 3 resolve | not yet sized | same parental line | MHC-I expansion dataset | VERIFIED, allele subset not yet chosen |
| thp1-rna | THP-1 (AML) RNA-seq, 3 reps | ENA PRJNA686824 (Nesvizhskii/Purcell immunopeptidogenomics, MCP 2021) | FASTQ | OPEN, ENA-verified | 3.4 GB | same line as thp1-imp | model RNA input | VERIFIED, not yet pulled |
| thp1-imp | THP-1 HLA-A*02:01 immunopeptidome, 3 reps | PRIDE PXD015039 (public 2020-05-13, Orbitrap Fusion) | .raw | OPEN, PRIDE API confirms | not yet sized | same line as thp1-rna | MHC-I expansion (cross-cancer, AML, different lab) | VERIFIED, not yet pulled |

### Verification notes (2026-08-05)

- **THP-1 answers the open FIGURES_PLAN question ("verify raw MS + RNA-seq open"): YES, both sides
  open.** The MCP 2021 immunopeptidogenomics paper generated the RNA-seq (PRJNA686824) and RE-ANALYSED
  a previously published immunopeptidome (PXD015039); their own re-analysis outputs are MSV000086922 /
  RMSV000000338.1. Use PXD015039 as the primary MS source.
- **thp1-rna poly(A) status: LIKELY but NOT STATED (checked 2026-08-05).** SRA says
  `library_selection = unspecified`. The paper's methods give only: "The library was prepared using an
  MGIEasy-stranded mRNA chemistry V2 kit, and sequencing used MGITech MGISEQ2000RS hardware". That is
  an mRNA-capture chemistry (oligo-dT based), so poly(A) is the strong inference -- but the authors do
  NOT state the selection mechanism, and no rRNA-depletion method is mentioned either.
  DO NOT treat this as confirmed. Per `feedback_rnaseq_must_be_polya`, verify EMPIRICALLY after
  download: run salmon and check (a) mapping rate, (b) expressed-universe size against the other human
  lines, (c) presence of tRNA / 7SL / other structural RNA that poly(A) selection should have removed.
  A Ribo-Zero total-RNA arm crushed the GSE243134 liver universe to 3,321 tx and was only caught that
  way. The dataset is 3.4 GB, so the empirical check is cheap.
- **CAVEAT (b721-rna): Smart-seq2 + Nextera XT.** Oligo-dT primed, so poly(A)-selective and NOT a
  total-RNA risk. But it is a full-length low-input protocol with markedly less uniform per-nucleotide
  coverage than the TruSeq libraries used for HBL-1 / DoHH2 / SU-DHL-4 -- and per-nt coverage is exactly
  what the model consumes. Check coverage uniformity + salmon mapping rate against the other lines
  before trusting the universe.
- **CAVEAT (b721 matching): the RNA-seq and the Ribo-seq/immunopeptidome come from DIFFERENT STUDIES.**
  RNA = Sarkizova 2020 HLA-C-transduced lines; Ribo-seq = Ouspenskaia 2022 HLA-A/B mono-allelic lines.
  Same B721.221 parental line, different HLA transductions. Acceptable for an expressed universe (HLA
  transduction should not materially move the transcriptome) but it is NOT same-sample matching, unlike
  HBL-1/DoHH2/SU-DHL-4 where RNA and MS come from one study.
- Group quota at verification time: 14.08 / 15 TB (0.92 TB free). Sized before pulling for that reason.

### Ribo-seq handling notes (b721-ribo)

Before ANY alignment, per the standing Ribo-seq rules:
1. **Adapter check on the FULL downloaded file**, not a streamed slice --
   `zcat X.fastq.gz | head -400000 | awk 'NR%4==2' | grep -c AGATCGGAAGAGC`. >50% positive means pass
   `-a AGATCGGAAGAGC` to cutadapt. A byte-range slice once reported 0% where the real file was 93%.
   Sanity-check against insert biology: footprints are ~30 nt, so a "trimmed" 75 nt read is a
   contradiction.
2. **STAR `--outFilterMultimapNmax 1`** -- single-mappers only. RiboCode's `process_bam.py` does NOT
   drop multi-mappers, so they silently inflate P-site counts and break periodicity.
3. **Drop rRNA / tRNA / miRNA / Mt_rRNA / Mt_tRNA loci** before the caller sees the BAM.
4. **Check the unique-mapping rate before trusting the library.** HBL-1's Ribo-seq was abandoned at
   ~2% unique (~98% rRNA); a large BAM does not mean usable reads. This dataset is the replacement
   for that abandoned arm, so the same check decides whether it is usable.

### BLOCKER: MassIVE FTP is unreachable from prism (2026-08-07)

The B721.221 immunopeptidome MS lives only on MassIVE (MSV000084172 Sarkizova mono-allelic,
MSV000080527 Abelin 2017). Both datasets resolve via the MassIVE proxi API over HTTPS and are
confirmed OPEN, but the FILES are only served over FTP and **port 21 is filtered from this cluster**:

  massive.ucsd.edu:21   no connection (filtered)
  massive.ucsd.edu:443  open  (metadata API works; no file-listing endpoint)
  ENA / PRIDE :443      200   (control -- both fine)

Tried and failed: `proxi/v0.1/files` (404), `QueryDatasets` (400), `DownloadResult` (405),
ftp v01/v02/v06 paths (all timeout), ProteomeCentral PX mirror (no record for either MSV).

CONSEQUENCE: B721 stays a Ribo-seq VALIDATION dataset (327M footprints, 26,920 measured ORF calls)
and does NOT join the immunopeptidome panel. The panel remains HBL-1 / DoHH2 / SU-DHL-4 + THP-1
once pulled, i.e. 4 of the 6 in locked decision D2.

TO UNBLOCK, any one of: (a) a network exception for massive.ucsd.edu:21, (b) download on a host with
FTP egress and copy to the group fs, (c) ask the authors for an HTTPS/S3 copy. Not actionable from
inside this session.

## THP-1 PULLED (2026-08-07) -- 4th immunopeptidome, both sides on disk

| id | what | source | tier | size | integrity |
|---|---|---|---|---|---|
| thp1-rna | THP-1 RNA-seq, 3 paired-end runs (SRR13279451/52/53) | ENA PRJNA686824 | OPEN | 6.7 GB (6 files) | **md5 verified against ENA `fastq_md5`, 6/6 OK, 0 failed** |
| thp1-imp | THP-1 BB7.2 immunopeptidome, 15 raw files | PRIDE PXD015039 | OPEN | 15.6 GB | see below -- PRIDE publishes no md5 |

Both from the same THP-1 (AML) line, so RNA-seq and immunopeptidome are matched at the cell-line
level. Downloaded SERIALLY on the head node per the standing rule, single stream per file.

**BB7.2 was chosen over the 21 W6/32 files in the same submission, deliberately.** W6/32 is a
pan-HLA class-I antibody and BB7.2 is HLA-A*02:01-specific. The narrower antibody gives a cleaner,
single-allele peptide population, which matches how the other three immunopeptidomes in the panel
were generated and avoids mixing allele specificities inside one search. The cost is fewer files
(15 vs 21) and therefore fewer spectra.

**Integrity: PRIDE publishes no per-file checksum**, so the download could only be size-checked, and
a size match is not sufficient (the ENA incident on SRR10846529 produced a correct-size, wrong-md5
14.5 GB file). The compensating check is the raw -> mzML conversion: a truncated `.raw` fails to
parse rather than silently yielding a short mzML. All **15/15 converted cleanly** (ThermoRawFileParser
2.0, `-f 2` indexed, default centroiding -- same settings as the macrophage and HBL-1 arms so search
results stay comparable), 4.2 GB of mzML, smallest 188 MB, 36,864 spectra in the first file. The
conversion job (`proteogenomics/scripts/thp1_raw2mzml.sbatch`) fails any file whose mzML is under
1 MB, so a silent partial conversion cannot pass.

**Poly(A) status is being verified empirically, not assumed** (standing rule: SRA
`library_selection=cDNA` does not prove poly(A) selection, and a Ribo-Zero total-RNA arm crushed the
GSE243134 liver universe to 3,321 transcripts). `process_rnaseq_thp1.sbatch` reports the fraction of
salmon TPM landing in rRNA / Mt / misc / sn(o)RNA biotypes; poly(A)-selected libraries sit in the low
single digits and total-RNA libraries run much higher.

## The poly(A) gate is now a real check, and it found something (2026-08-07)

`proteogenomics/scripts/polya_check.py`, run as step 0 of every pilot prep. **Two earlier versions
were silently wrong**, which is why this is a standalone script with calibrated thresholds rather
than an inline heredoc:

- **v1** looked `quant.sf`'s `Name` up in `tx2biotype.tsv`. Salmon carries the full pipe-delimited
  GENCODE header (`ENST00000641515.2|ENSG...|...`); the table is keyed on the bare versioned id.
  Every lookup missed, the junk total stayed 0.0, and it printed **"0.00%"** -- indistinguishable
  from a flawless library. A total-RNA arm would have printed the same thing.
- **v2** fixed the key split and added a "was enough TPM assigned a biotype" guard (97.8% assigned),
  and STILL reported 0.00%. Also meaningless: **the decoy-aware v49/vM38 salmon indexes contain
  pc + lncRNA transcripts ONLY.** There is no rRNA / snRNA / snoRNA / misc_RNA / Mt_rRNA / Mt_tRNA in
  the index, so no library can put TPM there. The biotype approach cannot detect total RNA through
  this index and would have stamped POLY(A)-SELECTED on a Ribo-Zero library.

**What actually discriminates** (calibrated on this project's known-good and known-bad libraries):
salmon mapping rate against the pc+lncRNA index; universe size at TPM >= 1; and above all the
**single-transcript TPM share**, because the structural contaminants that DO reach a pc+lncRNA index
are typed `lncRNA` and survive every biotype filter.

| library | mapped | tx TPM>=1 | max single-tx share | verdict |
|---|--:|--:|--:|---|
| GSE243134 liver totalRNA (known bad) | 21.5% | 4,451 | **79.8%** (`n-TKctt14`) | REJECT |
| Janich liver (in use) | 37-44% | 9,613-13,147 | **44-56%** (`Gm59647`) | **AMBIGUOUS** |
| Wang liver | 92.8% | 14,185 | 11.1% (`Alb`) | poly(A) |
| GSE302188 liver | 92.5% | 18,677 | 3.2% | poly(A) |
| THP-1 | 81.4% | 41,702 | 1.0% (`MT-CO3`) | poly(A) |

### The Janich flag is REAL, and chasing it down found a pipeline property worth knowing

All 7 Janich liver RNA replicates put **52-70% of salmon TPM into the same two transcripts that got
GSE243134 rejected**: `Gm59647` (a 7SL/SRP duplicate, 44-56%) and `n-TKctt14` (tRNA-derived, 7-15%).
Both are typed `lncRNA`, so no biotype filter removes them. Mapping rate is 37-44% throughout.

**But the model input is CLEAN.** In `data/packed_heldout_mouse_janich_liver`, both transcripts are
present and both are empty: `Gm59647` mean depth 7, `n-TKctt14` mean depth 0, together **0.000% of
the pack's 5.26e9 coverage mass**. The pack's top transcripts are Alb (4.5%), Trf, Apoe, Serpina3k --
ordinary liver mRNA.

The reason is the difference between the two quantification paths:

- **salmon** has no multimapper filter and uses EM, so it piles the reads from massively multi-copy
  structural loci (7SL, tRNA) onto these few transcripts, and TPM is length-normalised, which
  further amplifies short ones.
- **the pack coverage** comes from STAR with `--outFilterMultimapNmax 1` plus the ncRNA-locus filter,
  which removes essentially all of those reads before they reach the model.

**Consequence, stated precisely:** the contamination distorts Janich's UNIVERSE SELECTION (those two
transcripts made the TPM >= 1 cut carrying no real unique-mapping signal, and the TPM ranking of
everything else is depressed by ~50%), but it does NOT reach the per-nucleotide coverage the model
consumes. No trained model or prediction is affected. The flag is worth recording because **salmon
TPM and pack coverage can disagree by orders of magnitude on multi-copy loci**, so a library that
looks contaminated by TPM may still be fine as a model input, and the two must be checked separately
rather than one being taken as a proxy for the other.

### FIXED by read-level decontamination + re-quantification (2026-08-07)

`scripts/janich_decontam_quant.sbatch`: bowtie2 `--very-sensitive-local` against a 2-sequence index
(Gm59647 302 nt, n-TKctt14 100 nt), keep the non-aligning reads, re-run salmon against the
**unchanged shared** decoy-aware vM38 index. Output: `data/salmon_quant_janich_liver_decontam/`.

Reads were filtered rather than the index edited, because CLAUDE.md keeps one canonical index per
(tool, species, annotation) and a Janich-specific index would make its TPMs non-comparable with Wang
and GSE302188, the other two arms of the three-way liver ORF-call comparison. `--very-sensitive-local`
rather than end-to-end because a read may only partially overlap a 100 nt transcript.

| quantity | before | after |
|---|--:|--:|
| reads removed | -- | 6.2-11.0% |
| contaminant share of TPM | 52.4-70.4% | **0.0-4.9 TPM of 1e6** (~0.0003%) |
| top transcript | `Gm59647` (7SL) 44-56% | `Alb` / `Mup7` **6.2-7.6%** |
| transcripts at TPM >= 1 (mean) | 11,532 | **16,379** |

The top transcripts are now Alb, Mup7 and Mup14, i.e. ordinary highly-expressed mouse liver mRNA,
and the single-transcript share (6.2-7.6%) is in family with Wang liver's Alb at 11.1%. **The
corrected universe (16,379) now sits between Wang (14,185) and GSE302188 (18,677)**, where before it
was below both.

Only **6-11% of READS** carried **52-70% of TPM**, because TPM is length-normalised and these are
100-302 nt transcripts: a modest read count over a tiny effective length produces an enormous TPM.

**Cross-check against pure renormalisation.** Dropping the two rows from the original `quant.sf` and
rescaling to 1e6 predicts 16,469 transcripts at TPM >= 1; the read-level fix gives 16,379, a
difference of 39-144 transcripts per replicate (0.5-0.9%). So reads ambiguously mapping between a
contaminant and a real transcript are negligible here, and the cheap renormalisation would have been
a sound estimate. The read-level fix is still the one to keep, because it is exact and it also
removes those reads from the library-size denominator.

**What this did NOT fix, and could not.** Salmon mapping rate is 41.8% -> 40.0% and 37.2% -> 34.1%
across the fix, i.e. essentially unchanged (slightly lower, since the removed reads were ones that
DID map). The 34-42% rate is an independent property of the Janich library -- most plausibly
rRNA-depleted total RNA retaining intronic pre-mRNA that cannot map to a transcriptome index -- and
`polya_check.py` still flags the decontaminated quants on that axis alone. That flag is correct and
should stay; it is a statement about library chemistry, not about the contamination.

**PENDING DECISION (not actioned):** using the corrected universe means rebuilding
`data/packed_heldout_mouse_janich_liver` (10,431 tx -> ~15-16k), re-predicting both released models
on it, and re-running the three-way Wang/Janich/GSE243134 liver ORF-call comparison. That cascade
touches published comparisons and was left for an explicit go-ahead.

## Mouse immunopeptidome (D2, 6th dataset): candidates CHECKED and REJECTED (2026-08-07)

Recorded so these are not re-checked. The binding requirement is **matched RNA-seq of the same
sample**, because the model's input IS the query sample's own per-nucleotide RNA coverage; importing
another study's RNA-seq would reproduce the Söllner coverage-bias artefact that
`feedback_never_reuse_universe_across_datasets` exists to prevent.

| candidate | MS raw | matched RNA-seq | verdict |
|---|---|---|---|
| **PXD008733** (Schuster 2018, tissue-based murine MHC-I map, PMC6080492) | YES, 39 `.raw` + 60 mzXML, public 2018-06-19, C57BL/6 H2-Db/Kb, 19 tissues + 4 tumour lines | **NO** -- MS only, authors generated no transcriptome | **REJECTED**: no RNA-seq. Would require importing another study's mouse tissue RNA-seq. |
| **Rospo 2023** (MMRd CT26 colorectal, PMC10797964; H-2Kd/Dd/Ld) | **NO accession** -- data availability points only to "Additional file 2" (a processed peptide table) | YES, ENA PRJEB58630 | **REJECTED**: no raw MS. Scientifically the best fit (non-canonical antigens are its subject); worth an author request if the mouse point becomes required. |
| PXD020620 | n/a | n/a | **NOT MOUSE** (human DLBCL; it is the SU-DHL-4 / DoHH2 source already in the panel). |

CONSEQUENCE: the panel stands at 4 of the 6 in D2 (A549 tryptic + HBL-1 / SU-DHL-4 / DoHH2 HLA-I,
plus THP-1 landing). B721 is blocked on MassIVE FTP (above) and the mouse point has no verified
candidate. The FIGURES_PLAN lists the cross-species immunopeptidome as a **stretch** goal, so this
does not gate Fig 2; it should be stated as a limitation rather than left looking un-attempted.

---

## Three new human Ribo-seq datasets (2026-08-09, tasks 72-74)

All 22 FASTQs downloaded serially on the head node, single stream, **every one md5-verified against
ENA's published `fastq_md5`, 0 failures**. Manifests with per-run md5/bytes/layout at
`data/external/<dataset>/manifest.tsv`; fetch script `data/external/fetch_human_2026_08.sh`.

| dataset | GEO | runs | raw | arm composition |
|---|---|--:|--:|---|
| `GSE208041_thp1` | GSE208041 | 9 | 43.1 GB | 5 RPF + 4 RNA (paired) |
| `GSE39561_thp1` | GSE39561 | 3 | 11.5 GB | 3 RPF, **no RNA** |
| `GSE304796_cart` | GSE304796 | 6 | 16.6 GB | 3 RPF + 3 RNA |

### Library chemistry, established empirically -- each dataset needed a DIFFERENT recipe

Every one of these would have produced a plausible-looking BAM under the project's default
(`cutadapt -a AGATCGGAAGAGC`), and two of them would have been silently wrong.

**GSE208041 RPF -- N-PADDED, not pre-trimmed.** 99.54% of reads are exactly 35 nt with 0% TruSeq
adapter, which reads as "submitter already trimmed". They did trim, then **padded back to a fixed 35 nt
with N**: 34% of reads end in 8 Ns, 24% in 7. Aligning as-is charges every pad N as a mismatch.
Stripping the trailing N run recovers a textbook footprint distribution -- **84.4% in 26-34 nt, sharp
mode at 28 nt**. Recipe: strip trailing N, then `--minimum-length 20 --maximum-length 40`.

**GSE39561 -- POLY(A)-TAILED, not adapter-ligated.** 50 nt reads, 0% TruSeq adapter, so a naive check
calls them "not RPF-like". Wrong probe: the dominant 12-mer at read position 29 is `AAAAAAAAAAAA`
(59% of reads). This is the 2012-era Ingolia protocol, which poly(A)-tails the footprint before
reverse transcription. Cutting at the first poly-A run gives **69.6% in 26-34 nt, mode 26-27** --
genuinely Ribo-seq. Recipe: poly-A trim (`-a "A{6}"`), then the same length window. Standard TruSeq
trimming would have left the poly-A attached.

**GSE304796 -- the GEO labels ARE wrong, confirmed empirically.** GEO titles all six samples
"Ribo-seq", three as "RNA 1/2/3". After TruSeq trimming:

| samples | adapter read-through | 20-40 nt inserts | full-length 72 nt | verdict |
|---|--:|--:|--:|---|
| GSM9157460-62 ("RPF") | 72-75% | **51.7%** | 26.1% | ribosome footprints |
| GSM9157454-56 ("RNA") | 7-8% | 0.7% | **92.6%** | RNA-seq |

Decisive, and it matches the user's correction.

### Open issue: GSE304796 RPF inserts are 38-39 nt, not 28-31

The submitter's own GEO data-processing note says "**Unique molecular identifiers (UMIs) were pruned
from read sequences using UMI-tools**" (Eclipsebio pipeline v1) -- but SRA carries the RAW reads, so
the UMIs are still attached, and GEO does not state their length or position. 38-39 minus a ~10 nt UMI
lands on 28-29, which fits.

Per-position base composition and Shannon entropy over the first/last 12 nt of the insert do NOT
localise the UMI: entropy is 1.87-1.99 bits everywhere, in the RNA library as well as the RPF, so
composition cannot separate a random UMI block from aggregate genomic sequence. **To be resolved
empirically by 3-nt periodicity**, which is ground truth for Ribo-seq: trim a small grid of UMI
lengths and keep the one that maximises periodicity in the RiboCode metaplot. Do not guess.

### RNA chemistry: one arm is poly(A), one is NOT

- **GSE208041 RNA: poly(A)-selected.** Protocol states "Poly(A)-purified mRNA-seq libraries" via
  TruSeq Stranded mRNA. Paired 2x101. This is the only clean RNA of the three.
- **GSE304796 RNA: rRNA-DEPLETED total RNA**, not poly(A) -- "Both sample types were taken through
  rRNA depletion and RNA-seq samples were heat fragmented". This is the chemistry that crushed the
  GSE243134 liver universe to 3,321 tx (feedback_rnaseq_must_be_polya): Ribo-Zero total RNA retains
  tRNA and 7SL, which are typed `lncRNA` and so survive the biotype filter. **Run the poly(A) gate
  before packing and expect it may fail.**
- **GSE39561 has no RNA at all.** Its other three samples (GSM971691-93) are `THP-1_puro`, puromycin
  TIS-mapping, not RNA. The only route is pairing its cycloheximide RPFs with GSE208041's THP-1 RNA:
  same cell line, but a **CROSS-STUDY RNA/Ribo pairing** that must be labelled wherever its numbers
  appear, exactly like the GSE243134 Ribo + GSE302188 RNA arm.

### Disk

Group ceph was at 92.4% of its 15 TB quota (1.15 TB free) when the download started and 10.51 TB used
/ 4.49 TB free after -- other users freed space concurrently. Quota overflow kills SLURM jobs silently
with exit 120 (feedback_group_ceph_15tb_quota), so delete these 71 GB of FASTQ once the packs exist.
