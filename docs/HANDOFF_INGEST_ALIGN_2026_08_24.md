# Handoff: ingest, align, and archive public data (2026-08-24)

Supersedes `HANDOFF_DOWNLOADS_2026_08_24.md` (download-only). Scope is now:
**fetch -> align (4-way matrix) -> quantify -> archive.** No pack building, no model runs, no analysis.

---

## 0. Orientation — paths, environments, scripts

**Project root (all relative paths below are from here):**
```
/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
```
Shared references (NOT under the project root — one canonical copy, reused across projects):
```
/private/groups/carpenterlab/emalekos/genomes/          genome FASTA + GTF
/private/groups/carpenterlab/emalekos/STAR_indexes/     star_index_grch38_v49, star_index_grcm39_vM38
/private/groups/carpenterlab/emalekos/Salmon_indexes/   salmon_index_decoy_v49, salmon_index_decoy_vM38
/private/groups/carpenterlab/emalekos/conda_envs/       interpreters (below)
/private/groups/carpenterlab/emalekos/singularity_cache/  SIFs (below)
/private/warm-archive/carpenterlab/RNAZoo_meta/         warm archive, mirrors the project tree
```

**Environments — verified contents. There is NO single env that does everything.**

| env / image | has | LACKS | use for |
|---|---|---|---|
| `conda_envs/riboseq` | STAR, salmon, samtools, cutadapt, h5py, pysam, numpy | **torch** | alignment, quant, coverage hd5, BAM work |
| `conda_envs/ribocode` | RiboCode, h5py, pysam, numpy | torch, STAR | `ribocode_dropin.py`, ORF calling |
| `conda_envs/cas12a` | torch, numpy | **h5py** | torch-only utilities |
| SIF `...-orthrus-latest.img` | torch + **mamba_ssm** | — | any `--mixer mamba` train/predict (**CUDA-only**) |
| SIF `...-rinalmo-latest.img` | torch | mamba_ssm | transformer-mixer train/predict; runs on CPU |

Gotchas that have already cost time:
- `scripts/prepare/*` hardcode `cas12a`, which has no h5py — they die instantly on `import packlib`.
  Use `riboseq`.
- Always `singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1`. Without
  `--no-home --cleanenv` the host `~/.local` site-packages shadow the SIF's own.
- mamba_ssm has **no CPU fallback** (`causal_conv1d_fwd` asserts `x.is_cuda`), so even a 35-tx smoke
  test needs a GPU.

**Reusable scripts — prefer these over hand-rolling; each encodes a rule learned the hard way.**

| script | does |
|---|---|
| `scripts/bamlib.py` | BAM primitives. `multimap_cost()` compares postures by TOTAL RECORDS PER TX (never an NH split); `quickcheck()` catches truncation a header check misses |
| `scripts/archive_to_warm.sh` | copy → **checksum-verify** → index → stub → delete. Never hand-roll rsync |
| `scripts/fetch_heldout_rna.sh` | ENA fetch: serial, single-connection, md5-verified, 3 retries |
| `scripts/heldout/filter_tx_heldout.py` | transcriptome-BAM cleaner: ncRNA drop + cross-gene read drop |
| `scripts/rnaseq_coverage.py` | BAM → per-nt coverage hd5, auto strand/libtype |
| `scripts/prepare/prepare_pack.py` | build a pack; `--reuse-target` swaps ONLY the RNA channel |
| `scripts/dump_pred_profiles.py` | model → per-nt predictions npz (has the n_skip guard) |
| `scripts/eval_union.sbatch` | full standing-rule ORF-call eval; `PSUF`/`RUN`/`TAG`/`LABEL` overrides |
| `scripts/compare_mm10cov_arms.py` | test-metric table judged against the seed-noise floor |
| `scripts/check_figure_provenance.py` | enforces the figure `model`/`source` keys |

**Standing rules live in** `/private/groups/carpenterlab/emalekos/RNAZoo_meta/CLAUDE.md` (project)
and `~/.claude/CLAUDE.md` (cluster-wide). Read both before starting. Prior results and their
caveats are in `results.md`; design decisions in `methods.md`.

## 0b. DO NOT DUPLICATE

Two processes are live. Check both before starting.

| what | pid file | destination |
|---|---|---|
| Ruiz-Orera human iPSC-CM RNA (5 runs, 18.8 GB) | `/tmp/fetch_ruiz_pid.txt` | `data/external/heldout_rna_mm10/ruizorera_hsCM_rna/` |
| Reference genomes x6 (chained behind it) | `/tmp/fetch_refs_pid.txt` | `genomes/xspecies_refs/<label>/` |

Wait for `ALL_REFS_DONE` in `logs/fetch_refs_driver.log`.

**Already local -- verify at FULL DEPTH before fetching anything** (`find <root> -name "<ACC>*"`, no
`-maxdepth`; a 3.4 GB redundant Wang download happened on 2026-08-24 from a `-maxdepth 1` survey):
Drosophila r6.67 reference (`annotations/fly/`), Wang mouse liver RNA
(`data/external/wang2021_mouse/fastq/mouse_liver/`), Janich liver RNA, PRJEB34766, GSE208041,
GSE304796, GSE39561. Also check `data/ARCHIVE_INDEX.tsv` and any `*.ARCHIVED.md` stubs.

---

## 1. Data to ingest: PRJEB65856 non-human primates

45 runs, 221 GB. Full run table: `data/external/PRJEB65856_runs.tsv`.
Accession lists are in `HANDOFF_DOWNLOADS_2026_08_24.md` section 1 (still accurate).

| species | RNA | Ribo | total |
|---|---|---|---|
| gorilla (gg) | 6 runs / 38.0 GB | 3 runs / 18.1 GB | 56 GB |
| macaque (rm) | 10 / 49.3 GB | 7 / 27.7 GB | 77 GB |
| chimp (pt) | 11 / 58.9 GB | 8 / 29.2 GB | 88 GB |

Order: gorilla -> macaque -> chimp (smallest first). RNA is paired-end, Ribo single-end.
Fetch with `scripts/fetch_heldout_rna.sh <label> <ACC>...` (serial, head node, md5-verified, 3 retries).

---

## 2. The alignment matrix

**Per RNA library: 3 outputs.** **Per Ribo library: 2 outputs.**

| assay | tool | setting | why |
|---|---|---|---|
| RNA | STAR | `--outFilterMultimapNmax 1` | matches the historical posture; needed for like-for-like with existing packs |
| RNA | STAR | `--outFilterMultimapNmax 10` | the corrected input posture (mm10 recovers +70.6% coverage on short-insert libraries) |
| RNA | salmon | decoy-aware index, `-l A` | TPM / universe construction |
| Ribo | STAR | `--outFilterMultimapNmax 1` | **the only posture usable as a model TARGET** |
| Ribo | STAR | `--outFilterMultimapNmax 25` | diagnostic only: quantifies what mm1 discards |

### Non-negotiable per-assay flags

**Ribo-seq** (standing rule -- applies to every Ribo dataset on this cluster):
- `--alignEndsType EndToEnd`
- `--outFilterMultimapNmax 1` for the target arm. RiboCode does NOT drop multimappers itself, so
  multi-mappers silently inflate P-site counts and destroy periodicity.
- Drop rRNA/tRNA/miRNA/Mt_rRNA/Mt_tRNA loci before any caller sees the BAM.
- **Check raw FASTQs for 3' adapter before aligning.** Many submissions ship full TruSeq
  (`AGATCGGAAGAGC`) attached; without `cutadapt -a` you get 0% unique mapping and a ~5 KB BAM stub.
  Diagnose on the COMPLETE file (a streamed byte-range slice once reported 0% where the real file was
  93%): `zcat s.fq.gz | head -400000 | awk 'NR%4==2' | grep -c AGATCGGAAGAGC`; >50% means trim.
  Footprints also want `--minimum-length 20 --maximum-length 40`.
- The mm25 arm is **diagnostic only** -- never use it as a training target.

**RNA-seq:**
- No `--maximum-length` (inserts are not ~30 nt footprints).
- No `--discard-untrimmed` -- correct for RPF, wrong for RNA. Applying it to a CAR-T RNA arm once
  kept 16.4% of reads vs 99.0% without.
- **Verify the library is poly(A)-selected.** `library_selection=cDNA` does NOT prove it. Ribo-Zero
  total RNA retains tRNA/7SL, which are typed lncRNA and survive the biotype filter -- this crushed
  the GSE243134 liver universe to 3,321 tx. Read the GEO/ENA protocol text and record the answer.
- salmon: **decoy-aware index only** (`Salmon_indexes/salmon_index_decoy_*`). Transcriptome-only
  misattributes intronic/repeat reads to lncRNAs and biases TPM.

### Reference templates to copy, not re-derive

Per-dataset trimming is the part that breaks silently. Derive each new script from the closest
existing one, changing only the multimap flag and output paths:

**PAIRED-END templates already exist — do not write new ones.** The primate RNA is all PE (two
`fastq_bytes` per accession); the Ribo is SE.

| need | template | note |
|---|---|---|
| **RNA PE STAR, mm10** | `scripts/rebuild_coverage_mm10.sbatch` | already the mm10 posture; start here |
| **RNA PE STAR, mm1** | same file, flip `--outFilterMultimapNmax 10` -> `1` | |
| **RNA PE + salmon in one script** | `scripts/process_rnaseq_prjeb34766.sbatch` | STAR PE at `:73`, salmon PE at `:55` |
| **salmon PE/SE auto-branch** | `scripts/salmon_human_rna_polya_gate.sbatch` | PE at `:52`, SE at `:56` — cleanest branch pattern |
| salmon PE only | `scripts/salmon_chothani_decoy.sbatch:42` | decoy-aware, `-l A` |
| salmon SE only | `scripts/process_cart_gse304796.sbatch:159`, `scripts/janich_decontam_quant.sbatch:76` | likely unneeded: primate RNA is PE |
| RNA PE across multimap settings | `scripts/hep_rna_multimap_sweep.sbatch` | |
| **Ribo SE STAR** | `scripts/align_chothani_refetch.sbatch` | EndToEnd, `--trim-n -m 20 -M 40` |
| RNA SE, short-insert + poly-A trim | `scripts/heldout/align_janich_rna_mm10.sbatch` | SE only |
| RNA SE, standard TruSeq | `scripts/heldout/align_wang_rna_mm10.sbatch` | SE only |

New species need STAR + salmon indexes built first, into the SHARED trees
`STAR_indexes/star_index_<species>_<annver>/` and `Salmon_indexes/salmon_index_decoy_<annver>/`.
Never build per-project.

### Post-hoc transcriptome-BAM filtering: run the CROSS-GENE filter, not a best-score filter

`scripts/heldout/filter_tx_heldout.py <bam> <ncrna_tx.txt> <tx_to_gene.tsv>`

Run this on every Ribo transcriptome BAM (and on RNA BAMs if the universe is built from them).
It does two things, in order:

1. Drops records pointing at ncRNA transcripts (the rRNA/tRNA/miRNA/Mt drop required for Ribo).
2. Per read: **keep all records if they map within a single gene** (isoform multimapping is
   legitimate), **drop the whole read if it spans multiple genes** (cross-gene paralog ambiguity).

**Why that split, and why NOT a "keep the best alignment" filter.** These are two different
phenomena that both show up as a high `NH`:

- *Isoform multiplicity* -- `--quantMode TranscriptomeSAM` projects ONE genomic alignment onto every
  compatible transcript of the gene, so the copies have IDENTICAL alignment scores. "Best" cannot
  discriminate between them; picking one is a coin flip, not a criterion. It also barely matters:
  the profile head is a log-softmax over positions, so a uniform per-transcript inflation cancels,
  and `rnaseq_coverage.py` deliberately uses posture A "to match the P-site target" -- input and
  label are counted the same way, which is the property worth preserving.
- *Cross-gene paralog ambiguity* -- a repeat- or paralog-derived read landing on transcripts of
  DIFFERENT genes, donating signal to loci that do not have it. This is the real risk at mm10/mm25
  and is what the filter removes.

**Two operational constraints:**
- It **overwrites the BAM in place**. Archive the unfiltered BAM FIRST (intermediates are not
  deleted until the downstream result is verified).
- It relies on STAR's TranscriptomeSAM being **grouped by query name**, so it must run BEFORE any
  coordinate sort.

Not covered by this filter: choosing among DISTINCT genomic loci whose alignment scores actually
differ. That is genuinely unimplemented. Do not build it speculatively -- the mm1-vs-mm10
comparison now on disk for Janich is where evidence for whether it matters would come from.

### DO NOT bother with "best alignment only" flags

`--outSAMmultNmax 1`, `--outSAMprimaryFlag OneBestScore`, and both `--outMultimapperOrder` tie
orders are **no-ops on the transcriptome BAM**. Verified 2026-08-23: best-only output was
byte-identical to mm25 (BAM 11,094,927,983 / 11,094,938,192 vs 11,094,839,041 -- compression noise
only; GTF2I coverage identical at 1,797,673 across mm10/mm25/both tie orders).
Cause: `--outSAMmultNmax` governs the *genomic* BAM; `--quantMode TranscriptomeSAM` is a separate
path and emits all alignments up to `--outFilterMultimapNmax` regardless. One-alignment-per-read
must be done downstream by picking the best-scoring record per read.

---

## 3. Validation gates -- fail loudly, do not proceed on a bad BAM

1. `samtools quickcheck` every BAM you did not just write. A BAM truncated mid-write keeps a
   VALID HEADER, so `view -H | grep SO:coordinate` accepts it. This cost 92 min of RiboCode once.
2. Record STAR `Log.final.out` per run; check `Uniquely mapped reads %` and
   `% of reads mapped to multiple loci`. A large BAM does NOT mean usable reads -- a Janich library
   once mapped 0.03% because of untrimmed adapter.
3. After sorting, assert the record count is unchanged (`samtools view -c` before/after).
4. **NH in a transcriptome BAM is NOT genomic multimapping.** `--quantMode TranscriptomeSAM`
   expands one genomic alignment across every compatible isoform, so a genomically-unique read
   routinely carries `NH:i:17`. To compare mm1 vs mm25, compare TOTAL RECORDS PER TRANSCRIPT
   (`scripts/bamlib.py::multimap_cost`), never an NH split. An NH split once reported "99.8% of
   reads discarded" against STAR's own 30%.

---

## 4. Archive

Everything bulky gets archived, not deleted. Disk is the cheap resource; ENA served 0.30 MB/s for a
whole morning once.

```
scripts/archive_to_warm.sh [--dry-run] <local_path> "<why / how to restore>"
```
COPY -> **checksum-verify** -> index -> stub -> delete source. It deletes only after
`rsync --checksum --dry-run` reports zero differing paths. Leaves a row in `data/ARCHIVE_INDEX.tsv`
and a `<name>.ARCHIVED.md` stub with the literal restore command.

Archive: all FASTQs, all BAMs (both multimap postures), salmon quant dirs.
**Do NOT archive** anything an in-flight job is reading, or anything the next queued step needs --
check `squeue` first. Keep coverage hd5 and quant outputs local.

Quota: 12.76 TB / 15 TB (85.1%), 2.24 TB free. Overflow causes EDQUOT that **silently kills SLURM
jobs** (exit 120, no traceback, low MaxRSS). Check:
`getfattr -n ceph.dir.rbytes --only-values /private/groups/carpenterlab`

---

## 5. Compute placement

- **Downloads: head node, serial, single connection.** Compute nodes throttle one stream to
  ~0.4 MB/s and parallel array downloads trip ENA's per-IP limit.
- **Alignment: SLURM.** Never multi-core work on `mustard`. `medium` = 12 h cap, `long` = 14 d.
  GPU nodes often have hundreds of idle CPUs -- `--partition=gpu` with NO `--gres` starts fast.
- **Scripts need an explicit `--partition`**; the default `short` has a 1 h cap and will reject a
  12 h request.
- `/data/tmp` is **node-local**. Anything a SLURM job must read, or that must outlive the job, goes
  on the group fs. Staging into node-local scratch is only correct when the same job does it.
- Interpreters: `conda_envs/riboseq` has STAR/salmon/samtools/h5py/numpy but **no torch**;
  `conda_envs/cas12a` has torch but **no h5py**. `scripts/prepare/*` hardcode cas12a and will die on
  `import h5py` -- use riboseq.

---

## 6. Process hygiene

- **Kill by recorded PID**, never `pkill -f` / pattern match. Record `echo $! > /tmp/x.pid` at
  launch. Pattern-matching killed a monitoring process three times in one session (exit 144) because
  any watcher's command line contains the pattern too.
- **Never edit a shell script while it is running** -- bash reads by byte offset and will resume at
  a stale offset in new content.
- **One writer per path.** Two jobs sorting to the same `.sorted.bam` truncated a file mid-read and
  killed a 92-minute job.
- Read job state from `sacct -j <id> -o State`, not from whether an output file exists. Treat every
  terminal state (FAILED/TIMEOUT/CANCELLED/OUT_OF_MEMORY/NODE_FAIL), not just COMPLETED.

---

## 7. Definition of done

Per library: FASTQ md5-verified; STAR BAMs at both postures passing `quickcheck` with
`Log.final.out` retained; **unfiltered BAM archived, then `filter_tx_heldout.py` applied
pre-sort**; salmon quant dir (RNA only); adapter + poly(A) status recorded; bulky artifacts
archived with an `ARCHIVE_INDEX.tsv` row. Report any accession that exhausted its retries
and any library whose unique-mapping rate looks anomalous -- do not silently drop either.

## 8. Out of scope

Pack building, model training, evaluation. Also the zebrafish / Drosophila / yeast datasets --
their matched Ribo arms are not yet located (`docs/SPECIES_EXPANSION_PLAN.md` section D).
