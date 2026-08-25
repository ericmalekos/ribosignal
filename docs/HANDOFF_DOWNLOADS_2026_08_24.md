# Handoff: cross-species downloads (2026-08-24)

Scope of this handoff: **downloading FASTQs and reference genomes only.** No alignment, no pack
building, no analysis. Everything below is verified, not inferred -- where something is unverified
it says so explicitly.

---

## 0. DO NOT DUPLICATE: what is already running or already on disk

**Two download processes are live right now.** Check before starting anything.

| what | pid file | destination | check |
|---|---|---|---|
| Ruiz-Orera human iPSC-CM RNA (5 runs, 18.8 GB) | `/tmp/fetch_ruiz_pid.txt` | `data/external/heldout_rna_mm10/ruizorera_hsCM_rna/` | `kill -0 $(cat /tmp/fetch_ruiz_pid.txt)` |
| Reference genomes x6 (chained to start when the above finishes) | `/tmp/fetch_refs_pid.txt` | `genomes/xspecies_refs/<label>/` | `tail logs/fetch_refs_driver.log` |

The reference chain fetches, in order: gorilla, chimp, macaque, celegans, zebrafish, yeast.
Wait for `ALL_REFS_DONE` in `logs/fetch_refs_driver.log`.

**Already local, do NOT re-download:**
- Drosophila genome + GTF: `annotations/fly/dmel-all-chromosome-r6.67.fasta` + `dmel-all-r6.67.gtf`
  (FlyBase r6.67 -- NEWER than Ensembl BDGP6.54 and NCBI's FlyBase 6.54 mirror)
- Wang mouse liver RNA: `data/external/wang2021_mouse/fastq/mouse_liver/SRR526287{4,5}.fastq.gz`
  (md5-verified against ENA)
- Janich mouse liver RNA: `data/external/janich_liver_rnaseq/` (7 runs)
- PRJEB34766 liver RNA (7 ERR), GSE208041 THP-1, GSE304796 CAR-T, GSE39561 THP-1: all under
  `data/external/`

**BEFORE fetching anything, search the local tree at FULL DEPTH.** On 2026-08-24 a 3.4 GB Wang
re-download happened because the availability survey used `find -maxdepth 1` and missed the files at
depth 3. Use `find <root> -name "<ACC>*"` with no depth cap.

---

## 1. The job: PRJEB65856 non-human primates

Ribo-seq + RNA-seq from primate hearts and iPSC-CMs. The human arm is done/in flight; these three
species are the ask. Full run table saved at `data/external/PRJEB65856_runs.tsv`.

| species | assay | runs | size | accessions |
|---|---|---|---|---|
| gorilla (gg) | RNA-Seq | 6 | 38.0 GB | ERR12549902 ERR12549903 ERR12549904 ERR12549914 ERR12549915 ERR12549916 |
| gorilla (gg) | Ribo-Seq | 3 | 18.1 GB | ERR12549937 ERR12549938 ERR12549939 |
| chimp (pt) | RNA-Seq | 11 | 58.9 GB | ERR12549882 ERR12549883 ERR12549884 ERR12549885 ERR12549886 ERR12549896 ERR12549897 ERR12549898 ERR12549908 ERR12549909 ERR12549910 |
| chimp (pt) | Ribo-Seq | 8 | 29.2 GB | ERR12549917 ERR12549918 ERR12549919 ERR12549920 ERR12549921 ERR12549931 ERR12549932 ERR12549933 |
| macaque (rm) | RNA-Seq | 10 | 49.3 GB | ERR12549887 ERR12549888 ERR12549889 ERR12549890 ERR12549899 ERR12549900 ERR12549901 ERR12549911 ERR12549912 ERR12549913 |
| macaque (rm) | Ribo-Seq | 7 | 27.7 GB | ERR12549922 ERR12549923 ERR12549924 ERR12549925 ERR12549934 ERR12549935 ERR12549936 |

**Total: 45 runs, 221 GB.** At the observed 1.5-2.7 MB/s that is roughly 23-41 hours, serial.
Suggested order: gorilla (smallest, 56 GB) -> macaque (77 GB) -> chimp (88 GB).

RNA runs are paired-end (two files per accession); Ribo runs are single-end.

---

## 2. How to download -- use the existing script, do not hand-roll

```
bash scripts/fetch_heldout_rna.sh <label> <ACC> [ACC ...]
```

Writes to `data/external/heldout_rna_mm10/<label>/`, logs to `logs/fetch_<label>.log`.
It already implements every rule below. Suggested labels: `primate_gg_rna`, `primate_gg_ribo`, etc.

Launch pattern (record the PID -- see the kill rule):
```
nohup bash -c "
  bash scripts/fetch_heldout_rna.sh primate_gg_rna  ERR12549902 ERR12549903 ...
  bash scripts/fetch_heldout_rna.sh primate_gg_ribo ERR12549937 ERR12549938 ERR12549939
  ...
" > logs/fetch_primates_driver.log 2>&1 &
echo $! > /tmp/fetch_primates_pid.txt
```

### The rules the script encodes, and why each exists

- **HEAD NODE (`mustard`), SERIAL, ONE CONNECTION.** Compute nodes throttle a single stream to
  ~0.4 MB/s and parallel sbatch-array downloads trip ENA's per-IP limit (it starts refusing
  connections). Downloads are I/O-bound so they are a light head-node task even when long.
- **`wget -nv -c`, never `aria2c -x4`.** Multi-connection silently corrupts large files, producing
  the CORRECT SIZE with a bad md5.
- **`-nv` not `-q`.** `-q` suppresses wget's ERRORS as well as progress. On 2026-08-24 a transfer
  died at 13 MB of 1,594 MB and left an empty log because of `-q`.
- **VERIFY BY md5 AGAINST ENA, NEVER BY SIZE.** `fastq_md5` from the filereport API. A raced or
  truncated FASTQ routinely lands at exactly the expected length.
- **On md5 failure, DELETE and refetch; never `wget -c` onto bad bytes** -- resume keeps them and
  appends. The script retries 3x and exits nonzero with a failure summary.
- Two transient truncations occurred on 2026-08-24 and both recovered on retry, so a single failure
  is not evidence of a bad mirror copy.

---

## 3. Traps that cost time today -- please do not re-derive these

1. **`library_strategy` is unreliable for identifying Ribo vs RNA.** For PRJEB65856 it IS correct
   (`Ribo-Seq` / `RNA-Seq`, used in the table above). But for worm/fish/fly/yeast, ZERO runs are
   labelled `Ribo-Seq` -- that vocabulary term postdates their submissions, so ribosome profiling
   hides under `RNA-Seq` or `OTHER` and must be identified from **sample_title**.
2. **Never infer accessions from BAM filenames.** ERR12549926-930 were assumed to be the
   Ruiz-Orera RNA runs because they were the BAMs on disk; they are the **Ribo** runs. 3.8 GB was
   fetched before the mistake was caught. Always confirm with the ENA `library_strategy` field
   before fetching.
3. **Do not edit a shell script while it is running.** Bash reads scripts by byte offset; rewriting
   the file mid-execution makes it resume at a stale offset in new content.
4. **Kill by RECORDED PID, never by `pkill -f` / pattern-matching.** A pattern matching a directory
   name killed a monitoring process (exit 144) three separate times in one session, because any
   watcher's command line contains the pattern too. `echo $! > /tmp/x.pid` at launch, then
   `kill "$(cat /tmp/x.pid)"`.

---

## 4. Where things go, and the disk situation

- FASTQs: `data/external/heldout_rna_mm10/<label>/` (project fs, NOT `/data/tmp` -- that is
  node-local and invisible to SLURM jobs)
- References: `genomes/xspecies_refs/<label>/{genome.fna,genomic.gtf,PROVENANCE.json}` -- the SHARED
  tree, one canonical copy per (species, annotation version), reused across projects
- **Quota: 12.76 TB of 15 TB used (85.1%), 2.24 TB free.** 221 GB of primate FASTQ is fine. The
  group ceph quota is shared and hitting it causes EDQUOT that SILENTLY kills SLURM jobs (exit 120,
  no traceback). Check with:
  `getfattr -n ceph.dir.rbytes --only-values /private/groups/carpenterlab`
- Archive anything bulky when done rather than deleting:
  `scripts/archive_to_warm.sh <path> "<why / how to restore>"` (copy -> checksum-verify -> index ->
  stub -> delete). Check `data/ARCHIVE_INDEX.tsv` and `*.ARCHIVED.md` stubs BEFORE any re-download.

---

## 5. Definition of done

- Every accession in section 1 present under its label dir with `[md5 OK]` in the log.
- `logs/fetch_<label>.log` ends with `all md5 verified` (the script exits nonzero otherwise).
- `logs/fetch_refs_driver.log` contains `ALL_REFS_DONE`, and each
  `genomes/xspecies_refs/<label>/` holds `genome.fna`, `genomic.gtf`, `PROVENANCE.json`.
- Report any accession that exhausted its 3 retries -- do not silently drop it.

## 6. Explicitly NOT in scope

STAR/Salmon index building, alignment, pack building, evaluation. Also **not** the zebrafish /
Drosophila / yeast datasets: their matched Ribo arms are not yet located, and
`docs/SPECIES_EXPANSION_PLAN.md` section D lists what must be resolved first. Do not queue those.
