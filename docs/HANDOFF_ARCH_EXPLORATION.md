# Handoff: architecture exploration (CNN-free mixer stacks)

Self-contained. Scope: **train and evaluate alternative architectures for the per-nt Ribo-seq
signal model.** No data ingest, no pack building, no held-out rebuilds -- those are owned elsewhere
(`HANDOFF_INGEST_ALIGN_2026_08_24.md`). Do not touch `data/packed*`.

---

## 0. Orientation -- paths, environments, scripts

**Project root (all relative paths below are from here):**
```
/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model
```
Shared references (NOT under the project root -- one canonical copy, reused across projects):
```
/private/groups/carpenterlab/emalekos/genomes/          genome FASTA + GTF
/private/groups/carpenterlab/emalekos/STAR_indexes/     star_index_grch38_v49, star_index_grcm39_vM38
/private/groups/carpenterlab/emalekos/Salmon_indexes/   salmon_index_decoy_v49, salmon_index_decoy_vM38
/private/groups/carpenterlab/emalekos/conda_envs/       interpreters (below)
/private/groups/carpenterlab/emalekos/singularity_cache/  SIFs (below)
/private/warm-archive/carpenterlab/RNAZoo_meta/         warm archive, mirrors the project tree
```

**Environments -- verified contents. There is NO single env that does everything.**

| env / image | has | LACKS | use for |
|---|---|---|---|
| `conda_envs/riboseq` | STAR, salmon, samtools, cutadapt, h5py, pysam, numpy | **torch** | alignment, quant, coverage hd5, BAM work |
| `conda_envs/ribocode` | RiboCode, h5py, pysam, numpy | torch, STAR | `ribocode_dropin.py`, ORF calling |
| `conda_envs/cas12a` | torch, numpy | **h5py** | torch-only utilities |
| SIF `...-orthrus-latest.img` | torch + **mamba_ssm** | -- | `--mixer mamba` train/predict (GPU for training; CPU predict via `RIBO_MAMBA_CPU=1`) |
| SIF `...-rinalmo-latest.img` | torch | mamba_ssm | transformer-mixer train/predict; runs on CPU |

Gotchas that have already cost time:
- `scripts/prepare/*` hardcode `cas12a`, which has no h5py -- they die instantly on `import packlib`.
  Use `riboseq`.
- Always `singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1`. Without
  `--no-home --cleanenv` the host `~/.local` site-packages shadow the SIF's own.
- mamba_ssm has **no CPU fallback** (`causal_conv1d_fwd` asserts `x.is_cuda`), so even a 35-tx smoke
  test needs a GPU.

**Reusable scripts -- prefer these over hand-rolling; each encodes a rule learned the hard way.**

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

## 1. The question

The deployed model is a 10-block dilated CNN (receptive field 4,093 nt) followed by 2-4 global
mixer layers. **Can the CNN be removed entirely, leaving a pure Mamba (or transformer) stack?**

The CNN is **78% of parameters** (3.93M of 5.07M -- each `ResidualDilatedBlock` holds two
`Conv1d(256,256,3)`), so removing it frees real capacity for the mixer.

---

## 2. What is already known -- read before designing anything

### 2a. Mixer depth SATURATES at 4, and attention gets worse beyond it

Union universe, one-hot, nokozak, mm1, holdout Hepatocytes, seed 0. All on disk under
`results/loto/`.

| arch | mixer | layers | test Pearson | val | frame0 | period |
|---|---|---|---|---|---|---|
| mamba4 | mamba | 4 | **0.6851** | 0.6190 | 0.8196 | 0.2632 |
| mamba6 | mamba | 6 | 0.6843 | 0.6251 | 0.8217 | 0.2564 |
| mamba_ds32 | mamba | 2 (d_state 32) | 0.6725 | 0.5982 | 0.8276 | 0.2845 |
| mamba | mamba | 2 | 0.6711 | 0.6020 | 0.8234 | 0.2842 |
| attn | transformer | 2 | 0.6585 | 0.6099 | -- | 0.3145 |
| attn4 | transformer | 4 | 0.6575 | 0.6095 | 0.8251 | 0.3043 |
| attn6 | transformer | 6 | -- (stopped e10) | **0.5808** | -- | -- |
| mamba_ds64 | mamba | 2 (d_state 64) | -- (stopped e17) | 0.5904 | -- | -- |

Conclusions, already paid for: **4 mixer layers is the plateau** (4→6 gains nothing);
**Mamba > transformer at every matched depth** (+0.028 test at 4 layers); **larger d_state does not
help** (32 and 64 both worse than 16).

### 2b. Cost: an attention layer is ~3,840x a conv layer on this data

Union: 84,472 tx, mean 2,517 nt, p50 2,070, p90 5,069, p99 8,435, max 9,996.

| quantity | value |
|---|---|
| Σ L -- cost of ONE conv layer, O(L) | 212.6 M |
| Σ L² -- cost of ONE attention layer, O(L²) | **816.2 G** |
| longest 1% of tx share of attention cost | 9% |

Attention already dominates ~99.9% of mixer FLOPs. **Mamba is O(L)** -- same per-layer cost as a
conv -- which is why the Mamba variant is the tractable one.

Memory, one materialised 8-head attention matrix: **3.20 GB fp32 at L=10,000** (1.60 GB fp16).
A5500 is 24 GB. Note `forward()` calls `self.attn(xt, src_key_padding_mask=~mask)`, and in several
PyTorch versions that mask forces the non-flash math path -- verify before assuming O(L) memory.

### 2c. THE critical fact: there is no positional encoding anywhere

Grep confirms it: no PE, no RoPE, no ALiBi in `scripts/model.py`. **The convolutions supply order.**

Consequence: self-attention is permutation-equivariant, so a CNN-free *transformer*
**cannot represent reading frame (position mod 3) at all** -- and 3-nt periodicity is the signal.
A CNN-free transformer therefore REQUIRES adding a positional encoding, and one under which
`position mod 3` is recoverable. Mamba is inherently ordered and needs none.

### 2d. "Pure Mamba" is not actually CNN-free

Every Mamba block contains an internal depthwise `conv1d` with `mamba_d_conv=4`. So the real
experiment is *"replace an explicit multi-scale dilated CNN (RF 4,093 nt) with many short 4-nt convs
plus SSM state"* -- not *"remove convolution"*. Say it that way in any write-up.

### 2e. The long-range job is already precomputed

The ORF-candidate track (input channels) marks occupancy from each ATG through its in-frame stop.
The transcript-scale structure a mixer would have to learn is handed to the model as input. This is
the leading explanation for why extra mixer layers buy nothing, and it predicts a CNN-free model
will not gain much from added global context either.

---

## 3. The experiment

`n_blocks=0` is **already expressible and verified to work**: `dilations` becomes an empty tuple,
`self.blocks` an empty `ModuleList`, and the `for blk in self.blocks` loop is a no-op. `in_proj` is
`Conv1d(..., 1)` -- kernel size 1, pure per-position projection, no spatial mixing. So
`--n_blocks 0` gives a genuinely conv-free body.

**Primary arm (run this first):**
```
--n_blocks 0 --mixer mamba --n_attn_layers 12
```
Sweep `n_attn_layers` in {8, 12, 16}. Keep `--channels 256 --mamba_d_state 16` (d_state 32/64 are
already known losses). 12 Mamba layers ≈ the parameter budget freed by dropping the CNN.

**Secondary arm (only if the primary is interesting):** CNN-free transformer. Requires implementing
RoPE or ALiBi first -- see 2c. Expect it to be slow and to need fp16 + gradient checkpointing.
Given attn6 already regressed with the CNN still present, this is low prior.

**Intermediate arm worth including:** `--n_blocks 4` (RF ~120 nt) with 8-12 Mamba layers. If the
CNN's value is purely local codon/frame detection, a short CNN plus a deep SSM may beat both ends.
This is the most likely place to find a real win.

### Everything else must be held fixed for comparability
```
--emb_backend onehot --cov_norm global_mean --min_train_signal 50 --val_fold 0
--max_tx_per_tissue 6000 --epochs 30 --patience 6 --warmup 2 --dropout 0.1
--lr 3e-4 --weight_decay 1e-2 --budget 16000 --count_weight 0.1 --eval_cap 2500
--holdout Hepatocytes --train_tissues Fibroblast,VSMC,ES,Fat,HA_EC,HCAEC,HUVEC --seed 0
```
Env: `RIBO_PACK_SUFFIX=union`, `RIBO_ORF_TRACK=data/packed_union/orf_track_v2.npy`
(kozak:"none" -- verified), `RIBO_ONEHOT_FASTA=data/union_universe.fa`.
Template to copy: `scripts/train_loto_union.sbatch` (`MIXER=mamba N_ATTN=4`).

**`--kozak none` always** -- the heuristic gate is redundant/harmful, standing rule.

---

## 4. How to judge the result -- fixed BEFORE running

Test split only (n=70,883), identical across all arms. Seed-noise floor measured from the three
baseline seeds:

| metric | baseline mean | seed range | sd |
|---|---|---|---|
| test pearson_median | 0.6799 | **0.0098** | 0.0049 |
| test frame0_pred_median | 0.8199 | **0.0013** | 0.0007 |
| test period_pred_median | 0.2817 | 0.0345 | 0.0174 |

A pearson difference under ~0.010 is noise. **`frame0_pred_median` is the sensitive discriminator**
(spread an order of magnitude tighter). `period_pred_median` is noisiest -- never let it carry a
conclusion alone.

### Pre-registered failure criterion

**Periodicity is where a CNN-free model should break first.** The dilated CNN with k=3, dilation 1
is the component most responsible for codon/frame structure, and the sweep already shows mixers
trading periodicity for precision (attn period 0.3145 vs mamba4 0.2632 while mamba4 wins overall).
Register up front: *if `--n_blocks 0` loses frame0 beyond the 0.0013 floor while pearson is flat,
the CNN's contribution is frame structure and the answer is a short CNN, not none.*

Comparison script (test-only, judged against the floor):
`scripts/compare_mm10cov_arms.py` -- adapt its run list.

**Any promising arm must then pass the ORF-call evaluation** (`scripts/eval_union.sbatch`, which
takes `PSUF`/`RUN`/`TAG`/`LABEL` overrides): RiboCode drop-in on real / pred_obsdepth /
pred_preddepth **plus** the CDS-anchored Poisson sweep. Aggregate Pearson has already proven a poor
proxy -- arm B of the mm10 experiment improved frame0 (+0.0054) and periodicity (+0.0448) yet lost
0.072 non-canonical precision.

**When comparing Poisson arms across models, the ACHIEVED `cds_recall` must match, not just the
target.** Different operating points invert conclusions; this has already happened twice.

---

## 5. Infrastructure

- **A5500, 24 GB.** `--partition=gpu --gres=gpu:A5500:1`. Scripts need an **explicit
  `--partition`** -- the default `short` caps at 1 h and rejects longer requests.
- **mamba_ssm's FAST path is CUDA-only** (`causal_conv1d_fwd` asserts `x.is_cuda`), but there IS a
  CPU fallback: set **`RIBO_MAMBA_CPU=1`** and `dump_pred_profiles.py` swaps the CUDA kernels for
  mamba_ssm's own pure-PyTorch reference twins (`selective_scan_ref`, `mamba_inner_ref`,
  `causal_conv1d_ref`). Verified finite on the real model (mamba4, d_state 16) at L=1000 and 3000;
  ~7.9 s per 3,000 nt transcript single-threaded, so shard wide (`--nshards`). The subtlety if you
  re-implement it: those names were bound at import time in SEVERAL modules, so patching one is not
  enough -- `mamba_inner_ref` calls the `causal_conv1d_fn` inside `selective_scan_interface`. See
  `_enable_mamba_cpu()`. **GPUs here run 48/48 allocated for long stretches**, so a CPU smoke test
  often finishes before a queued GPU job starts. TRAINING still wants a GPU. SIF: `singularity_cache/ghcr.io-ericmalekos-rnazoo-orthrus-latest.img`.
  Transformer arms can use the rinalmo SIF and run on CPU.
- Always `singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1` -- without
  `--no-home --cleanenv` the host `~/.local` site-packages shadow the SIF's own.
- Baseline epoch time ~1,900 s at 42,000 train tx; runs stop by e20-e26 (patience 6). Budget
  ~12-14 h per arm; `--time=24:00:00` on `gpu`. `--resume` reads `<out>/last.pt`.
- **Gate any eval on `test_metrics.json`, not `best.pt`** -- `best.pt` appears after epoch 0, and
  gating on it once submitted an eval of a 1-epoch model.
- Interpreters: `conda_envs/riboseq` = h5py/numpy, **no torch**; `conda_envs/cas12a` = torch,
  **no h5py**. Anything touching `model.py` runs in a SIF.
- Quota 12.79/15 TB (85.2%), 2.21 TB free. Checkpoints are small (~20 MB), so training is not a
  disk risk -- but EDQUOT silently kills SLURM jobs (exit 120, no traceback).

---

## 6. Do not

- Re-run depth sweeps at 2/4/6 mixer layers with the CNN present -- done, table in 2a.
- Try `mamba_d_state` 32 or 64 -- both already lost.
- Build a CNN-free transformer without adding positional encoding first; it cannot encode frame.
- Set the Kozak heuristic.
- Touch `data/packed*` or any held-out pack.
- Report a Poisson-arm comparison without both models' achieved `cds_recall`.

## 7. Deliverable

Per arm: run dir under `results/loto/`, `args.json`, `history.json`, `test_metrics.json`, and a row
in a comparison table with all three test metrics against the seed-noise floor. Plus a one-paragraph
verdict on the pre-registered criterion in section 4: **did removing the CNN cost frame structure?**
A clean negative is a publishable result -- the point is to settle whether the dilated CNN is load-
bearing, not to find a winner.
