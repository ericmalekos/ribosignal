# riboseq_signal_model

Predict per-nucleotide Ribo-seq P-site signal along a transcript from sequence (+ optional matched
RNA-seq coverage). A dilated-CNN body with a transformer or Mamba mixer feeds two heads: a profile head
(multinomial over the transcript, the periodic P-site SHAPE) and a count head (log total P-sites, the
DEPTH). Trained on Chothani matched primary human tissues (GENCODE v49 / GRCh38), transcriptome
coordinates.

## Status (2026-08-04)

Model TRAINED + FULLY EVALUATED. It predicts translation from SEQUENCE, so it needs no Ribo-seq at
inference.

### TWO models ship (decision 2026-07-31, Task 26)

| | `orf_v2_mamba4` | `orf_v2_attn` |
|---|---|---|
| role | **primary / headline** | supplemental, but a maintained release |
| mixer | dilated CNN + 4 Bi-Mamba blocks | dilated CNN + 2 transformer layers |
| params | 7,521,026 | 5,071,106 |
| held-out Pearson (3 seeds) | **0.6799** (0.6753-0.6851) | 0.6595 (0.6585-0.6603) |
| device | **GPU only** (mamba-ssm CUDA kernels) | **CPU or GPU** |

mamba4 wins by +0.0204 with non-overlapping seed ranges (its worst seed beats attn's best), so the
separation is real and not a seed artefact. It is ~5x noisier across seeds though (spread 0.0098 vs
0.0018), so never quote a single-seed mamba4 number without the spread.

**attn is not deprecated.** It is the only one of the two that runs without a GPU, so it stays the
released inference path for CPU-bound users. Keep both sets of docs, figures, and checkpoints current.

Encoding: one-hot. FM embeddings (RiNALMo / Orthrus / HydraRNA) give no lift over one-hot (Task 15),
so the released models carry no foundation-model dependency.

### Landed work

Tasks 1-26 plus 46/52/53 (see `results.md`): held-out-chromosome baselines, FM + input/architecture
ablations, uORF/dORF-aware eval, 9-fold leave-one-tissue-out transfer, ORF-localization + RiboCode
drop-in calling, cross-study (human Ruiz-Orera) + cross-species (mouse Wang) transfer, depth crossover
(predicting beats measuring below ~29M P-sites), non-AUG ATG+CTG drop-in, the Kozak ablation, the
posture-B multimap sensitivity check, the 3-seed architecture decision, and the macrophage
proteogenomics application across 12 populations.

Key decision from Task 20: the hand-picked Kozak start-context factor is redundant with what the
sequence backbone learns, so the ORF-track default is `--kozak none`. Both shipped models are trained
on the no-Kozak track; only the OLD deployed checkpoint retains the heuristic (train/inference match).

## Documents (read in this order)

- **`results.md`** -- findings, Tasks 1-21, headline numbers and tables. Start here.
- **`methods.md`** -- what was done and why: target definition, inputs, splits, model, per-task
  methodology (sections 1-7 + subsections 5B..5N).
- **`KOZAK_PLAN.md`**, **`LOTO_PLAN.md`** -- plan + outcome for the Kozak ablation and the leave-one-
  tissue-out pilot.
- **`design_count_magnitude_transferability.md`** -- design rationale for the dual-head + transferability.
- **`proteogenomics/`** -- the MS application: its own `methods.md` + `results.md`, the `pgx` pipeline,
  and `DATA_PROVENANCE.md` (log every external dataset there before use).
- **`HANDOFF.md`** -- HISTORICAL. The 2026-07-08 kickoff brief, written before any code existed. Kept
  for provenance; it describes the project as not yet started and is not a guide to the current state.
- Auto-memory `project_riboseq_signal_model.md` (in the Claude memory dir) -- the running cross-session
  STATE block; read its tail first when resuming.

## Layout

- `scripts/` -- all build / train / eval code. Training entrypoints: `train.py` (fold-based held-out
  chromosome) and `train_loto.py` (leave-one-tissue-out). Model in `model.py`, data in `dataset.py`.
- `data/` -- `packed/` (Fibroblast target + coverage + ORF tracks + tx_order), `packed_<Tissue>/` (other
  tissues), `target/`, universe FASTA/TSV, external held-out sets under `data/external/`.
- `results/` -- per-run output dirs (`args.json`, `best.pt`, metrics JSON, `pertx.tsv`, `dropin*/`).
- `figures/` -- curated figures, each folder with a `FIGURE_DATA_INPUTS.md` provenance file.
- `logs/` -- captured train / eval / dropin / heldout / input-build logs.
- `manuscript/`, `tutorial/` -- writeup scaffolding.
- `release/` -- **which checkpoints actually ship.** Config, metrics, checksums and companion-artifact
  paths for the two released models. Read this before attributing any result to "the model".
- `env/` -- pinned conda specs + the Singularity images. NOTE: training runs in CONTAINERS, and attn
  and mamba4 use DIFFERENT images (only one carries the Mamba CUDA kernels).

## Run a released model (predict from sequence)

Inference uses `dump_pred_profiles.py` (per-nt predicted profile per transcript) or `eval_localization.py`
(ORF-level scoring). A trained run dir carries its full config in `args.json`; point the eval scripts at
`--run <results/.../run_dir>`. GPU via the `rnazoo-rinalmo` SIF, or CPU with the `--device cpu` variants.
See methods.md 5C (training protocol) and 5G/5H (eval + drop-in) for exact commands.

Pick by device: `orf_v2_mamba4` if a GPU is available (better, and the headline numbers), `orf_v2_attn`
if not. `--device cpu` on a mamba4 run dir will fail in the mamba-ssm CUDA kernels -- that is a hard
constraint of the dependency, not a configuration problem.

## Conventions

SLURM cluster `prism`: heavy work via sbatch (partitions short/medium/long/gpu), head node `mustard` for
light inspection + serial downloads only; outputs on the group filesystem, scratch under
`/data/tmp/emalekos`. Mitochondrial (chrM) genes are excluded everywhere (methods.md 7.1). Ribo-seq STAR
alignment keeps single-mappers only and drops rRNA/tRNA/miRNA loci. ASCII punctuation only.
