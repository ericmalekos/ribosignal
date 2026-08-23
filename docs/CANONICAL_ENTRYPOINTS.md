# Canonical entry points -- check here BEFORE writing a new script

Status: POLICY. Adopted 2026-08-13 after `train_canon.sbatch` was written from scratch alongside a
working `train_loto_union.sbatch` and reintroduced two bugs that script had already solved.

**The rule: if a task appears in this table, you extend that script. You do not write a new one.**
Extending means adding an env var with a default that preserves current behaviour, so every existing
caller is unaffected. If you believe no entry fits, the new script's header must name the entry you
considered and state why it did not -- in one sentence, in the file.

`scripts/check_entrypoints.py` lists every `*.sbatch` not registered here. Run it before committing.

## Why a registry and not just "be careful"

Three failures in one day, all the same shape: the copy took the visible part, the divergence was in
the part that has a default.

| what happened | what was copied | what diverged | cost |
|---|---|---|---|
| Janich/Wang aligned off-recipe | the recipe from `align_gse243134.sbatch`, "for consistency" | missing `--alignEndsType EndToEnd` | full 40-alignment rebuild |
| leukocyte RNA glob | the `*.bam | grep -v bai` idiom | that directory was NOT RNA-only | caught 1 min in |
| `train_canon.sbatch` | the ARGUMENTS of `train_loto_union.sbatch` | the `--env` block | 3 dead jobs |

Environment is the recurring blind spot because it does not appear in `args.json`, so a
"hyperparameters match the deployed run" check cannot see it.

## The registry

| task | THE script | parameterized by | notes |
|---|---|---|---|
| Ribo-seq alignment (both species) | `scripts/riboseq_align.sbatch` | samplesheet TSV: `dataset, run, fastq, adapter, species, [umi], [cutadapt_extra]`; `SHEET`, `OUTROOT` | replaced 5 one-off `align_*`. Species selects index / ncRNA list / tx2gene ONLY |
| LOTO training, union universe | `scripts/train_loto_union.sbatch` | `BACKEND`, `MIXER`, `N_ATTN`, `D_STATE`, `SEED`; env `RIBO_PACK_SUFFIX`, `RIBO_ORF_TRACK`, `RIBO_ONEHOT_FASTA` | selects the SIF from `MIXER` (mamba needs orthrus for `mamba_ssm`) |
| LOTO training, final-recipe packs | `scripts/train_canon.sbatch` | `MODEL`, `ORFTRACK`, `TAG` | **should have been `train_loto_union.sbatch` + a `SUFFIX` var.** Retained only because runs 36769694/95/96 are live against it; fold back afterwards (see below) |
| pack assembly from per-nt hd5 | `scripts/prepare/prepare_pack.py` | `--ribo-psites`, `--rna-coverage`, `--universe-tx`/`--ref-pack`, `--group`, `--species`, `--out` | |
| pack assembly from BAMs | `scripts/prepare/prepare_from_bams.py` | as above + `--ribo-bam`, `--rna-bam`, `--annot`, `--ncrna-tx`, `--tx-to-gene` | filter args REQUIRED unless `--no-ncrna-filter` is declared |
| canonical Chothani training packs | `scripts/build_canonical_chothani_packs.sbatch` | array index 0-7 = tissue | reuses stored RNA coverage; asserts universe identity |
| rebuild a held-out pack canonically | `scripts/rebuild_heldout_canon.sbatch` | array index = arm | explicit run lists, Ribo/RNA disjointness assert |
| per-tissue regeneration divergence | `scripts/chothani_canonical_diverge.sbatch` | array index 0-7 = tissue | |
| ablated ORF track | `scripts/make_ablated_orf_track.py` | `--in`, `--zero`, `--out` | writes `_meta.json` recording `zeroed_channels` |
| ORF-call drop-in scoring | `scripts/compare_dropin_calls.build_loader` | imported, never reimplemented | reimplementing it once produced a table off by 0.14 F1 |
| channel-ablation scoring | `scripts/score_channel_ablation.py` | `--ablation-dir`, `--base-dir` | |
| 3x3 factorial scoring | `scripts/score_liver3x3.py` | `--results`, `--out` | |
| pooled multi-dataset ORF reference | `scripts/merged_liver_ribocode.sbatch` | BAM set is in-script (28 mouse-liver) | joint RiboCode over all libraries, NOT a union of call sets: one periodicity test on pooled signal, per-library offsets. Deliberately not folded into `ribocode_dropin.sbatch`, which calls per-dataset. Species/tissue is hardcoded because it exists to arbitrate ONE question (mouse-liver over-calls); generalize it only when a second pooled reference is actually needed |
| over-call validation vs a pooled reference | `scripts/score_merged_liver_overcall.py` | `--merged`, `--out`, `--tx2gene` | must report the experiment-unique CONTROL and the supported ANCHOR, never a bare recovery rate |

## Outstanding consolidation

- **`train_canon.sbatch` -> fold into `train_loto_union.sbatch`.** Add `SUFFIX=${SUFFIX:-union}` and
  `ORFTRACK=${ORFTRACK:-...}` there, and the final-recipe retrain becomes
  `BACKEND=onehot MIXER=mamba SUFFIX=canon sbatch scripts/train_loto_union.sbatch`. Blocked only
  until the live runs finish, so their provenance keeps pointing at the script that produced them.
- The 6 scripts carrying `SUPERSEDED` headers stay on disk for provenance; they are deliberately NOT
  in this registry, so `check_entrypoints.py` flags them and the flag is the point.

## What "extend" means in practice

```bash
# WRONG -- a new file that must re-derive every --env line from memory
cp train_loto_union.sbatch train_canon.sbatch && $EDITOR train_canon.sbatch

# RIGHT -- one new variable, default preserves every existing caller
SUFFIX=${SUFFIX:-union}
  --env RIBO_PACK_SUFFIX="$SUFFIX" \
```
