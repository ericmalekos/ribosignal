
## Task #84 -- the deleted-annotation blast radius, scoped 2026-08-16

Found while checking a poster caption, not by looking for it. **Both** trees under
`biotype_probe/expression_context_human/data/` are gone:

- `ribocode_annot/`      (the RiboCode annotation: transcripts_cds.txt, transcripts.pickle, ...)
- `ribocode_per_tissue/` (per-tissue RiboCode CALLS, e.g. Hepatocytes/Hepatocytes_collapsed.txt)

**18 files reference them.** The previously logged `eval_union.sbatch` was one symptom, not the
problem. Enumerated by artifact rather than by memory, per
`feedback_blast_radius_scan_by_artifact`:

| references | files |
|---|---|
| `ribocode_per_tissue` (15) | `eval_localization.py`, `build_psite_target.py`, `verify_target_inputs.py`, `build_fibroblast_psite_target.py`, `plot_prediction_examples.py`, `replicate_concordance.py`, `ribocode_dropin_compare.sbatch`, `heldout/ribocode_heldout.sh`, `prepare/repack_union_{fib,others}.sbatch`, `prepare/validate_prepare_pack.sbatch`, `figures/S_riboseq_qc/{make_riboseq_qc.py,FIGURE_DATA_INPUTS.md}`, `figures/{prediction_examples,depth_crossover}/FIGURE_DATA_INPUTS.md` |
| `expression_context_human/data/ribocode_annot` (5) | `eval_canon_orfcalls.sbatch`, `heldout/dropin_sweep_{gse_poisson,hep_poisson,hep_pval}.sbatch`, `figures/prediction_examples/FIGURE_DATA_INPUTS.md` |

**The two halves are NOT equally recoverable, and that is the useful part:**

- **`ribocode_annot` is a repoint.** `data/human_ribocode_annot/` survives and carries
  `transcripts_cds.txt`, `transcripts.pickle`, `transcripts_sequence.fa`. The 4 sbatch files can
  be repointed at it. Verify the universe matches before trusting a result.
- **`ribocode_per_tissue` is GONE with no replacement.** Searched the whole experiments tree: no
  surviving `*_collapsed.txt` under any per-tissue path. These are per-tissue RiboCode calls and
  only re-running RiboCode restores them.

**Live consequence already observed:** `plot_prediction_examples.py --auto` raises
`FileNotFoundError`, so `figures/prediction_examples/*` cannot be regenerated at all. That is why
the AGGF1-vs-POLR1D discrepancy in its `FIGURE_DATA_INPUTS.md` could not be adjudicated by
re-running; the shipped PNG had to be taken as authoritative.

Not fixed. Repointing the 4 sbatch is cheap; deciding what to do about the 15 per-tissue
references needs a call on whether those calls get regenerated or the dependents retired.

---

## Task #93 RESOLVED 2026-08-18 -- why the data went, and the real root cause

**The loss is bigger than the scoping above recorded.** It is not two annotation directories. The
**entire `experiments/biotype_probe/` tree is gone**, including its own `methods.md`, `results.md`,
`scripts/`, `figures/` and all three sub-experiments (`short_orf_subset`, `embedding_geometry`,
`expression_context`). `RNAZoo/experiments/` now holds only `orthrus_probe` and
`riboseq_signal_model`.

### What removed it

Nothing in this project did, and that was checked rather than assumed:

- **No Claude session issued it.** All **130** transcripts under `~/.claude/projects/` were searched
  for an `rm` naming `biotype_probe` or `expression_context`. Zero matches.
- **No script could have.** 70 scripts contain `rm -rf "$VAR"`. Every one whose variable was traced
  resolves inside its own output tree (`$RUN/dropin_sweep_*`, `$SW/theta_*`,
  `$G/pred_*/dropin_poisson_sweep/*`). None can reach another experiment's directory.
- `experiments/` directory mtime is **2026-08-15 13:38:59**, which pins when the entry was removed.
- Not recoverable: `RNAZoo/.gitignore` line 52 is `/experiments/`, so it was **never in git**; the
  riboseq repo has no history of it; no tarball exists under `RNAZoo_meta/` or `carpenterlab/emalekos/`.

So it was removed externally, and most likely as legitimate housekeeping.

### The root cause is not the deletion. It is the cross-experiment dependency.

`riboseq_signal_model` was keeping **load-bearing inputs inside a different experiment's data
directory**. Deleting `biotype_probe` was somebody's reasonable call about their own experiment; it
broke this one because this one had reached across the boundary.

**This had already happened once and was already diagnosed here.** From
`scripts/build_human_ribocode_annot.sbatch`, written before the loss:

> "The mouse annotation lives at `data/mouse_ribocode_annot/` (it was relocated out of
> `expression_context_human` at some point; **17 scripts still pointed at the old path and were
> repointed**). This one is written alongside it, in-project, under the matching name -- **NOT back
> into another experiment's data dir, which is how the mouse one got lost track of in the first
> place.**"

The mouse annotation was rescued by that relocation. The human `ribocode_annot` and the whole of
`ribocode_per_tissue` were left behind and went with the tree.

### Standing consequence

**No `riboseq_signal_model` artifact may live outside `riboseq_signal_model/`.** Anything this
project needs to re-run gets written in-project. Task #97 therefore writes the regenerated
per-tissue calls into `riboseq_signal_model/data/`, never back into `biotype_probe`.

---

## Task #96 DONE 2026-08-18 -- annotation repoint

The count of "4 sbatch" was wrong: `eval_canon_orfcalls.sbatch` had **already** been fixed (its
remaining mention of the dead path is an explanatory comment on line 52; line 55 points at
`data/human_ribocode_annot_primary`). Three needed the change.

| file | was | now | why |
|---|---|---|---|
| `heldout/dropin_sweep_hep_poisson.sbatch` | dead `ribocode_annot` | `data/human_ribocode_annot_primary` | matches `eval_canon_orfcalls.sbatch`, documented there as the live annotation (GENCODE v49 primary, matching `star_index_grch38_v49`) |
| `heldout/dropin_sweep_hep_pval.sbatch` | dead `ribocode_annot` | `data/human_ribocode_annot_primary` | same |
| `heldout/dropin_sweep_gse_poisson.sbatch` | dead `ribocode_annot_mouse` | `data/mouse_ribocode_annot` | GSE120762 is **mouse**; precedent set by `liver3x3_*.sbatch` and `liver_released_dropin.sbatch` |

Verified rather than assumed: all three pass `bash -n`; both target annotations exist and carry
`transcripts.pickle`; the mouse annotation is confirmed mouse (66,618 coding tx, `ENSMUST` ids).

**Note for whoever uses these next.** Two human annotations survive and they are NOT identical:
`human_ribocode_annot` (from `RNAZoo_meta/annotations/gencode.v49.annotation.gtf`, 233,995 coding tx)
and `human_ribocode_annot_primary` (from `genomes/gencode.v49.primary_assembly.annotation.gtf`,
234,024 coding tx). The primary-assembly build has 29 MORE coding transcripts than the full
annotation, which is backwards from what the names imply and is not explained anywhere. It does not
affect this repoint, since both cover the transcript space either sweep needs, but it should be run
down before either annotation is used for a headline number.

---

## Future work: collapse redundant isoforms in the training set (noted 2026-08-23)

**Not yet actioned.** Recorded now so it is not rediscovered.

**The problem, measured.** Of 81,611 trainable transcripts over 17,368 genes (mean 4.70 isoforms
per gene, max 69), **72% of training nucleotides are near-duplicate isoform copies** -- 209 Mnt
total against 59 Mnt for longest-isoform-per-gene.

The per-transcript equal-weighted multinomial loss makes this worse than the nucleotide count
suggests, because a gene's gradient share scales with its isoform count:

| isoforms/gene | genes | % of training gradient |
|---|--:|--:|
| 1 | 4,585 | 5.6% |
| 2-3 | 4,632 | 13.8% |
| 4-7 | 4,954 | 31.7% |
| 8-15 | 2,595 | 32.7% |
| 16+ | **602** | **16.1%** |

**602 genes (3.5%) consume 16% of the gradient**, and half of it goes to genes with 8+ isoforms.
The model repeatedly sees near-identical sequence with near-identical targets -- a memorisation
pressure entirely separate from multimapping.

**The proposed collapse, with its constraint.** Keep full transcripts; keep EVERY distinct 5'UTR
(uORFs are a scored class and live there); 3'UTR diversity is not worth preserving (dORFs are the
weakest class, recall 0.098). So: **collapse transcripts that share a 5'UTR + CDS and differ only
in their 3' end**, retaining one representative per (5'UTR, CDS) group.

**Do NOT confuse this with posture B.** Posture B removed 66% of transcripts by EXPRESSION RANK and
cost -0.016 profile r and -0.043 to -0.068 count r on a common test set (task #100). This selects
by SEQUENCE REDUNDANCY and preserves the 5'UTR space posture B discarded. Different axis, and the
count head's sensitivity to losing low-expressed isoforms is the specific risk to watch.

**Any evaluation must score on a COMMON test set including the collapsed-away transcripts.**
Scoring only on retained transcripts manufactures a win -- exactly the artifact that made Task 21
report a gain that reversed on correction.
