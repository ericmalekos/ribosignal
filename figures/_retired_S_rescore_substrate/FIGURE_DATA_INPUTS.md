# SUPPLEMENTAL S_rescore_substrate -- semi-supervised rescoring is substrate-dependent

**What it shows:** why MS2Rescore helps the immunopeptidome searches and hurts the tryptic one, so
the pipeline's choice to report raw-hyperscore class-specific FDR everywhere reads as a decision
rather than an omission. Closes FIGURES_PLAN rigor point 4.

- **(a)** the controlling variable, measured on the current pgx searches: how much *trainable* signal
  the non-canonical class contains.
- **(b)** the consequence, from the recorded rescoring runs.

---

## `make_rescore_substrate.py` -> `rescore_substrate.{png,pdf}` + `rescore_substrate_values.json`

### Data inputs

**Panel (a)** -- computed live:

    proteogenomics/data/<DS>_pilot/pgx_mamba4/search/null_atg/rank1.tsv.gz
    DS in {A549, HBL1, SUDHL4, DoHH2}

The `null_atg` arm is used because it has the largest non-canonical space, which makes it the most
generous estimate of how much novel signal is available. Classification (`load_rank1`) and the
class-specific FDR cut (`threshold`) are imported from `pgx.report` and `S_fdr_rigor` respectively,
so there is exactly one definition of each in the repository.

**Panel (b)** -- hard-coded constants, deliberately. The rescoring experiments were run in the OLD
era (f0-threshold databases, pre-pgx, pre-frozen search parameters) and have not been repeated on
the released models. Recomputing them from today's search directories would silently mix eras, so
the numbers are transcribed from `proteogenomics/results.md` (sections "MS2Rescore on A549" and
"MS2Rescore on the HBL-1 immunopeptidome") and the panel is annotated as old-era.

### Panel (a)'s metric was chosen after the obvious one was measured and failed

The intuitive quantity is "what share of rank-1 PSMs touch a non-canonical ORF". It points the wrong
way:

| dataset | rank-1 PSMs | novel-touching | share | confident novel | **confident / 1,000 novel PSMs** |
|---|--:|--:|--:|--:|--:|
| A549 | 425,559 | 64,159 | 15.08% | 10 | **0.16** |
| HBL-1 | 21,460 | 912 | 4.25% | 20 | **21.93** |
| SU-DHL-4 | 5,925 | 278 | 4.69% | 28 | **100.72** |
| DoHH2 | 7,970 | 320 | 4.02% | 13 | **40.62** |

By raw share, A549 looks like the *most* non-canonical search in the panel, which would predict that
rescoring works best there. It does the opposite. The share is dominated by spurious low-scoring
matches into a 369k-sequence database: it measures noise, not signal.

What mokapot needs is confident examples to train on relative to the noise it must reject. On that
metric the substrates separate by 137x to 630x, in the direction the rescoring outcomes require.
Panel (a) plots that metric, and this paragraph exists so nobody re-derives the misleading one.

---

## What the panel shows (2026-08-07)

Rescoring is not a free improvement; whether it helps is set by the substrate.

- **A549, tryptic whole proteome (0.16 confident novel per 1,000 novel PSMs):** global rescoring
  optimises the canonical-dominated target/decoy separation and discards the raw hyperscore signal
  that was distinguishing the ~10 real novel matches. Class-specific FDR collapses 10 -> 1.
  Class-aware rescoring (the methodologically correct fix) cannot bootstrap at all: the SVM folds
  contain a single class at every `train_fdr` tried (0.05 / 0.1 / 0.25) for both the model and the
  null arm.
- **HBL-1, HLA-I immunopeptidome (21.93):** mokapot can learn the class. Rescoring cleans the naive
  null (18 -> 4, ms2pip demotes spurious matches in the large background) and recovers the model
  (10 -> 26). The model's advantage over the naive null widens roughly tenfold.

The three immunopeptidomes all sit in the learnable regime (21.9 / 100.7 / 40.6), so the HBL-1
outcome is expected to generalise across them, but only HBL-1 was actually rescored.

### Rebuild

    /private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 \
        figures/S_rescore_substrate/make_rescore_substrate.py

Deterministic. Panel (a) recomputes from disk; panel (b) is constant until the rescoring experiments
are repeated on the released models, at which point the `RESCORE` dict must be replaced by a live
read and this note removed.

## A549 is a KNOWN Rule 5 exception, frozen at 7 aa (ARCHIVED 2026-08-15)

A549 is **tryptic** (`trypsin`, `num_enzyme_termini = 2`), so Rule 5 puts its correct ORF floor at
**30 aa**. It sits at **7**, and it will not be rebuilt: the A549 work was archived on 2026-08-15
(`proteogenomics/data/_archive_a549_2026_08_15/README.md`).

It is the project's ONLY Rule 5 exception. The 12 macrophage populations are tryptic and now at
30 aa; HBL-1, SU-DHL-4, DoHH2, THP-1 and B721.221 are `nonspecific` / `termini = 0` HLA-I, where
7 aa is correct and permanent (HLA-I peptides are 8-11 residues).

**This does not invalidate the panel.** The tryptic-vs-HLA contrast is about digestion chemistry, not
ORF-length policy, and each substrate is internally consistent. What it does forbid is placing
A549's database size or discovery density beside an ORF-call panel and implying a single ORF set --
the ORF-call track filters at 90 nt = exactly 30 aa.

**A549 is deliberately NOT dropped**, because it carries the results least flattering to the model
(CPAT/CPC2 win the tryptic proteome; the one Poisson failure at -35; the largest global-vs-class FDR
discrepancy). Archiving the work is a decision to stop spending compute, not to delete the
counter-examples.
