# SUPPLEMENTAL S_fdr_rigor -- FDR control for non-canonical identifications

**What it shows:** what the project's class-specific FDR policy is actually worth, measured on the
shipped searches. Closes FIGURES_PLAN rigor point 3.

- **(a)** novel peptides at global vs class-specific 1% FDR, every arm x dataset (paired).
- **(b)** the mechanism: hyperscore distributions of canonical targets, novel targets and novel
  decoys for the worst case, with both cut-offs drawn on.
- **(c)** what predicts the error: size of the canonical class, not size of the database.

---

## `make_fdr_rigor.py` -> `fdr_rigor.{png,pdf}` + `fdr_rigor_values.json`

### Data inputs

    proteogenomics/data/<DS>_pilot/pgx_mamba4/search/<ARM>/rank1.tsv.gz
    proteogenomics/data/<DS>_pilot/pgx_mamba4/db/db_<ARM>.fasta      (novel-sequence count only)

    DS  in {A549, HBL1, SUDHL4, DoHH2}
    ARM in {null_atg, cpat, cpc2, model_standard, model_poisson}

Released **mamba4** model, frozen search parameters. The peptide/decoy classification is imported
from `pgx.report` (`load_rank1`), so this figure and the results tables can never drift apart in how
they define canonical / novel / decoy. Search arms with no `rank1.tsv.gz` are **printed as MISSING**
and recorded in the JSON, never silently treated as zero.

### The two thresholds, defined precisely

Both use the same routine (`threshold`): walk PSMs from highest hyperscore downward, tracking the
running decoy/target ratio, and keep the lowest score at which that ratio is still <= 0.01.

- **global**: targets = `canon_t` + `novel_t`, decoys = `novel_d` + `other_d`. This is what the
  search itself reports.
- **class-specific**: targets = `novel_t` only, decoys = `novel_d` (`REV_nuORF|`) only. This is what
  every novel count in this study uses.

`novel_global` / `novel_class` are DISTINCT novel peptide sequences at or above the respective cut.

---

## What the panel actually shows (2026-08-07)

**A global 1% FDR fails to control the non-canonical class in BOTH directions**, and which way it
errs is set by how much of the search the canonical class occupies:

| dataset | digestion | canonical peptides | worst global/class ratio |
|---|---|--:|--:|
| A549 | tryptic, TMT | 65,405 | **37.5x over-report** (null_atg: 375 vs 10) |
| HBL-1 | nonspecific, HLA-I | 2,429 | 0.6x to 2.3x |
| SU-DHL-4 | nonspecific, HLA-I | 813 | 0.4x to 1.0x |
| DoHH2 | nonspecific, HLA-I | 1,507 | **0.2x under-report** (cpc2: 2 vs 8) |

This is a stronger claim than the usual "a global cut inflates non-canonical discovery", because it
cannot be waved away as merely conservative: in the immunopeptidomes the global cut is stricter than
the class-specific one and **discards real identifications**. Either way, the number a global cut
produces is not a 1% FDR for that class.

Panel (b) shows why, on A549 against the naive AUG database. The novel-target and novel-decoy
hyperscore distributions are nearly superimposed below ~20 -- that class simply has little signal.
Canonical targets carry a long high-scoring tail, and because they outnumber novel matches by orders
of magnitude they set the global cut at **19.5**, inside the novel decoy bulk. The class-specific cut
lands at **30.3**. Everything between those two values is the 37.5x.

Panel (c) rules out the tempting alternative explanation. If database size drove the error, the four
`null_atg` points (242k-369k novel sequences) would cluster together; they do not (37.5x, 1.9x, 0.8x,
2.2x). Plotted against canonical-class size instead, the points order cleanly.

### Figure-construction note

Panel (a) sets `ylim` to 9x the largest global count before drawing the legend. Without that
headroom the `upper left` legend sat directly on top of the A549 `null_atg` global point (375) --
the single marker the panel exists to show.

### Rebuild

    /private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 \
        figures/S_fdr_rigor/make_fdr_rigor.py

Deterministic; no sampling or fitting.
