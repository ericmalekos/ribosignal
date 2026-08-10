#!/usr/bin/env python3
"""Figure B4: how much of the Ribo-seq profile is even PREDICTABLE, and how much the model captures.

Split 32 Fibroblast samples into two 16-sample half-pools, take the per-transcript Pearson between
them, and Spearman-Brown-correct to full depth. That is the reproducible CEILING -- the score a
perfect model would get, given that the assay does not perfectly reproduce itself. Bars are ceiling
vs the released model's achieved per-transcript Pearson, over four matched strata.

TWO BUGS FIXED 2026-08-08, both of the same kind (a number that could not be recomputed from a file):

1. The third bar was `pc_cds_codon`: a per-codon CDS ceiling of 0.9501 next to a model value of
   0.526. Those are different quantities. `scripts/replicate_concordance.py` records the model side
   of that stratum as `None` with the comment "our per-codon TBD" -- the model's per-codon CDS Pearson
   was never computed, so 0.526 was a periodicity score standing in for an elongation-shape score, and
   the "55% of ceiling" it produced was meaningless. That bar is replaced by the two strata the
   concordance script explicitly built to match the eval (`pc_uorf5` / `pc_dorf3` vs eval_extra's
   `uorf_5utr` / `dorf_3utr`). The per-codon CDS comparison can come back if and when the model side
   is actually computed.

2. MODEL was a hardcoded dict from the pre-union eval (pc 0.638, lncRNA 0.364). It is now READ from
   the released run's `extra_metrics.json`, and the run name is printed and written into
   `B4_values.json`, so a stale figure shows up as a stale run name rather than as a plausible number.

cas12a env. Regenerate after any retrain: `python3 make_replicate_ceiling.py`.

MODEL (2026-08-08). Defaults to **mamba4**, per locked decision D1b (2026-07-31): main Figure 1 is
`orf_v2_mamba4`, with `orf_v2_attn` moving to supplemental but staying shipped as the CPU inference
path. This panel had been on attn -- D1b was applied to the plan but never to the generators. Switch
with `FIG_MODEL=attn` for the supplemental variant; the chosen model is written into the values JSON.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/B4_replicate_ceiling"
REP = NEW / "results/replicate_concordance.json"
MODEL = os.environ.get("FIG_MODEL", "mamba4")
if MODEL not in ("mamba4", "attn"):
    raise SystemExit(f"FIG_MODEL={MODEL!r}; want 'mamba4' (main, per D1b) or 'attn' (supplemental)")
RUN = f"orf_v2_{MODEL}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes"
EXTRA = NEW / "results/loto" / RUN / "extra_metrics.json"
# The MAIN figure (mamba4, per D1b) keeps the plain filename; the supplemental attn variant gets
# an "_attn" suffix. Without this the two share one path and whichever ran last silently wins --
# which happened once already, with the attn build overwriting the mamba4 PDF.
SUF = "" if MODEL == "mamba4" else f"_{MODEL}"


# (ceiling key in replicate_concordance.json, model key in extra_metrics.json, display label).
# The pairing is not incidental: replicate_concordance.py builds pc_uorf5 / pc_dorf3 with the comment
# "matching eval_extra uorf_5utr / dorf_3utr", i.e. the same window over the same transcripts.
STRATA = [
    ("pc_whole",      "protein_coding", "PC\nwhole-tx"),
    ("lncRNA_whole",  "lncRNA",         "lncRNA\nwhole-tx"),
    ("pc_uorf5",      "uorf_5utr",      "PC 5'UTR\n(uORF)"),
    ("pc_dorf3",      "dorf_3utr",      "PC 3'UTR\n(dORF)"),
]
CEIL_C, MODEL_C = "#B0B0B0", "#2C6FBB"


def spearman_brown(r_half):
    return 2 * r_half / (1 + r_half)


def main():
    rep = json.load(open(REP))
    ceil = {}
    for ck, _, _ in STRATA:
        vals = [spearman_brown(s["median"][ck]) for s in rep["per_seed"] if ck in s["median"]]
        if not vals:
            raise SystemExit(f"no ceiling seeds for stratum {ck} in {REP}")
        ceil[ck] = float(np.mean(vals))

    ex = json.load(open(EXTRA))
    model = {}
    for _, mk, _ in STRATA:
        if mk not in ex or "profile_pearson_median" not in ex[mk]:
            raise SystemExit(f"no profile_pearson_median for {mk} in {EXTRA}")
        model[mk] = float(ex[mk]["profile_pearson_median"])

    x = np.arange(len(STRATA)); w = 0.38
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    c = [ceil[ck] for ck, _, _ in STRATA]
    m = [model[mk] for _, mk, _ in STRATA]
    ax.bar(x - w / 2, c, w, label="Reproducible ceiling\n(replicate concordance)", color=CEIL_C,
           edgecolor="white", lw=0.5)
    ax.bar(x + w / 2, m, w, label="Model achieved", color=MODEL_C, edgecolor="white", lw=0.5)
    for xi, (ck, mk, _) in enumerate(STRATA):
        ax.text(xi - w / 2, ceil[ck] + 0.008, f"{ceil[ck]:.2f}", ha="center", va="bottom",
                fontsize=7.5, color="#555")
        ax.text(xi + w / 2, model[mk] + 0.008, f"{model[mk]:.2f}", ha="center", va="bottom",
                fontsize=7.5, color=MODEL_C)
        ax.text(xi + w / 2, model[mk] / 2, f"{100 * model[mk] / ceil[ck]:.0f}%\nof ceiling",
                ha="center", va="center", fontsize=6.5, color="white", fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels([l for _, _, l in STRATA], fontsize=8.5)
    ax.set_ylabel("Per-transcript Pearson", fontsize=9); ax.set_ylim(0, 1.0)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fracs = [100 * model[mk] / ceil[ck] for ck, mk, _ in STRATA]
    ax.set_title(f"Model reaches {min(fracs):.0f}-{max(fracs):.0f}% of what the assay reproduces\n"
                 "(ceiling = Spearman-Brown-corrected half-pool concordance)", fontsize=9.5)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"B4_replicate_ceiling{SUF}.{ext}", dpi=300, bbox_inches="tight")

    out = {"run": RUN, "ceiling_source": str(REP.relative_to(NEW)),
           "model_source": str(EXTRA.relative_to(NEW)),
           "strata": [{"label": l.replace("\n", " "), "ceiling_key": ck, "model_key": mk,
                       "ceiling": ceil[ck], "model": model[mk],
                       "pct_of_ceiling": 100 * model[mk] / ceil[ck]}
                      for ck, mk, l in STRATA]}
    (HERE / f"B4_values{SUF}.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"run: {RUN}")
    for ck, mk, l in STRATA:
        print(f"  {l.replace(chr(10), ' '):<18} ceiling {ceil[ck]:.3f}  model {model[mk]:.3f}  "
              f"({100 * model[mk] / ceil[ck]:.0f}%)")


if __name__ == "__main__":
    main()
