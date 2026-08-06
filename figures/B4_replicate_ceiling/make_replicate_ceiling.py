#!/usr/bin/env python3
"""Figure B4: how much of the Ribo-seq profile is even PREDICTABLE, and how much the model captures. Split
32 Fibroblast samples into two 16-sample half-pools, per-transcript Pearson between them, Spearman-Brown-
corrected to full depth = the reproducible CEILING (what any model could reach). Bars: ceiling vs the model's
achieved per-transcript Pearson, for pc whole-tx, lncRNA whole-tx, and pc CDS per-codon (elongation). Honest
'how good / how much headroom' panel. Reusable: ceiling read from replicate_concordance.json; MODEL dict is
the current model's eval -- update on retrain (mm1 ablation). cas12a env.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/B4_replicate_ceiling"
REP = NEW / "results/replicate_concordance.json"

# strata -> (json key, display label). MODEL = current orf_v2_attn per-tx Pearson (results.md Task 22).
STRATA = [("pc_whole", "PC\nwhole-tx"), ("lncRNA_whole", "lncRNA\nwhole-tx"), ("pc_cds_codon", "PC CDS\nper-codon")]
MODEL = {"pc_whole": 0.638, "lncRNA_whole": 0.364, "pc_cds_codon": 0.526}
CEIL_C, MODEL_C = "#B0B0B0", "#2C6FBB"


def spearman_brown(r_half):
    return 2 * r_half / (1 + r_half)


def main():
    rep = json.load(open(REP))
    # ceiling = mean over seeds of Spearman-Brown-corrected half-pool median
    ceil = {}
    for key, _ in STRATA:
        vals = [spearman_brown(s["median"][key]) for s in rep["per_seed"] if key in s["median"]]
        ceil[key] = float(np.mean(vals))

    x = np.arange(len(STRATA)); w = 0.38
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    c = [ceil[k] for k, _ in STRATA]; m = [MODEL[k] for k, _ in STRATA]
    ax.bar(x - w/2, c, w, label="Reproducible ceiling\n(replicate concordance)", color=CEIL_C,
           edgecolor="white", lw=0.5)
    ax.bar(x + w/2, m, w, label="Model achieved", color=MODEL_C, edgecolor="white", lw=0.5)
    for xi, (k, _) in enumerate(STRATA):
        ax.text(xi - w/2, ceil[k] + 0.008, f"{ceil[k]:.2f}", ha="center", va="bottom", fontsize=7.5, color="#555")
        ax.text(xi + w/2, MODEL[k] + 0.008, f"{MODEL[k]:.2f}", ha="center", va="bottom", fontsize=7.5, color=MODEL_C)
        frac = 100 * MODEL[k] / ceil[k]
        ax.text(xi + w/2, MODEL[k] / 2, f"{frac:.0f}%\nof ceiling", ha="center", va="center", fontsize=6.5,
                color="white", fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels([l for _, l in STRATA], fontsize=8.5)
    ax.set_ylabel("Per-transcript Pearson", fontsize=9); ax.set_ylim(0, 1.0)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Profile is highly reproducible; the model captures ~60% of it\n(honest ceiling + headroom)",
                 fontsize=9.5)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"B4_replicate_ceiling.{ext}", dpi=300, bbox_inches="tight")
    (HERE / "B4_values.json").write_text(json.dumps({"ceiling": ceil, "model": MODEL}, indent=2) + "\n")
    for k, l in STRATA:
        print(f"  {l.replace(chr(10),' ')}: ceiling {ceil[k]:.3f}  model {MODEL[k]:.3f}  ({100*MODEL[k]/ceil[k]:.0f}%)")


if __name__ == "__main__":
    main()
