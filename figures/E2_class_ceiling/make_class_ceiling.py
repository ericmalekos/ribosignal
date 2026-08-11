#!/usr/bin/env python3
"""Figure E2: model performance per ORF class, against the between-experiment ceiling.

The single most misread number in this project is non-canonical F1. Judged against 1.0 it looks like
failure; judged against what two real Ribo-seq experiments achieve against EACH OTHER it is a
different story per class. This figure puts the two on the same axis so the comparison cannot be
skipped.

For each ORF class: a shaded band spanning the six ordered experiment-vs-experiment pairs (the
ceiling), with model cells overlaid as points, split by arm. Where a point sits inside the band, the
model is at experimental agreement. Where it sits below, the gap is the honest deficit.

The `internal` row is the reason this figure exists in this form: the band sits at 0.539-0.620 and
the model points sit at essentially zero. An overall-F1 plot hides that completely.

Reads results/mouse_liver_3x3/scored/. cas12a env.
"""
from __future__ import annotations

import collections
import csv
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

ROOT = pathlib.Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model")
S = ROOT / "results/mouse_liver_3x3/scored"
OUT = ROOT / "figures/E2_class_ceiling"
CEIL_C = "#1F4E79"
ARM_C = {"pred_obsdepth": "#C0392B", "pred_preddepth": "#E8963C"}
CLASSES = ["annotated", "uORF", "novel", "Overlap_uORF", "dORF", "internal"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ceil = collections.defaultdict(list)
    for r in csv.DictReader(open(S / "ribo_vs_ribo_by_class.tsv"), delimiter="\t"):
        ceil[r["orf_class"]].append(float(r["f1"]))
    mod = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in csv.DictReader(open(S / "model_vs_ribo_by_class.tsv"), delimiter="\t"):
        mod[r["orf_class"]][r["arm"]].append(float(r["f1"]))

    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    rows, vals = [], []
    for i, c in enumerate(CLASSES):
        y = len(CLASSES) - 1 - i
        cv = ceil.get(c, [])
        if cv:
            ax.barh(y, max(cv) - min(cv), left=min(cv), height=0.52,
                    color=CEIL_C, alpha=0.22, zorder=1,
                    edgecolor=CEIL_C, linewidth=1.0)
            ax.plot([np.mean(cv)] * 2, [y - 0.26, y + 0.26], color=CEIL_C, lw=2.0, zorder=2)
        for j, (arm, col) in enumerate(ARM_C.items()):
            mv = mod[c].get(arm, [])
            if not mv:
                continue
            off = (j - 0.5) * 0.20
            ax.scatter(mv, [y + off] * len(mv), s=17, c=col, alpha=0.55, zorder=3,
                       edgecolor="none")
            ax.scatter([np.mean(mv)], [y + off], s=95, c=col, zorder=4,
                       edgecolor="white", linewidth=1.3)
        rows.append(c)
        vals.append(dict(orf_class=c,
                         ceiling_min=min(cv) if cv else None, ceiling_max=max(cv) if cv else None,
                         **{f"model_{a}_mean": (float(np.mean(mod[c][a])) if mod[c].get(a) else None)
                            for a in ARM_C}))
    ax.set_yticks(range(len(CLASSES)))
    ax.set_yticklabels(list(reversed(CLASSES)), fontsize=9.5)
    ax.set_xlabel("F1 against the reference dataset's observed calls")
    ax.set_xlim(-0.02, 1.02)
    ax.grid(axis="x", alpha=0.25, linewidth=0.6)
    ax.set_title("Model vs the between-experiment ceiling, by ORF class",
                 loc="left", fontsize=11, fontweight="bold")
    # Legend must name EVERY encoding on the panel. Three of them were unlabelled in the first
    # version (small vs large dot, and the vertical line inside the band), which is exactly the kind
    # of thing a reader should not have to reverse-engineer from the source.
    handles = [
        Patch(facecolor=CEIL_C, alpha=0.22, edgecolor=CEIL_C, linewidth=1.0,
              label="experiment vs experiment: range over the 6 ordered pairs"),
        Line2D([0], [0], color=CEIL_C, lw=2.0, label="mean of those 6 pairs"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="0.45", markersize=4.5,
               alpha=0.55, label="small dot: one model cell (2 models x 9 RNA/Ribo = 18 per arm)"),
    ] + [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=col, markersize=9,
               markeredgecolor="white", markeredgewidth=1.3,
               label=f"large dot: mean of those 18 ({arm})")
        for arm, col in ARM_C.items()
    ]
    ax.legend(handles=handles, frameon=False, fontsize=7.6, loc="upper left",
              bbox_to_anchor=(0.0, -0.13), ncol=2, handletextpad=0.6,
              columnspacing=1.4, labelspacing=0.5)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"E2_class_ceiling.{ext}", dpi=200, bbox_inches="tight")

    (OUT / "E2_values.json").write_text(json.dumps(vals, indent=2))
    with open(OUT / "E2_values.tsv", "w") as fh:
        fh.write("orf_class\tceiling_min\tceiling_max\tmodel_obsdepth_mean\tmodel_preddepth_mean\n")
        for v in vals:
            fh.write(f"{v['orf_class']}\t{v['ceiling_min']:.4f}\t{v['ceiling_max']:.4f}\t"
                     f"{v['model_pred_obsdepth_mean']:.4f}\t{v['model_pred_preddepth_mean']:.4f}\n")
    print(f"  wrote {OUT}/E2_class_ceiling.{{pdf,png}}")
    for v in vals:
        print(f"    {v['orf_class']:<14} ceiling {v['ceiling_min']:.3f}-{v['ceiling_max']:.3f}   "
              f"model obs={v['model_pred_obsdepth_mean']:.3f} pred={v['model_pred_preddepth_mean']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
