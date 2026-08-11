#!/usr/bin/env python3
"""Figure E1: the RNA-seq input does not have to come from the same experiment as the Ribo-seq.

The 3 Ribo x 3 RNA factorial as a heatmap, one panel per prediction arm. Rows are the RNA input,
columns the Ribo-seq reference, cells are F1. The DIAGONAL is the matched condition.

The result is the absence of a pattern: if provenance mattered the diagonal would be the brightest
cell in every column, and it is not. Matched-minus-mismatched averages +0.0030 F1 across 12 model x
reference x arm combinations, range -0.0057 to +0.0079, and goes the wrong way twice. The diagonal is
outlined so a reader can check that for themselves rather than take the summary number on trust.

What DOES vary is the row: Wang RNA (2 samples) is consistently the weakest input regardless of which
Ribo-seq it is predicting. Library quality matters; provenance does not.

Reads results/mouse_liver_3x3/scored/. cas12a env.
"""
from __future__ import annotations

import csv
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

ROOT = pathlib.Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model")
S = ROOT / "results/mouse_liver_3x3/scored"
OUT = ROOT / "figures/E1_rna_provenance"
DS = ["janich", "gse243134", "wang"]
NSAMP = {"janich": 7, "gse243134": 19, "wang": 2}
MODEL = "mamba4"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [r for r in csv.DictReader(open(S / "model_vs_ribo.tsv"), delimiter="\t")
            if r["model"] == MODEL]
    deltas = [float(r["delta"]) for r in
              csv.DictReader(open(S / "matched_vs_mismatched_rna.tsv"), delimiter="\t")]

    arms = [("pred_obsdepth", "A  Predicted shape at observed depth"),
            ("pred_preddepth", "B  Standalone (no Ribo-seq at inference)")]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.1))
    grids = {}
    for ax, (arm, title) in zip(axes, arms):
        M = np.full((3, 3), np.nan)
        for i, rna in enumerate(DS):
            for j, ref in enumerate(DS):
                m = [r for r in rows if r["arm"] == arm and r["rna_input"] == rna
                     and r["ribo_reference"] == ref]
                if m:
                    M[i, j] = float(m[0]["f1"])
        grids[arm] = M.tolist()
        # Per-arm colour scale: the whole point is that within-panel variation is tiny, and a shared
        # 0-1 scale would render both panels flat and unreadable.
        im = ax.imshow(M, cmap="viridis", vmin=np.nanmin(M), vmax=np.nanmax(M))
        for i in range(3):
            for j in range(3):
                if np.isnan(M[i, j]):
                    continue
                mid = (np.nanmin(M) + np.nanmax(M)) / 2
                ax.text(j, i, f"{M[i, j]:.3f}", ha="center", va="center", fontsize=9.5,
                        color="white" if M[i, j] < mid else "black",
                        fontweight="bold" if i == j else "normal")
                if i == j:
                    ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                           edgecolor="#C0392B", linewidth=2.6, zorder=5))
        ax.set_xticks(range(3)); ax.set_xticklabels(DS, fontsize=8.5)
        ax.set_yticks(range(3))
        ax.set_yticklabels([f"{d}\n({NSAMP[d]} samples)" for d in DS], fontsize=8.5)
        ax.set_xlabel("Ribo-seq reference")
        ax.set_ylabel("RNA-seq input")
        ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).ax.tick_params(labelsize=7.5)

    fig.suptitle(f"Red outline = matched RNA. Matched minus mismatched: "
                 f"{np.mean(deltas):+.4f} F1 mean "
                 f"({min(deltas):+.4f} to {max(deltas):+.4f}, negative in "
                 f"{sum(1 for d in deltas if d < 0)} of {len(deltas)})",
                 fontsize=9, y=0.02)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"E1_rna_provenance.{ext}", dpi=200, bbox_inches="tight")

    (OUT / "E1_values.json").write_text(json.dumps(
        dict(model=MODEL, datasets=DS, n_rna_samples=NSAMP, grids=grids,
             delta_mean=float(np.mean(deltas)), delta_min=min(deltas), delta_max=max(deltas)),
        indent=2))
    with open(OUT / "E1_values.tsv", "w") as fh:
        fh.write("arm\trna_input\tribo_reference\tf1\tmatched\n")
        for arm, _ in arms:
            for i, rna in enumerate(DS):
                for j, ref in enumerate(DS):
                    fh.write(f"{arm}\t{rna}\t{ref}\t{grids[arm][i][j]:.4f}\t{int(i == j)}\n")
    print(f"  wrote {OUT}/E1_rna_provenance.{{pdf,png}}")
    print(f"    matched-minus-mismatched: {np.mean(deltas):+.4f} mean, "
          f"{min(deltas):+.4f} to {max(deltas):+.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
