#!/usr/bin/env python3
"""Depth-crossover figure: ORF-calling F1 (and non-canonical recall) as a function of Ribo-seq
sequencing depth, with the depth-independent "predict from RNA-seq (no Ribo-seq)" line overlaid.

Reads depth_curve.tsv + depth_curve_summary.json (subsample_depth_curve.py). Two stacked panels vs
absolute test-set P-sites (log x): (top) drop-in F1, (bottom) non-canonical ORF recall. The
measurement curve (RiboCode on the real profile thinned to depth f) is markers + line with per-seed
std error bars; the model's pred_preddepth (standalone) and pred_obsdepth (predicted shape, full
depth) are horizontal reference lines. The crossover -- where the measurement F1 drops to the
pred_preddepth line -- is marked. cas12a env (numpy + matplotlib).

Usage: plot_depth_crossover.py --curve <depth_curve.tsv> --summary <depth_curve_summary.json> [--out png]
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def load(curve):
    rows = list(csv.DictReader(open(curve), delimiter="\t"))
    pred = {}
    meas = defaultdict(lambda: defaultdict(list))  # depth_frac -> metric -> [seed vals]
    for r in rows:
        if r["kind"] in ("pred_preddepth", "pred_obsdepth"):
            pred[r["kind"]] = {"f1": float(r["f1"]), "noncanon": float(r["noncanon_recall"])}
        elif r["kind"] == "measure":
            f = float(r["depth_frac"]); a = float(r["abs_psites"])
            meas[(f, a)]["f1"].append(float(r["f1"]))
            meas[(f, a)]["noncanon"].append(float(r["noncanon_recall"]))
    pts = sorted(meas.items())  # [((f,abs), {metric:[...]}), ...]
    return pred, pts


def series(pts, metric):
    x = np.array([a for (_f, a), _ in pts])
    m = np.array([np.mean(d[metric]) for _, d in pts])
    s = np.array([np.std(d[metric]) for _, d in pts])
    return x, m, s


def draw(ax, x, m, s, pred, key, title, ylab):
    ax.errorbar(x, m, yerr=s, marker="o", ms=5, lw=1.6, color="#1f77b4", capsize=3,
                label="measure at depth f (RiboCode on thinned real)")
    if "pred_preddepth" in pred:
        y = pred["pred_preddepth"][key]
        ax.axhline(y, ls="--", lw=1.8, color="#d62728",
                   label=f"predict from RNA-seq, NO Ribo-seq ({y:.3f})")
    if "pred_obsdepth" in pred:
        y = pred["pred_obsdepth"][key]
        ax.axhline(y, ls=":", lw=1.4, color="#9467bd",
                   label=f"predicted shape, full depth ({y:.3f})")
    ax.set_xscale("log")
    ax.set_ylabel(ylab)
    ax.set_title(title, fontsize=10)
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend(fontsize=7.5, loc="lower right", framealpha=0.9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curve", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    pred, pts = load(args.curve)
    summ = json.loads(Path(args.summary).read_text())
    full = summ["full_depth_psites"]
    cross = summ.get("crossover")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 8), sharex=True)
    x, m, s = series(pts, "f1")
    draw(ax1, x, m, s, pred, "f1", "ORF-calling drop-in F1 vs deep official calls", "drop-in F1")
    xn, mn, sn = series(pts, "noncanon")
    draw(ax2, xn, mn, sn, pred, "noncanon", "Non-canonical ORF recall (uORF / dORF / novel)",
         "non-canonical recall")
    ax2.set_xlabel("Ribo-seq sequencing depth (test-set P-sites; log scale)")

    xmin = min(a for (_f, a), _ in pts) / 1.6
    xmax = max(a for (_f, a), _ in pts) * 1.6
    for ax in (ax1, ax2):
        ax.set_xlim(xmin, xmax)
        if cross:
            ax.axvline(cross["abs_psites"], color="#333", lw=1.2, ls="-.")
    if cross:
        ax1.annotate(f"crossover ~{cross['abs_psites']:.1e} P-sites\n(f~{cross['depth_frac']:.3g}): below "
                     f"this, predicting\nfrom RNA-seq beats measuring",
                     xy=(cross["abs_psites"], pred["pred_preddepth"]["f1"]),
                     xytext=(0.05, 0.62), textcoords="axes fraction", fontsize=8,
                     arrowprops=dict(arrowstyle="->", color="#333", lw=1.2),
                     bbox=dict(boxstyle="round", fc="#fff3cd", ec="#856404", alpha=0.95))
    fig.suptitle(f"Predict vs measure: Hepatocytes held out, full test-set depth {full:.2e} P-sites",
                 fontsize=11, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    out = Path(args.out) if args.out else Path(args.curve).with_name("depth_crossover.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
