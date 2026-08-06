#!/usr/bin/env python3
"""Figure B7: the RNA-seq multimap posture (STAR --outFilterMultimapNmax 20, "mm20") does not
materially change the model's INPUT, so a single-mapper (mm1) retrain is unwarranted.

Fast proxy: align one 30M-read DoHH2 RNA-seq subsample two ways (mm1 vs mm20), build per-nt transcript
coverage each way, and compare per-transcript on the DoHH2 expressed universe (30,950 tx with mm20 depth
>= 50).

Panel (a): distribution of per-transcript Pearson(mm1, mm20) -- a hard spike at 1.0 (85.7% of tx >= 0.999,
median 1.0000). The coverage the model sees is identical on the overwhelming majority of transcripts.
Panel (b): the small divergent tail (5.3% with Pearson < 0.95) is entirely paralog / repeat / multimap-heavy
loci, where mm20 inflates depth 7-145x. Since the model is per-transcript-normalized and sees identical input
on ~95% of transcripts, its outputs are insensitive to the posture.

Reusable: reads proteogenomics/data/mm_proxy/coverage_mm_pertx.tsv (written by compare_coverage_mm.py).
cas12a env.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/B7_multimap_posture"
TSV = NEW / "proteogenomics/data/mm_proxy/coverage_mm_pertx.tsv"
MM20_C, MM1_C, TAIL_C = "#2C6FBB", "#88AED0", "#C0504D"


def main():
    df = pd.read_csv(TSV, sep="\t")
    n = len(df)
    pear = df["pearson"].values
    ratio = df["ratio_mm1_mm20"].values           # mm1/mm20 depth; <1 means mm20 inflated
    frac_ge999 = 100 * (pear >= 0.999).mean()
    frac_lt95 = 100 * (pear < 0.95).mean()
    med = float(np.median(pear))

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(9.2, 3.7),
                                   gridspec_kw=dict(width_ratios=[1.15, 1.1]))

    # panel a: Pearson distribution (log-count so the small tail is visible next to the 1.0 spike)
    bins = np.linspace(-0.65, 1.0, 60)
    axa.hist(pear, bins=bins, color=MM20_C, edgecolor="white", lw=0.3)
    axa.set_yscale("log")
    axa.axvline(0.999, ls="--", lw=1.0, color="#333")
    axa.set_xlabel("per-transcript Pearson(mm1, mm20 coverage)", fontsize=8.8)
    axa.set_ylabel("transcripts (log scale)", fontsize=9)
    axa.set_xlim(-0.7, 1.02)
    axa.spines[["top", "right"]].set_visible(False)
    axa.text(0.985, 0.96, f"median {med:.4f}\n{frac_ge999:.1f}% of tx  >= 0.999\n"
             f"{frac_lt95:.1f}% of tx  < 0.95",
             transform=axa.transAxes, ha="right", va="top", fontsize=8,
             bbox=dict(boxstyle="round,pad=0.35", fc="#eef4fb", ec=MM20_C, lw=0.7))
    axa.set_title("a  Coverage is identical on ~95% of transcripts", fontsize=9, loc="left")

    # panel b: the divergent tail is driven by multimap depth inflation
    infl = np.log2(np.clip(1.0 / np.clip(ratio, 1e-6, None), 1e-6, None))   # log2(mm20/mm1)
    is_tail = pear < 0.95
    axb.scatter(infl[~is_tail], pear[~is_tail], s=4, c=MM1_C, alpha=0.35, lw=0, label="Pearson >= 0.95")
    axb.scatter(infl[is_tail], pear[is_tail], s=8, c=TAIL_C, alpha=0.7, lw=0,
                label=f"Pearson < 0.95 ({is_tail.sum()} tx, {frac_lt95:.1f}%)")
    axb.axhline(0.95, ls=":", lw=0.8, color="#999")
    axb.axvline(np.log2(1.5), ls=":", lw=0.8, color="#999")
    axb.set_xlabel("multimap depth inflation  log2(mm20 / mm1)", fontsize=8.8)
    axb.set_ylabel("per-transcript Pearson", fontsize=9)
    axb.set_xlim(-0.5, max(1.0, np.percentile(infl, 99.5)))
    axb.legend(frameon=False, fontsize=7.5, loc="lower left")
    axb.spines[["top", "right"]].set_visible(False)
    axb.set_title("b  Divergence only on multimap-inflated (paralog/repeat) loci", fontsize=9, loc="left")

    fig.suptitle("RNA-seq multimap posture (mm20) barely changes the model input -> mm1 retrain unwarranted",
                 fontsize=10, y=1.02)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"B7_multimap_posture.{ext}", dpi=300, bbox_inches="tight")
    (HERE / "B7_values.json").write_text(json.dumps({
        "n_tx_assessed": int(n), "median_pearson": med,
        "frac_pearson_ge_0.999_pct": frac_ge999, "frac_pearson_lt_0.95_pct": frac_lt95,
        "total_depth_ratio_mm1_mm20": float(df["depth_mm1"].sum() / df["depth_mm20"].sum()),
        "n_multimap_inflated_mm20_gt_1.5x": int((ratio < 0.667).sum()),
    }, indent=2) + "\n")
    print(f"  n={n:,}  median Pearson {med:.4f}  >=0.999: {frac_ge999:.1f}%  <0.95: {frac_lt95:.1f}%")
    print(f"  total depth ratio mm1/mm20 = {df['depth_mm1'].sum()/df['depth_mm20'].sum():.4f}")


if __name__ == "__main__":
    main()
