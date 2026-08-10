#!/usr/bin/env python3
"""Figure A3: the model generalizes EVENLY across 9 leave-one-tissue-out folds; the apparent spread in
whole-transcript profile Pearson is held-out TARGET quality, not a generalization gradient.

For each of 9 Chothani tissues, a model is trained on the other 8 and evaluated on the held-out tissue
(results/loto_9fold/onehot_orf_v2_attn_holdout_<tissue>/). Raw pc profile Pearson ranges 0.23 (Brain) to
0.64 (Hepatocytes) -- which looks like uneven generalization until you notice it is 0.81-correlated with
the held-out tissue's OWN observed 3-nt periodicity (period_obs): Brain's Ribo-seq has almost no
periodicity (period_obs 0.044), so there is little periodic signal in the target to match. The model's
predicted periodicity (period_pred) stays high everywhere, and the abundance/count transfer is even
(pc count Pearson 0.59-0.88, CV 0.11 vs 0.24 for profile).

Panel (a): per-fold pc profile Pearson vs held-out period_obs, with the regression line and r.
Panel (b): pc count Pearson (abundance transfer) across the 9 tissues -- high and even.

CHECKPOINT (audited 2026-08-08, DELIBERATELY NOT REPOINTED). This is the one figure in the project
still built on the pre-union recipe, and that is a decision rather than an oversight. The 9-fold LOTO
sweep exists only as `results/loto_9fold/onehot_orf_v2_attn_holdout_<tissue>`; there is no union
equivalent, and producing one means retraining nine models (12-17 GPU-hours each). Every other panel
was moved to the released checkpoints, so this one must be LABELLED wherever it appears -- the caption
should say "9-fold LOTO, pre-union Fibroblast-universe recipe" rather than silently sitting beside
released-model panels.

The claim it makes is about the SHAPE of the spread across folds (that it tracks held-out target
periodicity at r=0.81, rather than being a generalization gradient), which is a property of the LOTO
design and the data, not of the checkpoint. The absolute Pearson values ARE checkpoint-dependent and
should not be quoted next to released-model numbers.

Reusable: reads results/loto_9fold/*/{test_metrics,extra_metrics}.json. cas12a env.
"""
from __future__ import annotations

import glob
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/A3_loto_spread"
BASE = NEW / "results/loto_9fold"
PROF_C, COUNT_C = "#2C6FBB", "#3f9b6d"


def load():
    rows = []
    for d in sorted(glob.glob(str(BASE / "*/"))):
        t = os.path.basename(d.rstrip("/")).split("holdout_")[-1]
        tmf, emf = os.path.join(d, "test_metrics.json"), os.path.join(d, "extra_metrics.json")
        if not (os.path.exists(tmf) and os.path.exists(emf)):
            continue
        tm, em = json.load(open(tmf)), json.load(open(emf))
        rows.append(dict(
            tissue=t,
            pc_prof=tm["protein_coding"]["pearson_median"],
            period_obs=tm["protein_coding"]["period_obs_median"],
            period_pred=tm["protein_coding"]["period_pred_median"],
            pc_count=em["protein_coding"]["count_pearson"]))
    return rows


def main():
    rows = load()
    rows.sort(key=lambda r: r["period_obs"])
    tis = [r["tissue"] for r in rows]
    prof = np.array([r["pc_prof"] for r in rows])
    po = np.array([r["period_obs"] for r in rows])
    count = np.array([r["pc_count"] for r in rows])
    r_prof = float(np.corrcoef(po, prof)[0, 1])

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(9.4, 3.8),
                                   gridspec_kw=dict(width_ratios=[1.1, 1.15]))

    # panel a: profile Pearson vs held-out target periodicity
    axa.scatter(po, prof, s=48, c=PROF_C, edgecolor="white", lw=0.7, zorder=3)
    m, b = np.polyfit(po, prof, 1)
    xs = np.linspace(po.min(), po.max(), 50)
    axa.plot(xs, m * xs + b, "--", color="#888", lw=1.2, zorder=2)
    for r in rows:
        dy = 0.018 if r["tissue"] not in ("HCAEC", "ES") else -0.028
        axa.annotate(r["tissue"], (r["period_obs"], r["pc_prof"]),
                     textcoords="offset points", xytext=(4, 3 if dy > 0 else -9),
                     fontsize=6.8, color="#333")
    axa.set_xlabel("held-out tissue observed 3-nt periodicity\n(period_obs = target quality)", fontsize=8.6)
    axa.set_ylabel("pc whole-tx profile Pearson", fontsize=9)
    axa.text(0.04, 0.96, f"r = {r_prof:.2f}\n(spread is target quality,\nnot a generalization gradient)",
             transform=axa.transAxes, ha="left", va="top", fontsize=8,
             bbox=dict(boxstyle="round,pad=0.35", fc="#eef4fb", ec=PROF_C, lw=0.7))
    axa.spines[["top", "right"]].set_visible(False)
    axa.set_title("a  Profile-Pearson spread = held-out data quality", fontsize=9, loc="left")

    # panel b: abundance/count transfer, even across tissues
    order = np.argsort(count)
    y = np.arange(len(tis))
    axb.barh(y, count[order], color=COUNT_C, edgecolor="white", lw=0.5, height=0.7)
    for yi, v, ti in zip(y, count[order], np.array(tis)[order]):
        axb.text(v + 0.008, yi, f"{v:.2f}", va="center", fontsize=7.2, color="#2c6b4a")
    axb.set_yticks(y)
    axb.set_yticklabels(np.array(tis)[order], fontsize=7.8)
    axb.set_xlim(0, 1.0)
    axb.axvline(np.median(count), ls=":", lw=1.0, color="#666")
    axb.text(np.median(count), len(tis) - 0.4, f"median {np.median(count):.2f}",
             fontsize=7, color="#666", ha="center", va="bottom")
    axb.set_xlabel("pc count (abundance) Pearson", fontsize=9)
    axb.spines[["top", "right"]].set_visible(False)
    axb.set_title("b  Abundance transfer is even across tissues", fontsize=9, loc="left")

    fig.suptitle("Even cross-tissue generalization (9-fold leave-one-tissue-out)", fontsize=10.5, y=1.02)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"A3_loto_spread.{ext}", dpi=300, bbox_inches="tight")
    (HERE / "A3_values.json").write_text(json.dumps({
        "per_fold": rows,
        "corr_periodobs_profilePearson": r_prof,
        "pc_profile_pearson_range": [float(prof.min()), float(prof.max())],
        "pc_profile_pearson_CV": float(prof.std() / prof.mean()),
        "pc_count_pearson_range": [float(count.min()), float(count.max())],
        "pc_count_pearson_CV": float(count.std() / count.mean()),
        "pc_count_pearson_median": float(np.median(count)),
    }, indent=2) + "\n")
    print(f"  9 folds; profile Pearson {prof.min():.3f}-{prof.max():.3f} (CV {prof.std()/prof.mean():.2f})")
    print(f"  count   Pearson {count.min():.3f}-{count.max():.3f} (CV {count.std()/count.mean():.2f}) "
          f"median {np.median(count):.3f}")
    print(f"  corr(period_obs, profile Pearson) = {r_prof:.3f}")


if __name__ == "__main__":
    main()
