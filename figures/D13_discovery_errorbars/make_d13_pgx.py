#!/usr/bin/env python3
"""Figure D13 (pgx rebuild): the database-design tradeoff, on the RELEASED models.

Supersedes make_discovery_errorbars.py, which read the f0-threshold searches (old checkpoint,
pred_frame0 >= 0.5 selection, unfrozen calibrate_mass = 2). Same three-panel story, current data:

(a) GAIN   -- discovery efficiency: novel peptides per 1,000 novel DB sequences.
(b) COST   -- canonical displacement: GENCODE PSMs lost vs the gencode-only baseline (dPSM).
(c) DRIVER -- database size: novel sequences added (= added decoy load), the thing that causes both.

Reads the pgx report JSONs directly, so the figure tracks whatever the pipeline last produced:
  <dataset>_pilot/pgx_{mamba4,attn}/report/table.json

Both released models are drawn, because the MS application is model-agnostic (49 vs 48 novel peptides
across the four datasets) and the figure should not imply the GPU-only model is required.

cas12a env.  Run: python make_d13_pgx.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt      # noqa: E402
import numpy as np                   # noqa: E402

NEW = Path(__file__).resolve().parents[2]
PD = NEW / "proteogenomics/data"
HERE = NEW / "figures/D13_discovery_errorbars"

DATASETS = [("HBL1", "HBL-1"), ("DoHH2", "DoHH2"), ("SUDHL4", "SU-DHL-4"), ("A549", "A549")]
# (arm key, display, colour). null_nc exists only for A549 (tryptic); nonspecific digestion of a
# ~3M-sequence near-cognate null exceeds MSFragger's ~2e9 peptide cap.
ARMS = [("model_poisson", "model (Poisson)", "#2C6FBB"),
        ("model_standard", "model (theta=1)", "#7FA8D8"),
        ("null_atg", "null AUG", "#B0B0B0"),
        ("null_nc", "null near-cognate", "#8A8A8A")]


def load(tag):
    out = {}
    for key, _ in DATASETS:
        p = PD / f"{key}_pilot" / f"pgx_{tag}" / "report" / "table.json"
        if p.exists():
            out[key] = json.load(open(p))["arms"]
    return out


def main():
    data = {t: load(t) for t in ("mamba4", "attn")}
    missing = [f"{t}/{k}" for t in data for k, _ in DATASETS if k not in data[t]]
    if missing:
        print(f"[warn] absent (drawn as gaps, NOT zeros): {missing}")

    x = np.arange(len(DATASETS))
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 5.0))
    # Panel (b) plots dPSM SIGNED. A zero must be visibly "measured and zero" -- an absent bar would
    # read as missing data, and zero canonical cost is the model's headline result on 6 of 8 arms.
    panels = [
        ("a", "GAIN", "novel peptides per 1,000 DB seqs", "log",
         lambda a: (a["novel_peptides"] / a["db_novel_seqs"] * 1000) if a["db_novel_seqs"] else np.nan),
        ("b", "COST", "change in GENCODE PSMs vs baseline", "symlog",
         lambda a: float(a["d_gencode_psms"])),
        ("c", "DRIVER", "novel DB sequences (decoy load)", "log",
         lambda a: float(a["db_novel_seqs"])),
    ]
    n_arm = len(ARMS)
    w = 0.8 / (n_arm * 2)

    for ax, (tag, title, ylab, scale, fn) in zip(axes, panels):
        zero_xs = []
        for i, (arm, disp, col) in enumerate(ARMS):
            for j, (mdl, hatch) in enumerate((("mamba4", None), ("attn", "///"))):
                vals, xs = [], []
                for k, (key, _) in enumerate(DATASETS):
                    a = data[mdl].get(key, {}).get(arm)
                    if a is None:
                        continue                       # genuinely absent -> no bar, not a zero
                    v = fn(a)
                    xs.append(k + (i * 2 + j) * w - 0.4 + w / 2)
                    vals.append(v)
                    if tag == "b" and v == 0:
                        zero_xs.append(xs[-1])
                if vals:
                    ax.bar(xs, vals, w, color=col, hatch=hatch, edgecolor="white", linewidth=0.4,
                           label=f"{disp} [{mdl}]" if ax is axes[0] else None)
        ax.set_yscale(scale)
        if tag == "b":
            ax.axhline(0, color="0.25", lw=0.9)
            # explicit marker + label so a zero cannot be mistaken for a missing measurement
            if zero_xs:
                ax.plot(zero_xs, [0] * len(zero_xs), marker="D", ms=3.4, ls="none",
                        color="#2C6FBB", zorder=5, clip_on=False)
                # corner text, not an arrow into the plot area -- an annotation placed near the
                # bars ran behind them and became unreadable
                ax.text(0.02, 0.97, "◆  0 = measured, no canonical cost",
                        transform=ax.transAxes, ha="left", va="top",
                        fontsize=7.5, color="#2C6FBB")
            ax.set_ylabel(ylab + "\n(negative = canonical IDs destroyed)")
        else:
            ax.set_ylabel(ylab)
        ax.set_xticks(x)
        ax.set_xticklabels([d for _, d in DATASETS])
        ax.set_title(f"({tag}) {title}", loc="left", fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=7.5, ncol=4, frameon=False,
               loc="lower center", bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Database-design tradeoff: predicting which ORFs to search keeps the DB compact "
                 "(hatched = attn, solid = mamba4)", fontsize=10)
    fig.tight_layout(rect=(0, 0.07, 1, 0.94))
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"D13_discovery_tradeoff_pgx.{ext}", dpi=200, bbox_inches="tight")

    vals = {mdl: {k: {arm: {"density_per_1k": (a["novel_peptides"] / a["db_novel_seqs"] * 1000)
                           if a["db_novel_seqs"] else None,
                           "novel_peptides": a["novel_peptides"],
                           "db_novel_seqs": a["db_novel_seqs"],
                           "d_gencode_psms": a["d_gencode_psms"]}
                    for arm, a in arms.items() if arm != "gencode"}
                for k, arms in data[mdl].items()} for mdl in data}
    (HERE / "D13_pgx_values.json").write_text(json.dumps(vals, indent=2) + "\n")
    print(f"wrote {HERE}/D13_discovery_tradeoff_pgx.png (+ .pdf, + D13_pgx_values.json)")
    for mdl in ("mamba4", "attn"):
        for k, _ in DATASETS:
            a = data[mdl].get(k, {}).get("model_poisson")
            if a:
                print(f"  {mdl:7s} {k:8s} poisson density "
                      f"{a['novel_peptides'] / a['db_novel_seqs'] * 1000:6.3f}/1k  dPSM {a['d_gencode_psms']:+}")


if __name__ == "__main__":
    main()
