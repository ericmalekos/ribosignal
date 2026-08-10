#!/usr/bin/env python3
"""Supplemental FDR-rigor panel (FIGURES_PLAN rigor point 3).

THE QUESTION A REVIEWER WILL ASK: "are these novel peptides real, or are they what falls out of
searching a bigger database at a threshold that was never valid for them?"

The answer this project uses everywhere is a CLASS-SPECIFIC cut: novel targets competed against
`REV_nuORF|` decoys only, never against the global target-decoy pool. This panel shows what that
choice is worth, measured on the shipped searches rather than argued from principle.

WHAT THE MEASUREMENT ACTUALLY SHOWS, and it is not the simple story. A global 1% FDR does not
merely over-report the non-canonical class -- it fails to control it in EITHER direction, and which
way it errs depends on how much of the search the canonical class occupies:

  * A549 (tryptic, 65,405 canonical peptides): the global cut is set almost entirely by canonical
    matches, which are numerous and high-scoring. Against the naive AUG database it admits 375 novel
    peptides where a class-specific cut admits 10 -- a 37.5x over-report.
  * The three immunopeptidomes (nonspecific, 813-2,429 canonical peptides): the canonical class is
    too sparse to dominate, and the global cut lands STRICTER than the class-specific one, ratios
    down to 0.2. Real identifications are discarded.

Either way the reported number is not a 1% FDR for that class. That is the point, and it is a
stronger point than "global inflates by 10x" because it cannot be dismissed as conservative.

  (a) The discrepancy: novel peptides at global vs class-specific 1% FDR, every arm x dataset.
  (b) Why it happens: score distributions of canonical targets, novel targets and novel decoys for
      the worst case. The global cut sits well inside the novel class's decoy bulk.
  (c) What predicts the error: the size of the canonical class, not the size of the database.

Usage: make_fdr_rigor.py [--out <dir>]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "proteogenomics" / "scripts"))
from pgx.report import load_rank1  # noqa: E402

DATASETS = ["A549", "HBL1", "SUDHL4", "DoHH2"]
ARMS = ["null_atg", "cpat", "cpc2", "model_standard", "model_poisson"]
ARM_LABEL = {"null_atg": "null (naive AUG)", "cpat": "CPAT", "cpc2": "CPC2",
             "model_standard": "model (theta=1)", "model_poisson": "model (Poisson)"}
ARM_COLOR = {"null_atg": "#9a9a9a", "cpat": "#7b5ea7", "cpc2": "#c07ab8",
             "model_standard": "#2f6f9f", "model_poisson": "#1b4f72"}
MODEL = "mamba4"
FDR = 0.01
DEMO = ("A549", "null_atg")    # panel (b): the arm where the discrepancy is largest (37.5x)


def threshold(rows, tset, dset, fdr=FDR):
    """Lowest-score cut whose running decoy/target ratio stays <= fdr, walking scores downward."""
    ranked = sorted(((h, k) for (_, h, k, _) in rows if k in tset | dset), key=lambda x: -x[0])
    t = d = 0
    last = None
    for h, k in ranked:
        if k in tset:
            t += 1
        else:
            d += 1
        if t and d / t <= fdr:
            last = h
    return last


def npept(rows, cut):
    """Distinct novel peptide sequences at or above `cut`."""
    if cut is None:
        return 0
    return len({p for (p, h, k, _) in rows if k == "novel_t" and h >= cut})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent))
    a = ap.parse_args()
    out = Path(a.out)

    recs, missing, demo_rows = [], [], None
    for ds in DATASETS:
        for arm in ARMS:
            d = ROOT / f"proteogenomics/data/{ds}_pilot/pgx_{MODEL}/search/{arm}"
            if not (d / "rank1.tsv.gz").exists():
                missing.append(f"{ds}/{arm}")
                continue
            rows = load_rank1(d)
            if (ds, arm) == DEMO:
                demo_rows = rows
            # GLOBAL: all targets vs all decoys, the threshold the search itself reports.
            g = threshold(rows, {"canon_t", "novel_t"}, {"novel_d", "other_d"})
            # CLASS-SPECIFIC: novel targets vs novel decoys only.
            c = threshold(rows, {"novel_t"}, {"novel_d"})
            fa = ROOT / f"proteogenomics/data/{ds}_pilot/pgx_{MODEL}/db/db_{arm}.fasta"
            n_db = 0
            if fa.exists():
                with open(fa) as fh:
                    n_db = sum(1 for ln in fh if ln.startswith(">") and "nuORF|" in ln
                               and not ln.startswith(">REV_"))
            n_canon = len({p for (p, h, k, _) in rows if k == "canon_t" and g is not None and h >= g})
            recs.append({"dataset": ds, "arm": arm, "db_novel_seqs": n_db,
                         "global_cut": g, "class_cut": c, "canon_peptides_global": n_canon,
                         "novel_global": npept(rows, g), "novel_class": npept(rows, c)})
    if missing:
        print("MISSING search arms (absent from the panel, not silently zeroed):")
        for m in missing:
            print(f"  {m}")

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))

    # (a) paired global vs class-specific novel counts
    ax = axes[0]
    xs, labels = [], []
    for i, ds in enumerate(DATASETS):
        rs = [r for r in recs if r["dataset"] == ds]
        for j, r in enumerate(rs):
            x = i * (len(ARMS) + 1.2) + j
            xs.append(x)
            ax.plot([x, x], [r["novel_class"], r["novel_global"]], "-", color="#bbbbbb", lw=1,
                    zorder=1)
            ax.plot(x, r["novel_global"], "v", ms=6, color="#c0392b", zorder=3)
            ax.plot(x, r["novel_class"], "o", ms=5.5, color=ARM_COLOR[r["arm"]], zorder=3)
        if rs:
            labels.append((i * (len(ARMS) + 1.2) + (len(rs) - 1) / 2, ds))
    ax.set_yscale("symlog", linthresh=10)
    ax.set_xticks([p for p, _ in labels]); ax.set_xticklabels([t for _, t in labels], fontsize=9)
    ax.set_ylabel("novel peptides at 1% FDR", fontsize=9)
    ax.set_title("(a) a global cut does not control the non-canonical class\n"
                 "red = global 1% cut, coloured = class-specific cut (used throughout)", fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    # Headroom so the legend cannot occlude the A549 null_atg global point (375, the whole point of
    # the panel). An earlier version put the legend upper-left and hid exactly that marker.
    top = max(r["novel_global"] for r in recs) if recs else 10
    ax.set_ylim(0, top * 9)
    handles = [plt.Line2D([], [], marker="o", ls="", color=ARM_COLOR[k], label=ARM_LABEL[k])
               for k in ARMS]
    handles.append(plt.Line2D([], [], marker="v", ls="", color="#c0392b", label="global 1% FDR"))
    ax.legend(handles=handles, fontsize=7, loc="upper center", ncol=3, framealpha=0.95)

    # (b) why: score distributions for one search
    ax = axes[1]
    if demo_rows is not None:
        canon = np.array([h for (_, h, k, _) in demo_rows if k == "canon_t"])
        nt = np.array([h for (_, h, k, _) in demo_rows if k == "novel_t"])
        nd = np.array([h for (_, h, k, _) in demo_rows if k == "novel_d"])
        hi = float(np.percentile(np.concatenate([canon, nt, nd]), 99.5))
        bins = np.linspace(0, hi, 60)
        ax.hist(canon, bins=bins, density=True, color="#2f6f9f", alpha=0.5, label="canonical targets")
        ax.hist(nt, bins=bins, density=True, histtype="step", lw=1.8, color="#c0392b",
                label="novel targets")
        ax.hist(nd, bins=bins, density=True, histtype="step", lw=1.8, color="#999999", ls="--",
                label="novel decoys (REV_nuORF|)")
        r = next(r for r in recs if (r["dataset"], r["arm"]) == DEMO)
        for v, c, lab in ((r["global_cut"], "#c0392b", "global 1% cut"),
                          (r["class_cut"], "#1b4f72", "class-specific 1% cut")):
            if v is not None:
                ax.axvline(v, color=c, lw=1.6, ls=":")
                ax.text(v, ax.get_ylim()[1] * 0.97, f" {lab}\n {v:.1f}", fontsize=7.5, color=c,
                        va="top", rotation=0)
        ax.set_xlabel("MSFragger hyperscore", fontsize=9)
        ax.set_ylabel("density", fontsize=9)
        ax.set_title(f"(b) why: {DEMO[0]}, {ARM_LABEL[DEMO[1]]}\n"
                     "the global cut sits inside the novel decoy bulk", fontsize=9)
        ax.legend(fontsize=7.5, loc="center right", framealpha=0.95)

    # (c) what predicts the error: the size of the canonical class, not of the database
    ax = axes[2]
    for r in recs:
        if r["novel_class"] <= 0 or r["canon_peptides_global"] <= 0:
            continue
        ax.plot(r["canon_peptides_global"], r["novel_global"] / r["novel_class"], "o", ms=7,
                color=ARM_COLOR[r["arm"]], alpha=0.85)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.axhline(1.0, color="k", ls="--", lw=1)
    ax.set_xlabel("canonical peptides in the search (global 1% FDR)", fontsize=9)
    ax.set_ylabel("global / class-specific novel peptides", fontsize=9)
    lo, hi = ax.get_xlim()
    ax.text(lo * 1.15, 1.08, "global cut happens to agree", fontsize=7.5)
    ax.text(lo * 1.15, ax.get_ylim()[1] * 0.55, "global OVER-reports\n(A549, tryptic)", fontsize=7.5,
            color="#c0392b")
    ax.text(lo * 1.15, ax.get_ylim()[0] * 1.35, "global UNDER-reports\n(immunopeptidomes)",
            fontsize=7.5, color="#1b4f72")
    ax.set_title("(c) the error is driven by the canonical class,\nnot by database size", fontsize=9)
    ax.grid(alpha=0.25, which="both")

    fig.suptitle("Supplemental: FDR control for non-canonical identifications", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    for ext in ("png", "pdf"):
        fig.savefig(out / f"fdr_rigor.{ext}", dpi=200)
    (out / "fdr_rigor_values.json").write_text(
        json.dumps({"model": MODEL, "fdr": FDR, "records": recs, "missing": missing}, indent=2))

    print(f"\n{'dataset':<9}{'arm':<16}{'DB novel':>10}{'global':>8}{'class':>7}{'inflation':>11}")
    for r in recs:
        infl = (r["novel_global"] / r["novel_class"]) if r["novel_class"] else float("nan")
        print(f"{r['dataset']:<9}{r['arm']:<16}{r['db_novel_seqs']:>10,}"
              f"{r['novel_global']:>8}{r['novel_class']:>7}{infl:>11.1f}")
    print(f"\nwrote {out}/fdr_rigor.png|.pdf|_values.json")


if __name__ == "__main__":
    main()
