#!/usr/bin/env python3
"""Supplemental: when semi-supervised rescoring helps, and when it does not (FIGURES_PLAN rigor 4).

MS2Rescore was run on both substrates in this project and gave opposite answers: nothing on the A549
tryptic whole proteome, a large gain on the HBL-1 immunopeptidome. That is not a contradiction, and
showing why is worth a small panel -- it demonstrates that the pipeline's choice to report RAW
hyperscore class-specific FDR everywhere is a considered decision rather than an omission.

THE MECHANISM. mokapot (inside MS2Rescore) is semi-supervised: it trains a classifier to separate
targets from decoys and needs a populated, LEARNABLE target class. The controlling variable is
therefore how much confident signal the non-canonical class contains, not how much of it there is.

MEASURING THIS CORRECTED THE OBVIOUS GUESS, so the metric in panel (a) is chosen deliberately. The
intuitive quantity -- what share of rank-1 PSMs touch a non-canonical ORF -- points the WRONG way:
A549 is 15.1% against 4.0-4.7% for the three immunopeptidomes. That share is dominated by spurious
low-scoring matches into a 369k-sequence database and measures noise, not learnable signal. The
quantity that matters to a semi-supervised classifier is the density of CONFIDENT novel targets
inside that novel-touching population -- how many real examples exist to train on per unit of noise.
On that metric the substrates separate by two to three orders of magnitude, in the direction the
rescoring outcomes require.

  * Tryptic whole proteome: ~10 confident novel targets among ~64k novel-touching PSMs. Global
    rescoring optimises the canonical-dominated separation and discards the raw hyperscore signal
    that was distinguishing the handful of real matches, so class-specific FDR collapses 10 -> 1.
    Class-aware rescoring cannot bootstrap at all (single-class SVM folds at every train_fdr tried).
  * Immunopeptidome: nonspecific digestion leaves a far higher confident-target density, mokapot can
    learn the class, and rescoring both cleans the naive null (18 -> 4) and recovers the model
    (10 -> 26).

  (a) The controlling variable, measured on the CURRENT pgx searches.
  (b) The consequence, from the recorded rescoring runs.

PROVENANCE WARNING, and it is why panel (b) is drawn from constants. The rescoring experiments were
run in the OLD era (f0-threshold databases, pre-pgx, pre-frozen search parameters) and have NOT been
repeated on the released models. Recomputing panel (b) from today's search directories would
silently mix eras. The numbers are therefore transcribed from results.md, and the panel says so.
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

sys.path.insert(0, str(ROOT / "figures" / "S_fdr_rigor"))
from make_fdr_rigor import threshold  # noqa: E402  -- one definition of the FDR cut, not two

MODEL = "mamba4"
ARM = "null_atg"     # the arm with the largest non-canonical space -> the most generous density estimate
DATASETS = [("A549", "tryptic\nwhole proteome"), ("HBL1", "HLA-I\nimmunopeptidome"),
            ("SUDHL4", "HLA-I\nimmunopeptidome"), ("DoHH2", "HLA-I\nimmunopeptidome")]
SUB_COLOR = {"tryptic\nwhole proteome": "#c0392b", "HLA-I\nimmunopeptidome": "#2f6f9f"}

# Transcribed from results.md. OLD-era runs (f0-threshold DBs, pre-pgx); see the module docstring.
RESCORE = {
    "A549 (tryptic)": {"model_raw": 10, "model_rescored": 1, "null_raw": 13, "null_rescored": 6,
                       "note": "class-aware rescoring untrainable"},
    "HBL-1 (HLA-I)": {"model_raw": 10, "model_rescored": 26, "null_raw": 18, "null_rescored": 4,
                      "note": "mokapot learns the novel class"},
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent))
    a = ap.parse_args()
    out = Path(a.out)

    dens, missing = [], []
    for ds, sub in DATASETS:
        d = ROOT / f"proteogenomics/data/{ds}_pilot/pgx_{MODEL}/search/{ARM}"
        if not (d / "rank1.tsv.gz").exists():
            missing.append(f"{ds}/{ARM}")
            continue
        rows = load_rank1(d)
        n = len(rows)
        n_nov = sum(1 for (_, _, k, _) in rows if k == "novel_t")
        cut = threshold(rows, {"novel_t"}, {"novel_d"})
        n_conf = len({p for (p, h, k, _) in rows
                      if k == "novel_t" and cut is not None and h >= cut})
        dens.append({"dataset": ds, "substrate": sub, "n_psm": n, "n_novel_psm": n_nov,
                     "novel_share": n_nov / n if n else float("nan"),
                     "n_confident_novel": n_conf,
                     "confident_per_1k_novel_psm": 1000.0 * n_conf / n_nov if n_nov else float("nan")})
    if missing:
        print("MISSING (absent from the panel):", ", ".join(missing))

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4))

    # (a) how dense is the non-canonical class in the search?
    ax = axes[0]
    x = np.arange(len(dens))
    ax.bar(x, [d["confident_per_1k_novel_psm"] for d in dens],
           color=[SUB_COLOR[d["substrate"]] for d in dens], width=0.62)
    for i, d in enumerate(dens):
        ax.text(i, d["confident_per_1k_novel_psm"] * 1.12,
                f"{d['n_confident_novel']} confident\nof {d['n_novel_psm']:,} novel PSMs",
                ha="center", va="bottom", fontsize=7)
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels([d["dataset"] for d in dens], fontsize=9)
    ax.set_ylabel("confident novel targets / 1,000 novel PSMs", fontsize=9)
    ax.set_ylim(top=max(d["confident_per_1k_novel_psm"] for d in dens) * 9)
    ax.set_title(f"(a) the controlling variable: trainable signal in the\n"
                 f"non-canonical class (current pgx, {MODEL}, {ARM} arm)", fontsize=9)
    ax.grid(axis="y", alpha=0.25, which="both")
    handles = [plt.Line2D([], [], marker="s", ls="", color=c, label=s.replace("\n", " "))
               for s, c in SUB_COLOR.items()]
    ax.legend(handles=handles, fontsize=8, loc="upper left", framealpha=0.95)

    # (b) what rescoring then does
    ax = axes[1]
    labels = list(RESCORE)
    w = 0.19
    pos = np.arange(len(labels))
    series = [("model, raw", "model_raw", "#2f6f9f", ""),
              ("model, rescored", "model_rescored", "#2f6f9f", "//"),
              ("null, raw", "null_raw", "#9a9a9a", ""),
              ("null, rescored", "null_rescored", "#9a9a9a", "//")]
    for j, (lab, key, col, hatch) in enumerate(series):
        vals = [RESCORE[k][key] for k in labels]
        ax.bar(pos + (j - 1.5) * w, vals, width=w, color=col, hatch=hatch, edgecolor="white",
               label=lab)
        for p, v in zip(pos + (j - 1.5) * w, vals):
            ax.text(p, v + 0.5, str(v), ha="center", fontsize=7.5)
    ax.set_xticks(pos); ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("novel peptides at 1% class-specific FDR", fontsize=9)
    ax.set_ylim(0, 33)
    ax.set_title("(b) the consequence: rescoring collapses the tryptic case\n"
                 "and sharpens the immunopeptidome one", fontsize=9)
    ax.legend(fontsize=7.5, ncol=2, loc="upper left", framealpha=0.95)
    ax.grid(axis="y", alpha=0.25)
    ax.annotate("OLD-era runs (f0-threshold DBs, pre-pgx);\nnot repeated on the released models",
                xy=(0.5, -0.235), xycoords="axes fraction", ha="center", fontsize=7.5,
                color="#666666")

    fig.suptitle("Supplemental: semi-supervised rescoring is substrate-dependent", fontsize=11)
    fig.tight_layout(rect=[0, 0.035, 1, 0.93])
    for ext in ("png", "pdf"):
        fig.savefig(out / f"rescore_substrate.{ext}", dpi=200)
    (out / "rescore_substrate_values.json").write_text(json.dumps(
        {"model": MODEL, "arm": ARM, "density": dens, "rescore_recorded": RESCORE,
         "missing": missing}, indent=2))

    print(f"\n{'dataset':<9}{'substrate':<22}{'rank1':>10}{'novel PSM':>11}{'share':>8}"
          f"{'confident':>11}{'conf/1k':>10}")
    for d in dens:
        print(f"{d['dataset']:<9}{d['substrate'].replace(chr(10), ' '):<22}{d['n_psm']:>10,}"
              f"{d['n_novel_psm']:>11,}{d['novel_share'] * 100:>7.2f}%"
              f"{d['n_confident_novel']:>11}{d['confident_per_1k_novel_psm']:>10.2f}")
    print(f"\nwrote {out}/rescore_substrate.png|.pdf|_values.json")


if __name__ == "__main__":
    main()
