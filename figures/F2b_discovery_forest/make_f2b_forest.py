#!/usr/bin/env python3
"""Main Fig 2b -- is the model's discovery advantage a REPRODUCIBLE PATTERN or a one-off?

FIGURES_PLAN names this the key weakness of Fig 2: the immunopeptidome case originally rested on a
single dataset (HBL-1) with small counts, and a brief communication needs a pattern. The fix
specified there is a forest plot of the model/null discovery-rate ratio across datasets. This is
that panel, on 5 datasets x 2 released models x 2 threshold arms.

THE METRIC IS A RATE RATIO, NOT A COUNT, and that choice is the whole point. Raw novel-peptide count
rewards whichever arm shipped the bigger database, because a bigger database absorbs more spectra;
`null_atg` wins on raw count in 2 of 5 datasets while carrying 25-180x more sequences. Discovery
density (novel peptides per 1,000 database sequences) asks the question that actually matters for
search-space selection: per sequence you commit to searching, how much do you find? The ratio
model/null then makes datasets with very different absolute yields comparable on one axis.

WHY A LOG AXIS AND WHY NO ERROR BARS. The ratios span roughly 20x to 2,700x, so linear is unreadable.
Error bars are deliberately omitted rather than invented: these are single-search point estimates,
not replicated measurements, and the honest uncertainty statement is the peptide count printed
beside each point (single digits to low double digits in every case). A confidence interval computed
from a Poisson assumption on those counts would imply a replication structure that does not exist.

Usage: make_f2b_forest.py [--out <dir>]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
# (label, substrate). Order = tryptic first, then the four HLA-I sets.
DATASETS = [("A549", "tryptic"), ("HBL1", "HLA-I"), ("SUDHL4", "HLA-I"),
            ("DoHH2", "HLA-I"), ("THP1", "HLA-I")]
MODELS = ["mamba4", "attn"]
ARMS = ["model_poisson", "model_standard"]
STYLE = {("mamba4", "model_poisson"): ("#1b4f72", "o", "mamba4 Poisson"),
         ("attn", "model_poisson"): ("#2f6f9f", "o", "attn Poisson"),
         ("mamba4", "model_standard"): ("#8e6b3f", "s", "mamba4 theta=1"),
         ("attn", "model_standard"): ("#c1934a", "s", "attn theta=1")}


def density(rec):
    n = rec.get("db_novel_seqs", 0)
    return (1000.0 * rec["novel_peptides"] / n) if n else float("nan")


def load(ds, model):
    for nm in ("table_all.json", "table.json"):
        f = ROOT / f"proteogenomics/data/{ds}_pilot/pgx_{model}/report/{nm}"
        if f.exists():
            return json.loads(f.read_text())["arms"]
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent))
    a = ap.parse_args()
    out = Path(a.out)

    recs, missing = [], []
    for ds, sub in DATASETS:
        for model in MODELS:
            arms = load(ds, model)
            if arms is None or "null_atg" not in arms:
                missing.append(f"{ds}/{model}")
                continue
            dn = density(arms["null_atg"])
            for arm in ARMS:
                if arm not in arms:
                    missing.append(f"{ds}/{model}/{arm}")
                    continue
                dm = density(arms[arm])
                if not (dn > 0 and dm == dm):
                    continue
                recs.append({"dataset": ds, "substrate": sub, "model": model, "arm": arm,
                             "ratio": dm / dn, "model_density": dm, "null_density": dn,
                             "novel_peptides": arms[arm]["novel_peptides"],
                             "db_novel_seqs": arms[arm]["db_novel_seqs"],
                             "null_peptides": arms["null_atg"]["novel_peptides"],
                             "null_seqs": arms["null_atg"]["db_novel_seqs"]})
    if missing:
        print("MISSING (absent from the figure, not imputed):")
        for m in missing:
            print(f"  {m}")

    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    ypos, ylabels = [], []
    y = 0.0
    for ds, sub in DATASETS:
        block = [r for r in recs if r["dataset"] == ds]
        if not block:
            continue
        for j, r in enumerate(block):
            col, mark, _ = STYLE[(r["model"], r["arm"])]
            yy = y + j * 0.17
            ax.plot(r["ratio"], yy, mark, ms=7, color=col, zorder=3)
            if r["arm"] == "model_poisson" and r["model"] == "attn":
                ax.annotate(f"{r['novel_peptides']} pept / {r['db_novel_seqs']:,} seqs",
                            xy=(r["ratio"], yy), xytext=(6, 0), textcoords="offset points",
                            fontsize=7, va="center", color="#444444")
        ypos.append(y + 0.17 * (len(block) - 1) / 2)
        ylabels.append(f"{ds}\n({sub})")
        y += 1.0

    ax.axvline(1.0, color="k", ls="--", lw=1.2)
    ax.text(1.0, len(ylabels) - 0.55, "  no advantage over\n  the naive null", fontsize=8, va="bottom")
    ax.set_xscale("log")
    ax.set_yticks(ypos); ax.set_yticklabels(ylabels, fontsize=9)
    ax.set_ylim(-0.4, len(ylabels) - 0.35)
    ax.invert_yaxis()
    ax.set_xlabel("discovery-density ratio, model / naive AUG null\n"
                  "(novel peptides per 1,000 database sequences, 1% class-specific FDR)", fontsize=9)
    ax.grid(axis="x", alpha=0.3, which="both")
    # Headroom on the right so the annotations are not clipped, and the legend on the LEFT where the
    # panel is empty -- every point sits above 10x, so nothing is occluded there. A lower-right
    # legend covered the THP-1 row and truncated its label.
    # The x=1 reference line MUST stay in view: an axis starting at 3.6 hid the line while
    # keeping its "no advantage" label, which reads as if the line were somewhere off-panel.
    lo = 0.6
    hi = max(r["ratio"] for r in recs) * 14.0
    ax.set_xlim(lo, hi)
    handles = [plt.Line2D([], [], marker=m, ls="", color=c, label=lab)
               for (c, m, lab) in STYLE.values()]
    ax.legend(handles=handles, fontsize=8, loc="upper left", framealpha=0.95)
    ax.set_title("Fig 2b -- the model selects a denser search space than naive enumeration,\n"
                 "reproducibly across 5 datasets and both released models", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out / f"F2b_discovery_forest.{ext}", dpi=200)
    (out / "F2b_values.json").write_text(json.dumps({"records": recs, "missing": missing}, indent=2))

    print(f"\n{'dataset':<8}{'model':<8}{'arm':<16}{'model dens':>11}{'null dens':>11}{'ratio':>10}")
    for r in recs:
        print(f"{r['dataset']:<8}{r['model']:<8}{r['arm']:<16}{r['model_density']:>11.3f}"
              f"{r['null_density']:>11.3f}{r['ratio']:>10.1f}x")
    rr = [r["ratio"] for r in recs]
    print(f"\n  n = {len(rr)} points, ALL above 1.0: {all(x > 1 for x in rr)}")
    print(f"  median ratio {np.median(rr):.1f}x, range {min(rr):.1f}x to {max(rr):.1f}x")
    print(f"\nwrote {out}/F2b_discovery_forest.png|.pdf|F2b_values.json")


if __name__ == "__main__":
    main()
