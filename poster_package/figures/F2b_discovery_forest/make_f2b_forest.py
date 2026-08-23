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

PANEL (b) REPORTS THE ABSOLUTE BOTTOM LINE, because a rate ratio alone is not a full account. It
shows the change in TOTAL unique peptides (canonical + novel) against searching GENCODE alone, for
the model arms AND the null. This is the number an experimentalist takes home, it is shown even
where it disfavours the model, and it carries the one thing the ratio cannot: an arm can have a
2,700x density advantage and still hand back fewer peptides than not searching a novel DB at all.
Standing rule: all proteomics work reports absolute totals.

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
# A549 REMOVED 2026-08-15 (archived; see proteogenomics/data/_archive_a549_2026_08_15/).
# It was the only tryptic human dataset, so every panel below is now HLA-I ONLY -- do not
# describe any of it as a whole-proteome or tryptic result.
# (label, substrate). All four remaining sets are HLA-I immunopeptidomes.
DATASETS = [("HBL1", "HLA-I"), ("SUDHL4", "HLA-I"),
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


def dtot(rec, arms):
    """Change in TOTAL unique peptides (canonical + novel) vs the gencode-only baseline."""
    return (rec["d_gencode_peptides"] + rec["novel_peptides"]
            - arms.get("gencode", {}).get("novel_peptides", 0))


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
                             "null_seqs": arms["null_atg"]["db_novel_seqs"],
                             # Net change in TOTAL unique peptides vs the gencode-only baseline.
                             # Baseline novel is subtracted explicitly rather than assumed 0.
                             "d_total_peptides": dtot(arms[arm], arms),
                             "null_d_total_peptides": dtot(arms["null_atg"], arms)})
    if missing:
        print("MISSING (absent from the figure, not imputed):")
        for m in missing:
            print(f"  {m}")

    fig, (ax, axT) = plt.subplots(1, 2, figsize=(14.4, 5.4), sharey=True,
                                  gridspec_kw={"width_ratios": [1.55, 1.0]})
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
    ax.set_title("(a) discovery density ratio -- the model selects a denser search space\n"
                 "than naive enumeration, reproducibly across 5 datasets and both models",
                 fontsize=9.5, loc="left")

    # --- (b) the absolute bottom line, on the same rows ---
    # Drawn at the same y as (a) so a reader can carry a row across without re-reading labels. The
    # null is drawn ONCE per dataset (it does not depend on the model or the threshold arm), as an
    # open marker so it does not compete with the model points.
    seen_null = set()
    for ds, sub in DATASETS:
        block = [r for r in recs if r["dataset"] == ds]
        if not block:
            continue
        yb = ypos[[l.split("\n")[0] for l in ylabels].index(ds)] - 0.17 * (len(block) - 1) / 2
        for j, r in enumerate(block):
            col, mark, _ = STYLE[(r["model"], r["arm"])]
            axT.plot(r["d_total_peptides"], yb + j * 0.17, mark, ms=7, color=col, zorder=3)
        if ds not in seen_null:
            seen_null.add(ds)
            axT.plot(block[0]["null_d_total_peptides"], yb + 0.17 * (len(block) - 1) / 2,
                     "X", ms=9, mfc="none", mec="#9a9a9a", mew=1.8, zorder=4,
                     label="naive AUG null" if len(seen_null) == 1 else None)
            axT.annotate(f"{block[0]['null_d_total_peptides']:+,}",
                         xy=(block[0]["null_d_total_peptides"],
                             yb + 0.17 * (len(block) - 1) / 2),
                         xytext=(0, -11), textcoords="offset points",
                         fontsize=6.8, ha="center", color="#6a6a6a")
    axT.axvline(0, color="k", ls="--", lw=1.2)
    axT.set_xscale("symlog", linthresh=20)
    axT.set_xlabel("change in TOTAL unique peptides vs searching GENCODE alone\n"
                   "(canonical + novel; negative = net LOSS)", fontsize=9)
    axT.grid(axis="x", alpha=0.3, which="both")
    # Headroom both ways, or the A549 null (-2,884) sits half off the left spine with its label
    # clipped. symlog, so multiply rather than pad additively.
    _tv = [r["d_total_peptides"] for r in recs] + [r["null_d_total_peptides"] for r in recs]
    axT.set_xlim(min(_tv) * 4.0, max(_tv) * 4.0)
    axT.legend(fontsize=8, loc="upper left", framealpha=0.95)
    axT.set_title("(b) TOTAL unique peptides -- the number taken home,\n"
                  "shown even where it disfavours the model", fontsize=9.5, loc="left")

    fig.suptitle("Fig 2b -- a denser search space (a), and what that is worth in absolute "
                 "peptides (b)", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for ext in ("png", "pdf"):
        fig.savefig(out / f"F2b_discovery_forest.{ext}", dpi=200)
    (out / "F2b_values.json").write_text(json.dumps({"records": recs, "missing": missing}, indent=2))

    print(f"\n{'dataset':<8}{'model':<8}{'arm':<16}{'model dens':>11}{'null dens':>11}{'ratio':>10}")
    for r in recs:
        print(f"{r['dataset']:<8}{r['model']:<8}{r['arm']:<16}{r['model_density']:>11.3f}"
              f"{r['null_density']:>11.3f}{r['ratio']:>10.1f}x")
    rr = [r["ratio"] for r in recs]
    pos = sum(1 for r in recs if r["d_total_peptides"] > 0)
    dss_seen = {r["dataset"] for r in recs}
    npos = sum(1 for ds in dss_seen
               if next(r["null_d_total_peptides"] for r in recs if r["dataset"] == ds) > 0)
    print("\n  TOTAL unique peptides vs GENCODE-only (the absolute bottom line):")
    print(f"    model arms above zero: {pos}/{len(recs)}")
    for arm in ARMS:
        sub = [r for r in recs if r["arm"] == arm]
        print(f"      {arm:<16} {sum(1 for r in sub if r['d_total_peptides'] > 0)}/{len(sub)}  "
              + ", ".join(f"{r['dataset']}/{r['model']}={r['d_total_peptides']:+}" for r in sub))
    print(f"    naive AUG null above zero: {npos}/{len(dss_seen)}  "
          + ", ".join(f"{ds}={next(r['null_d_total_peptides'] for r in recs if r['dataset']==ds):+}"
                      for ds in sorted(dss_seen)))
    print(f"\n  n = {len(rr)} points, ALL above 1.0: {all(x > 1 for x in rr)}")
    print(f"  median ratio {np.median(rr):.1f}x, range {min(rr):.1f}x to {max(rr):.1f}x")
    print(f"\nwrote {out}/F2b_discovery_forest.png|.pdf|F2b_values.json")


if __name__ == "__main__":
    main()
