#!/usr/bin/env python3
"""D12 -- model-selected ORFs vs the CPAT / CPC2 coding-potential baseline (locked decision D3).

WHY THIS IS THE FIGURE A REVIEWER ASKS FOR. Every other comparison in this study scores the model
against naive enumeration (`null_atg`), which only establishes "better than listing every ORF".
CPAT and CPC2 are the standard sequence-only coding-potential selectors and they shrink the search
space by the same order of magnitude as the model, using sequence composition alone -- no RNA-seq,
no translation model. They are the honest competitor.

All four arms enumerate from EXACTLY the same candidate pool as `null_atg` (enforced by the
regression test `test_coding_potential_pool_matches_null_arm`), so the only thing that differs
between them is which ORFs are kept. That is what makes this a controlled comparison of selection
rules rather than of pipelines.

  (a) Novel peptides found, per dataset. The raw deliverable.
  (b) Discovery efficiency: novel peptides per 1,000 database sequences. Normalises out the fact
      that the arms build databases of different sizes.
  (c) Canonical cost: change in GENCODE peptides against the gencode-only baseline. A selection rule
      that finds novel peptides by displacing canonical ones has not gained anything.

THE RESULT SPLITS BY ASSAY, and that split is the finding. Comparing the four SELECTION RULES on raw
novel-peptide count: CPAT/CPC2 win the tryptic whole proteome (25 / 23 vs the model's 11), the model
wins all three HLA-I immunopeptidomes (14 vs 10, 22 vs 14, 32 vs 8). CPAT and CPC2 score a TRANSCRIPT's
coding potential from sequence composition, so they keep long, codon-biased, ORF-like sequences --
the population a tryptic digest samples well and the least novel population there is. The model
scores PER-NUCLEOTIDE translation from the sample's own RNA-seq, so it keeps short, non-canonical,
cell-type-specific ORFs, which is what HLA-I presentation actually samples.

STATED PLAINLY BECAUSE THE FIGURE SHOWS IT: on RAW count the unselected `null_atg` arm beats every
selection rule on HBL-1 (20) and SU-DHL-4 (28). Keeping every candidate does find more peptides -- at
25x to 30x the database size, and, on A549, at a cost of 2,894 canonical peptides. Panel (b) is the
honest comparison for that reason, and there the model's Poisson arm leads on all four datasets
(4.6 / 7.6 / 6.2 / 7.2 per 1,000 sequences vs 0.03-0.11 for the null).

Usage: make_cpat_cpc2.py [--model mamba4] [--out <dir>]
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
DATASETS = [("A549", "tryptic"), ("HBL1", "HLA-I"), ("SUDHL4", "HLA-I"), ("DoHH2", "HLA-I")]
ARMS = ["null_atg", "cpat", "cpc2", "model_standard", "model_poisson"]
LABEL = {"null_atg": "null (naive AUG)", "cpat": "CPAT", "cpc2": "CPC2",
         "model_standard": "model (theta=1)", "model_poisson": "model (Poisson)"}
COLOR = {"null_atg": "#9a9a9a", "cpat": "#7b5ea7", "cpc2": "#c07ab8",
         "model_standard": "#2f6f9f", "model_poisson": "#1b4f72"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="mamba4")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent))
    a = ap.parse_args()
    out = Path(a.out)

    tab, missing = {}, []
    for ds, sub in DATASETS:
        f = ROOT / f"proteogenomics/data/{ds}_pilot/pgx_{a.model}/report/table_all.json"
        if not f.exists():
            missing.append(f"{ds}: {f}")
            continue
        j = json.loads(f.read_text())
        row = {}
        for arm in ARMS:
            r = j["arms"].get(arm)
            if r is None:
                missing.append(f"{ds}/{arm}")
                continue
            n_db = r["db_novel_seqs"]
            row[arm] = {
                "novel_peptides": r["novel_peptides"],
                "db_novel_seqs": n_db,
                "per_1k": 1000.0 * r["novel_peptides"] / n_db if n_db else float("nan"),
                "d_gencode_peptides": r["d_gencode_peptides"],
            }
        tab[ds] = {"substrate": sub, "arms": row}
    if missing:
        print("MISSING (absent from the figure, not zeroed):")
        for m in missing:
            print(f"  {m}")

    dss = [d for d, _ in DATASETS if d in tab]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.5))
    w = 0.16
    pos = np.arange(len(dss))

    def grouped(ax, key, ylabel, title, logy=False):
        for j, arm in enumerate(ARMS):
            vals = [tab[d]["arms"].get(arm, {}).get(key, np.nan) for d in dss]
            ax.bar(pos + (j - 2) * w, vals, width=w, color=COLOR[arm], label=LABEL[arm])
        ax.set_xticks(pos)
        ax.set_xticklabels([f"{d}\n({tab[d]['substrate']})" for d in dss], fontsize=8.5)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(title, fontsize=9)
        ax.grid(axis="y", alpha=0.25)
        if logy:
            ax.set_yscale("log")

    grouped(axes[0], "novel_peptides", "novel peptides (1% class-specific FDR)",
            "(a) novel peptides found\n* = best SELECTIVE arm (null is the no-selection reference)")
    # The asterisk ranks the four SELECTION RULES only. null_atg keeps every candidate, so it is the
    # reference the rules are trying to beat per unit of database, not a competitor on raw count --
    # and on raw count it does win HBL-1 and SU-DHL-4, which panel (b) is there to put in context.
    for i, d in enumerate(dss):
        arms = tab[d]["arms"]
        sel = [k for k in arms if k != "null_atg"]
        best = max(sel, key=lambda k: arms[k]["novel_peptides"])
        v = arms[best]["novel_peptides"]
        axes[0].text(pos[i] + (ARMS.index(best) - 2) * w, v * 1.03, "*", ha="center", va="bottom",
                     fontsize=14, color=COLOR[best])
    # headroom for the asterisks and the legend, or both clip
    axes[0].set_ylim(0, max(r["novel_peptides"] for d in dss for r in tab[d]["arms"].values()) * 1.45)
    axes[0].legend(fontsize=7.5, ncol=2, loc="upper right", framealpha=0.95)

    grouped(axes[1], "per_1k", "novel peptides per 1,000 DB sequences",
            "(b) discovery efficiency\nnormalises out database size", logy=True)

    ax = axes[2]
    for j, arm in enumerate(ARMS):
        vals = [tab[d]["arms"].get(arm, {}).get("d_gencode_peptides", np.nan) for d in dss]
        ax.bar(pos + (j - 2) * w, vals, width=w, color=COLOR[arm], label=LABEL[arm])
    ax.axhline(0, color="k", lw=1)
    ax.set_xticks(pos)
    ax.set_xticklabels([f"{d}\n({tab[d]['substrate']})" for d in dss], fontsize=8.5)
    ax.set_ylabel("change in GENCODE peptides vs gencode-only", fontsize=9)
    # A549's null arm (-2,894) is ~40x every other bar; on a linear axis it flattens the panel to a
    # single spike and the selective arms become indistinguishable from zero. symlog keeps the
    # magnitude visible AND resolves the small differences that decide between the selective arms.
    ax.set_yscale("symlog", linthresh=10)
    ax.set_title("(c) canonical cost (symlog)\nnegative = canonical identifications displaced",
                 fontsize=9)
    ax.grid(axis="y", alpha=0.25, which="both")

    fig.suptitle(f"D12 -- model-selected ORFs vs CPAT / CPC2 coding potential "
                 f"(released {a.model}, frozen search params)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    for ext in ("png", "pdf"):
        fig.savefig(out / f"D12_cpat_cpc2.{ext}", dpi=200)
    (out / "D12_values.json").write_text(json.dumps(
        {"model": a.model, "datasets": tab, "missing": missing}, indent=2))

    print(f"\n{'dataset':<9}{'arm':<17}{'novel pept':>11}{'DB seqs':>10}{'per 1k':>9}{'dGENCODE':>10}")
    for d in dss:
        for arm in ARMS:
            r = tab[d]["arms"].get(arm)
            if r is None:
                continue
            print(f"{d:<9}{arm:<17}{r['novel_peptides']:>11}{r['db_novel_seqs']:>10,}"
                  f"{r['per_1k']:>9.2f}{r['d_gencode_peptides']:>+10}")
    print(f"\nwrote {out}/D12_cpat_cpc2.png|.pdf|D12_values.json")


if __name__ == "__main__":
    main()
