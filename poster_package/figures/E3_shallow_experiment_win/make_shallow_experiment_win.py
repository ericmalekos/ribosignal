#!/usr/bin/env python3
"""Figure E3: the standalone model out-recalls a real, shallow Ribo-seq experiment on novel ORFs and uORFs.

The claim this figure has to earn is a practical one, not a leaderboard one. Deep ribosome profiling
beats the model and the figure says so. The point is that the realistic alternative most labs face is
not a deep library, it is a SHALLOW one -- and against that, a model needing no ribosome profiling at
all recovers more non-canonical ORFs.

Panel A establishes that Wang is a typical shallow library rather than a strawman: pooled P-sites
against 3-nt periodicity for the three mouse-liver datasets, all pooled through one pipeline onto one
shared 22,974-tx universe. Wang sits lowest on both axes.

Panel B is the result. For novel ORFs and uORFs, against each of the two DEEPER references, the recall
achieved by the shallow experiment (Wang) and by the standalone model (`pred_preddepth`, no Ribo-seq
for the query sample at any point). The arrow is the margin.

Deliberately NOT shown as a win, and stated in the caption instead: the model loses to both deeper
experiments on recall, and to all three on precision and therefore F1. Reporting only the favourable
comparison would be the same error as quoting non-canonical F1 against 1.0 instead of against the
between-experiment ceiling.

All numbers read from results/mouse_liver_3x3_canon/scored/, so this regenerates if the factorial is re-run.
cas12a env.
"""
from __future__ import annotations

import csv
import collections
import glob
import json
import pathlib
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model")
S = ROOT / "results/mouse_liver_3x3_canon/scored"
OUT = ROOT / "figures/E3_shallow_experiment_win"
EXP_C, MOD_C = "#1F4E79", "#C0392B"          # observed / predicted, house colours
POOLS = {"janich": "pool_janich", "gse243134": "pool_gse243134", "wang": "pool_wang"}
BAMS = {"janich": "mouse_janich_liver_ribo", "gse243134": "mouse_gse243134_liver",
        "wang": "mouse_wang_liver_ribo"}


def library_stats():
    """Pooled P-sites and mean f0 periodicity per dataset, from the packs and metaplots configs."""
    ribo_only = set()
    rr = ROOT / "data/liver3x3/gse243134_ribo_runs.txt"
    if rr.exists():
        ribo_only = set(rr.read_text().split())
    out = {}
    for k, pool in POOLS.items():
        ps = int(np.load(ROOT / f"data/liver3x3/{pool}/target_counts.npy").sum())
        f0 = []
        for cfg in glob.glob(str(ROOT / f"data/liver3x3/{pool}/_work/psites/*_pre_config.txt")):
            for line in open(cfg):
                m = re.match(r"#\s+\d+\s+[\d.]+%\s+\d+\s+\d+\s+\d+\s+\d+\s+([\d.]+)%", line)
                if m:
                    f0.append(float(m.group(1)))
        logs = sorted((ROOT / f"data/heldout_bam/{BAMS[k]}").glob("*Log.final.out"))
        if k == "gse243134" and ribo_only:
            logs = [x for x in logs if x.name.split(".")[0] in ribo_only]
        out[k] = dict(psites=ps, f0=float(np.mean(f0)) if f0 else float("nan"), n=len(logs))
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lib = library_stats()

    ceil = collections.defaultdict(dict)
    for r in csv.DictReader(open(S / "ribo_vs_ribo_by_class.tsv"), delimiter="\t"):
        ceil[(r["reference"], r["orf_class"])][r["compared"]] = dict(
            precision=float(r["precision"]), recall=float(r["recall"]), f1=float(r["f1"]))
    mrows = list(csv.DictReader(open(S / "model_vs_ribo_by_class.tsv"), delimiter="\t"))

    # ONE FIXED CELL PER ROW, NOT AN ARGMAX.
    #
    # This was `max(cand, key=recall)` over a candidate list filtered only by class, reference and
    # arm -- so it maximised over 6 cells (3 RNA inputs x 2 architectures) while the shallow
    # experiment it is compared against is a single fixed value with no equivalent best-of. That is
    # asymmetric in the model's favour, and the `rna_input` it recorded was whichever cell won
    # rather than a stated choice.
    #
    # Fixed to attn / gse243134 RNA / pred_preddepth. Measured 2026-08-15, all 12 attn cells:
    #
    #   class  ref         gse243134   janich    wang
    #   novel  gse243134     0.6953    0.6953   0.6833
    #   novel  janich        0.7188    0.7121   0.6955
    #   uORF   gse243134     0.6844    0.6571   0.2869
    #   uORF   janich        0.7113    0.6933   0.3165
    #
    # THE NUMBERS DO NOT CHANGE. gse243134 RNA is the argmax cell in three rows, and in the fourth
    # (novel/gse243134) it TIES janich exactly at 0.6953 (tp 461 of 663 both ways) -- `max` had
    # simply returned whichever it saw first. So the rule becomes stateable at zero cost to the
    # result.
    #
    # THE SPREAD IS NOT NOISE FOR uORF, WHICH IS WHY THE ARGMAX LOOKED HARMLESS. Wang RNA collapses
    # uORF recall to 0.29-0.32 against ~0.68-0.71 for the other two inputs -- a >0.38 gap, two
    # orders of magnitude beyond E1's ~0.0041 mean RNA-provenance effect on F1. The argmax was not
    # selecting noise there; it was steering away from a genuinely bad cell. For `novel` the
    # spread IS noise (0.6833-0.6953).
    #
    # WHY THE WANG RNA CELL IS BAD IS *NOT* DEPTH. This comment previously read "Wang is the
    # shallowest RNA arm (2 samples), and uORFs sit in 5'UTRs where coverage is thinnest". The
    # sample count is right and the depth claim is BACKWARDS. Measured 2026-08-16 on the shared
    # 22,974-tx universe, so the three are directly comparable:
    #
    #   arm        n_lib   RNA coverage total   mean/nt   tx>0    top-100 share
    #   gse243134     19        3,242,315,370     57.54   22503           31.1%
    #   janich         7        6,853,032,580    121.62   22569           30.7%
    #   wang           2        6,884,416,652    122.17   22605           39.6%
    #
    # Wang RNA is the DEEPEST arm by total and by mean/nt, and covers the most transcripts. What
    # separates it is EVENNESS: its top 100 transcripts hold 39.6% of all coverage against ~31%
    # for the other two, from only 2 libraries. So the uORF collapse is a concentration/replicate
    # effect, not a depth one, and the precise mechanism is NOT established by these numbers --
    # state it as unexplained rather than substituting a second guess for the first.
    #
    # Note also that n_lib and coverage run in OPPOSITE directions here (19 libraries buys the
    # LEAST coverage). Library count is not a depth proxy on this dataset.
    FIX_MODEL, FIX_RNA = "attn", "gse243134"
    pairs = []
    # `annotated` added 2026-08-16 at the poster session's request. Same fixed cell and the same
    # prf class convention as the other rows -- it is the high-consensus class, and showing it beside
    # novel/uORF is what makes the non-canonical gap legible rather than assertable.
    for cls in ("annotated", "novel", "uORF"):
        for ref in ("gse243134", "janich"):
            cand = [r for r in mrows if r["orf_class"] == cls and r["ribo_reference"] == ref
                    and r["arm"] == "pred_preddepth"
                    and r["model"] == FIX_MODEL and r["rna_input"] == FIX_RNA]
            if not cand or "wang" not in ceil[(ref, cls)]:
                continue
            if len(cand) != 1:
                raise SystemExit(f"expected exactly 1 cell for {cls}/{ref}, got {len(cand)}")
            b = cand[0]
            pairs.append(dict(orf_class=cls, reference=ref,
                              exp_recall=ceil[(ref, cls)]["wang"]["recall"],
                              model_recall=float(b["recall"]),
                              model=b["model"], rna_input=b["rna_input"],
                              n_ref=int(b["n_ref"]), cell_selection="fixed (not argmax)"))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.2, 4.3),
                                   gridspec_kw=dict(width_ratios=[1.0, 1.35]))

    # ---- Panel A: the three libraries ----
    for k, v in lib.items():
        shallow = (k == "wang")
        axA.scatter(v["psites"] / 1e6, v["f0"], s=190 if shallow else 130,
                    c=MOD_C if shallow else EXP_C, zorder=3,
                    edgecolor="white", linewidth=1.4)
        axA.annotate(f"{k}\n{v['n']} samples", (v["psites"] / 1e6, v["f0"]),
                     textcoords="offset points", xytext=(0, -30 if shallow else 14),
                     ha="center", fontsize=8.5,
                     color=MOD_C if shallow else EXP_C, fontweight="bold" if shallow else "normal")
    axA.set_xlabel("pooled P-sites (millions)")
    axA.set_ylabel("3-nt periodicity, mean f0 (%)")
    axA.set_title("A  The three mouse-liver libraries", loc="left", fontsize=10, fontweight="bold")
    axA.grid(alpha=0.25, linewidth=0.6)
    axA.margins(0.22)

    # ---- Panel B: shallow experiment vs standalone model ----
    labels, y = [], []
    for i, p in enumerate(pairs):
        yy = len(pairs) - 1 - i
        y.append(yy)
        labels.append(f"{p['orf_class']}\nvs {p['reference']}")
        axB.plot([p["exp_recall"], p["model_recall"]], [yy, yy], color="0.75", lw=2.2, zorder=1)
        axB.scatter(p["exp_recall"], yy, s=95, c=EXP_C, zorder=3, edgecolor="white", linewidth=1.2)
        axB.scatter(p["model_recall"], yy, s=95, c=MOD_C, zorder=3, edgecolor="white", linewidth=1.2)
        axB.annotate(f"+{p['model_recall']-p['exp_recall']:.3f}",
                     ((p["exp_recall"] + p["model_recall"]) / 2, yy),
                     textcoords="offset points", xytext=(0, 9), ha="center",
                     fontsize=8.5, color=MOD_C, fontweight="bold")
    axB.set_yticks(y)
    axB.set_yticklabels(labels, fontsize=8.5)
    axB.set_xlabel("recall of the reference's calls")
    axB.set_title("B  Shallow experiment vs standalone model (no Ribo-seq)",
                  loc="left", fontsize=10, fontweight="bold")
    axB.grid(axis="x", alpha=0.25, linewidth=0.6)
    # Extra headroom BELOW the lowest row: an in-axes legend at lower right sat on top of the
    # `uORF vs janich` dumbbell, which reaches x=0.688.
    axB.margins(x=0.16)
    axB.set_ylim(-0.95, len(pairs) - 0.4)
    axB.scatter([], [], s=95, c=EXP_C, label="Wang (2 samples, shallow, 75.4% f0)")
    axB.scatter([], [], s=95, c=MOD_C, label="model, pred_preddepth (no Ribo-seq)")
    axB.legend(frameon=False, fontsize=8, loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, -0.02), handletextpad=0.4, columnspacing=1.2)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"E3_shallow_experiment_win.{ext}", dpi=200, bbox_inches="tight")

    vals = dict(library=lib, pairs=pairs,
                note=("Standalone model recall vs the shallow Wang experiment, on the two deeper "
                      "references. Model loses to both deeper experiments on recall and to all "
                      "three on precision/F1; see FIGURE_DATA_INPUTS.md."))
    (OUT / "E3_values.json").write_text(json.dumps(vals, indent=2))
    with open(OUT / "E3_values.tsv", "w") as fh:
        fh.write("orf_class\treference\texp_recall_wang\tmodel_recall\tmargin\tmodel\trna_input\tn_ref\n")
        for p in pairs:
            fh.write(f"{p['orf_class']}\t{p['reference']}\t{p['exp_recall']:.4f}\t"
                     f"{p['model_recall']:.4f}\t{p['model_recall']-p['exp_recall']:+.4f}\t"
                     f"{p['model']}\t{p['rna_input']}\t{p['n_ref']}\n")
    print(f"  wrote {OUT}/E3_shallow_experiment_win.{{pdf,png}} + values")
    for p in pairs:
        print(f"    {p['orf_class']:<6} vs {p['reference']:<10} wang={p['exp_recall']:.3f} "
              f"model={p['model_recall']:.3f} ({p['model']}) margin={p['model_recall']-p['exp_recall']:+.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
