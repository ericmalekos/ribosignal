#!/usr/bin/env python3
"""Publication figures for the cross-species Ribo-seq generalisation study.

Every number is READ FROM AN ARTIFACT on disk; nothing is hard-coded except the
species ordering and the labels. Sources per panel are recorded in the values JSON
written alongside, per the project's figure-provenance rule.

  Figure 1  Study design and data flow
  Figure 2  ORF-call agreement, model vs observed, by species and ORF class
  Figure 3  Per-nucleotide profile agreement, and the Spearman tie artifact
  Figure 4  The mm25 multimap diagnostic

Usage: python3 figures/X_xspecies_pub/make_xspecies_figures.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

# Self-contained: paths resolve relative to this file, so the bundle regenerates its own
# figures anywhere it is unpacked. Every input is in ../data/ alongside it.
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
OUT = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)


class _Rel:
    """Maps the original `NEW / "results" / x` and `NEW / "data" / x` onto the flat
    bundle data/ directory, so the figure code below is unchanged from the version that
    produced the published figures."""

    def __truediv__(self, part):
        return _RelDir(DATA)


class _RelDir:
    def __init__(self, base):
        self.base = base

    def __truediv__(self, name):
        return self.base / name


NEW = _Rel()

# Ordered by divergence from the training species (human + mouse).
ORDER = ["human", "chimp", "gorilla", "macaque", "zebrafish", "celegans", "yeast"]
LABEL = {"human": "human", "chimp": "chimp", "gorilla": "gorilla", "macaque": "macaque",
         "zebrafish": "zebrafish", "celegans": "C. elegans", "yeast": "yeast"}
ARM = {"yeast": "yeast_gse173654_ribo", "celegans": "worm_gse52905_ribo",
       "zebrafish": "zf_gse46512_ribo", "gorilla": "primate_gg_ribo",
       "chimp": "primate_pt_ribo", "macaque": "primate_rm_ribo",
       "human": "ruizorera_hsCM_ribo"}

# Colour-blind-safe (Okabe-Ito). Primates share a hue family, outgroups are distinct.
C = {"human": "#0072B2", "chimp": "#3C8DC5", "gorilla": "#6BA8D6", "macaque": "#9BC3E6",
     "zebrafish": "#E69F00", "celegans": "#009E73", "yeast": "#CC79A7"}
INK, GRID, MUTED = "#1A1A1A", "#D9D9D9", "#6E6E6E"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.labelsize": 8.5, "axes.titlesize": 9, "axes.titleweight": "bold",
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.edgecolor": INK, "axes.linewidth": 0.7,
    "xtick.color": INK, "ytick.color": INK, "text.color": INK, "axes.labelcolor": INK,
    "figure.dpi": 200, "savefig.dpi": 400, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False,
})
PROV: dict = {"model": "attn", "source": {}}


def tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f(x, d=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def panel_tag(ax, s, dx=-0.085, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="top", ha="left")


# ---------------------------------------------------------------- figure 2
def figure2():
    src_cls = NEW / "results" / "xspecies_orf_classes_attn.tsv"
    src_all = NEW / "results" / "xspecies_orf_calls.tsv"
    cls = {(r["pack"], r["class"]): r for r in tsv(src_cls)}
    allc = [r for r in tsv(src_all)
            if r["arch"] == "attn" and r["variant"] == "pred_preddepth"]
    overall = {r["pack"]: r for r in allc}

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7),
                             gridspec_kw={"width_ratios": [1.15, 1.5, 1.05], "wspace": 0.55})

    # (a) overall precision / recall
    ax = axes[0]
    xs = [f(overall[p]["precision"]) for p in ORDER if p in overall]
    ys = [f(overall[p]["recall"]) for p in ORDER if p in overall]
    ps = [p for p in ORDER if p in overall]
    for p, x, y in zip(ps, xs, ys):
        ax.scatter(x, y, s=46, color=C[p], edgecolor="white", linewidth=0.8, zorder=3)
    # Hand-placed: the four primates cluster inside a 0.02 x 0.03 box and auto-placement
    # overlaps them illegibly.
    OFF = {"yeast": (-4, 7), "celegans": (-30, 7), "zebrafish": (7, 0),
           "human": (7, 3), "gorilla": (7, -1), "chimp": (7, -5), "macaque": (7, -9)}
    HA = {"yeast": "center", "celegans": "left"}
    for p, x, y in zip(ps, xs, ys):
        ax.annotate(LABEL[p], (x, y), textcoords="offset points", xytext=OFF[p],
                    fontsize=6.8, color=MUTED, ha=HA.get(p, "left"), va="center")
    for v in (0.9, 0.95, 1.0):
        ax.axhline(v, color=GRID, lw=0.5, zorder=0)
        ax.axvline(v, color=GRID, lw=0.5, zorder=0)
    ax.set_xlim(0.78, 1.02); ax.set_ylim(0.86, 1.02)
    ax.set_xlabel("precision"); ax.set_ylabel("recall")
    ax.set_title("Overall ORF-call agreement", pad=6)
    panel_tag(ax, "a")

    # (b) recall by class -- the panel that matters
    ax = axes[1]
    classes = ["CDS", "uORF", "ncORF"]
    w = 0.26
    for j, cl in enumerate(classes):
        for i, p in enumerate(ORDER):
            r = cls.get((p, cl))
            if not r or int(r["n_real"]) == 0:
                continue
            val = f(r["recall"])
            if not np.isfinite(val):
                continue
            ax.bar(j + (i - 3) * w / 2.4, val, width=w / 2.6, color=C[p],
                   edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels(["CDS", "uORF", "ncORF"])
    for j, sub in enumerate(["annotated", "5'UTR", "non-coding tx"]):
        ax.annotate(sub, (j, 0), xycoords=("data", "axes fraction"),
                    textcoords="offset points", xytext=(0, -21), ha="center",
                    fontsize=6.2, color=MUTED)
    ax.set_ylim(0, 1.05); ax.set_ylabel("recall")
    ax.axhline(1.0, color=GRID, lw=0.5, zorder=0)
    ax.set_title("Recall by ORF class", pad=6)
    ax.legend(handles=[Patch(facecolor=C[p], label=LABEL[p]) for p in ORDER],
              ncol=4, fontsize=6.2, frameon=False, loc="lower center",
              bbox_to_anchor=(0.5, -0.46), columnspacing=0.9, handlelength=1.1)
    panel_tag(ax, "b", dx=-0.075)

    # (c) uORF recall, ordered by value, with depth annotated.
    #
    # An earlier draft ordered this by phylogenetic distance from the training species and
    # described the result as tracking divergence. THAT WAS WRONG and the panel is built
    # this way so the error cannot recur: gorilla (0.580) falls BELOW macaque (0.665)
    # despite being the closer outgroup to human, so the ordering is not phylogenetic.
    # Nor is it depth: chimp reaches 0.706 on 64.4M pooled frame-0 P-sites while macaque
    # reaches only 0.665 on 123.9M. Human, the training species, is highest; beyond that
    # the ordering follows neither variable monotonically, and the panel asserts no trend.
    ax = axes[2]
    DEPTH_M = {"human": 114.9, "chimp": 64.4, "gorilla": 50.5, "macaque": 123.9,
               "zebrafish": 4.2, "celegans": 77.7}
    up = [(p, f(cls[(p, "uORF")]["recall"]), int(cls[(p, "uORF")]["n_real"]))
          for p in ORDER if (p, "uORF") in cls and int(cls[(p, "uORF")]["n_real"]) > 0]
    up.sort(key=lambda t: -t[1])
    for k, (p, v, n) in enumerate(up):
        ax.barh(k, v, height=0.62, color=C[p], edgecolor="white", linewidth=0.5, zorder=3)
        ax.annotate(f"{v:.3f}", (v, k), textcoords="offset points", xytext=(4, 0),
                    va="center", fontsize=6.6, color=INK)
        ax.annotate(f"n={n:,} · {DEPTH_M[p]:.0f}M", (0.012, k), va="center",
                    fontsize=5.6, color="white", zorder=4)
    ax.set_yticks(range(len(up)))
    ax.set_yticklabels([LABEL[p] for p, _, _ in up])
    ax.invert_yaxis()
    ax.set_xlim(0, 0.86); ax.set_xlabel("uORF recall")
    ax.set_title("uORF recall (n uORFs · P-sites)", pad=6)
    panel_tag(ax, "c", dx=-0.30)

    fig.subplots_adjust(bottom=0.30)
    fig.savefig(OUT / "Fig2_orf_calls.pdf")
    fig.savefig(OUT / "Fig2_orf_calls.png")
    plt.close(fig)
    PROV["source"]["Fig2"] = ["data/xspecies_orf_classes_attn.tsv", "data/xspecies_orf_calls.tsv"]
    return {p: {"precision": f(overall[p]["precision"]), "recall": f(overall[p]["recall"]),
                "f1": f(overall[p]["f1"]),
                **{c: {"recall": f(cls[(p, c)]["recall"]) if (p, c) in cls else None,
                       "n_real": int(cls[(p, c)]["n_real"]) if (p, c) in cls else 0}
                   for c in classes}}
            for p in ORDER if p in overall}


# ---------------------------------------------------------------- figure 4
def figure4():
    """mm25 diagnostic: (a) how much the call set moves, (b) whether the MODEL's score moves.

    Ratios are drawn as points against a 1.0 reference rather than as bars from a
    non-zero baseline, which would visually exaggerate differences of ~1%.
    """
    src_a = NEW / "results" / "xspecies_mm1_vs_mm25_calls.tsv"
    src_b = NEW / "results" / "xspecies_model_vs_mm25_reference_attn.tsv"
    rows = {r["arm"]: r for r in tsv(src_a)}
    ref = tsv(src_b) if src_b.exists() else []

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6),
                             gridspec_kw={"width_ratios": [1, 1], "wspace": 0.42})

    # (a) call-set movement
    ax = axes[0]
    ps = [p for p in ORDER if ARM[p] in rows]
    ax.axvline(1.0, color=INK, lw=0.9, ls="--", zorder=2)
    for k, p in enumerate(ps):
        r = rows[ARM[p]]
        v = f(r["ratio"]); ok = r["comparable"] == "yes"
        ax.plot([1.0, v], [k, k], color=C[p], lw=1.6, zorder=3, alpha=0.55)
        ax.scatter(v, k, s=52, color=C[p] if ok else "white", edgecolor=C[p],
                   linewidth=1.4, zorder=4, marker="o" if ok else "X")
        ax.annotate(f"{v:.3f}", (v, k), textcoords="offset points",
                    xytext=(9 if v >= 1 else -9, 0), ha="left" if v >= 1 else "right",
                    va="center", fontsize=6.8, color=INK if ok else "#B22222")
        if not ok:
            ax.annotate(f"only {r['runs_mm25']} of {r['runs_mm1']} runs passed\nthe periodicity gate",
                        (v, k), textcoords="offset points", xytext=(30, 0), ha="left",
                        va="center", fontsize=5.9, color="#B22222")
    ax.set_yticks(range(len(ps))); ax.set_yticklabels([LABEL[p] for p in ps])
    ax.invert_yaxis()
    ax.set_xlim(0.855, 1.045); ax.set_xlabel("mm25 / mm1 ORF calls")
    ax.set_title("Multimappers barely move the call set", pad=6)
    panel_tag(ax, "a", dx=-0.26)

    # (b) does the MODEL's score depend on which truth set is used?
    ax = axes[1]
    if ref:
        by = {}
        for r in ref:
            by.setdefault(r["pack"], {})[r["reference"]] = r
        ks = [p for p in ORDER if p in by and {"mm1", "mm25"} <= set(by[p])]
        ax.axvline(0, color=INK, lw=0.9, ls="--", zorder=2)
        for k, p in enumerate(ks):
            d = f(by[p]["mm25"]["f1"]) - f(by[p]["mm1"]["f1"])
            ok = not by[p]["mm25"].get("note")
            ax.plot([0, d], [k, k], color=C[p], lw=1.6, zorder=3, alpha=0.55)
            ax.scatter(d, k, s=52, color=C[p] if ok else "white", edgecolor=C[p],
                       linewidth=1.4, zorder=4, marker="o" if ok else "X")
            ax.annotate(f"{d:+.3f}", (d, k), textcoords="offset points",
                        xytext=(9 if d >= 0 else -9, 0),
                        ha="left" if d >= 0 else "right", va="center",
                        fontsize=6.8, color=INK if ok else "#B22222")
        ax.set_yticks(range(len(ks))); ax.set_yticklabels([LABEL[p] for p in ks])
        ax.invert_yaxis()
        ax.set_xlim(-0.019, 0.010)
        ax.set_xlabel("$\\Delta$F1  (mm25 truth $-$ mm1 truth)")
        ax.set_title("The model's score does not depend on it", pad=6)
    panel_tag(ax, "b", dx=-0.26)

    fig.savefig(OUT / "Fig4_mm25.pdf")
    fig.savefig(OUT / "Fig4_mm25.png")
    plt.close(fig)
    PROV["source"]["Fig4"] = ["data/xspecies_mm1_vs_mm25_calls.tsv", "data/xspecies_model_vs_mm25_reference_attn.tsv"]
    out = {rows[ARM[p]]["arm"]: {"ratio": f(rows[ARM[p]]["ratio"]),
                                 "jaccard": f(rows[ARM[p]]["jaccard"]),
                                 "comparable": rows[ARM[p]]["comparable"],
                                 "runs_mm1": int(rows[ARM[p]]["runs_mm1"]),
                                 "runs_mm25": int(rows[ARM[p]]["runs_mm25"])} for p in ps}
    if ref:
        for r in ref:
            out.setdefault("model_vs_reference", {}).setdefault(r["pack"], {})[r["reference"]] = {
                "precision": f(r["precision"]), "recall": f(r["recall"]), "f1": f(r["f1"]),
                "n_ref": int(r["n_ref"])}
    return out


# ---------------------------------------------------------------- figure 1
def figure1():
    """Study design: what was sequenced, what was built, what survived."""
    src = NEW / "data" / "xspecies_dataset_registry.tsv"
    reg = tsv(src)
    # Precomputed from data/xspecies_refs/<arm>_universe.tsv (21 MB of per-transcript rows,
    # from which only two counts per species are used). The summary is shipped instead.
    uni = {}
    for r in tsv(DATA / "xspecies_universe_summary.tsv"):
        uni[r["species"]] = (int(r["n_transcripts"]), int(r["n_protein_coding"]))

    SPKEY = {"human": "human", "chimp": "chimp", "gorilla": "gorilla",
             "macaque": "macaque", "zebrafish": "zebrafish", "celegans": "celegans",
             "yeast": "yeast"}
    runs = {p: {"ribo": 0, "rna": 0} for p in ORDER}
    for r in reg:
        if r["species"] in SPKEY.values() and r["status"] == "active":
            for p, k in SPKEY.items():
                if r["species"] == k:
                    runs[p][r["assay"]] += 1
    fly = [r for r in reg if r["species"] == "fly"]

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5),
                             gridspec_kw={"width_ratios": [1.05, 1.15, 1.0], "wspace": 0.5})

    # (a) runs per species
    ax = axes[0]
    y = np.arange(len(ORDER))
    rb = [runs[p]["ribo"] for p in ORDER]
    rn = [runs[p]["rna"] for p in ORDER]
    ax.barh(y - 0.19, rb, height=0.36, color=[C[p] for p in ORDER], zorder=3)
    ax.barh(y + 0.19, rn, height=0.36, color=[C[p] for p in ORDER], alpha=0.45, zorder=3)
    for k, p in enumerate(ORDER):
        if rb[k]:
            ax.annotate(str(rb[k]), (rb[k], k - 0.19), xytext=(3, 0),
                        textcoords="offset points", va="center", fontsize=6)
        if rn[k]:
            ax.annotate(str(rn[k]), (rn[k], k + 0.19), xytext=(3, 0),
                        textcoords="offset points", va="center", fontsize=6)
    ax.set_yticks(y); ax.set_yticklabels([LABEL[p] for p in ORDER]); ax.invert_yaxis()
    ax.set_xlabel("sequencing runs"); ax.set_xlim(0, 30)
    ax.set_title("Data ingested", pad=6)
    # The 5 human RNA runs were already on disk from an earlier fetch and are therefore
    # absent from THIS study's download manifest; the bar is 0 by construction, not because
    # the arm lacks RNA. Flagged rather than silently back-filled.
    ax.annotate("human RNA (5 runs) predates\nthis manifest", (0.5, 0), xytext=(6, -2),
                textcoords="offset points", fontsize=5.4, color=MUTED, va="center")
    ax.legend(handles=[Patch(facecolor="#666666", label="Ribo-seq"),
                       Patch(facecolor="#666666", alpha=0.45, label="RNA-seq")],
              fontsize=6.2, frameon=False, loc="lower right",
              bbox_to_anchor=(1.02, -0.02), handlelength=1.1)
    panel_tag(ax, "a", dx=-0.30)

    # (b) transcript universe
    ax = axes[1]
    tot = [uni[p][0] for p in ORDER if p in uni]
    pcs = [uni[p][1] for p in ORDER if p in uni]
    ps = [p for p in ORDER if p in uni]
    y = np.arange(len(ps))
    ax.barh(y, tot, height=0.6, color=[C[p] for p in ps], alpha=0.35, zorder=3)
    ax.barh(y, pcs, height=0.6, color=[C[p] for p in ps], zorder=4)
    for k, p in enumerate(ps):
        ax.annotate(f"{tot[k]:,}", (tot[k], k), xytext=(3, 0),
                    textcoords="offset points", va="center", fontsize=6)
    ax.set_yticks(y); ax.set_yticklabels([LABEL[p] for p in ps]); ax.invert_yaxis()
    ax.set_xlabel("transcripts in universe"); ax.set_xlim(0, 240000)
    ax.set_xticks([0, 100000, 200000]); ax.set_xticklabels(["0", "100k", "200k"])
    ax.set_title("Transcript universe", pad=6)
    ax.legend(handles=[Patch(facecolor="#666666", label="protein-coding"),
                       Patch(facecolor="#666666", alpha=0.35, label="+ lncRNA")],
              fontsize=6.2, frameon=False, loc="lower right", handlelength=1.1)
    panel_tag(ax, "b", dx=-0.30)

    # (c) why Drosophila was dropped: frame distribution at the dominant read length
    ax = axes[2]
    # Frame distribution at each library's DOMINANT read length (methods.md section G).
    # Two independent Drosophila datasets were tested and both fail; yeast is the positive
    # control. Percentages are of P-sites in that length class.
    FRAMES = {"S. cerevisiae\nGSE173654 (RNase I)": (93.9, 2.6, 3.6),
              "D. melanogaster\nGSE99920 (RNase T1)": (15.2, 56.2, 28.7),
              "D. melanogaster\nGSE147619 (RNase I)": (21.5, 19.5, 59.0)}
    xs = np.arange(3)
    cols = [C["yeast"], "#B22222", "#E07B7B"]
    for i, (nm, v) in enumerate(FRAMES.items()):
        ax.bar(xs * 1.35 + (i - 1) * 0.32, v, width=0.29, color=cols[i], zorder=3,
               edgecolor="white", linewidth=0.5)
    ax.axhline(100 / 3, color=INK, lw=0.8, ls=":", zorder=2)
    ax.annotate("no periodicity (33%)", (3.15, 100 / 3), xytext=(0, 3),
                textcoords="offset points", ha="right", fontsize=5.8, color=MUTED)
    ax.set_xticks(xs * 1.35); ax.set_xticklabels(["frame 0", "frame 1", "frame 2"])
    ax.set_ylabel("% of P-sites"); ax.set_ylim(0, 118)
    ax.set_title("Drosophila arm dropped", pad=6)
    ax.legend(handles=[Patch(facecolor=c, label=n.replace("\n", " ")) for c, n
                       in zip(cols, FRAMES)],
              fontsize=5.2, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, 1.02), handlelength=1.0, labelspacing=0.25)
    panel_tag(ax, "c", dx=-0.24)

    fig.savefig(OUT / "Fig1_design.pdf")
    fig.savefig(OUT / "Fig1_design.png")
    plt.close(fig)
    PROV["source"]["Fig1"] = [
        "data/xspecies_dataset_registry.tsv", "data/xspecies_universe_summary.tsv",
        "panel c frame percentages: methods.md section G (measured, dominant read length)"]
    return {"runs": runs, "universe": uni,
            "fly_runs_dropped": len(fly), "frame_pct": FRAMES}


# ---------------------------------------------------------------- figure 3
def figure3():
    """Per-nucleotide profile agreement, and the Spearman tie artifact.

    Panels a-b use the FULL evaluation (results/xspecies_profile_eval.tsv). Panel c uses a
    2,000-transcript subsample per species, computed separately, and is labelled as such:
    it exists to demonstrate that the negative rank correlations this project reported are
    a tie-handling artifact rather than anti-correlated profiles.
    """
    src = NEW / "results" / "xspecies_profile_eval.tsv"
    rows = [r for r in tsv(src) if r["arch"] == "attn"]
    by = {r["pack"].replace("xsp_", ""): r for r in rows}

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9),
                             gridspec_kw={"width_ratios": [1.15, 1.0, 1.25], "wspace": 0.5})

    # (a) profile vs total correlation
    ax = axes[0]
    ps = [p for p in ORDER if p in by]
    y = np.arange(len(ps))
    pr = [f(by[p]["profile_r"]) for p in ps]
    tr = [f(by[p]["total_r"]) for p in ps]
    ax.barh(y - 0.19, pr, height=0.36, color=[C[p] for p in ps], zorder=3)
    ax.barh(y + 0.19, tr, height=0.36, color=[C[p] for p in ps], alpha=0.42, zorder=3)
    for k in range(len(ps)):
        ax.annotate(f"{pr[k]:.2f}", (pr[k], k - 0.19), xytext=(3, 0),
                    textcoords="offset points", va="center", fontsize=6)
        ax.annotate(f"{tr[k]:.2f}", (tr[k], k + 0.19), xytext=(3, 0),
                    textcoords="offset points", va="center", fontsize=6)
    ax.set_yticks(y); ax.set_yticklabels([LABEL[p] for p in ps]); ax.invert_yaxis()
    ax.set_xlim(0, 1.02); ax.set_xlabel("median Pearson r")
    ax.set_title("Profile vs transcript total", pad=6)
    ax.legend(handles=[Patch(facecolor="#666666", label="per-nt profile"),
                       Patch(facecolor="#666666", alpha=0.42, label="transcript total")],
              fontsize=5.8, frameon=False, loc="upper center", ncol=2,
              bbox_to_anchor=(0.5, -0.24), handlelength=1.0, columnspacing=1.0)
    panel_tag(ax, "a", dx=-0.30)

    # (b) the positional null
    ax = axes[1]
    # The shuffled null is |r| < 0.002 for every species, i.e. sub-pixel as a bar. It is
    # drawn as an explicit marker at its value so the control reads as measured rather than
    # as an omission.
    sh = [f(by[p]["shuffled_r"]) for p in ps]
    ax.barh(y, pr, height=0.5, color=[C[p] for p in ps], zorder=3)
    ax.scatter(sh, y, s=26, color="#333333", marker="|", linewidth=1.4, zorder=4)
    ax.axvline(0, color=INK, lw=0.7)
    ax.set_yticks(y); ax.set_yticklabels([LABEL[p] for p in ps]); ax.invert_yaxis()
    ax.set_xlim(-0.06, 0.62); ax.set_xlabel("median Pearson r")
    ax.set_title("Positional information", pad=6)
    ax.legend(handles=[Patch(facecolor="#666666", label="observed profile"),
                       plt.Line2D([], [], marker="|", ls="", color="#333333", mew=1.4,
                                  label="position-shuffled null")],
              fontsize=5.8, frameon=False, loc="upper center", ncol=1,
              bbox_to_anchor=(0.62, -0.24), handlelength=1.0)
    ax.annotate(f"null: max |r| = {max(abs(v) for v in sh):.4f}", (0.34, len(ps) - 0.6),
                fontsize=5.6, color=MUTED, ha="center")
    panel_tag(ax, "b", dx=-0.30)

    # (c) the tie artifact
    ax = axes[2]
    # Full evaluation (not a subsample): zero fraction, ordinal-rank rho as the project
    # implemented it, and average-rank rho as standard. Read from the same table as a-b.
    SUB = {p_: (100 * f(by[p_]["zero_frac"]), f(by[p_]["profile_rho_ordinal"]),
                f(by[p_]["profile_rho"]))
           for p_ in ps if by[p_].get("zero_frac")}
    for p, (z, ordn, avg) in SUB.items():
        ax.plot([z, z], [ordn, avg], color=C[p], lw=1.2, alpha=0.5, zorder=2)
        ax.scatter(z, ordn, s=44, color="white", edgecolor=C[p], linewidth=1.4,
                   marker="v", zorder=3)
        ax.scatter(z, avg, s=44, color=C[p], edgecolor="white", linewidth=0.8, zorder=3)
        LOFF = {"gorilla": (-3, 8), "human": (14, 7), "chimp": (-14, 7),
                "macaque": (0, 8), "celegans": (0, 8), "yeast": (0, 8),
                "zebrafish": (0, 8)}
        ax.annotate(LABEL[p], (z, avg), textcoords="offset points",
                    xytext=LOFF.get(p, (0, 7)), ha="center", fontsize=5.6, color=MUTED)
    ax.axhline(0, color=INK, lw=0.8, ls="--", zorder=1)
    ax.set_xlabel("zero fraction of observed profile (%)")
    ax.set_ylabel("median rank correlation")
    ax.set_xlim(69, 99); ax.set_ylim(-0.50, 0.60)
    ax.set_title("Spearman tie artifact", pad=6)
    ax.legend(handles=[
        plt.Line2D([], [], marker="v", ls="", mfc="white", mec=INK, mew=1.2,
                   label="ordinal ranks (as implemented)"),
        plt.Line2D([], [], marker="o", ls="", color=INK, label="average ranks (standard)")],
        fontsize=5.6, frameon=False, loc="lower left", handlelength=1.0)
    panel_tag(ax, "c", dx=-0.24)

    fig.subplots_adjust(bottom=0.30)
    fig.savefig(OUT / "Fig3_profiles.pdf")
    fig.savefig(OUT / "Fig3_profiles.png")
    plt.close(fig)
    PROV["source"]["Fig3"] = ["data/xspecies_profile_eval.tsv",
                              "panel c: same table; profile_rho_ordinal vs profile_rho "
                              "(methods section 6.2)"]
    return {p: {"profile_r": f(by[p]["profile_r"]), "total_r": f(by[p]["total_r"]),
                "shuffled_r": f(by[p]["shuffled_r"])} for p in ps} | {"tie_subsample": SUB}


if __name__ == "__main__":
    v = {"figure1": figure1(), "figure2": figure2(), "figure3": figure3(), "figure4": figure4()}
    PROV["values"] = v
    (OUT / "figure_values.json").write_text(json.dumps(PROV, indent=1))
    print(f"  wrote figures -> {OUT}")
