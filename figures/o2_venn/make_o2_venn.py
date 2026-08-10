#!/usr/bin/env python3
"""Three-way Venn: Wang observed / Janich observed / model predicted ORF calls.

What it shows that an F1 table cannot: the Wang-Janich overlap IS the between-experiment reproducibility
ceiling, made visual, and it exposes the region the pairwise F1s hide -- ORFs the model recovers that only
ONE real experiment saw. Those are not false positives; they are the model siding with one experiment
where the other lacked depth, which is the "as informative as a second experiment" claim made concrete.

Layout: rows = transcript space (INTERSECTION of TPM>=1, then UNION), cols = ORF class (CDS, uORF, novel).
One figure per (model arm, operating point), per the standing two-arm rule -- an uncalibrated panel
inflates the model-only lobe with calibration slack rather than biology.

Transcript space (standing rule, feedback_orf_call_transcript_space): every panel states its own n, and
all three call sets are restricted to the same space. The INTERSECTION is the primary comparison -- all
three parties had a fair shot at every transcript. The UNION is shown to expose the coverage asymmetry:
a transcript expressed only in Wang cannot be called by Janich, so its absence there is not evidence.

Circles are AREA-PROPORTIONAL (matplotlib_venn 1.1.2 solves radii and offsets from the set sizes),
so region area tracks the count rather than being schematic.
cas12a env. Run: python make_o2_venn.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn3 as mpl_venn3

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
ECH = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/biotype_probe/expression_context_human")
HERE = NEW / "figures/o2_venn"

WANG_OBS = NEW / "data/heldout_psites/mouse_wang_liver/mouse_wang_liver_collapsed.txt"
JANICH_OBS = ECH / "data/ribocode_mouse_liver/Liver_5samp/Liver_5samp_collapsed.txt"
QW = NEW / "data/salmon_quant_wang_liver/wang_merged_quant.sf"
QJ = NEW / "data/salmon_quant_janich_liver/janich_merged_quant.sf"

ARMS = {
    "wang":   {"theta1": NEW / "results/o2_liver/wang_nokozak/dropin/pred_preddepth_collapsed.txt",
               "calib":  NEW / "results/o2_liver/wang_nokozak/dropin_sweep_poisson/theta_0.2/"
                               "pred_preddepth_collapsed.txt",
               "label": "model on WANG RNA"},
    "janich": {"theta1": NEW / "results/o2_liver/janich_wanguni/dropin/pred_preddepth_collapsed.txt",
               "calib":  NEW / "results/o2_liver/janich_wanguni/dropin_sweep_poisson/theta_0.2/"
                               "pred_preddepth_collapsed.txt",
               "label": "model on JANICH RNA"},
}
CLASSES = [("annotated", "CDS"), ("uORF", "uORF"), ("novel", "novel")]
C_W, C_J, C_M = "#2C6FBB", "#B5714F", "#6E8B74"       # Wang, Janich, model


def coord_key(o):
    p = o.rsplit("_", 3)
    return "_".join(p[-3:]) if len(p) >= 4 else o


def load(path):
    rows = []
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        oi, ti, ty = h.index("ORF_ID"), h.index("transcript_id"), h.index("ORF_type")
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            rows.append((coord_key(p[oi]), p[ty], p[ti]))
    return rows


def tpm(path):
    d = {}
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.split("\t")
            d[p[0].split("|")[0]] = float(p[3])
    return d


def sel(rows, cls, keep):
    return {k for k, t, tx in rows if t == cls and (keep is None or tx in keep)}


def venn3(ax, A, B, C, labels, colors, title, sub):
    """AREA-PROPORTIONAL 3-set Venn (matplotlib_venn solves circle radii + offsets from the set
    sizes, so region area tracks count). Returns the 7 region counts."""
    v = {
        "A": len(A - B - C), "B": len(B - A - C), "C": len(C - A - B),
        "AB": len((A & B) - C), "AC": len((A & C) - B), "BC": len((B & C) - A),
        "ABC": len(A & B & C),
    }
    if not (A or B or C):
        ax.text(0, 0, "no calls", ha="center", va="center", fontsize=7, color="#888888")
        ax.set_title(title, fontsize=8.6, fontweight="bold", pad=6)
        ax.axis("off")
        return v
    d = mpl_venn3([set(A), set(B), set(C)], set_labels=labels, ax=ax,
                  set_colors=colors, alpha=0.42)
    for t in (d.set_labels or []):
        if t is not None:
            t.set_fontsize(7.2)
            t.set_fontweight("bold")
    for t in (d.subset_labels or []):
        if t is not None:
            t.set_fontsize(6.9)
    # colour the set labels to match their circles
    for t, c in zip(d.set_labels or [], colors, strict=False):
        if t is not None:
            t.set_color(c)
    ax.set_title(title, fontsize=8.6, fontweight="bold", pad=6)
    ax.text(0.5, -0.13, sub, ha="center", va="top", fontsize=6.2, color="#5b6270",
            transform=ax.transAxes)
    return v


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    wo, jo = load(WANG_OBS), load(JANICH_OBS)
    tw, tj = tpm(QW), tpm(QJ)
    ew = {t for t, v in tw.items() if v >= 1}
    ej = {t for t, v in tj.items() if v >= 1}
    spaces = [("INTERSECTION of TPM>=1", ew & ej), ("UNION of TPM>=1", ew | ej)]
    print(f"expressed: Wang {len(ew):,}  Janich {len(ej):,}  "
          f"intersection {len(ew & ej):,}  union {len(ew | ej):,}", file=sys.stderr)

    for arm, cfg in ARMS.items():
        for op, key in (("theta1", "theta = 1 (uncalibrated)"),
                        ("calib", "CDS-anchored Poisson (theta = 0.2)")):
            path = cfg[op]
            if not path.exists():
                print(f"  skip {arm}/{op}: {path.name} not present yet", file=sys.stderr)
                continue
            mo = load(path)
            # Per-region set sizes for the whole 2x3 grid, written beside the figure. A Venn is the
            # least readable place to recover a number from, so the seven regions are tabulated.
            venn_rows = []
            fig, axes = plt.subplots(2, 3, figsize=(11.4, 8.0), facecolor="white")
            fig.suptitle(f"Wang vs Janich vs {cfg['label']}   |   {key}",
                         fontsize=11.5, fontweight="bold", y=0.985)
            for r, (sp_lab, keep) in enumerate(spaces):
                for c, (cls, cname) in enumerate(CLASSES):
                    A, B, C = sel(wo, cls, keep), sel(jo, cls, keep), sel(mo, cls, keep)
                    only1 = len((A ^ B) & C)            # model agrees with exactly one experiment
                    tot1 = len(A ^ B)
                    sub = (f"model recovers {only1:,}/{tot1:,} "
                           f"({only1/tot1*100:.0f}%) of single-experiment ORFs" if tot1 else "")
                    venn3(axes[r][c], A, B, C,
                          ["Wang obs", "Janich obs", "model"], [C_W, C_J, C_M],
                          f"{cname}  --  {sp_lab}", sub)
                    venn_rows.append({
                        "arm": arm, "threshold": op, "orf_class": cls, "tx_space": sp_lab,
                        "n_tx_in_space": len(keep),
                        "wang_only": len(A - B - C), "janich_only": len(B - A - C),
                        "model_only": len(C - A - B), "wang_janich": len((A & B) - C),
                        "wang_model": len((A & C) - B), "janich_model": len((B & C) - A),
                        "all_three": len(A & B & C),
                        "wang_total": len(A), "janich_total": len(B), "model_total": len(C),
                        "single_experiment_orfs": tot1, "model_recovers_of_those": only1})
                axes[r][0].text(-0.13, 0.5, f"{sp_lab}\n({len(keep):,} tx)", rotation=90,
                                ha="center", va="center", fontsize=8, fontweight="bold",
                                color="#33415C", transform=axes[r][0].transAxes)
            fig.text(0.5, 0.012,
                     "Wang n Janich = the between-experiment reproducibility ceiling. "
                     "Regions where the model overlaps exactly one experiment are ORFs it recovers that "
                     "only that experiment detected -- not false positives.",
                     ha="center", fontsize=7.2, color="#5b6270")
            fig.tight_layout(rect=[0.03, 0.03, 1, 0.96])
            name = f"o2_venn_{arm}_{op}"
            for ext in ("pdf", "png"):
                fig.savefig(HERE / f"{name}.{ext}", dpi=300, bbox_inches="tight")
            plt.close(fig)
            cols = list(venn_rows[0]) if venn_rows else []
            with open(HERE / f"{name}.tsv", "w") as fh:
                fh.write("\t".join(cols) + "\n")
                for row in venn_rows:
                    fh.write("\t".join(str(row[c]) for c in cols) + "\n")
            print(f"  wrote {name}.pdf / .png / .tsv ({len(venn_rows)} regions)", file=sys.stderr)


if __name__ == "__main__":
    main()
