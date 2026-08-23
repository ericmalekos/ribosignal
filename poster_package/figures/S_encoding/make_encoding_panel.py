#!/usr/bin/env python3
"""Supplemental figure: one-hot sequence encoding vs RNA foundation-model embeddings.

THE POINT, which is accessibility and not a horse race: the simplest, dependency-free encoding is as
good as a 650M-parameter RNA language model on this task. That is why the released models carry no
foundation-model dependency and run on CPU.

All four arms are the SAME recipe, the SAME holdout and the SAME 33,918 scored transcripts, differing
only in `--emb_backend`, so the comparison is internally clean.

RECIPE CAVEAT, stated on the figure itself. These are `noBrain_nokozak_mm1` runs on the FIBROBLAST
universe, not the union universe the released models use. Union FM embeddings were never extracted, so
a union version of this panel would mean three more 12-17 h trainings. The ENCODING conclusion does not
depend on the universe -- the four arms move together -- but the absolute Pearson values are
recipe-specific and must not be quoted beside released-model numbers.

Reads `test_metrics.json` from each run; writes the plotted values to `S_encoding.tsv` alongside the
figure. cas12a env.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/S_encoding"
RECIPE = "noBrain_nokozak_mm1"

# (display label, run dir, colour). one-hot first: it is the reference, not a competitor.
ARMS = [
    ("one-hot\n(shipped)", f"orf_v2_attn_onehot_{RECIPE}_holdout_Hepatocytes", "#2C6FBB"),
    ("RiNALMo\n(650M)", f"orf_v2_attn_rinalmo_{RECIPE}_holdout_Hepatocytes", "#7D8CA3"),
    ("Orthrus\n(Mamba)", f"orf_v2_attn_orthrus_{RECIPE}_holdout_Hepatocytes", "#9AA7B4"),
    ("HydraRNA\n(full-length)", f"orf_v2_attn_hydrarna_{RECIPE}_holdout_Hepatocytes", "#B6C0CA"),
]
# (test_metrics path, display label, axis label)
METRICS = [
    (("protein_coding", "pearson_median"), "protein-coding\nprofile Pearson", "per-transcript Pearson"),
    (("lncRNA", "pearson_median"), "lncRNA\nprofile Pearson", "per-transcript Pearson"),
    (("protein_coding", "period_pred_median"), "predicted 3-nt\nperiodicity (pc)", "periodicity"),
]


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    data, ns = {}, set()
    for lab, run, _ in ARMS:
        f = NEW / "results/loto" / run / "test_metrics.json"
        if not f.exists():
            raise SystemExit(f"missing {f}")
        t = json.loads(f.read_text())
        data[lab] = t
        ns.add(t["n"])
    if len(ns) != 1:
        raise SystemExit(f"arms scored different transcript counts {ns} -- not comparable")
    n = ns.pop()

    fig, axes = plt.subplots(1, len(METRICS), figsize=(9.4, 3.3))
    for ax, (path, title, ylab) in zip(axes, METRICS):
        vals = [data[lab][path[0]][path[1]] for lab, _, _ in ARMS]
        x = np.arange(len(ARMS))
        ax.bar(x, vals, 0.62, color=[c for _, _, c in ARMS], edgecolor="white", lw=0.5)
        base = vals[0]
        ax.axhline(base, color="#2C6FBB", lw=0.8, ls=":")
        for xi, v in enumerate(vals):
            ax.text(xi, v + max(vals) * 0.012, f"{v:.4f}", ha="center", va="bottom", fontsize=7)
            if xi:
                # 3 dp, not 4: at 4 dp the delta label is wider than the 0.62-wide bar.
                ax.text(xi, v / 2, f"{v - base:+.3f}", ha="center", va="center", fontsize=6.4,
                        color="white", fontweight="bold")
        ax.set_xticks(x); ax.set_xticklabels([lab for lab, _, _ in ARMS], fontsize=7.2)
        ax.set_ylim(0, max(vals) * 1.16)
        ax.set_ylabel(ylab, fontsize=8)
        ax.set_title(title, fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("The dependency-free encoding matches RNA foundation-model embeddings\n"
                 f"held-out Hepatocytes, n={n:,} transcripts; identical recipe, only --emb_backend "
                 f"differs (dotted line = one-hot)", fontsize=9.5)
    fig.text(0.5, 0.005,
             f"Recipe: {RECIPE} on the Fibroblast universe, NOT the union universe of the released "
             "models -- union FM embeddings were never extracted. Encoding ordering is the finding; "
             "absolute values are recipe-specific.",
             ha="center", fontsize=6.6, color="#666")
    fig.tight_layout(rect=(0, 0.035, 1, 0.90))
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"S_encoding.{ext}", dpi=300, bbox_inches="tight")

    cols = ["encoding", "run", "n_scored", "pc_profile_pearson_median", "lncRNA_profile_pearson_median",
            "all_profile_pearson_median", "pc_period_pred_median", "pc_period_obs_median",
            "delta_pc_vs_onehot"]
    base_pc = data[ARMS[0][0]]["protein_coding"]["pearson_median"]
    with open(HERE / "S_encoding.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for lab, run, _ in ARMS:
            t = data[lab]
            fh.write("\t".join([
                lab.replace("\n", " "), run, str(t["n"]),
                f"{t['protein_coding']['pearson_median']:.6f}",
                f"{t['lncRNA']['pearson_median']:.6f}",
                f"{t['all']['pearson_median']:.6f}",
                f"{t['protein_coding']['period_pred_median']:.6f}",
                f"{t['protein_coding']['period_obs_median']:.6f}",
                f"{t['protein_coding']['pearson_median'] - base_pc:+.6f}"]) + "\n")

    print(f"n={n:,} scored transcripts, recipe {RECIPE}")
    for lab, _, _ in ARMS:
        t = data[lab]
        print(f"  {lab.replace(chr(10), ' '):<22}pc {t['protein_coding']['pearson_median']:.4f}   "
              f"lnc {t['lncRNA']['pearson_median']:.4f}   "
              f"periodPRED {t['protein_coding']['period_pred_median']:.4f}   "
              f"(pc {t['protein_coding']['pearson_median'] - base_pc:+.4f} vs one-hot)")
    print(f"wrote {HERE}/S_encoding.pdf/.png/.tsv")


if __name__ == "__main__":
    main()
