#!/usr/bin/env python3
"""Figure B6: what each input modality actually contributes, on the DEPLOYED recipe.

The panel exists to answer one question the proteogenomics argument depends on: is the model
cell-type-specific, or is it a sequence-only ORF prior that would emit the same thing for every cell
type? Zeroing the RNA-seq coverage channel and retraining answers it directly.

The measured answer is a double dissociation, and it is NOT the answer the claim assumed:

  * SHAPE is sequence-borne. Sequence-only recovers 98.5% of the deployed pc profile Pearson, and the
    3-nt periodicity actually SHARPENS without RNA-seq (the coverage channel adds a smooth envelope
    that dilutes fine periodicity).
  * MAGNITUDE is RNA-seq-borne. The count head loses 0.193 Pearson without it, and RNA-seq alone beats
    sequence alone at per-transcript depth.

So the cell-type specificity runs through the COUNT HEAD -- which ORFs clear the caller's depth
threshold -- not through the profile shape, which is nearly cell-type-invariant. The figure is drawn
to make that split legible: two panels, shape on the left, magnitude on the right, with the
sequence-only bar annotated as a percentage of the deployed model on the shape side only.

Every value is read from the runs' extra_metrics.json / test_metrics.json at build time, and the run
directory names are written into B6_values.json. Nothing is hardcoded -- five Fig 1 panels were found
on 2026-08-08 to have been silently built from a superseded checkpoint, so figures in this project
either read their provenance or assert it.

cas12a env. Regenerate: `python3 make_input_ablation.py`.
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
HERE = NEW / "figures/B6_input_ablation"

# (label, run dir, bar colour). `both` is the deployed run itself, not a re-train, so the contrast
# carries no extra seed noise.
ARMS = [
    ("both\n(deployed)", "results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes",
     "#2C6FBB"),
    ("sequence only\n(RNA-seq zeroed)",
     "results/ablation/orf_v2_attn_onehot_union_nokozak_input_emb_holdout_Hepatocytes", "#6FA8DC"),
    ("RNA-seq only\n(sequence zeroed)",
     "results/ablation/orf_v2_attn_onehot_union_nokozak_input_cov_holdout_Hepatocytes", "#C8A45B"),
]


def load(run):
    d = NEW / run
    e = json.loads((d / "extra_metrics.json").read_text())
    t = json.loads((d / "test_metrics.json").read_text())
    return {
        "pc_profile": e["protein_coding"]["profile_pearson_median"],
        "lnc_profile": e["lncRNA"]["profile_pearson_median"],
        "pc_count": e["protein_coding"]["count_pearson"],
        "pc_period": t["protein_coding"]["period_pred_median"],
        "n": t["n"],
    }


def main():
    vals = [(lab, load(run), c) for lab, run, c in ARMS]
    ns = {v["n"] for _, v, _ in vals}
    if len(ns) != 1:
        raise SystemExit(f"arms scored different transcript counts {ns} -- not comparable")
    n = ns.pop()
    base = vals[0][1]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.6, 3.7))
    x = np.arange(len(vals)); w = 0.36

    # -- left: SHAPE (profile Pearson pc + lncRNA, and periodicity as a line) --
    axL.bar(x - w / 2, [v["pc_profile"] for _, v, _ in vals], w, label="protein-coding",
            color=[c for _, _, c in vals], edgecolor="white", lw=0.5)
    axL.bar(x + w / 2, [v["lnc_profile"] for _, v, _ in vals], w, label="lncRNA",
            color=[c for _, _, c in vals], edgecolor="white", lw=0.5, alpha=0.55)
    for i, (_, v, _) in enumerate(vals):
        axL.text(i - w / 2, v["pc_profile"] + 0.012, f"{v['pc_profile']:.3f}", ha="center",
                 va="bottom", fontsize=7.5)
        axL.text(i + w / 2, v["lnc_profile"] + 0.012, f"{v['lnc_profile']:.3f}", ha="center",
                 va="bottom", fontsize=7.5, color="#666")
    # The headline number: sequence-only as a fraction of the deployed model, on shape. Two lines and
    # inside the bar, so it cannot run past the bar edge into the lncRNA bar beside it.
    # Percentage only. "of deployed" was tried underneath at fontsize 6.5 and still overflowed the
    # 0.36-wide bar into the lncRNA bar beside it; the subtitle already establishes the comparison.
    frac = 100 * vals[1][1]["pc_profile"] / base["pc_profile"]
    axL.text(1 - w / 2, vals[1][1]["pc_profile"] / 2, f"{frac:.1f}%", ha="center", va="center",
             fontsize=10, color="white", fontweight="bold")
    axL.set_xlabel("bar label = % of the deployed model's protein-coding profile Pearson",
                   fontsize=6.8, color="#666", labelpad=2)
    axL.set_ylabel("Per-transcript profile Pearson (median)", fontsize=9)
    axL.set_ylim(0, 0.80)
    axL.set_title("SHAPE is sequence-borne", fontsize=10)
    axL.legend(frameon=False, fontsize=7.5, loc="upper right")

    # -- right: MAGNITUDE (count Pearson) --
    axR.bar(x, [v["pc_count"] for _, v, _ in vals], w * 1.6,
            color=[c for _, _, c in vals], edgecolor="white", lw=0.5)
    for i, (_, v, _) in enumerate(vals):
        # The middle bar's value label sits BELOW its top: the drop arrow occupies the space above it,
        # and an above-bar label collided with the arrow's magnitude text.
        below = (i == 1)
        axR.text(i, v["pc_count"] + (-0.030 if below else 0.012), f"{v['pc_count']:.3f}",
                 ha="center", va="top" if below else "bottom", fontsize=8,
                 color="white" if below else "black")
    drop = vals[1][1]["pc_count"] - base["pc_count"]
    axR.annotate("", xy=(1, vals[1][1]["pc_count"]), xytext=(1, base["pc_count"]),
                 arrowprops=dict(arrowstyle="<->", color="#B03A2E", lw=1.2))
    axR.text(1.30, (base["pc_count"] + vals[1][1]["pc_count"]) / 2, f"{drop:+.3f}",
             fontsize=9, color="#B03A2E", fontweight="bold", va="center", ha="left")
    axR.set_ylabel("Protein-coding count Pearson", fontsize=9)
    axR.set_ylim(0, 1.0)
    axR.set_title("MAGNITUDE is RNA-seq-borne", fontsize=10)

    for ax in (axL, axR):
        ax.set_xticks(x)
        ax.set_xticklabels([lab for lab, _, _ in vals], fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("Cell-type specificity enters through the count head, not the profile shape\n"
                 f"(held-out Hepatocytes, n={n:,} transcripts, deployed union recipe)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"B6_input_ablation.{ext}", dpi=300, bbox_inches="tight")

    out = {"n_scored": n,
           "arms": [{"label": lab.replace("\n", " "), "run": run, **v}
                    for (lab, run, _), (_, v, _) in zip(ARMS, vals)],
           "seq_only_pct_of_deployed_shape": frac,
           "count_pearson_loss_without_rnaseq": drop}
    (HERE / "B6_values.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"n={n:,}")
    for lab, v, _ in vals:
        print(f"  {lab.replace(chr(10), ' '):<32} pc_profile={v['pc_profile']:.4f} "
              f"lnc={v['lnc_profile']:.4f} count={v['pc_count']:.4f} period={v['pc_period']:.4f}")
    print(f"  sequence-only keeps {frac:.1f}% of shape; count Pearson {drop:+.3f} without RNA-seq")


if __name__ == "__main__":
    main()
