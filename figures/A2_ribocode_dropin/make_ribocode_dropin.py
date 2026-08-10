#!/usr/bin/env python3
"""Figure A2: the model's PREDICTED per-nt profile, dropped into RiboCode in place of real Ribo-seq,
reproduces RiboCode's OWN ORF calls at F1 ~0.91 -- i.e. the prediction is good enough to call ORFs with no
Ribo-seq for the query. Panel a: precision / recall / F1 for the fully-predicted drop-in (predicted profile
AND predicted depth = 'pred_preddepth_vs_real'), with the observed-depth variant for reference. Panel b:
recall by ORF type. Reusable: reads dropin_metrics.json, regenerates if the model is retrained. cas12a env.

CHECKPOINT (repointed 2026-08-08). This read `orf_v2_attn_onehot_holdout_Hepatocytes` -- the
pre-nokozak / pre-mm1 / pre-union checkpoint -- while the project shipped the union models. Run
directories are named for the DATASET, not the CHECKPOINT, so the staleness was invisible from the
path. It surfaced only in a scan of every figure generator, which found five of the six Fig 1 panels
in the same state. The panel now reads the released run, and asserts the checkpoint below rather than
trusting the path: on release the number moved 0.923 -> 0.909.

MODEL (2026-08-08). Defaults to **mamba4**, per locked decision D1b (2026-07-31): main Figure 1 is
`orf_v2_mamba4`, with `orf_v2_attn` moving to supplemental but staying shipped as the CPU inference
path. This panel had been on attn -- D1b was applied to the plan but never to the generators. Switch
with `FIG_MODEL=attn` for the supplemental variant; the chosen model is written into the values JSON.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/A2_ribocode_dropin"
MODEL = os.environ.get("FIG_MODEL", "mamba4")
if MODEL not in ("mamba4", "attn"):
    raise SystemExit(f"FIG_MODEL={MODEL!r}; want 'mamba4' (main, per D1b) or 'attn' (supplemental)")
RUN = NEW / f"results/loto/orf_v2_{MODEL}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes"
SRC = RUN / "dropin/dropin_metrics.json"
# The npz `meta` field is the only place the checkpoint identity actually lives. Assert it, so a
# repointed-but-stale source fails loudly instead of producing a plausible number.
RELEASED_TAG = "union_noBrain_nokozak_mm1"
# The MAIN figure (mamba4, per D1b) keeps the plain filename; the supplemental attn variant gets
# an "_attn" suffix. Without this the two share one path and whichever ran last silently wins --
# which happened once already, with the attn build overwriting the mamba4 PDF.
SUF = "" if MODEL == "mamba4" else f"_{MODEL}"

MODEL_C, OBS_C = "#2C6FBB", "#88AED0"
TYPE_ORDER = ["annotated", "uORF", "Overlap_uORF", "dORF", "Overlap_dORF", "internal", "novel"]


def check_checkpoint():
    npz = RUN / "dropin/pred_profiles.npz"
    if not npz.exists():
        raise SystemExit(f"no {npz} -- cannot verify the checkpoint this figure is built from")
    meta = str(np.atleast_1d(np.load(npz, allow_pickle=True)["meta"])[0])
    if RELEASED_TAG not in meta:
        raise SystemExit(f"ABORT: {npz} came from '{meta}', which is not a released "
                         f"({RELEASED_TAG}) checkpoint. Re-dump before rebuilding this figure.")
    return meta


def main():
    print(f"checkpoint: {check_checkpoint()}")
    d = json.load(open(SRC))
    full = d["pred_preddepth_vs_real"]      # predicted profile + predicted depth (no Ribo-seq at all)
    obsd = d["pred_obsdepth_vs_real"]       # predicted profile, observed depth (reference)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.2, 3.5), gridspec_kw=dict(width_ratios=[1, 1.5]))

    # panel a: precision / recall / F1
    metrics = ["precision", "recall", "f1"]
    x = np.arange(len(metrics)); w = 0.38
    ax1.bar(x - w/2, [full[m] for m in metrics], w, label="fully predicted", color=MODEL_C,
            edgecolor="white", lw=0.5)
    ax1.bar(x + w/2, [obsd[m] for m in metrics], w, label="pred. profile,\nobs. depth", color=OBS_C,
            edgecolor="white", lw=0.5)
    for xi, m in enumerate(metrics):
        ax1.text(xi - w/2, full[m] + 0.005, f"{full[m]:.3f}", ha="center", va="bottom", fontsize=7, color=MODEL_C)
    ax1.set_xticks(x); ax1.set_xticklabels(["Precision", "Recall", "F1"], fontsize=8.5)
    ax1.set_ylim(0.5, 1.0); ax1.set_ylabel("vs RiboCode's own calls", fontsize=9)
    ax1.legend(frameon=False, fontsize=7, loc="lower right")
    ax1.spines[["top", "right"]].set_visible(False)
    ax1.set_title(f"Drop-in ORF calling\nF1 = {full['f1']:.3f}", fontsize=9.5)

    # panel b: recall by ORF type (fully predicted)
    rpt = full["recall_per_type"]
    types = [t for t in TYPE_ORDER if t in rpt]
    rec = [rpt[t]["recall"] for t in types]
    nrl = [rpt[t]["n_real"] for t in types]
    cols = ["#444" if t == "annotated" else "#2C6FBB" for t in types]
    y = np.arange(len(types))[::-1]
    ax2.barh(y, rec, color=cols, edgecolor="white", lw=0.5, height=0.66)
    for yi, (r, n) in zip(y, zip(rec, nrl)):
        ax2.text(r + 0.01, yi, f"{r:.2f}  (n={n})", va="center", fontsize=7, color="#333")
    ax2.set_yticks(y); ax2.set_yticklabels(types, fontsize=8)
    ax2.set_xlim(0, 1.12); ax2.set_xlabel("Recall of RiboCode's calls", fontsize=9)
    ax2.spines[["top", "right"]].set_visible(False)
    ax2.set_title("Recall by ORF type (fully predicted)", fontsize=9.5)

    fig.suptitle("Predicted profile reproduces the ORF caller's own calls (held-out Hepatocytes)",
                 fontsize=10, y=1.02)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"A2_ribocode_dropin{SUF}.{ext}", dpi=300, bbox_inches="tight")
    (HERE / f"A2_values{SUF}.json").write_text(json.dumps(
        {"fully_predicted": {k: full[k] for k in ("precision", "recall", "f1")},
         "recall_per_type": {t: rpt[t] for t in types}}, indent=2) + "\n")
    print(f"wrote A2_ribocode_dropin.pdf/.png; fully-predicted F1={full['f1']:.3f} "
          f"P={full['precision']:.3f} R={full['recall']:.3f}")
    for t in types:
        print(f"  {t}: recall {rpt[t]['recall']:.3f} (n={rpt[t]['n_real']})")


if __name__ == "__main__":
    main()
