#!/usr/bin/env python3
"""Poster panel P9: UpSet of ORF calls -- three mouse-liver Ribo-seq experiments and the model.

The F1 tables say how well the model reproduces each experiment. They cannot say WHERE the
disagreement lives, and in particular they cannot separate two very different things that both cost
F1: an ORF the model missed that all three experiments agree on, and an ORF the model called that no
experiment saw. This panel separates them.

Four sets, all on one shared universe (22,974 transcripts) and one filtering path:

  Janich / GSE243134 / Wang   each experiment's OWN observed RiboCode calls (`real`)
  model                       the STANDALONE arm (`pred_preddepth`) -- no Ribo-seq at inference

The bars to read first are the extremes: all four (model agrees with the consensus of three
independent experiments), model only (called by the model and by NO experiment -- the candidate false
positives, and the honest upper bound on over-calling), and all-three-but-not-the-model (a systematic
blind spot, which is a different failure needing a different fix).

FOUR VIEWS ARE EMITTED, because the classes behave nothing alike and pooling them hides that:

  (no suffix)  all ORF classes
  _cds         annotated CDS -- the large, easy, high-consensus class
  _lncrna      lncRNA ORFs (RiboCode's `novel`) -- small, and where the disagreement lives
  _uorf        uORFs in the 5'UTR, non-overlapping -- a separate class, never inside `novel`

In the pooled view the CDS bars set the y-scale and everything else is invisible, which is the whole
reason for the split. All four are filtered from ONE load, so they share a transcript space exactly
rather than being separate re-loads that could drift.

`novel` AND `lncRNA` ARE THE SAME SET HERE, not nested: 100% of novel calls sit on lncRNA
transcripts in every set checked (954/954 observed GSE243134, 2,774/2,774 model, 649/649 Wang). That
follows from the definitions -- RiboCode calls an ORF `novel` when its transcript has no annotated
CDS, and in a protein-coding + lncRNA universe those transcripts are lncRNAs. The view is labelled
lncRNA because that is the informative name; there is no further split to make.

WHICH MODEL SET. `pred_preddepth` depends on (model, RNA input) and is identical down each RNA
column, so there is no single "the model" set -- there are three, one per RNA. These panels use the
GSE243134 RNA (the deepest, 19 libraries). The alternatives are written to the values JSON so the
choice is checkable rather than implicit.

Filtering is `compare_dropin_calls.build_loader`, the same path behind every other drop-in number
here, and keying is genomic (gene_id, ORF_gstop) so a different representative isoform of one ORF is
not counted as a disagreement (feedback_orf_call_transcript_space).

cas12a env. Regenerate: `python3 make_upset_mouse3x3.py`.
"""
from __future__ import annotations

import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NEW = pathlib.Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
                   "RNAZoo/experiments/riboseq_signal_model")
sys.path.insert(0, str(NEW / "scripts"))
from compare_dropin_calls import build_loader  # noqa: E402

HERE = NEW / "figures/P9_upset_mouse3x3"
MODEL = "mamba4"                 # the single-model family, and the loader/universe reference
MODEL_RNA = "gse243134"
DATASETS = ["janich", "gse243134", "wang"]
TX2GENE = NEW / "data/tx2biotype_mouse.tsv"
OBSNAMES = ["Janich", "GSE243134", "Wang"]
# Two families. The 5-set one asks a question the 4-set one cannot: when the model calls an ORF no
# experiment saw, does the OTHER architecture call it too? attn is a transformer, mamba4 an SSM.
#
# THE ORIGINAL READING OF THAT BAR WAS WRONG AND IS CORRECTED HERE. It was built expecting that two
# independent architectures agreeing on the same "false positive" would be evidence of a real ORF the
# experiments were too shallow to see. The 28-library merged reference falsifies that: calls made by
# BOTH models and no experiment are corroborated at 8.7%, indistinguishable from either model alone
# (7.7% / 8.6%) and far below the 74.5% experiment-unique control (results.md 2026-08-14). Two
# architectures agreeing is a SHARED SYSTEMATIC ERROR, not independent evidence, so the shared bar
# must not be read as a confidence filter.
FAMILIES = [("", ["mamba4"]), ("_2models", ["attn", "mamba4"])]
# Measured false-positive rates for the "model only" bar, from the merged 28-library reference.
# Read at draw time so the annotation cannot drift from the scorer's output.
OVERCALL = NEW / "results/merged_liver_ribocode/overcall/overcall_meta.json"
OVERCALL_CLASS = NEW / "results/merged_liver_ribocode/overcall/overcall_recovery_by_class.tsv"
MAX_BARS = 20                    # 5 sets -> up to 31 intersections; the tail is reported, not hidden
MODELCOL = "#c05621"
OBSCOL = "#2b6cb0"
# One view per ORF class that behaves differently. RiboCode's `novel` is VERIFIED here to be 100%
# lncRNA-hosted in every set (954/954 observed GSE243134, 2,774/2,774 model, 649/649 Wang), so
# "novel" and "lncRNA ORF" are the same set in this universe, not nested -- the view is labelled
# lncRNA because that is the informative name. uORF is a SEPARATE RiboCode class, never inside
# `novel`, and gets its own view. Overlap_uORF is deliberately excluded from the uORF view: it
# overlaps the annotated CDS and behaves differently; its counts are printed in the summary.
VARIANTS = [("", None, "all ORF classes"),
            ("_cds", {"annotated"}, "annotated CDS"),
            ("_lncrna", {"novel"}, "lncRNA ORFs (RiboCode 'novel')"),
            ("_uorf", {"uORF"}, "uORFs (5'UTR, non-overlapping)")]


def cell(model, ribo, rna):
    return NEW / f"results/mouse_liver_3x3_canon/{model}_ribo-{ribo}_rna-{rna}"


# View -> the ORF class its rate is looked up under in the by-class table. The pooled view has no
# single class and uses the overall row instead.
VIEW_CLASS = {"": None, "_cds": "annotated", "_lncrna": "novel", "_uorf": "uORF"}


def overcall_rate(suffix, n_models):
    """Measured corroboration rate for the 'model only' bar, or None if the scorer has not run.

    Returned as (fp_percent, recovered_percent, control_percent) so the annotation states the
    measurement AND the control it is judged against -- a bare 'X% false positives' would be
    unreadable without knowing that a genuinely real, singly-observed ORF only corroborates 74.5%
    of the time either.
    """
    if not OVERCALL.exists() or not OVERCALL_CLASS.exists():
        return None
    meta = json.loads(OVERCALL.read_text())
    ctrl = [r for r in meta["rows"] if r["kind"].startswith("experiment-unique")]
    control = 100 * sum(r["recovered"] for r in ctrl) / max(1, sum(r["n"] for r in ctrl))
    grp = "model_unique_both" if n_models > 1 else f"model_unique_{MODEL}"
    cls = VIEW_CLASS.get(suffix)
    if cls is None:
        d = meta["detail"].get(grp)
        rec = 100 * d["rate"] if d else None
    else:
        rec = None
        for line in OVERCALL_CLASS.read_text().splitlines()[1:]:
            f = line.split("\t")
            if f[0] == grp and f[1] == cls:
                rec = 100 * float(f[4])
    return None if rec is None else (100 - rec, rec, control)


def draw(sets, keys, n_obs_sets, setnames, suffix, vlabel, alt, n_test_tx,
         var_sfx=""):
    # var_sfx is the VARIANT suffix alone ("", "_cds", ...). `suffix` is family+variant and
    # names the output files; keying the class lookup on it silently fell back to the pooled
    # rate for every 2-model view, printing one number for all four classes.
    universe = set().union(*sets.values())
    if not universe:
        print(f"  SKIP {vlabel}: no calls in any set")
        return None
    combos = {}
    for k in universe:
        pat = tuple(k in sets[s] for s in keys)
        combos[pat] = combos.get(pat, 0) + 1
    order_all = sorted([p for p in combos if any(p)], key=lambda p: -combos[p])
    order = order_all[:MAX_BARS]
    dropped = order_all[MAX_BARS:]
    n_drop = sum(combos[p] for p in dropped)

    fig = plt.figure(figsize=(9.8, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.5, 1.25], width_ratios=[1, 3.6],
                          hspace=0.06, wspace=0.06)
    axBar = fig.add_subplot(gs[0, 1])
    axDot = fig.add_subplot(gs[1, 1], sharex=axBar)
    axSet = fig.add_subplot(gs[1, 0], sharey=axDot)

    xs = np.arange(len(order))
    vals = [combos[p] for p in order]
    n_m = len(keys) - n_obs_sets
    bcol = []
    for p in order:
        n_obs = sum(p[:n_obs_sets])
        n_mod = sum(p[n_obs_sets:])
        if n_mod == n_m and n_obs == 0:
            bcol.append(MODELCOL)                      # every model, no experiment
        elif n_mod == 0 and n_obs == n_obs_sets:
            bcol.append("#718096")                     # every experiment, no model
        elif n_mod == n_m and n_obs == n_obs_sets:
            bcol.append("#276749")                     # full consensus
        else:
            bcol.append("#a0aec0")
    axBar.bar(xs, vals, color=bcol, zorder=3)
    for x, v in zip(xs, vals):
        axBar.text(x, v + max(vals) * 0.015, f"{v:,}", ha="center", fontsize=7)
    axBar.set_ylabel("ORFs in this exact intersection", fontsize=8.5)
    axBar.grid(axis="y", alpha=0.25, zorder=0)
    axBar.tick_params(labelbottom=False, labelsize=8)
    axBar.set_ylim(0, max(vals) * 1.10)
    axBar.set_xlim(-0.7, len(order) + 2.6)
    axBar.spines[["top", "right"]].set_visible(False)

    for i, p in enumerate(order):
        for j, on in enumerate(p):
            axDot.plot(i, j, "o", ms=8,
                       color=(MODELCOL if j >= n_obs_sets else OBSCOL) if on else "#e2e8f0",
                       zorder=3)
        on_idx = [j for j, o in enumerate(p) if o]
        if len(on_idx) > 1:
            axDot.plot([i, i], [min(on_idx), max(on_idx)], "-", color="0.35", lw=1.4, zorder=2)
    # Set names as TEXT to the RIGHT. y-ticks do not survive here: axSet shares this axis, so
    # clearing its ticks also cleared these labels, and wspace leaves no room on the left anyway.
    axDot.set_yticks([])
    for j, nm in enumerate(setnames):
        is_model = j >= n_obs_sets
        axDot.text(len(order) - 0.05, j, nm, ha="left", va="center", fontsize=8,
                   color=(MODELCOL if is_model else "0.2"),
                   fontweight="bold" if is_model else "normal")
    axDot.set_ylim(-0.6, len(keys) - 0.4)
    axDot.set_xlim(-0.7, len(order) + 2.6)
    axDot.invert_yaxis()
    axDot.tick_params(labelbottom=False, length=0)
    for sp in axDot.spines.values():
        sp.set_visible(False)
    axDot.grid(False)

    tot = [len(sets[k]) for k in keys]
    axSet.barh(range(len(keys)), tot,
               color=[OBSCOL] * n_obs_sets + [MODELCOL] * (len(keys) - n_obs_sets),
               height=0.5, zorder=3)
    for j, t in enumerate(tot):
        axSet.text(t * 1.04, j, f"{t:,}", va="center", ha="right", fontsize=7, color="0.2")
    axSet.set_xlim(0, max(tot) * 1.38)
    axSet.invert_xaxis()
    axSet.invert_yaxis()
    axSet.tick_params(left=False, labelleft=False, labelsize=7)
    axSet.set_xlabel("set size", fontsize=8)
    axSet.spines[["top", "right", "left"]].set_visible(False)

    mlab = "model" if n_m == 1 else "both models"
    # The "model only" bar is no longer labelled "candidate FPs": a 28-library merged reference
    # measured how many are corroborated, so the panel states the measurement rather than a guess.
    oc = overcall_rate(var_sfx, n_m)
    fplab = (f"{mlab} only\n{oc[0]:.0f}% are FALSE POSITIVES\n"
             f"({oc[1]:.0f}% corroborated by a 28-library\nmerge vs {oc[2]:.0f}% control)"
             if oc else f"{mlab} only\n(candidate FPs)")
    # The FP label grew from 2 lines to 4 when it gained the measurement, and at the old yfrac it
    # sat on top of the green annotation and the bar counts. Headroom above the tallest bar plus a
    # wider vertical spread keeps all three legible.
    axBar.set_ylim(0, max(vals) * 1.30)
    ann = [(tuple([False] * n_obs_sets + [True] * n_m), fplab, MODELCOL, 1.02, 1.3),
           (tuple([True] * len(keys)), f"all 3 experiments + {mlab}", "#276749", 0.58, 2.6),
           (tuple([True] * n_obs_sets + [False] * n_m), f"all 3 experiments,\n{mlab} missed",
            "#718096", 0.30, 4.2)]
    for pat, lab, c, yfrac, dx in ann:
        if pat not in combos:
            continue
        i = order.index(pat)
        axBar.annotate(lab, xy=(i, combos[pat]), xytext=(i + dx, max(vals) * yfrac),
                       ha="left", va="center", fontsize=7.5, color=c, fontweight="bold",
                       arrowprops=dict(arrowstyle="->", color=c, lw=1.0,
                                       connectionstyle="arc3,rad=-0.15"))

    # Name the models actually drawn. The global MODEL is the loader/universe reference, not the
    # set list, so hardcoding it labelled the 5-set panels "mamba4" while plotting attn as well.
    drawn = " + ".join(setnames[n_obs_sets:])
    fig.suptitle(f"Mouse liver, {vlabel}: model vs three independent Ribo-seq experiments"
                 f"\n({drawn}, standalone arm, {MODEL_RNA} RNA; genomic keying)",
                 fontsize=9.5, y=0.995)
    HERE.mkdir(parents=True, exist_ok=True)
    fig.savefig(HERE / f"P9_upset_mouse3x3{suffix}.pdf", bbox_inches="tight")
    fig.savefig(HERE / f"P9_upset_mouse3x3{suffix}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    out = {"model": MODEL, "model_rna": MODEL_RNA, "arm": "pred_preddepth", "view": vlabel,
           "keying": "genomic (gene_id, ORF_gstop)", "n_test_tx": n_test_tx,
           "set_sizes": {k: len(sets[k]) for k in keys},
           "alt_model_set_sizes_by_rna": alt,
           "intersections": [{"pattern": dict(zip(keys, p)), "n": combos[p]} for p in order_all],
           "bars_drawn": len(order), "bars_omitted": len(dropped), "orfs_in_omitted_bars": n_drop}
    (HERE / f"P9_values{suffix}.json").write_text(json.dumps(out, indent=2))
    if dropped:
        axBar.text(0.995, 0.97, f"+{len(dropped)} smaller intersections ({n_drop:,} ORFs) not drawn",
                   transform=axBar.transAxes, ha="right", va="top", fontsize=6.5, color="0.45")
    named = {lab.replace("\n", " "): combos.get(pat, 0) for pat, lab, *_rest in ann}
    return {"sizes": {k: len(sets[k]) for k in keys}, "named": named,
            "omitted": (len(dropped), n_drop)}


def main():
    base = cell(MODEL, MODEL_RNA, MODEL_RNA)
    lc, keep_tx, _ = build_loader(str(base / "pred_profiles.npz"), key="genomic",
                                  tx2gene=str(TX2GENE))
    # Load ONCE, keeping the full dicts: each call carries its ORF `type`, and every view below is a
    # filter on that. Re-loading per view would risk the three panels drifting apart.
    raw = {}
    for ds in DATASETS:
        raw[ds] = lc(str(cell(MODEL, ds, ds) / "real_collapsed.txt"), is_pred=False)
    for m in ("attn", "mamba4"):
        f = cell(m, MODEL_RNA, MODEL_RNA) / "pred_preddepth_collapsed.txt"
        if f.exists():
            raw[m] = lc(str(f), is_pred=True)
        else:
            print(f"  WARNING: no {m} dump at {f}")

    alt = {}
    for rna in DATASETS:
        if rna == MODEL_RNA:
            continue
        d = cell(MODEL, rna, rna)
        if (d / "pred_preddepth_collapsed.txt").exists():
            alt[rna] = len(lc(str(d / "pred_preddepth_collapsed.txt"), is_pred=True))

    for fam_sfx, models in FAMILIES:
        models = [m for m in models if m in raw]
        if not models:
            continue
        keys = DATASETS + models
        setnames = OBSNAMES + [m for m in models]
        print(f"  === family: {' + '.join(models)} ===")
        for suffix, classes, vlabel in VARIANTS:
            sets = {k: {kk for kk, vv in raw[k].items()
                        if classes is None or vv["type"] in classes} for k in keys}
            res = draw(sets, keys, len(DATASETS), setnames, fam_sfx + suffix, vlabel,
                       alt, len(keep_tx), var_sfx=suffix)
            if res is None:
                continue
            print(f"    {vlabel}:")
            print("        sizes  " + ", ".join(f"{k}={n:,}" for k, n in res["sizes"].items()))
            for lab, n in res["named"].items():
                print(f"        {lab:<36} {n:,}")
            nd, no = res["omitted"]
            if nd:
                print(f"        (+{nd} smaller intersections omitted from the plot, {no:,} ORFs)")
    # Overlap_uORF is not plotted; report it so the omission is visible rather than silent.
    ov = {k: sum(1 for v in raw[k].values() if v["type"] == "Overlap_uORF") for k in raw}
    print("  Overlap_uORF (NOT plotted, overlaps the CDS): "
          + ", ".join(f"{k}={n:,}" for k, n in ov.items()))


if __name__ == "__main__":
    main()
