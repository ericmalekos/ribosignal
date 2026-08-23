#!/usr/bin/env python3
"""Poster panel P5: the model reproduces ORF calls on four independent human held-out datasets.

This is the poster's validation claim, and nothing in the project showed all four arms together
before. Each arm is a different KIND of held-out, which is the point -- they fail differently:

  hepatocyte   LOTO holdout. Same study, same protocol, tissue withheld from training. The easiest
               arm, and the one where ~99.7% of transcripts were still seen in another tissue
               (quantified separately: memorisation is worth ~0.10 profile r at matched depth).
  iPSC-CM      Ruiz-Orera 2024 (ENA PRJEB65856). Different lab, different protocol, and
               cardiomyocyte is absent from the 9-tissue training panel -- a genuinely
               out-of-distribution CELL TYPE, not just a withheld one.
  THP-1        GSE208041. Different lab, monocytic leukaemia line.
  CAR-T        GSE304796. Different lab, primary engineered T cells.

Two arms are always reported (feedback_two_arm_orf_calling):
  pred_obsdepth   predicted SHAPE scaled by the dataset's observed depth
  pred_preddepth  fully standalone -- no Ribo-seq at inference at all
The standalone arm is the one that matters for the application, since it is what you get from
sequence + RNA-seq alone.

Panel B is the honest half. Overall F1 is dominated by annotated CDS (the large majority of every
reference set), so a single number reads as ~0.9 everywhere and hides the real behaviour. Splitting
recall into annotated vs non-canonical shows the model at 0.99 on canonical CDS and 0.55-0.65 on
non-canonical -- which is where the remaining work is, and where the proteogenomics application
actually operates.

Every arm states n_ref (feedback_orf_call_transcript_space), and the panel records which pack each
number came from, because CAR-T's original pack carried ~29% PCR duplicates (dedup ran on the genome
BAM while the pack was built from the undeduplicated transcriptome BAM). The generator PREFERS a
canonical CAR-T rescore when one exists and says which it used, so a stale number cannot be plotted
silently.

cas12a env. Regenerate: `python3 make_heldout_human_panel.py`.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/P5_heldout_human"
# DEFAULT REMAINS mamba4 (locked decision D1b: mamba4 is the main-figure model). `--model attn`
# renders the transformer instead, which is what the POSTER package uses -- the poster emphasises the
# transformer because it is the architecture most readers know and the only CPU-runnable one
# (mamba_ssm dispatches to causal_conv1d_cuda). Both are available for all four arms and the two are
# statistically indistinguishable, so this is a presentation choice, not a result change.
#
# THE MODEL NAME IS WRITTEN INTO P5_values.json AND ONTO THE FIGURE. An earlier poster README printed
# this panel's mamba4 numbers under an "attn" heading; the label is now carried by the artifact.
MODEL = "mamba4"
ARMS = ["pred_obsdepth", "pred_preddepth"]
NONCANON = ["novel", "uORF", "internal", "Overlap_uORF", "dORF", "Overlap_dORF"]

# Two source shapes. `dropin` arms carry recall_grouped already; `scored` arms are per-class TSV rows
# that must be pooled into annotated / non-canonical here. Built as functions of MODEL so --model
# cannot leave a path pointing at the other architecture.
def _dropin_paths():
    """Per-arm dropin metrics, PREFERRING the final-recipe pack where one exists.

    The hepatocyte arm has two: its original `dropin/` (built on data/packed_union, the mm1 packs)
    and `dropin_canonpack/` (the same deployed checkpoint re-dumped against packed_canon_Hepatocytes,
    produced by eval_canon_orfcalls for task #86). The latter is on the final recipe and is used when
    present, exactly as the scored arms prefer human_orf_calls_canon.

    MEASURED DIFFERENCE, attn standalone: F1 0.909 -> 0.911, n_ref 16,329 -> 16,236 (-0.6%). That is
    NOT the ~2-point drop THP-1 and CAR-T showed, and the reason is that packed_union was already
    ncRNA + cross-gene filtered -- it differs from packed_canon only by the 2026-08-13 realignment,
    worth 0.8% of reference calls on this tissue. Do not assume the THP-1/CAR-T delta transfers to
    an arm whose pack was already filtered.
    """
    hep_canon = (f"results/loto/orf_v2_{MODEL}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/"
                 "dropin_canonpack/dropin_metrics.json")
    ipsc_canon = (f"results/heldout/human_ruizorera/released_{MODEL}/"
                  "dropin_canonpack/dropin_metrics.json")
    hep = hep_canon if (NEW / hep_canon).exists() else (
        f"results/loto/orf_v2_{MODEL}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/"
        "dropin/dropin_metrics.json")
    return {
        "hepatocyte\n(LOTO holdout)": hep,
        "iPSC-CM\n(Ruiz-Orera, cross-study)": (
            ipsc_canon if (NEW / ipsc_canon).exists()
            else f"results/heldout/human_ruizorera/released_{MODEL}/dropin/dropin_metrics.json"),
    }
SCORED_TSV = "results/human_orf_calls/scored/human_by_class.tsv"
SCORED_TSV_CANON = "results/human_orf_calls_canon/scored/human_by_class.tsv"
SCORED = {"THP-1\n(GSE208041)": "gse208041", "CAR-T\n(GSE304796)": "cart"}

# OVERALL F1 COMES FROM human_overall.tsv, NOT FROM SUMMING THE CLASS ROWS.
# `score_human_orf_calls.prf(pred, ref, subset=(c,))` filters BOTH sides by class, so an ORF that
# both callers found at the same genomic key but LABELLED DIFFERENTLY (uORF vs Overlap_uORF, internal
# vs uORF -- adjacent-boundary confusions) lands in neither class's intersection and disappears from
# the per-class table, while `human_overall.tsv` still counts it. The gap is small but present in
# every row measured (1-8 calls; THP-1/attn is 7, which is exactly F1 0.8501 vs 0.8506).
# Summing class rows therefore UNDERSTATES F1 slightly. results.md quotes the overall file, so
# reading it here also makes the panel and the text agree.
OVERALL_TSV = "results/human_orf_calls/scored/human_overall.tsv"
OVERALL_TSV_CANON = "results/human_orf_calls_canon/scored/human_overall.tsv"

# THIS PANEL MIXED SUBSTRATES, AND USED TO REPORT IT FOR ONLY HALF ITS ARMS.
# `prov` records a source PATH for the two SCORED arms (THP-1, CAR-T) but the two DROPIN arms
# carried no pack information at all -- and those were exactly the ones not on the final recipe.
# A reader saw provenance printed for the good arms and nothing for the questionable ones.
# Every arm now reports its substrate. State as of 2026-08-15, after task #86 and #92:
#
#   arm          pack                                     on final recipe?
#   hepatocyte   data/packed_canon_Hepatocytes            YES (re-dumped, task #86)
#   THP-1        data/packed_heldout_..._gse208041_canon  YES
#   CAR-T        data/packed_heldout_human_cart_canon     YES
#   iPSC-CM      data/packed_heldout_human_ruizorera      NO -- rebuild in flight
#
# THE SIZE OF THE CORRECTION IS NOT THE SAME FOR EVERY ARM, so do not apply one blanket adjustment:
#   THP-1       0.851 -> 0.831, n_ref -7.5%   (pack was NOT ncRNA/cross-gene filtered)
#   CAR-T       0.904 -> 0.886, n_ref -6.8%   (same, plus ~29% PCR duplicates)
#   hepatocyte  0.909 -> 0.911, n_ref -0.6%   (packed_union WAS already filtered; only the
#                                              2026-08-13 realignment differs)
# The THP-1/CAR-T delta transfers only to packs that skipped the filter entirely.
#
# iPSC-CM's status is UNKNOWN rather than confirmed bad: its BAMs, source hd5 and align script are
# all gone, so nothing survives to inspect. Its FASTQs were re-downloaded from ENA (PRJEB65856,
# 15 files, md5-verified, 2026-08-15) and are being realigned on the final recipe, which settles it
# by construction rather than by inference.
SUBSTRATE = {
    "hepatocyte\n(LOTO holdout)": ("data/packed_canon_Hepatocytes", True,
                                   "re-dumped on the final-recipe pack (dropin_canonpack, task #86). "
                                   "Its previous packed_union/mm1 dropin gave F1 0.909 / n_ref 16,329; "
                                   "the change is +0.002 / -0.6%, far smaller than THP-1's -7.5% "
                                   "because packed_union was ALREADY ncRNA + cross-gene filtered"),
    "iPSC-CM\n(Ruiz-Orera, cross-study)": ("data/packed_heldout_human_ruizorera_canon", True,
                                           "REBUILT on the final recipe 2026-08-16 (task #92): "
                                           "FASTQs re-fetched from ENA PRJEB65856 (md5-verified), "
                                           "5 Ribo runs realigned EndToEnd+mm1+filter, universe and "
                                           "RNA coverage reused verbatim. P-sites -2.5%, n_ref -0.5% "
                                           "-- so the old pack HAD been filtered after all"),
    "THP-1\n(GSE208041)": ("data/packed_heldout_human_gse208041_canon", True, ""),
    "CAR-T\n(GSE304796)": ("data/packed_heldout_human_cart_canon", True, ""),
}


def from_dropin(path):
    j = json.loads((NEW / path).read_text())
    out = {}
    for arm in ARMS:
        d = j[f"{arm}_vs_real"]
        g = d["recall_grouped"]
        out[arm] = dict(f1=d["f1"], precision=d["precision"], recall=d["recall"],
                        n_ref=d["n_real"],
                        annot_recall=g["annotated"]["recall"],
                        noncanon_recall=g["noncanonical"]["recall"],
                        n_annot=g["annotated"]["n_real"], n_noncanon=g["noncanonical"]["n_real"])
    return out, str(path)


def from_scored(tsv, dataset):
    rows = [r for r in csv.DictReader(open(NEW / tsv), delimiter="\t")
            if r["dataset"] == dataset and r["model"] == MODEL]
    if not rows:
        return None, None
    # The overall file sits beside the by-class one and is authoritative for F1 (see OVERALL_TSV).
    ov_path = NEW / str(tsv).replace("human_by_class.tsv", "human_overall.tsv")
    ov = {}
    if ov_path.exists():
        ov = {r["arm"]: r for r in csv.DictReader(open(ov_path), delimiter="\t")
              if r["dataset"] == dataset and r["model"] == MODEL}
    out = {}
    for arm in ARMS:
        a = [r for r in rows if r["arm"] == arm]
        ann = next((r for r in a if r["orf_class"] == "annotated"), None)
        nc = [r for r in a if r["orf_class"] in NONCANON]
        if not ann or not nc:
            continue
        tp_nc = sum(int(r["tp"]) for r in nc)
        nref_nc = sum(int(r["n_ref"]) for r in nc)
        tp_all = tp_nc + int(ann["tp"])
        nref_all = nref_nc + int(ann["n_ref"])
        npred_all = sum(int(r["n_pred"]) for r in a)
        prec = tp_all / npred_all if npred_all else 0.0
        rec = tp_all / nref_all if nref_all else 0.0
        d = dict(f1=2 * prec * rec / (prec + rec) if (prec + rec) else 0.0,
                 precision=prec, recall=rec, n_ref=nref_all,
                 annot_recall=float(ann["recall"]),
                 noncanon_recall=tp_nc / nref_nc if nref_nc else 0.0,
                 n_annot=int(ann["n_ref"]), n_noncanon=nref_nc,
                 f1_source="summed_by_class")
        if arm in ov:
            d.update(f1=float(ov[arm]["f1"]), precision=float(ov[arm]["precision"]),
                     recall=float(ov[arm]["recall"]), n_ref=int(ov[arm]["n_ref"]),
                     f1_source="human_overall.tsv",
                     tp_class_disagreement=int(ov[arm]["tp"]) - tp_all)
        out[arm] = d
    return out, str(tsv)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=MODEL, choices=["mamba4", "attn"],
                    help="mamba4 = locked main-figure default (D1b); attn = the poster's transformer")
    ap.add_argument("--suffix", default=None,
                    help="output filename suffix. Defaults to '' for mamba4 and '_attn' for attn, so "
                         "rendering one never overwrites the other.")
    a_ = ap.parse_args()
    globals()["MODEL"] = a_.model
    sfx = a_.suffix if a_.suffix is not None else ("" if a_.model == "mamba4" else f"_{a_.model}")
    DROPIN = _dropin_paths()

    data, prov = {}, {}
    for lab, p in DROPIN.items():
        if (NEW / p).exists():
            data[lab], prov[lab] = from_dropin(p)
        else:
            print(f"  MISSING {lab}: {p}")
    # PER-DATASET preference, not per-file. The final-recipe rescore lands one dataset at a time, so a
    # whole-file switch silently DROPS any dataset not yet rescored -- it did exactly that to THP-1,
    # turning a 4-arm panel into a 3-arm one with only a log line. Each dataset now takes the
    # canonical source if present and falls back to the original, and `prov` records which, so a
    # mixed panel is visible rather than implied.
    canon_ok = (NEW / SCORED_TSV_CANON).exists()
    for lab, ds in SCORED.items():
        d = src = None
        if canon_ok:
            d, src = from_scored(SCORED_TSV_CANON, ds)
        if not d:
            d, src = from_scored(SCORED_TSV, ds)
            if d and canon_ok:
                print(f"  NOTE {lab}: not yet in the final-recipe rescore -- using the ORIGINAL pack")
        if d:
            data[lab], prov[lab] = d, src
        else:
            print(f"  MISSING {lab}: dataset {ds} in neither scored table")
    if not data:
        raise SystemExit("no arms available")
    labs = [l for l in list(DROPIN) + list(SCORED) if l in data]
    cart_lab = next((l for l in SCORED if "CAR-T" in l), None)
    cart_canon = bool(cart_lab and "canon" in str(prov.get(cart_lab, "")))
    # Substrate for EVERY arm, not just the scored ones.
    n_canon = sum(1 for l in labs if SUBSTRATE.get(l, ("", False, ""))[1])
    print(f"  substrate: {n_canon}/{len(labs)} arms on-recipe")
    for l in labs:
        pack, canon, note = SUBSTRATE.get(l, ("(unknown)", False, "not in SUBSTRATE map"))
        print(f"    {l.replace(chr(10),' '):<34} {'ON-RECIPE    ' if canon else 'OFF-RECIPE   '} "
              f"{pack}{('  -- ' + note) if note else ''}")
    if n_canon != len(labs):
        print(f"  WARNING: this panel MIXES substrates. Moving an arm onto the final recipe lowers "
              f"its score (THP-1 0.851->0.831, CAR-T 0.904->0.886) because the REFERENCE shrinks, "
              f"not because the model changes -- 2x2 decomposition on THP-1: reference -0.0740, "
              f"prediction -0.0001. THP-1 has NO UMIs, so its loss is the ncRNA + cross-gene "
              f"filter alone (~9.7% of records); CAR-T additionally carried ~29% PCR duplicates. "
              f"The {len(labs)-n_canon} off-recipe arm(s) are therefore likely OVERSTATED.")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.0, 3.8), constrained_layout=True,
                                   gridspec_kw={"width_ratios": [1.0, 1.15]})
    x = np.arange(len(labs))
    w = 0.36
    ac = {"pred_obsdepth": "#2b6cb0", "pred_preddepth": "#c05621"}
    an = {"pred_obsdepth": "predicted shape x observed depth",
          "pred_preddepth": "standalone (no Ribo-seq at inference)"}

    for i, arm in enumerate(ARMS):
        v = [data[l][arm]["f1"] for l in labs]
        axA.bar(x + (i - 0.5) * w, v, w, label=an[arm], color=ac[arm], zorder=3)
        for xi, vi in zip(x + (i - 0.5) * w, v):
            axA.text(xi, vi + 0.012, f"{vi:.3f}", ha="center", fontsize=7)
    axA.set_xticks(x)
    axA.set_xticklabels(labs, fontsize=7.5)
    axA.set_ylabel("ORF-call F1 vs the dataset's own\nobserved RiboCode calls", fontsize=8.5)
    axA.set_ylim(0, 1.08)
    axA.set_title("A  Held-out ORF calling, 4 independent human datasets", fontsize=9,
                  loc="left", pad=26)
    # Legend ABOVE the axes: inside the axes it covered the n= annotations along the bottom.
    axA.legend(fontsize=7, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2,
               frameon=False)
    axA.grid(axis="y", alpha=0.25, zorder=0)
    for xi, l in zip(x, labs):
        axA.text(xi, 0.03, f"n={data[l][ARMS[0]]['n_ref']:,}", ha="center", fontsize=6.5, color="0.35")

    for i, arm in enumerate(ARMS):
        a = [data[l][arm]["annot_recall"] for l in labs]
        n = [data[l][arm]["noncanon_recall"] for l in labs]
        axB.bar(x + (i - 0.5) * w, a, w, color=ac[arm], alpha=0.95, zorder=3,
                label=f"annotated CDS -- {an[arm]}" if i == 0 else None)
        axB.bar(x + (i - 0.5) * w, n, w, color=ac[arm], alpha=0.45, zorder=4,
                hatch="///", edgecolor="white", linewidth=0.4)
    axB.set_xticks(x)
    axB.set_xticklabels(labs, fontsize=7.5)
    axB.set_ylabel("recall", fontsize=8.5)
    axB.set_ylim(0, 1.08)
    axB.set_title("B  Solid = annotated CDS, hatched = non-canonical", fontsize=9, loc="left")
    axB.grid(axis="y", alpha=0.25, zorder=0)
    med_a = np.median([data[l][a2]["annot_recall"] for l in labs for a2 in ARMS])
    med_n = np.median([data[l][a2]["noncanon_recall"] for l in labs for a2 in ARMS])
    axB.axhline(med_a, color="0.3", lw=0.8, ls=":", zorder=2)
    axB.axhline(med_n, color="0.3", lw=0.8, ls=":", zorder=2)
    axB.text(len(labs) - 0.42, med_a + 0.015, f"annotated median {med_a:.2f}", fontsize=7, ha="right")
    axB.text(len(labs) - 0.42, med_n + 0.015, f"non-canonical median {med_n:.2f}", fontsize=7, ha="right")

    arch = "transformer (attn)" if MODEL == "attn" else "Mamba SSM (mamba4)"
    fig.suptitle(f"Held-out human validation, {arch}: near-ceiling on annotated CDS, and the "
                 "non-canonical classes are where the work remains", fontsize=9.5)
    HERE.mkdir(parents=True, exist_ok=True)
    fig.savefig(HERE / f"P5_heldout_human{sfx}.pdf")
    fig.savefig(HERE / f"P5_heldout_human{sfx}.png", dpi=200)
    plt.close(fig)

    # `tsv` no longer exists: the source is per-dataset now, so record the mapping, not one path.
    out = {"model": MODEL, "arms": ARMS,
           "source_per_dataset": {l: str(prov[l]) for l in SCORED if l in prov},
           "cart_source_is_canonical": cart_canon,
           # Substrate for every arm, so a consumer can see the mix without reading this script.
           "substrate": {l: {"pack": SUBSTRATE.get(l, ("(unknown)", False, ""))[0],
                             "canonical": SUBSTRATE.get(l, ("", False, ""))[1],
                             "note": SUBSTRATE.get(l, ("", False, ""))[2]} for l in labs},
           "n_arms_canonical": n_canon, "n_arms": len(labs),
           "mixed_substrate": n_canon != len(labs),
           "mixed_substrate_warning": (
               "The final-recipe alignment redo LOWERS scores where applied (THP-1 0.851->0.831, "
               "CAR-T 0.904->0.886). Arms marked canonical=false predate it and are likely "
               "overstated by a similar margin. iPSC-CM cannot be rebuilt from local data: its "
               "source hd5 and its 5 ENA BAMs are both deleted."),
           "provenance": prov, "data": {l: data[l] for l in labs}}
    (HERE / f"P5_values{sfx}.json").write_text(json.dumps(out, indent=2))
    print(f"  wrote {HERE/('P5_heldout_human'+sfx+'.pdf')}   model={MODEL}")
    for l in labs:
        d = data[l][ARMS[1]]
        print(f"    {l.replace(chr(10),' '):<36} standalone F1 {d['f1']:.3f}  "
              f"annot {d['annot_recall']:.3f}  non-canon {d['noncanon_recall']:.3f}  n={d['n_ref']:,}")


if __name__ == "__main__":
    main()
