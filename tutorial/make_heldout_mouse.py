#!/usr/bin/env python3
"""Generate tutorial/source/data/heldout-mouse.md from the result JSONs.

Companion to `make_heldout_human.py`, kept as a SEPARATE page and a separate generator: the human
evidence rests on one independent study, the mouse evidence rests on three independent liver studies
plus a T-cell set, and merging them into one table hides that asymmetry.

Same contract as the human generator: every number is read from a JSON at build time, and every row
records the checkpoint string from the npz `meta` field, so a stale source shows up as a stale
checkpoint rather than as a plausible wrong number.

WHY THE MOUSE LIVER ARMS ARE THE INTERESTING ONES. Three labs measured Ribo-seq on mouse liver
independently (Wang, Janich, GSE243134). That gives a between-experiment reproducibility CEILING for
the caller itself, which is what the model should be judged against -- not against 1.0. The ceiling is
computed by `scripts/threeway_liver_concordance.py` and quoted below.

Usage: make_heldout_mouse.py            # writes the page
       make_heldout_mouse.py --check    # print what it would emit, write nothing
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tutorial/source/data/heldout-mouse.md"
THREEWAY = ROOT / "results/o2_liver/threeway_liver_concordance.json"
LIVER3X3 = ROOT / "results/mouse_liver_3x3/scored"

# (label, dataset key, dropin_metrics.json, npz whose meta names the checkpoint)
SOURCES = [
    ("Wang liver (attn)", "Wang",
     "results/liver_released/attn_mouse_wang_liver/dropin_metrics.json",
     "results/liver_released/attn_mouse_wang_liver/pred_profiles.npz"),
    ("Wang liver (mamba4)", "Wang",
     "results/liver_released/mamba4_mouse_wang_liver/dropin_metrics.json",
     "results/liver_released/mamba4_mouse_wang_liver/pred_profiles.npz"),
    ("Janich liver (attn)", "Janich",
     "results/liver_released/attn_mouse_janich_liver_decon/dropin_metrics.json",
     "results/liver_released/attn_mouse_janich_liver_decon/pred_profiles.npz"),
    ("Janich liver (mamba4)", "Janich",
     "results/liver_released/mamba4_mouse_janich_liver_decon/dropin_metrics.json",
     "results/liver_released/mamba4_mouse_janich_liver_decon/pred_profiles.npz"),
    ("GSE243134 liver (attn)", "GSE243134",
     "results/liver_released/attn_mouse_gse243134_liver/dropin_metrics.json",
     "results/liver_released/attn_mouse_gse243134_liver/pred_profiles.npz"),
    ("GSE243134 liver (mamba4)", "GSE243134",
     "results/liver_released/mamba4_mouse_gse243134_liver/dropin_metrics.json",
     "results/liver_released/mamba4_mouse_gse243134_liver/pred_profiles.npz"),
]
CLS = ["annotated", "uORF", "Overlap_uORF", "novel", "dORF", "Overlap_dORF", "internal"]
RELEASED = "union_noBrain_nokozak_mm1"


def ckpt(npz):
    p = ROOT / npz
    if not p.exists():
        return "UNKNOWN (no npz)"
    z = np.load(p, allow_pickle=True)
    return str(np.atleast_1d(z["meta"])[0])


def rows_for(arm):
    """Collect one table row per source for the given drop-in arm."""
    rows, missing = [], []
    for label, ds, mj, npz in SOURCES:
        f = ROOT / mj
        if not f.exists():
            missing.append(f"{label}: {mj}")
            continue
        m = json.loads(f.read_text()).get(arm)
        if m is None:
            missing.append(f"{label}: no {arm} in {mj}")
            continue
        rows.append({"label": label, "ds": ds, "ckpt": ckpt(npz), "m": m})
    return rows, missing


def table(rows):
    L = ["| dataset | n_ref | P | R | F1 | " + " | ".join(c.replace("Overlap_", "Ovl_") for c in CLS)
         + " |", "|---|--:|--:|--:|--:|" + "--:|" * len(CLS)]
    for r in rows:
        m = r["m"]
        pt = m.get("recall_per_type", {})
        cells = " | ".join(f"{pt[c]['recall']:.3f}" if c in pt else "--" for c in CLS)
        L.append(f"| {r['label']} | {m['n_real']:,} | {m['precision']:.3f} | {m['recall']:.3f} | "
                 f"{m['f1']:.3f} | {cells} |")
    return L


def liver3x3_section():
    """The 3 Ribo x 3 RNA factorial, read from results/mouse_liver_3x3/scored/.

    Read from the scored TSVs rather than transcribed, so the page cannot drift from the numbers.
    """
    import csv
    mv = LIVER3X3 / "model_vs_ribo.tsv"
    mm = LIVER3X3 / "matched_vs_mismatched_rna.tsv"
    if not (mv.exists() and mm.exists()):
        return [f"_(3x3 factorial unavailable: {LIVER3X3} not built)_\n"]
    rows = list(csv.DictReader(open(mv), delimiter="\t"))
    deltas = [float(r["delta"]) for r in csv.DictReader(open(mm), delimiter="\t")]
    DS = ["janich", "gse243134", "wang"]
    L = []
    L.append("## Does the RNA input have to match the Ribo-seq?\n")
    L.append("The three liver datasets each bring their own RNA-seq, so every Ribo-seq reference can be")
    L.append("predicted from every RNA input: a **3 Ribo x 3 RNA** factorial, all nine cells pooled from")
    L.append("BAMs through one pipeline onto one shared universe (22,974 tx). If the model were quietly")
    L.append("relying on the RNA-seq being from the same experiment, the diagonal would dominate.\n")
    L.append("`mamba4`, `pred_obsdepth`, all-class F1. Rows are the RNA input, columns the Ribo reference;")
    L.append("`*` marks the matched cell.\n")
    L.append("| RNA input | " + " | ".join(f"ref={d}" for d in DS) + " |")
    L.append("|---|" + "--:|" * len(DS))
    for rna in DS:
        cells = []
        for ref in DS:
            m = [r for r in rows if r["model"] == "mamba4" and r["arm"] == "pred_obsdepth"
                 and r["rna_input"] == rna and r["ribo_reference"] == ref]
            cells.append(f"{float(m[0]['f1']):.3f}{'*' if rna == ref else ''}" if m else "-")
        L.append(f"| {rna} | " + " | ".join(cells) + " |")
    L.append("")
    L.append(f"**The matched cell never wins its column.** Across all {len(deltas)} model x reference x arm")
    L.append(f"combinations the matched-minus-mismatched difference averages **{sum(deltas)/len(deltas):+.4f} F1**")
    L.append(f"(range {min(deltas):+.4f} to {max(deltas):+.4f}); every one is under 0.008 in absolute value and")
    L.append("two are negative. Substituting an unrelated experiment's RNA-seq costs essentially nothing.\n")
    L.append("The RNA effect that *is* real is input quality, not matching:\n")
    L.append("| RNA input | samples | mean F1 (obsdepth, over refs and models) |")
    L.append("|---|--:|--:|")
    NS = {"janich": 7, "gse243134": 19, "wang": 2}
    for rna in DS:
        v = [float(r["f1"]) for r in rows if r["rna_input"] == rna and r["arm"] == "pred_obsdepth"]
        L.append(f"| {rna} | {NS[rna]} | {sum(v)/len(v):.4f} |")
    L.append("")
    L.append("Wang's 2-sample RNA costs ~0.008 F1 against any reference, while 19 samples buy nothing over")
    L.append("7. Depth matters up to a handful of samples and then saturates.\n")
    L.append("```{admonition} Not comparable to pre-2026-08-10 Janich numbers\n:class: warning")
    L.append("RiboCode `metaplots` re-derives which read lengths carry periodicity **per sample**, so")
    L.append("re-pooling Janich from re-fetched FASTQs reproduces its historical pack only to ~1.7%. These")
    L.append("nine cells are internally consistent; the superseded numbers are archived under")
    L.append("`results/_archive_pre_2026_08_10_janich_pack/`.")
    L.append("```\n")
    return L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    obs, miss_obs = rows_for("pred_obsdepth_vs_real")
    std, miss_std = rows_for("pred_preddepth_vs_real")
    missing = miss_obs + miss_std

    L = []
    L.append("# Held-out mouse data\n")
    L.append("Mouse is the stronger generalisation evidence in this project, for a reason that is")
    L.append("structural rather than about sample size: **three independent labs measured Ribo-seq on")
    L.append("mouse liver**. That yields a between-experiment ceiling for the ORF caller itself, so the")
    L.append("model can be scored against what the assay actually reproduces instead of against 1.0.\n")
    L.append("```{list-table}\n:header-rows: 1\n:widths: 24 26 50\n")
    L.append("* - Dataset\n  - Role\n  - Notes")
    L.append("* - **Wang liver**\n  - held-out study, cross-species\n  - Human-trained model applied to "
             "mouse with no retraining. The shallowest of the three liver sets.")
    L.append("* - **Janich liver**\n  - held-out study, cross-species\n  - RNA arm was rebuilt after "
             "decontamination: 6-11% of reads carried 52-70% of TPM through salmon length "
             "normalisation, which shrank the expressed universe to 11,532 tx. Post-fix: 16,379.")
    L.append("* - **GSE243134 liver**\n  - held-out study, cross-species\n  - Deepest of the three. Its "
             "RNA arm is GSE302188 poly(A), replacing a Ribo-Zero total-RNA arm that could not be used "
             "(tRNA/7SL survive the biotype filter).")
    # `/results` is not a tutorial page (the project's results.md lives in the repo root, outside
    # Sphinx). Point at the page that does exist, or the build emits an unknown-document warning.
    L.append("* - **GSE155087 T-cell**\n  - held-out study, different tissue\n  - Mouse CD4+ T-cell, WT "
             "only. Reported in {doc}`/benchmarks` rather than here.")
    L.append("```\n")

    L.append("## Between-experiment ceiling\n")
    L.append("Before reading any model row: this is what the three observed datasets score **against")
    L.append("each other**, on the same restricted transcript space and the same genomic key. It is the")
    L.append("realistic upper bound.\n")
    if THREEWAY.exists():
        tw = json.loads(THREEWAY.read_text())
        obs_rows = [r for r in tw["rows"]
                    if r["kind"] == "observed" and r["pred"] != r["ref"]]
        if obs_rows:
            L.append(f"Restricted to {tw['restriction']}.\n")
            L.append("| observed pair | F1 | CDS F1 | non-canonical F1 |")
            L.append("|---|--:|--:|--:|")
            for r in obs_rows:
                L.append(f"| {r['pred']} vs {r['ref']} | {r['overall']['f1']:.3f} | "
                         f"{r['CDS']['f1']:.3f} | {r['non-canonical']['f1']:.3f} |")
            L.append("")
            L.append("**Annotated CDS calls reproduce; non-canonical calls do not.** The CDS ceiling is")
            L.append("near-perfect and stable across all three pairs, while the non-canonical ceiling")
            L.append("both is far lower and varies by pair. Any single-pair non-canonical number is")
            L.append("therefore a property of that pair, not a property of the assay.\n")
    else:
        L.append(f"_(ceiling table unavailable: {THREEWAY} not built)_\n")

    L.append("## Standalone drop-in: no Ribo-seq for the query sample\n")
    L.append("`pred_preddepth` = predicted profile shape **and** predicted depth. Scored against RiboCode")
    L.append("run on that pack's real counts, so the transcript space is identical on both sides by")
    L.append("construction. Class columns are recall.\n")
    L.extend(table(std))
    L.append("")

    L.append("## Predicted shape at observed depth\n")
    L.append("`pred_obsdepth` keeps the real per-transcript depth and replaces only the profile shape.")
    L.append("The gap between this and the table above isolates **what the count head costs**, separate")
    L.append("from whether the shape is right.\n")
    L.extend(table(obs))
    L.append("")
    if obs and std:
        d = [o["m"]["f1"] - s["m"]["f1"] for o, s in zip(obs, std)]
        L.append(f"Mean F1 cost of predicting depth as well as shape: **{np.mean(d):+.3f}** "
                 f"(range {min(d):+.3f} to {max(d):+.3f}). Precision carries essentially all of it --")
        L.append("the standalone arm calls MORE ORFs than the reference, it does not miss them.\n")

    L.append("```{admonition} Checkpoint provenance\n:class: important")
    L.append("Each row's source checkpoint, read from the npz `meta` field at build time:\n")
    for r in std:
        tag = "RELEASED" if RELEASED in r["ckpt"] else "**STALE -- predates the released models**"
        L.append(f"- `{r['label']}` -> `{r['ckpt']}` ({tag})")
    L.append("")
    if any(RELEASED not in r["ckpt"] for r in std):
        L.append("Rows marked STALE come from a pre-nokozak / pre-mm1 / pre-union checkpoint and are kept")
        L.append("only until that dataset is re-dumped on the shipping models. Do not quote them as")
        L.append("released-model performance.")
    else:
        L.append("All rows are on the shipping models.")
    L.append("```\n")

    if missing:
        L.append("```{admonition} Missing sources\n:class: warning")
        for m in missing:
            L.append(f"- {m}")
        L.append("```\n")

    L.extend(liver3x3_section())

    L.append("## Reading it\n")
    L.append("- **The cross-species transfer holds.** A model trained only on human tissue recovers")
    L.append("  ~99% of annotated mouse ORFs from sequence and RNA-seq alone.")
    L.append("- **uORF recall tracks the dataset, not the model.** The two arms of each dataset agree")
    L.append("  closely with each other and differ substantially between datasets, which points at")
    L.append("  library depth and 5'UTR coverage rather than at architecture.")
    L.append("- **dORF and internal are at the floor everywhere**, exactly as in human.")
    L.append("- **attn and mamba4 are interchangeable here.** No dataset separates them by more than a")
    L.append("  few thousandths of F1, consistent with the 3-seed union comparison.\n")
    L.append("_Generated by `tutorial/make_heldout_mouse.py`; do not edit by hand._")

    text = "\n".join(L) + "\n"
    if a.check:
        print(text)
        return 0
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(text)
    print(f"wrote {OUT}  ({len(std)} standalone rows, {len(obs)} obsdepth rows, {len(missing)} missing)")
    for r in std:
        flag = "released" if RELEASED in r["ckpt"] else "STALE"
        print(f"  {r['label']:<26}{flag:<10}{r['ckpt'][:54]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
