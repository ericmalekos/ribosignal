#!/usr/bin/env python3
"""Generate tutorial/source/data/heldout-human.md from the result JSONs.

GENERATED, NOT HAND-WRITTEN, on purpose. The previous version of this content lived as prose in
`data/test-sets.md` with values like "~0.93" and "n/a", and those placeholders silently outlived the
runs they described. Worse, the numbers that WERE concrete came from
`orf_v2_attn_onehot_holdout_Hepatocytes` -- a pre-nokozak, pre-mm1, pre-union checkpoint -- while the
project ships the union models. Directory names encode the DATASET, not the checkpoint, so that
staleness is invisible unless you read the npz `meta` field.

So: every number here is read from a JSON at build time, and every row records the checkpoint it came
from. A stale source shows up as a stale checkpoint string in the page rather than as a plausible
wrong number.

MOUSE IS A SEPARATE PAGE (`heldout-mouse.md`) and a separate generator, because the mouse held-out
sets are mid-refresh on the released models.

Usage: make_heldout_human.py            # writes the page
       make_heldout_human.py --check    # print what it would emit, write nothing
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tutorial/source/data/heldout-human.md"

# (label, what it tests, dropin_metrics.json, npz whose meta names the checkpoint)
SOURCES = [
    ("Hepatocytes (attn)", "held-out TISSUE, same study",
     "results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin/dropin_metrics.json",
     "results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin/pred_profiles.npz"),
    ("Hepatocytes (mamba4)", "held-out TISSUE, same study",
     "results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin/dropin_metrics.json",
     "results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin/pred_profiles.npz"),
    ("Ruiz-Orera (attn)", "held-out STUDY, independent lab",
     "results/heldout/human_ruizorera/released_attn/dropin/dropin_metrics.json",
     "results/heldout/human_ruizorera/released_attn/dropin/pred_profiles.npz"),
    ("Ruiz-Orera (mamba4)", "held-out STUDY, independent lab",
     "results/heldout/human_ruizorera/released_mamba4/dropin/dropin_metrics.json",
     "results/heldout/human_ruizorera/released_mamba4/dropin/pred_profiles.npz"),
]
CLS = ["annotated", "uORF", "Overlap_uORF", "novel", "dORF", "Overlap_dORF", "internal"]
ARM = "pred_preddepth_vs_real"      # standalone: predicted shape AND predicted depth, no Ribo-seq


def ckpt(npz):
    p = ROOT / npz
    if not p.exists():
        return "UNKNOWN (no npz)"
    z = np.load(p, allow_pickle=True)
    return str(np.atleast_1d(z["meta"])[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    rows, missing = [], []
    for label, tests, mj, npz in SOURCES:
        f = ROOT / mj
        if not f.exists():
            missing.append(f"{label}: {mj}")
            continue
        j = json.loads(f.read_text())
        b = j.get(ARM)
        if b is None:
            missing.append(f"{label}: no {ARM} in {mj}")
            continue
        rows.append({"label": label, "tests": tests, "ckpt": ckpt(npz), "m": b})

    released = "union_noBrain_nokozak_mm1"
    L = []
    L.append("# Held-out human data\n")
    L.append("Human generalisation rests on **one genuinely independent study**. Everything else human is")
    L.append("either the same study with a tissue withheld, or a proteogenomics cell line. That is worth")
    L.append("stating plainly: if a reviewer discounts Ruiz-Orera, there is no second human study behind it.\n")
    L.append("```{list-table}\n:header-rows: 1\n:widths: 26 30 44\n")
    L.append("* - Dataset\n  - Role\n  - What it does and does not test")
    L.append("* - **Hepatocytes**\n  - held-out tissue (LOTO)\n  - Same Chothani panel, same protocol, same lab, "
             "one tissue withheld from training. Tests cross-TISSUE generalisation only.")
    L.append("* - **Ruiz-Orera**\n  - held-out study\n  - iPSC-cardiomyocyte, different lab, different cell type. "
             "The only true cross-STUDY human test.")
    L.append("* - **B721.221**\n  - ground truth\n  - Deepest Ribo-seq in the project (18.8 M P-sites). Used as the "
             "measured-translation reference, not as a training holdout. See {doc}`/benchmarks`.")
    L.append("* - **HBL-1**\n  - proteogenomics\n  - Weakest library here (2 libs, 128 K P-sites, 34 nt mode, lowest "
             "f0). Weight its result below the others.")
    L.append("```\n")

    L.append("## Standalone drop-in: precision and recall by ORF class\n")
    L.append("`pred_preddepth` = predicted profile shape **and** predicted depth, i.e. no Ribo-seq for the query")
    L.append("sample at any point. Scored against RiboCode run on that pack's real counts, so the transcript space")
    L.append("is identical on both sides by construction.\n")
    hdr = "| dataset | n_ref | P | R | F1 | " + " | ".join(c.replace("Overlap_", "Ovl_") for c in CLS) + " |"
    L.append(hdr)
    L.append("|---|--:|--:|--:|--:|" + "--:|" * len(CLS))
    for r in rows:
        m = r["m"]
        pt = m.get("recall_per_type", {})
        cells = " | ".join(f"{pt[c]['recall']:.3f}" if c in pt else "--" for c in CLS)
        L.append(f"| {r['label']} | {m['n_real']:,} | {m['precision']:.3f} | {m['recall']:.3f} | "
                 f"{m['f1']:.3f} | {cells} |")
    L.append("")
    L.append("Class columns are **recall** (the class breakdown is recall-only in the source metrics).\n")

    L.append("```{admonition} Checkpoint provenance\n:class: important")
    L.append("Each row's source checkpoint, read from the npz `meta` field at build time:\n")
    for r in rows:
        tag = "RELEASED" if released in r["ckpt"] else "**STALE -- predates the released models**"
        L.append(f"- `{r['label']}` -> `{r['ckpt']}` ({tag})")
    L.append("")
    if any(released not in r["ckpt"] for r in rows):
        L.append("Rows marked STALE come from a pre-nokozak / pre-mm1 / pre-union checkpoint and are kept only")
        L.append("until that dataset is re-dumped on the shipping models. Do not quote them as released-model")
        L.append("performance.")
    else:
        L.append("All rows are on the shipping models.")
    L.append("```\n")

    if missing:
        L.append("```{admonition} Missing sources\n:class: warning")
        for m in missing:
            L.append(f"- {m}")
        L.append("```\n")

    L.append("## Reading it\n")
    L.append("- **Annotated ORFs transfer essentially intact** across tissue and study.")
    L.append("- **dORF and internal are at the floor everywhere.** They are the hardest classes for the caller,")
    L.append("  not just for the model.")
    L.append("- **uORF is the most dataset-dependent class**, and that spread is the largest unexplained gap in")
    L.append("  the table.\n")
    L.append("_Generated by `tutorial/make_heldout_human.py`; do not edit by hand._")

    text = "\n".join(L) + "\n"
    if a.check:
        print(text)
        return 0
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(text)
    print(f"wrote {OUT}  ({len(rows)} rows, {len(missing)} missing)")
    for r in rows:
        flag = "released" if released in r["ckpt"] else "STALE"
        print(f"  {r['label']:<24}{flag:<10}{r['ckpt'][:54]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
