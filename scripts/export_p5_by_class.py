#!/usr/bin/env python3
"""Per-ORF-class recall for all four human held-outs -> figures/P5_heldout_human/P5_by_class_<model>.tsv

WHY. P5 already reads per-class rows and pools them into annotated / non-canonical, discarding the
breakdown. The only human per-class file that ever reached the package was A2's hepatocyte one, so a
class-resolved panel could show exactly one dataset -- the easiest of the four (same study, ~99.7% of
its transcripts seen in another training tissue) and, until 2026-08-16, the only one still off-recipe.
No new computation is needed; the rows exist in two shapes and are simply not exported.

  scored arms (THP-1, CAR-T)      results/human_orf_calls*/scored/human_by_class.tsv
  dropin arms (hepatocyte, iPSC)  <run>/dropin*/dropin_metrics.json -> recall_per_type

RECALL ONLY, DELIBERATELY. `score_human_orf_calls.prf(pred, ref, subset=(c,))` filters BOTH sides by
class, so an ORF both callers found but labelled differently (uORF vs Overlap_uORF, internal vs uORF)
falls out of every class bucket. That makes per-class PRECISION and F1 unsound -- they divide by a
n_pred that excludes those calls. Recall is safe: its denominator is the reference's own class
assignment, which is unaffected. n_pred is emitted for context but must NOT be used to derive
precision.

EACH ROW CARRIES `on_final_recipe`, because these four arms reached the final recipe on different
dates and a class table with a mixed substrate invites exactly the comparison it should prevent.

  python3 scripts/export_p5_by_class.py --model attn
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib

NEW = pathlib.Path(__file__).resolve().parents[1]
CLASSES = ["annotated", "novel", "uORF", "internal", "Overlap_uORF", "dORF", "Overlap_dORF"]
ARM = "pred_preddepth"

# (label, kind, path template). Dropin paths prefer the final-recipe dump, mirroring P5 itself.
DROPIN = {
    "hepatocyte": ["results/loto/orf_v2_{m}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin_canonpack/dropin_metrics.json",
                   "results/loto/orf_v2_{m}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin/dropin_metrics.json"],
    "iPSC-CM": ["results/heldout/human_ruizorera/released_{m}/dropin_canonpack/dropin_metrics.json",
                "results/heldout/human_ruizorera/released_{m}/dropin/dropin_metrics.json"],
}
SCORED = {"THP-1": "gse208041", "CAR-T": "cart"}
SCORED_TSV = ["results/human_orf_calls_canon/scored/human_by_class.tsv",
              "results/human_orf_calls/scored/human_by_class.tsv"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="attn", choices=["attn", "mamba4"])
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = pathlib.Path(a.out) if a.out else (
        NEW / f"figures/P5_heldout_human/P5_by_class_{a.model}.tsv")

    rows = []
    # ---- dropin arms: recall_per_type ----------------------------------------------------------
    for lab, cands in DROPIN.items():
        path = next((c.format(m=a.model) for c in cands if (NEW / c.format(m=a.model)).exists()), None)
        if not path:
            print(f"  MISSING {lab}: no dropin_metrics.json"); continue
        j = json.loads((NEW / path).read_text())
        blk = j.get(f"{ARM}_vs_real", {})
        per = blk.get("recall_per_type") or {}
        if not per:
            print(f"  MISSING {lab}: no recall_per_type in {path}"); continue
        final = "dropin_canonpack" in path
        for c in CLASSES:
            d = per.get(c)
            if not d:
                continue
            # recall_per_type entries are {recovered, n_real, recall} (A2's schema).
            tp = d.get("recovered", d.get("tp"))
            nref = d.get("n_real", d.get("n_ref"))
            rec = d.get("recall")
            if rec is None and tp is not None and nref:
                rec = tp / nref
            rows.append(dict(dataset=lab, model=a.model, arm=ARM, on_final_recipe=final,
                             orf_class=c, tp=tp, n_pred="", n_ref=nref, recall=rec,
                             source=path))
        print(f"  {lab:<12} {len([r for r in rows if r['dataset']==lab])} classes  "
              f"(on_final_recipe={final})")

    # ---- scored arms: human_by_class.tsv -------------------------------------------------------
    tsv = next((t for t in SCORED_TSV if (NEW / t).exists()), None)
    if tsv:
        allrows = list(csv.DictReader(open(NEW / tsv), delimiter="\t"))
        final = "canon" in tsv
        for lab, ds in SCORED.items():
            sel = [r for r in allrows if r["dataset"] == ds and r["model"] == a.model
                   and r["arm"] == ARM]
            if not sel:
                print(f"  MISSING {lab}: not in {tsv}"); continue
            for r in sel:
                if r["orf_class"] not in CLASSES:
                    continue
                rows.append(dict(dataset=lab, model=a.model, arm=ARM, on_final_recipe=final,
                                 orf_class=r["orf_class"], tp=int(r["tp"]),
                                 n_pred=int(r["n_pred"]), n_ref=int(r["n_ref"]),
                                 recall=float(r["recall"]), source=tsv))
            print(f"  {lab:<12} {len(sel)} classes  (on_final_recipe={final})")

    cols = ["dataset", "model", "arm", "on_final_recipe", "orf_class", "tp", "n_pred", "n_ref",
            "recall", "source"]
    with open(out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            r = dict(r)
            if isinstance(r.get("recall"), float):
                r["recall"] = f"{r['recall']:.4f}"
            w.writerow(r)
    print(f"\n  wrote {out}  ({len(rows)} rows, "
          f"{len({r['dataset'] for r in rows})} datasets)")
    print("  NOTE: recall only. Per-class precision/F1 are unsound here (prf filters both sides by "
          "class); n_pred is context, not a precision denominator.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
