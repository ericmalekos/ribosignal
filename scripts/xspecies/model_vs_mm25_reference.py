#!/usr/bin/env python3
"""Score the MODEL's predicted ORF calls against two different ground truths: the mm1
observed call set and the mm25 one.

THE QUESTION. Every precision/recall/F1 this project reports uses mm1 observed calls as
truth. If the truth set changes (mm25 admits multimappers), does the model look better or
worse? A model whose apparent accuracy is stable across both references is being measured
on real signal; one that moves a lot is partly being scored on posture artifacts.

The model predictions are UNCHANGED between the two rows: the same
pred_preddepth_collapsed.txt from the same drop-in run. Only the reference moves. So any
difference is a property of the truth set, not of the model.

Both references pass through compare_dropin_calls.build_loader with identical settings
(test-tx restriction from the same pred_profiles.npz, pval <= 0.05, ORF >= 90 nt, and the
0.5x enrichment gate on the predicted side only), keyed (gene_id, ORF_gstop).

CHIMP IS NOT COMPARABLE and is labelled: five of its eight runs lost 3-nt periodicity under
mm25 and were dropped by RiboCode, so its mm25 reference rests on less data than its mm1 one.

Usage: python3 scripts/xspecies/model_vs_mm25_reference.py [--arch attn]
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
sys.path.insert(0, str(NEW / "scripts"))
from compare_dropin_calls import build_loader  # noqa: E402

RUN = "orf_v2_{a}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes"
# dropin dataset -> (species tx2biotype, full-data arm label)
ARMS = {
    "xsp_yeast":     ("yeast",        "yeast_gse173654_ribo"),
    "xsp_celegans":  ("celegans",     "worm_gse52905_ribo"),
    "xsp_zebrafish": ("zebrafish",    "zf_gse46512_ribo"),
    "xsp_gorilla":   ("gorilla",      "primate_gg_ribo"),
    "xsp_chimp":     ("chimp",        "primate_pt_ribo"),
    "xsp_macaque":   ("macaque",      "primate_rm_ribo"),
    "xsp_human":     ("human_refseq", "ruizorera_hsCM_ribo"),
}
CONFOUNDED = {"xsp_chimp"}   # mm25 reference built from 3/8 runs, mm1 from 8/8


def prf(pred: set, ref: set):
    tp = len(pred & ref)
    p = tp / len(pred) if pred else float("nan")
    r = tp / len(ref) if ref else float("nan")
    f = 2 * p * r / (p + r) if pred and ref and (p + r) else float("nan")
    return p, r, f, tp


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", default="attn")
    ap.add_argument("--variant", default="pred_preddepth")
    args = ap.parse_args()

    rows = []
    for ds, (sp, arm) in ARMS.items():
        dd = NEW / "results" / "xspecies_dropin" / f"{ds}_{args.arch}"
        prof = NEW / "results" / "loto" / RUN.format(a=args.arch) / f"dropin_{ds}" / "pred_profiles.npz"
        predf = dd / f"{args.variant}_collapsed.txt"
        ref1 = NEW / "data" / "xspecies_psites" / arm / f"{arm}_collapsed.txt"
        ref25 = NEW / "data" / "xspecies_psites_mm25" / arm / f"{arm}_collapsed.txt"
        if not (prof.exists() and predf.exists() and ref1.exists() and ref25.exists()):
            print(f"  skip {ds}", file=sys.stderr)
            continue
        lc, _, _ = build_loader(prof, key="genomic",
                                tx2gene=NEW / "data" / "annot" / sp / "tx2biotype.tsv")
        pred = set(lc(predf, is_pred=True))
        for tag, ref in (("mm1", ref1), ("mm25", ref25)):
            R = set(lc(ref))
            p, r, f, tp = prf(pred, R)
            rows.append({"pack": ds.replace("xsp_", ""), "arch": args.arch,
                         "reference": tag, "n_ref": len(R), "n_pred": len(pred), "tp": tp,
                         "precision": round(p, 3), "recall": round(r, 3), "f1": round(f, 3),
                         "note": "reference from fewer runs" if (tag == "mm25" and ds in CONFOUNDED) else ""})
    if not rows:
        print("nothing to score", file=sys.stderr); return 1

    h = f"  {'pack':<11}{'ref':<7}{'n_ref':>9}{'n_pred':>9}{'TP':>9}{'prec':>7}{'rec':>7}{'F1':>7}   delta F1"
    print("\n" + h); print("  " + "-" * (len(h) - 2))
    by = {}
    for r in rows:
        by.setdefault(r["pack"], {})[r["reference"]] = r
    for pack, d in by.items():
        for tag in ("mm1", "mm25"):
            r = d[tag]
            dl = ""
            if tag == "mm25" and "mm1" in d:
                dl = f"{r['f1'] - d['mm1']['f1']:+.3f}"
                if r["note"]:
                    dl += "  (not comparable)"
            print(f"  {pack:<11}{tag:<7}{r['n_ref']:>9,}{r['n_pred']:>9,}{r['tp']:>9,}"
                  f"{r['precision']:>7.3f}{r['recall']:>7.3f}{r['f1']:>7.3f}   {dl}")
    out = NEW / "results" / f"xspecies_model_vs_mm25_reference_{args.arch}.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"\n  -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
