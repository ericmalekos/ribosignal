#!/usr/bin/env python3
"""Task 68: does better RNA-seq change ORF-call precision and recall, and through which path?

RNA-seq enters the model twice and the two paths are separable:

  (1) UNIVERSE  -- salmon TPM >= 1 decides which transcripts exist at all
  (2) COVERAGE  -- the per-nt RNA depth track the model consumes as an input channel

Three arms, Ribo-seq held FIXED across all three (the same 5 Janich P-site hd5, read from arm a's
provenance rather than re-listed, so the arms cannot drift):

  arm a  universe Janich-decon (14,887 tx)   coverage Janich       results/liver_released/*_mouse_janich_liver_decon
  arm b  universe Janich-decon (14,887 tx)   coverage PRJEB34766   results/rna_quality_factorial/*_rnaq_b_*
  arm c  universe PRJEB34766   (16,474 tx)   coverage PRJEB34766   results/rna_quality_factorial/*_rnaq_c_*

  a vs b  = pure COVERAGE effect (universe held byte-identical, so the ORF track and transcript space
            are identical and the contrast is not confounded by what got called)
  b vs c  = pure UNIVERSE effect
  a vs c  = the full "swap the RNA-seq experiment" effect

READ THE COMPARISON CAREFULLY. Each arm is scored against ITS OWN pack's `real_collapsed.txt`, so the
reference moves between arms b and c (different universes) even though the Ribo-seq did not. Between a
and b the reference is identical, which is why a-vs-b is the clean contrast and b-vs-c is not a
like-for-like F1 comparison -- for b vs c the informative columns are n_ref and the per-class counts,
not F1. Both are printed; the summary calls out which is which.

Filtering is imported from `compare_dropin_calls.build_loader` -- the same convention as every other
drop-in number in the project (test-tx restriction, pval <= 0.05, ORF >= 90 nt, predicted enrichment
>= 0.5x uniform). See score_released_two_arm.py for why that is imported rather than copied.

Usage: score_rna_quality_factorial.py [--key genomic|transcript] [--out <json>]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_dropin_calls import build_loader  # noqa: E402

NONCANON = {"uORF", "Overlap_uORF", "dORF", "Overlap_dORF", "novel", "internal"}
TX2GENE = ROOT / "data/mouse_tx2biotype.tsv"

# (arm, model, dump dir)
ARMS = [
    ("a  uni=Janich  cov=Janich", "attn", "results/liver_released/attn_mouse_janich_liver_decon"),
    ("a  uni=Janich  cov=Janich", "mamba4", "results/liver_released/mamba4_mouse_janich_liver_decon"),
    ("b  uni=Janich  cov=34766", "attn",
     "results/rna_quality_factorial/attn_rnaq_b_uniJanich_cov34766"),
    ("b  uni=Janich  cov=34766", "mamba4",
     "results/rna_quality_factorial/mamba4_rnaq_b_uniJanich_cov34766"),
    ("c  uni=34766   cov=34766", "attn",
     "results/rna_quality_factorial/attn_rnaq_c_uni34766_cov34766"),
    ("c  uni=34766   cov=34766", "mamba4",
     "results/rna_quality_factorial/mamba4_rnaq_c_uni34766_cov34766"),
]


def prf(pred, ref, subset=None):
    if subset is not None:
        pred = {k: v for k, v in pred.items() if v["type"] in subset}
        ref = {k: v for k, v in ref.items() if v["type"] in subset}
    m = len(set(pred) & set(ref))
    p = m / len(pred) if pred else 0.0
    r = m / len(ref) if ref else 0.0
    return {"n_pred": len(pred), "n_ref": len(ref), "n_match": m, "precision": p, "recall": r,
            "f1": 2 * p * r / (p + r) if p + r else 0.0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--key", default="genomic", choices=["genomic", "transcript"])
    ap.add_argument("--out", default=str(ROOT / "results/rna_quality_factorial/factorial_scores.json"))
    a = ap.parse_args()

    rows, missing = [], []
    for arm, model, d in ARMS:
        dd = ROOT / d
        need = {n: dd / f"{n}_collapsed.txt" for n in ("real", "pred_obsdepth", "pred_preddepth")}
        prof = dd / "pred_profiles.npz"
        absent = [str(p) for p in list(need.values()) + [prof] if not p.exists()]
        if absent:
            missing.append(f"{arm} / {model}: missing {', '.join(absent)}")
            continue
        lc, keep_tx, _ = build_loader(prof, key=a.key, tx2gene=str(TX2GENE))
        R = lc(need["real"])
        for variant in ("pred_obsdepth", "pred_preddepth"):
            P = lc(need[variant], is_pred=True)
            rows.append({"arm": arm, "model": model, "variant": variant, "n_test_tx": len(keep_tx),
                         "overall": prf(P, R), "CDS": prf(P, R, {"annotated"}),
                         "non-canonical": prf(P, R, NONCANON)})

    hdr = (f"{'arm':<26}{'model':<8}{'variant':<16}{'nTx':>8}{'nPred':>8}{'nRef':>8}"
           f"{'P':>7}{'R':>7}{'F1':>7} | {'ncP':>6}{'ncR':>6}{'ncF1':>6}")
    print(hdr); print("-" * len(hdr))
    prev = None
    for r in rows:
        if prev and prev != r["arm"]:
            print()
        prev = r["arm"]
        o, n = r["overall"], r["non-canonical"]
        print(f"{r['arm']:<26}{r['model']:<8}{r['variant']:<16}{r['n_test_tx']:>8,}"
              f"{o['n_pred']:>8,}{o['n_ref']:>8,}{o['precision']:>7.3f}{o['recall']:>7.3f}"
              f"{o['f1']:>7.3f} | {n['precision']:>6.3f}{n['recall']:>6.3f}{n['f1']:>6.3f}")

    # --- the two contrasts, stated separately, because only one is like-for-like ---
    def get(arm_pre, model, variant):
        for r in rows:
            if r["arm"].startswith(arm_pre) and r["model"] == model and r["variant"] == variant:
                return r
        return None

    print("\n=== a vs b: PURE COVERAGE effect (universe held identical -> same reference, like-for-like) ===")
    for model in ("attn", "mamba4"):
        for variant in ("pred_obsdepth", "pred_preddepth"):
            A, B = get("a", model, variant), get("b", model, variant)
            if not (A and B):
                continue
            if A["overall"]["n_ref"] != B["overall"]["n_ref"]:
                print(f"  WARNING {model}/{variant}: n_ref differs ({A['overall']['n_ref']:,} vs "
                      f"{B['overall']['n_ref']:,}) -- the universe was NOT held")
            print(f"  {model:<8}{variant:<16} F1 {A['overall']['f1']:.3f} -> {B['overall']['f1']:.3f} "
                  f"({B['overall']['f1'] - A['overall']['f1']:+.3f})   "
                  f"P {A['overall']['precision']:.3f} -> {B['overall']['precision']:.3f} "
                  f"({B['overall']['precision'] - A['overall']['precision']:+.3f})   "
                  f"R {A['overall']['recall']:.3f} -> {B['overall']['recall']:.3f} "
                  f"({B['overall']['recall'] - A['overall']['recall']:+.3f})   "
                  f"ncF1 {A['non-canonical']['f1']:.3f} -> {B['non-canonical']['f1']:.3f} "
                  f"({B['non-canonical']['f1'] - A['non-canonical']['f1']:+.3f})")

    print("\n=== b vs c: UNIVERSE effect -- NOT like-for-like (the reference set changes with the "
          "universe) ===")
    print("    Read n_ref and n_pred, not the F1 delta: arm c's reference is RiboCode over a larger")
    print("    transcript space, so its F1 is computed against a different denominator.")
    for model in ("attn", "mamba4"):
        for variant in ("pred_obsdepth", "pred_preddepth"):
            B, C = get("b", model, variant), get("c", model, variant)
            if not (B and C):
                continue
            print(f"  {model:<8}{variant:<16} tx {B['n_test_tx']:,} -> {C['n_test_tx']:,}   "
                  f"nRef {B['overall']['n_ref']:,} -> {C['overall']['n_ref']:,}   "
                  f"nPred {B['overall']['n_pred']:,} -> {C['overall']['n_pred']:,}   "
                  f"(F1 {B['overall']['f1']:.3f} vs {C['overall']['f1']:.3f}, different denominators)")

    if missing:
        print("\nMISSING:")
        for m in missing:
            print(f"  {m}")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps({"key": a.key, "rows": rows, "missing": missing}, indent=2))
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
