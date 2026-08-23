#!/usr/bin/env python3
"""Task #86: did training on CANONICAL packs change ORF CALLING?

Consumes `eval_canon_orfcalls.sbatch`. Five arms, all dumped against the SAME canonical Hepatocytes
pack, so the only thing differing between the deployed and final-recipe arms is the TRAINING substrate:

    canon/attn          canon/attn_nof012          canon/mamba4
    deployed/attn       deployed/mamba4            (re-dumped on the final-recipe pack)

WHY THIS AND NOT PROFILE METRICS. Profile Pearson is a wash between the two training substrates
(0.659 canon vs comparable deployed), and the standing rule is that a retrain is judged on how ORF
CALLS change (feedback_orf_call_eval_standing). Nothing in test_metrics.json answers "do cleaner
targets produce better non-canonical ORF calls".

THE `real` ARM IS THE CONTROL, AND IT IS ASSERTED, NOT ASSUMED. It depends only on the pack, so all
five arms must call the SAME ORFs. If they diverge, the arms are not on the substrate this script
claims and every downstream number is void -- a hard failure, not a warning.

THE CONTROL IS ON THE CALL SET, NOT ON FILE BYTES. An md5 of real_collapsed.txt FAILS here, and the
failure is spurious: RiboCode writes `alt_ORF_type` by joining a Python set, whose iteration order is
not stable between processes, so 4 of 18,813 rows differ as `internal,uORF` vs `uORF,internal` while
describing the identical ORF. The ORF_ID sets are identical and `ORF_type` -- the only classification
anything here scores on -- agrees everywhere. So the check compares the ORF_ID set and ORF_type, and
reports alt_ORF_type ordering separately as known cosmetic non-determinism.

FILTERING IS THE PROJECT'S SINGLE PATH. `compare_dropin_calls.build_loader` (pval<=0.05, ORF>=90 nt
= exactly 30 aa, predicted enrichment>=0.5x uniform) and genomic keying on (gene_id, ORF_gstop), so
choosing a different representative isoform of the same ORF is not scored as a disagreement.
`prf` is imported from score_human_orf_calls rather than reimplemented.

READ THE PER-CLASS TABLE, NOT THE OVERALL F1. Overall F1 is dominated by annotated CDS and reads
~0.9 for everything; the question lives in the non-canonical classes.

  python3 scripts/compare_canon_orfcalls.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from compare_dropin_calls import build_loader          # noqa: E402
from score_human_orf_calls import prf                  # noqa: E402

NEW = HERE.parent
CLASSES = ["annotated", "novel", "uORF", "internal", "Overlap_uORF", "dORF", "Overlap_dORF"]
NONCANON = [c for c in CLASSES if c != "annotated"]
ARMS = ["pred_obsdepth", "pred_preddepth"]

RUNS = [
    ("canon/attn", "results/loto_canon/orf_v2_attn_onehot_canon_noBrain_nokozak_holdout_Hepatocytes/dropin"),
    ("canon/attn_nof012", "results/loto_canon/orf_v2_attn_onehot_canon_noBrain_nokozak_nof012_holdout_Hepatocytes/dropin"),
    ("canon/mamba4", "results/loto_canon/orf_v2_mamba4_onehot_canon_noBrain_nokozak_holdout_Hepatocytes/dropin"),
    ("deployed/attn", "results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin_canonpack"),
    ("deployed/mamba4", "results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes/dropin_canonpack"),
]


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tx2gene", default=str(NEW / "data" / "tx2biotype.tsv"))
    ap.add_argument("--out", default=str(NEW / "results/loto_canon/orfcall_comparison"))
    a = ap.parse_args()
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    present = [(lbl, NEW / rel) for lbl, rel in RUNS if (NEW / rel / "real_collapsed.txt").exists()]
    missing = [lbl for lbl, rel in RUNS if not (NEW / rel / "real_collapsed.txt").exists()]
    if missing:
        print(f"  NOT YET AVAILABLE: {missing}")
    if not present:
        raise SystemExit("no arms available yet -- is eval_canon_orfcalls.sbatch still running?")

    # ---- control: every arm must have called the SAME ORFs ------------------------------------
    def real_calls(d):
        """-> ({ORF_ID: ORF_type}, {ORF_ID: alt_ORF_type})"""
        import csv as _csv
        t, alt = {}, {}
        with open(d / "real_collapsed.txt") as fh:
            for row in _csv.DictReader(fh, delimiter="\t"):
                t[row["ORF_ID"]] = row["ORF_type"]
                alt[row["ORF_ID"]] = row.get("alt_ORF_type", "")
        return t, alt

    ref_lbl, ref_d = present[0]
    ref_t, ref_alt = real_calls(ref_d)
    print(f"  CONTROL: `real` call set across {len(present)} arms (reference = {ref_lbl}, "
          f"{len(ref_t):,} calls)")
    bad, alt_only = [], []
    for lbl, d in present[1:]:
        t, alt = real_calls(d)
        if set(t) != set(ref_t):
            bad.append(f"{lbl}: ORF_ID set differs "
                       f"(+{len(set(t)-set(ref_t))} / -{len(set(ref_t)-set(t))})")
            continue
        mism = [k for k in t if t[k] != ref_t[k]]
        if mism:
            bad.append(f"{lbl}: ORF_type differs on {len(mism)} calls (e.g. {mism[:2]})")
        n_alt = sum(1 for k in alt if alt[k] != ref_alt[k])
        # Set-ordering only, if the comma-split contents match.
        reorder = sum(1 for k in alt if alt[k] != ref_alt[k]
                      and set(alt[k].split(",")) == set(ref_alt[k].split(",")))
        if n_alt:
            alt_only.append(f"{lbl}: {n_alt} alt_ORF_type differ ({reorder} pure reordering)")
        print(f"    {lbl:<20} ORF_IDs identical, ORF_type identical"
              + (f", {n_alt} alt_ORF_type reordered" if n_alt else ""))
    if bad:
        for b in bad:
            print(f"    FAIL {b}")
        raise SystemExit("ABORT: the `real` arm differs between runs in a way that matters. They are "
                         "NOT on one substrate, so no comparison below would be valid.")
    if alt_only:
        print("    (alt_ORF_type ordering is RiboCode serialising a Python set; nothing here scores "
              "on it)")

    # ---- score every arm against the shared final-recipe reference ------------------------------
    rows = []
    for lbl, d in present:
        lc, _, _ = build_loader(d / "pred_profiles.npz", tx2gene=a.tx2gene)
        real = lc(d / "real_collapsed.txt", is_pred=False)
        for arm in ARMS:
            f = d / f"{arm}_collapsed.txt"
            if not f.exists():
                continue
            pred = lc(f, is_pred=True)
            o = prf(pred, real)
            rec = dict(run=lbl, arm=arm, scope="overall", **o)
            rows.append(rec)
            nc = prf(pred, real, subset=tuple(NONCANON))
            rows.append(dict(run=lbl, arm=arm, scope="non-canonical", **nc))
            for c in CLASSES:
                rows.append(dict(run=lbl, arm=arm, scope=c, **prf(pred, real, subset=(c,))))

    cols = ["run", "arm", "scope", "precision", "recall", "f1", "n_pred", "n_ref", "tp"]
    with open(out / "orfcall_comparison.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c])
                               for c in cols) + "\n")
    (out / "orfcall_comparison.json").write_text(json.dumps(
        {"reference": "shared canonical Hepatocytes pack (data/packed_canon_Hepatocytes)",
         "real_control": "ORF_ID set + ORF_type identical across arms", "n_arms": len(present), "missing": missing,
         "filter": "pval<=0.05, ORF>=90nt (30 aa), pred enrichment>=0.5x; key (gene_id, ORF_gstop)",
         "rows": rows}, indent=2))

    # ---- the answer ---------------------------------------------------------------------------
    def get(run, arm, scope):
        return next((r for r in rows if (r["run"], r["arm"], r["scope"]) == (run, arm, scope)), None)

    print(f"\n  STANDALONE ARM (pred_preddepth), scored on the shared final-recipe reference")
    print(f"    {'run':<20}{'overall F1':>11}{'annot F1':>10}{'non-canon F1':>14}"
          f"{'nc recall':>11}{'nc prec':>9}{'n_ref':>9}")
    for lbl, _ in present:
        o, an, nc = (get(lbl, "pred_preddepth", s) for s in ("overall", "annotated", "non-canonical"))
        if not o:
            continue
        print(f"    {lbl:<20}{o['f1']:>11.4f}{an['f1']:>10.4f}{nc['f1']:>14.4f}"
              f"{nc['recall']:>11.4f}{nc['precision']:>9.4f}{nc['n_ref']:>9,}")

    for arch in ("attn", "mamba4"):
        c, dpl = get(f"canon/{arch}", "pred_preddepth", "non-canonical"), \
                 get(f"deployed/{arch}", "pred_preddepth", "non-canonical")
        if c and dpl:
            d_f1 = c["f1"] - dpl["f1"]
            verdict = ("final-recipe training HELPS" if d_f1 > 0.005 else
                       "final-recipe training HURTS" if d_f1 < -0.005 else
                       "NO MATERIAL DIFFERENCE (|delta| <= 0.005)")
            print(f"\n  {arch}: non-canonical F1 canon {c['f1']:.4f} vs deployed {dpl['f1']:.4f} "
                  f"({d_f1:+.4f}) -> {verdict}")
    print(f"\n  wrote {out}/orfcall_comparison.{{tsv,json}}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
