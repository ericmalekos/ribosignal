#!/usr/bin/env python3
"""Compare predicted vs observed RiboCode ORF calls for every cross-species arm.

Reuses scripts/compare_dropin_calls.build_loader rather than re-parsing the collapsed
files, because the filtering is the analysis: test-transcript restriction, raw
pval_combined <= 0.05, ORF >= 90 nt, and on the PREDICTED side a >= 0.5x-uniform
enrichment gate that removes the diffuse softmax leak into 3'UTRs. Hand-loading the same
files once produced a 4x overcount of model-specific calls and flipped a conclusion; the
loader exists so that cannot happen twice.

Keying is (gene_id, ORF_gstop) -- the genomic ORF locus -- not ORF_ID. One genomic ORF
appears once per transcript it sits on, and which isoform RiboCode picks as the collapse
representative shifts when the transcript universe changes.

Per species, tx2gene MUST be that species' own table. The loader raises if none of the
test transcripts match, which is what catches a wrong-assembly map.

Three call sets per arm, all from the same caller and statistics:
  real            observed P-site profile      -> the reference set
  pred_obsdepth   predicted shape, observed depth
  pred_preddepth  predicted shape, predicted depth (fully de novo)

Usage: python3 scripts/xspecies/compare_xspecies_orf_calls.py [--arch attn]
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
sys.path.insert(0, str(NEW / "scripts"))
from compare_dropin_calls import build_loader  # noqa: E402

RUN = "orf_v2_{a}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes"
# dataset label -> the species whose annotation tables it was built against
SPECIES = {
    "xsp_yeast": "yeast", "xsp_celegans": "celegans", "xsp_zebrafish": "zebrafish",
    "xsp_gorilla": "gorilla", "xsp_chimp": "chimp", "xsp_macaque": "macaque",
    "xsp_human": "human_refseq",
}


def prf(pred: set, real: set) -> tuple[float, float, float]:
    tp = len(pred & real)
    p = tp / len(pred) if pred else float("nan")
    r = tp / len(real) if real else float("nan")
    f = 2 * p * r / (p + r) if pred and real and (p + r) else float("nan")
    return p, r, f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", default="attn")
    ap.add_argument("--tsv", type=Path, default=NEW / "results" / "xspecies_orf_calls.tsv")
    args = ap.parse_args()

    rows = []
    for ds, sp in SPECIES.items():
        dd = NEW / "results" / "xspecies_dropin" / f"{ds}_{args.arch}"
        prof = (NEW / "results" / "loto" / RUN.format(a=args.arch)
                / f"dropin_{ds}" / "pred_profiles.npz")
        need = [dd / f"{v}_collapsed.txt" for v in
                ("real", "pred_obsdepth", "pred_preddepth")]
        if not prof.exists() or not all(p.exists() for p in need):
            print(f"  skip {ds} (incomplete)", file=sys.stderr)
            continue
        t2g = NEW / "data" / "annot" / sp / "tx2biotype.tsv"
        lc, keep_tx, keep_genes = build_loader(prof, key="genomic", tx2gene=t2g)
        real = lc(need[0])
        po = lc(need[1], is_pred=True)
        pp = lc(need[2], is_pred=True)
        types = Counter(v["type"] for v in real.values())
        n_ann = types.get("annotated", 0)
        for name, pred in (("pred_obsdepth", po), ("pred_preddepth", pp)):
            p, r, f = prf(set(pred), set(real))
            # Overall F1 is dominated by annotated CDS, which is the easy majority
            # (65-100% of every real set). Recall broken out per ORF type is what
            # distinguishes "recovers the annotated proteome" from "finds the
            # non-canonical events", and a false-positive count on the types ABSENT
            # from the real set is where the known 3'UTR over-call shows up.
            per_type = {}
            for t in ("annotated", "uORF", "Overlap_uORF", "dORF", "Overlap_dORF",
                      "internal", "novel"):
                keys = {k for k, v in real.items() if v["type"] == t}
                per_type[f"recall_{t}"] = (round(len(keys & set(pred)) / len(keys), 3)
                                           if keys else "")
                per_type[f"n_{t}"] = len(keys)
            fp = set(pred) - set(real)
            fp_types = Counter(pred[k]["type"] for k in fp)
            rows.append({
                "pack": ds.replace("xsp_", ""), "arch": args.arch, "variant": name,
                "test_genes": len(keep_genes), "real_calls": len(real),
                "real_annotated": n_ann,
                "real_pct_annotated": round(100 * n_ann / len(real), 1) if real else 0,
                "pred_calls": len(pred), "tp": len(set(pred) & set(real)),
                "precision": round(p, 3), "recall": round(r, 3), "f1": round(f, 3),
                "fp_total": len(fp),
                "fp_dORF": fp_types.get("dORF", 0) + fp_types.get("Overlap_dORF", 0),
                "fp_uORF": fp_types.get("uORF", 0) + fp_types.get("Overlap_uORF", 0),
                "fp_novel": fp_types.get("novel", 0),
                **per_type,
            })
        print(f"  {ds:<14} test_genes={len(keep_genes):>7,} real={len(real):>7,} "
              f"({100*n_ann/len(real):.0f}% annotated)", file=sys.stderr)

    if not rows:
        print("nothing to compare", file=sys.stderr)
        return 1
    hdr = (f"  {'pack':<11}{'variant':<16}{'real':>8}{'ann%':>7}{'pred':>8}"
           f"{'TP':>8}{'prec':>7}{'rec':>7}{'F1':>7}")
    print("\n" + hdr)
    print("  " + "-" * (len(hdr) - 2))
    for r in rows:
        print(f"  {r['pack']:<11}{r['variant']:<16}{r['real_calls']:>8,}"
              f"{r['real_pct_annotated']:>7.1f}{r['pred_calls']:>8,}{r['tp']:>8,}"
              f"{r['precision']:>7.3f}{r['recall']:>7.3f}{r['f1']:>7.3f}")

    print("\n  Recall by ORF type (n in the real set), pred_preddepth (standalone):")
    h2 = (f"  {'pack':<11}{'annotated':>16}{'uORF':>14}{'dORF':>13}"
          f"{'internal':>13}{'novel':>14}{'FP dORF':>9}")
    print(h2); print("  " + "-" * (len(h2) - 2))
    for r in rows:
        if r["variant"] != "pred_preddepth":
            continue
        def cell(t):
            n = r.get(f"n_{t}", 0)
            v = r.get(f"recall_{t}", "")
            return f"{v:.3f} ({n:,})" if v != "" else f"-  ({n:,})"
        print(f"  {r['pack']:<11}{cell('annotated'):>16}{cell('uORF'):>14}"
              f"{cell('dORF'):>13}{cell('internal'):>13}{cell('novel'):>14}"
              f"{r['fp_dORF']:>9,}")
    args.tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"\n  -> {args.tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
