#!/usr/bin/env python3
"""Score the human ORF calls: model predictions against each dataset's own observed RiboCode calls.

Two datasets carry an observed arm: GSE208041 (THP-1) and GSE304796 (CAR-T). GSE39561 has none --
its library shows no 3-nt periodicity, so no P-sites can be placed -- and it is excluded here for a
second reason as well: because it borrows GSE208041's RNA and universe, its standalone predictions
are IDENTICAL to GSE208041's (verified: same call sets, differing only in the item order inside the
`alt_ORF_type` field). It is not an independent third measurement and must not be reported as one.

Filtering comes from `compare_dropin_calls.build_loader` and keying is genomic (gene_id, ORF_gstop),
matching every other drop-in number in this project. Per-class precision AND recall, because a
recall-only class table cannot separate "found them" from "called them everywhere".
"""
import argparse
import csv
import itertools
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from compare_dropin_calls import build_loader  # noqa: E402

MODELS = ["attn", "mamba4"]
# GSE39561 omitted deliberately: no observed arm, and its predictions duplicate GSE208041's.
DATASETS = ["gse208041", "cart"]
CLASSES = ["annotated", "novel", "uORF", "internal", "Overlap_uORF", "dORF", "Overlap_dORF"]


def prf(pred, ref, subset=None):
    if subset:
        pred = {k: v for k, v in pred.items() if v["type"] in subset}
        ref = {k: v for k, v in ref.items() if v["type"] in subset}
    if not ref or not pred:
        return dict(precision=0.0, recall=0.0, f1=0.0, n_pred=len(pred), n_ref=len(ref), tp=0)
    tp = len(set(pred) & set(ref))
    p, r = tp / len(pred), tp / len(ref)
    return dict(precision=p, recall=r, f1=2 * p * r / (p + r) if (p + r) else 0.0,
                n_pred=len(pred), n_ref=len(ref), tp=tp)


def write_tsv(path, rows):
    if not rows:
        print(f"  (no rows for {path.name})")
        return
    hdr = list(rows[0].keys())
    with open(path, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k])
                               for k in hdr) + "\n")
    print(f"  wrote {path.name}  ({len(rows)} rows)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tx2gene",
                    default=str(pathlib.Path(__file__).resolve().parents[1] / "data" / "tx2biotype.tsv"))
    a = ap.parse_args()
    res, out = pathlib.Path(a.results), pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    overall, byclass = [], []
    for model, ds in itertools.product(MODELS, DATASETS):
        d = res / f"{model}_{ds}"
        if not (d / "pred_profiles.npz").exists() or not (d / "real_collapsed.txt").exists():
            continue
        lc, _, _ = build_loader(d / "pred_profiles.npz", tx2gene=a.tx2gene)
        real = lc(d / "real_collapsed.txt", is_pred=False)
        for arm in ("pred_obsdepth", "pred_preddepth"):
            f = d / f"{arm}_collapsed.txt"
            if not f.exists():
                continue
            pc = lc(f, is_pred=True)
            s = prf(pc, real)
            overall.append(dict(dataset=ds, model=model, arm=arm, **s))
            for c in CLASSES:
                sc = prf(pc, real, subset=(c,))
                byclass.append(dict(dataset=ds, model=model, arm=arm, orf_class=c, **sc))
        print(f"  {model}/{ds}: real={len(real):,} calls")

    write_tsv(out / "human_overall.tsv", overall)
    write_tsv(out / "human_by_class.tsv", byclass)
    return 0


if __name__ == "__main__":
    sys.exit(main())
