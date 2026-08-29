#!/usr/bin/env python3
"""Break the cross-species ORF calls into the three classes of biological interest:
CDS, uORF and ncORF, and dump the calls themselves so they can be looked at.

CLASS DEFINITIONS, and why the RiboCode ORF_type alone is not enough:

  CDS    ORF_type == "annotated" on a protein_coding transcript. The annotated proteome.
  uORF   ORF_type in {uORF, Overlap_uORF} on a protein_coding transcript. An ORF in the
         5'UTR; Overlap_uORF starts in the 5'UTR and runs into the CDS out of frame.
  ncORF  ANY called ORF on a NON-coding transcript (lncRNA, antisense_RNA,
         ncRNA_pseudogene, ...). This is a biotype question, not an ORF_type question:
         RiboCode types these "novel" because there is no annotated CDS to compare to,
         and "novel" ALSO contains novel ORFs on protein_coding transcripts, which are a
         different claim entirely. Splitting on ORF_type alone conflates the two.

  (dORF, internal, and novel-on-coding are reported as `other` for completeness; they are
   not part of what was asked for but dropping them silently would make the class counts
   fail to reconcile against the totals in results/xspecies_orf_calls.tsv.)

MEMBERSHIP IS DEFINED BY THE SHARED LOADER, NOT BY THIS FILE. build_loader applies the
test-tx restriction, pval <= 0.05, ORF >= 90 nt and the predicted-side 0.5x enrichment
gate, keyed on (gene_id, ORF_gstop). This script re-reads the collapsed file ONLY to
attach the biotype and coordinates to keys the loader already accepted, so the classes
sum to the totals already reported. Hand-loading these files instead of using the loader
once produced a 4x overcount and flipped a conclusion.

Usage: python3 scripts/xspecies/orf_classes.py [--arch attn] [--dump-dir results/xspecies_orf_classes]
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
SPECIES = {
    "xsp_yeast": "yeast", "xsp_celegans": "celegans", "xsp_zebrafish": "zebrafish",
    "xsp_gorilla": "gorilla", "xsp_chimp": "chimp", "xsp_macaque": "macaque",
    "xsp_human": "human_refseq",
}
CODING = "protein_coding"
CLASSES = ("CDS", "uORF", "ncORF", "other")


def classify(orf_type: str, tx_type: str) -> str:
    if tx_type != CODING:
        return "ncORF"
    if orf_type == "annotated":
        return "CDS"
    if orf_type in ("uORF", "Overlap_uORF"):
        return "uORF"
    return "other"


def annotate(path: Path, keys: set) -> dict:
    """Attach biotype + coordinates to the keys the loader accepted."""
    out = {}
    with path.open() as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            k = (r["gene_id"], r["ORF_gstop"])
            if k not in keys or k in out:
                continue
            r["_class"] = classify(r["ORF_type"], r["transcript_type"])
            out[k] = r
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", default="attn")
    ap.add_argument("--variant", default="pred_preddepth")
    ap.add_argument("--dump-dir", type=Path,
                    default=NEW / "results" / "xspecies_orf_classes")
    args = ap.parse_args()
    args.dump_dir.mkdir(parents=True, exist_ok=True)

    rows, dumped = [], 0
    for ds, sp in SPECIES.items():
        dd = NEW / "results" / "xspecies_dropin" / f"{ds}_{args.arch}"
        prof = (NEW / "results" / "loto" / RUN.format(a=args.arch)
                / f"dropin_{ds}" / "pred_profiles.npz")
        rp, pp = dd / "real_collapsed.txt", dd / f"{args.variant}_collapsed.txt"
        if not (prof.exists() and rp.exists() and pp.exists()):
            print(f"  skip {ds} (incomplete)", file=sys.stderr)
            continue
        lc, _, _ = build_loader(prof, key="genomic",
                                tx2gene=NEW / "data" / "annot" / sp / "tx2biotype.tsv")
        real, pred = lc(rp), lc(pp, is_pred=True)
        ra, pa = annotate(rp, set(real)), annotate(pp, set(pred))

        # dump the calls themselves, one file per species per class
        for cls in CLASSES:
            sub = [r for r in ra.values() if r["_class"] == cls]
            if not sub:
                continue
            sub.sort(key=lambda r: -int(r["ORF_length"]))
            f = args.dump_dir / f"{ds.replace('xsp_','')}_{cls}_real.tsv"
            cols = ["ORF_type", "transcript_id", "transcript_type", "gene_id", "gene_name",
                    "chrom", "strand", "ORF_length", "ORF_gstart", "ORF_gstop",
                    "Psites_sum_frame0", "pval_combined", "recovered_by_model"]
            with f.open("w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                                   extrasaction="ignore", lineterminator="\n")
                w.writeheader()
                for r in sub:
                    r["recovered_by_model"] = "yes" if (r["gene_id"], r["ORF_gstop"]) in pred else "no"
                    w.writerow(r)
            dumped += 1

        for cls in CLASSES:
            R = {k for k, r in ra.items() if r["_class"] == cls}
            P = {k for k, r in pa.items() if r["_class"] == cls}
            tp = len(R & P)
            rows.append({
                "pack": ds.replace("xsp_", ""), "arch": args.arch, "variant": args.variant,
                "class": cls, "n_real": len(R), "n_pred": len(P), "tp": tp,
                "precision": round(tp / len(P), 3) if P else "",
                "recall": round(tp / len(R), 3) if R else "",
                "median_len_nt": (sorted(int(ra[k]["ORF_length"]) for k in R)[len(R) // 2]
                                  if R else ""),
            })

    if not rows:
        print("nothing to classify", file=sys.stderr)
        return 1
    tsv = NEW / "results" / f"xspecies_orf_classes_{args.arch}.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    for cls in CLASSES:
        print(f"\n  === {cls} ===")
        print(f"  {'pack':<11}{'real':>9}{'pred':>9}{'TP':>9}{'prec':>8}{'rec':>8}{'medLen':>9}")
        print("  " + "-" * 61)
        for r in rows:
            if r["class"] != cls:
                continue
            print(f"  {r['pack']:<11}{r['n_real']:>9,}{r['n_pred']:>9,}{r['tp']:>9,}"
                  f"{str(r['precision']):>8}{str(r['recall']):>8}{str(r['median_len_nt']):>9}")
    print(f"\n  -> {tsv}\n  -> {dumped} per-class call files in {args.dump_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
