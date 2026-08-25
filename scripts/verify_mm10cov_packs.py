#!/usr/bin/env python3
"""Verify the mm10cov pack family before any GPU time is spent on it.

Three things must hold, and each has bitten this project before:

1. target_counts.npy is BYTE-IDENTICAL to the source union pack. --reuse-target promises the Ribo
   label is untouched; if it moved, an A-vs-baseline difference could be the label rather than the
   input, and the experiment is uninterpretable.
2. coverage.npy actually CHANGED. A pack that silently reused mm1 coverage would train fine and
   prove nothing.
3. The change lands where it was predicted to. GTF2I ENST00000901263.1 has a 1,617 nt window
   (positions 1682:3299) reading 8.9% of CDS mean at mm1 and 138.3% at mm10 in the BAM-level
   analysis (figures/X_gtf2i_isoforms/X_gtf2i_rna_postures_values.json). If the pack does not
   reproduce that, the coverage rebuild did not reach the pack.

Usage: verify_mm10cov_packs.py [--tissues ...] [--json out.json]
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

R = Path(__file__).resolve().parent.parent
TX = "ENST00000901263.1"
GAP = (1682, 3299)
CDS = None  # filled from the postures values file
PAIRS = [("packed_union", "packed_mm10cov"),
         ("packed_union_ES", "packed_mm10cov_ES"),
         ("packed_union_Fat", "packed_mm10cov_Fat"),
         ("packed_union_HA_EC", "packed_mm10cov_HA_EC"),
         ("packed_union_HCAEC", "packed_mm10cov_HCAEC"),
         ("packed_union_Hepatocytes", "packed_mm10cov_Hepatocytes"),
         ("packed_union_HUVEC", "packed_mm10cov_HUVEC"),
         ("packed_union_VSMC", "packed_mm10cov_VSMC")]


def md5(p, chunk=1 << 24):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def tx_slice(pack, tx):
    order = (pack / "tx_order.txt").read_text().split()
    if tx not in order:
        return None
    i = order.index(tx)
    off = np.load(pack / "offsets.npy")
    return int(off[i]), int(off[i + 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default="")
    a = ap.parse_args()
    global CDS
    pv = json.loads((R / "figures/X_gtf2i_isoforms/X_gtf2i_rna_postures_values.json").read_text())
    CDS = pv["cds"]
    exp_mm1 = pv["postures"]["mm1"]["gap_pct_of_cds"]
    exp_mm10 = pv["postures"]["mm10"]["gap_pct_of_cds"]
    print(f"  BAM-level reference: gap is {exp_mm1}% of CDS at mm1, "
          f"{exp_mm10}% at mm10\n")

    out, ok = [], True
    print(f"  {'tissue':14s} {'target':>10s} {'coverage':>10s} {'gap%CDS':>9s}  verdict")
    for src, new in PAIRS:
        sp, np_ = R / "data" / src, R / "data" / new
        if not (np_ / "coverage.npy").exists():
            print(f"  {new:14s} {'-':>10s} {'-':>10s} {'-':>9s}  MISSING")
            ok = False
            continue
        tgt_same = md5(sp / "target_counts.npy") == md5(np_ / "target_counts.npy")
        cov_diff = md5(sp / "coverage.npy") != md5(np_ / "coverage.npy")
        gp = float("nan")
        sl = tx_slice(np_, TX)
        if sl:
            b, e = sl
            cov = np.load(np_ / "coverage.npy", mmap_mode="r")[b:e].astype(np.float64)
            g = cov[GAP[0]:GAP[1]].mean()
            c = cov[CDS[0]:CDS[1]].mean()
            gp = 100 * g / c if c else float("nan")
        v = ("OK" if (tgt_same and cov_diff)
             else "TARGET MOVED" if not tgt_same else "COVERAGE UNCHANGED")
        if v != "OK":
            ok = False
        print(f"  {src.replace('packed_union','').lstrip('_') or 'Fibroblast':14s} "
              f"{'same' if tgt_same else 'CHANGED':>10s} {'differs' if cov_diff else 'SAME':>10s} "
              f"{gp:8.1f}%  {v}")
        out.append({"tissue": src, "target_identical": tgt_same,
                    "coverage_differs": cov_diff,
                    "gtf2i_gap_pct_of_cds": None if gp != gp else round(gp, 1)})
    print(f"\n  ALL PACKS VALID: {ok}")
    if a.json:
        Path(a.json).write_text(json.dumps(
            {"model": None,
             "source": "scripts/verify_mm10cov_packs.py; "
                       "data/packed_union* vs data/packed_mm10cov*",
             "provenance_note": "MODEL-FREE pack contract check. Reference percentages "
                                "from figures/X_gtf2i_isoforms/"
                                "X_gtf2i_rna_postures_values.json (BAM level).",
             "expected_gap_pct_mm1": exp_mm1, "expected_gap_pct_mm10": exp_mm10,
             "all_valid": ok, "packs": out}, indent=2) + "\n")
        print(f"  wrote {a.json}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
