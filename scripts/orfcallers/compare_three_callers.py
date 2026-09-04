#!/usr/bin/env python3
"""Three-way ORF-caller comparison on chr1: RiboCode, Ribo-TISH, RiboTaper.

All three score the SAME predicted profiles (build_density is imported by both drop-ins, and the
RiboTaper bams round-trip to the same P-sites), so any difference is the caller.

TWO THINGS THAT HAD TO BE FIXED FIRST, both caught by a cross-check that returned zero overlap:

1. COORDINATE CONVENTIONS DIFFER BY CALLER AND BY STRAND. Keyed on the ORF's genomic START
   (the AUG), the offsets against RiboCode's ORF_gstart are, calibrated empirically here rather
   than assumed:
       Ribo-TISH   + strand +1,  - strand  0
       RiboTaper   + strand  0,  - strand +1
   Keying on the STOP coordinate does not work at all: the three disagree by -1 and +2 there,
   because they differ over whether the stop codon is inside the ORF.

2. CLASS VOCABULARIES DIFFER. They are mapped onto one set below. Classes with no counterpart in
   a caller are reported as absent, never as zero.

The offsets are re-derived on every run and printed, so a change in any caller's output format
shows up as a shifted offset instead of silently as lost overlap.
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path

NEW = Path(__file__).resolve().parents[2]
RUN = NEW / "results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes"
GP = re.compile(r"^(\S+):(\d+)-(\d+):([+-])$")
ARMS = ["real", "pred_obsdepth", "pred_preddepth"]

# one vocabulary. None means the caller has no counterpart for that class.
CANON = {
    "RiboCode":  {"annotated": "canonical", "uORF": "uORF", "dORF": "dORF", "novel": "novel",
                  "internal": "internal", "Overlap_uORF": "Overlap_uORF",
                  "Overlap_dORF": "Overlap_dORF"},
    "Ribo-TISH": {"Annotated": "canonical", "Truncated": "canonical", "Extended": "canonical",
                  "Truncated:Known": "canonical", "Extended:Known": "canonical",
                  "Extended:CDSFrameOverlap": "canonical",
                  "5'UTR": "uORF", "5'UTR:Known": "uORF",
                  "5'UTR:CDSFrameOverlap": "Overlap_uORF",
                  "3'UTR": "dORF", "3'UTR:Known": "dORF",
                  "3'UTR:CDSFrameOverlap": "Overlap_dORF",
                  "Internal": "internal", "Internal:Known": "internal",
                  "Internal:CDSFrameOverlap": "internal",
                  "Novel": "novel", "Novel:Known": "novel", "Novel:CDSFrameOverlap": "novel"},
    # RiboTaper's uORF/dORF are STRICTLY non-overlapping: CCDS_orf_finder.R:989 overwrites the
    # uORF label with Overl_uORF whenever the ORF runs past annotated_start, and :988 assigns
    # Overl_dORF separately. Verified on the observed arm: 0 of 192 uORF calls have
    # stop_pos > annotated_start, 0 of 21 dORF calls have start_pos < annotated_stop.
    # Overl_uORF is mapped even though chr1 produced none, so a future call cannot be dropped
    # silently by the `cc is None: continue` branch.
    "RiboTaper": {"ORFs_ccds": "canonical", "uORF": "uORF", "dORF": "dORF",
                  "ncORFS": "novel", "Overl_uORF": "Overlap_uORF",
                  "Overl_dORF": "Overlap_dORF", "nonccds_coding_ORFs": "novel"},
}
CLASSES = ["canonical", "uORF", "Overlap_uORF", "dORF", "Overlap_dORF", "internal", "novel"]


def strip(x):
    return x.split(".")[0]


def load_ribocode(arm):
    out = []
    with (RUN / "dropin" / f"{arm}_collapsed.txt").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["chrom"] != "chr1":
                continue
            out.append((strip(r["gene_id"]), r["strand"], int(r["ORF_gstart"]), r["ORF_type"]))
    return out


def load_ribotish(arm):
    out = []
    with (RUN / "ribotish/longest" / f"{arm}.txt").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m = GP.match(r["GenomePos"])
            if not m or m.group(1) != "chr1":
                continue
            a, b, st = int(m.group(2)), int(m.group(3)), m.group(4)
            out.append((strip(r["Gid"]), st, a if st == "+" else b, r["TisType"]))
    return out


def load_ribotaper(arm):
    out = []
    with (RUN / "ribotaper" / arm / "ORFs_max_filt").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c, a, b = r["ORF_id_gen"].rsplit("_", 2)
            a, b, st = int(a), int(b), r["strand"]
            out.append((strip(r["gene_id"]), st, a if st == "+" else b, r["category"]))
    return out


LOAD = {"RiboCode": load_ribocode, "Ribo-TISH": load_ribotish, "RiboTaper": load_ribotaper}


def calibrate(rows, ref_by_strand):
    """Best per-strand start offset against RiboCode canonical, re-derived every run."""
    off = {}
    for st in ("+", "-"):
        ref = ref_by_strand[st]
        cand = [(g, p) for (g, s, p, c) in rows if s == st and CANON_OF(rows, c) == "canonical"]
        best = (0, -1)
        for o in range(-4, 5):
            n = len({(g, p + o) for g, p in cand} & ref)
            if n > best[1]:
                best = (o, n)
        off[st] = best[0]
    return off


def keyed(caller, rows, off):
    """(gene, strand, calibrated start) -> canonical class; drops classes with no mapping."""
    m = CANON[caller]
    d = {}
    for g, s, p, c in rows:
        cc = m.get(c)
        if cc is None:
            raise SystemExit(f"unmapped {caller} category {c!r}: add it to CANON rather than "
                             f"letting it be dropped silently")
        d[(g, s, p + off[s])] = cc
    return d


def CANON_OF(rows, c):
    for caller, m in CANON.items():
        if c in m:
            return m[c]
    return None


def prf(pred, real):
    cls = {}
    cls.update(pred)
    cls.update(real)
    tp, fp, fn = Counter(), Counter(), Counter()
    for k in set(real) | set(pred):
        c = cls[k]
        if k in real and k in pred:
            tp[c] += 1
        elif k in pred:
            fp[c] += 1
        else:
            fn[c] += 1
    out = {}
    for c in set(tp) | set(fp) | set(fn):
        P = tp[c] / max(1, tp[c] + fp[c])
        R = tp[c] / max(1, tp[c] + fn[c])
        out[c] = {"P": P, "R": R, "F1": 0.0 if P + R == 0 else 2 * P * R / (P + R),
                  "TP": tp[c], "FP": fp[c], "FN": fn[c], "n_real": tp[c] + fn[c]}
    T, Fp, Fn = sum(tp.values()), sum(fp.values()), sum(fn.values())
    P = T / max(1, T + Fp)
    R = T / max(1, T + Fn)
    out["ALL"] = {"P": P, "R": R, "F1": 0.0 if P + R == 0 else 2 * P * R / (P + R),
                  "TP": T, "FP": Fp, "FN": Fn, "n_real": T + Fn}
    return out


def main():
    raw = {c: {a: LOAD[c](a) for a in ARMS} for c in LOAD}
    rc_real = raw["RiboCode"]["real"]
    ref = {st: {(g, p) for (g, s, p, c) in rc_real if s == st and c == "annotated"}
           for st in ("+", "-")}
    offs = {"RiboCode": {"+": 0, "-": 0}}
    for c in ("Ribo-TISH", "RiboTaper"):
        offs[c] = calibrate(raw[c]["real"], ref)
    print("calibrated start-coordinate offsets vs RiboCode ORF_gstart:")
    for c, o in offs.items():
        print(f"   {c:<10} + {o['+']:+d}   - {o['-']:+d}")

    K = {c: {a: keyed(c, raw[c][a], offs[c]) for a in ARMS} for c in LOAD}

    print("\n=== chr1 call counts, common class vocabulary ===")
    print(f"{'class':<14}" + "".join(f"{c[:9]:>11}" * 3 for c in LOAD))
    print(f"{'':<14}" + "".join(f"{a[:9]:>11}" for c in LOAD for a in ARMS))
    for cl in CLASSES:
        row = f"{cl:<14}"
        for c in LOAD:
            for a in ARMS:
                n = sum(1 for v in K[c][a].values() if v == cl)
                row += f"{(n if n else '-'):>11}"
        print(row)

    print("\n=== predicted vs observed, per caller, per class ===")
    res = {}
    for c in LOAD:
        for a in ("pred_obsdepth", "pred_preddepth"):
            m = prf(K[c][a], K[c]["real"])
            res[f"{c}|{a}"] = m
            print(f"\n-- {c}, {a}")
            print(f"{'class':<14}{'P':>8}{'R':>8}{'F1':>8}{'TP':>8}{'FP':>8}{'FN':>8}{'n_real':>8}")
            for cl in CLASSES + ["ALL"]:
                if cl not in m:
                    continue
                v = m[cl]
                print(f"{cl:<14}{v['P']:>8.4f}{v['R']:>8.4f}{v['F1']:>8.4f}"
                      f"{v['TP']:>8,}{v['FP']:>8,}{v['FN']:>8,}{v['n_real']:>8,}")

    print("\n=== caller vs caller on the OBSERVED profile (chr1) ===")
    names = list(LOAD)
    overlap = {}
    for cl in CLASSES:
        sets = {c: {k for k, v in K[c]["real"].items() if v == cl} for c in names}
        if not any(sets.values()):
            continue
        print(f"\n-- {cl}: " + "  ".join(f"{c} {len(sets[c]):,}" for c in names))
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = sets[names[i]], sets[names[j]]
                if not a or not b:
                    continue
                jac = len(a & b) / len(a | b)
                overlap[f"{names[i]}|{names[j]}|{cl}"] = {"shared": len(a & b),
                                                          "jaccard": round(jac, 4)}
                print(f"     {names[i]} & {names[j]}: shared {len(a & b):,}  jaccard {jac:.3f}")
        if all(sets.values()):
            allx = set.intersection(*sets.values())
            overlap[f"ALL3|{cl}"] = len(allx)
            print(f"     all three: {len(allx):,}")

    out = NEW / "results/orfcaller_comparison_chr1.json"
    out.write_text(json.dumps({
        "model": "mamba4", "restriction": "chr1",
        "run": RUN.name, "offsets": offs,
        "source": "RiboCode dropin/*_collapsed.txt; Ribo-TISH ribotish/longest/*.txt; "
                  "RiboTaper ribotaper/*/ORFs_max_filt. Same density arrays in all three.",
        "class_map": CANON,
        "counts": {c: {a: dict(Counter(K[c][a].values())) for a in ARMS} for c in LOAD},
        "observed_overlap": overlap,
        "pred_vs_real": res}, indent=1) + "\n")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    raise SystemExit(main())
