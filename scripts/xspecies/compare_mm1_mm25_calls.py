#!/usr/bin/env python3
"""Compare mm1 vs mm25 RiboCode ORF calls per cross-species arm.

WHAT THIS MEASURES AND WHY IT IS A DIAGNOSTIC, NOT A RESULT.
mm25 keeps reads aligning to up to 25 loci. RiboCode does not drop multi-mappers itself,
so those reads add P-site density that the frame test then reads as periodicity. The
question this answers is how much that inflates the call set, per genome.

The mouse-liver pair already on disk (data/mm25_diagnostic/psites_combined_mm{1,25}/)
shows the expected signature: total calls 23,476 -> 31,120 (1.33x), but annotated rises
only 5% while `novel` goes 4,355 -> 10,993 (2.5x). Extra calls concentrate in the
unannotated classes, which is what multimapper-manufactured periodicity looks like.

KEYING is (gene_id, ORF_gstop), the genomic ORF locus, NOT ORF_ID: one genomic ORF
appears once per transcript it sits on, and the collapse representative shifts when the
input changes. Filters match the project default: pval_combined <= 0.05, ORF >= 90 nt.

Usage: python3 scripts/xspecies/compare_mm1_mm25_calls.py
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
MM1 = NEW / "data" / "xspecies_psites"
MM25 = NEW / "data" / "xspecies_psites_mm25"
ARMS = ["yeast_gse173654_ribo", "worm_gse52905_ribo", "zf_gse46512_ribo",
        "primate_gg_ribo", "primate_pt_ribo", "primate_rm_ribo",
        "ruizorera_hsCM_ribo", "fly_gse99920_ribo"]
CLASSES = ("annotated", "uORF", "Overlap_uORF", "dORF", "Overlap_dORF", "internal", "novel")


def runs_used(pre_config: Path) -> int:
    """How many sequencing runs actually contributed to the pooled call.

    RiboCode's `<ds>_pre_config.txt` carries one UNCOMMENTED row per (run, read length)
    that passed the 3-nt periodicity test: fields are sample, bam, keep-flag, read_length,
    psite_offset. A run whose periodicity fails contributes no row at all and is silently
    absent from the pooled call.

    This is the confound that makes a raw call-count ratio unreadable. Chimp passed 8/8
    runs at mm1 but only 3/8 at mm25, so its 0.884x "fewer calls" is a DEPTH effect from
    losing five libraries, not evidence that mm25 calls fewer ORFs. Any arm where these
    counts differ is not a like-for-like comparison and is flagged accordingly.
    """
    if not pre_config.exists():
        return 0
    runs = set()
    for ln in pre_config.read_text().splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split("\t")
        if len(f) >= 2:
            runs.add(f[0])
    return len(runs)


def load(p: Path, max_pval=0.05, min_len=90) -> dict:
    d = {}
    if not p.exists():
        return d
    with p.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                if float(r["pval_combined"]) > max_pval:
                    continue
                if int(r["ORF_length"]) < min_len:
                    continue
            except (ValueError, KeyError):
                continue
            d.setdefault((r["gene_id"], r["ORF_gstop"]), r)
    return d


def main() -> int:
    rows = []
    for arm in ARMS:
        a = load(MM1 / arm / f"{arm}_collapsed.txt")
        b = load(MM25 / arm / f"{arm}_collapsed.txt")
        if not a or not b:
            print(f"  skip {arm} (mm1={len(a):,} mm25={len(b):,})", file=sys.stderr)
            continue
        ca, cb = Counter(r["ORF_type"] for r in a.values()), Counter(r["ORF_type"] for r in b.values())
        ra = runs_used(MM1 / arm / f"{arm}_pre_config.txt")
        rb = runs_used(MM25 / arm / f"{arm}_pre_config.txt")
        rows.append({
            "runs_mm1": ra, "runs_mm25": rb,
            "comparable": "yes" if ra == rb else "NO_depth_differs",
            "arm": arm, "mm1": len(a), "mm25": len(b), "ratio": round(len(b) / len(a), 3),
            "shared": len(set(a) & set(b)), "mm1_only": len(set(a) - set(b)),
            "mm25_only": len(set(b) - set(a)),
            "jaccard": round(len(set(a) & set(b)) / len(set(a) | set(b)), 3),
            **{f"mm1_{c}": ca.get(c, 0) for c in CLASSES},
            **{f"mm25_{c}": cb.get(c, 0) for c in CLASSES},
            "annotated_ratio": round(cb.get("annotated", 0) / ca["annotated"], 3) if ca.get("annotated") else "",
            "novel_ratio": round(cb.get("novel", 0) / ca["novel"], 3) if ca.get("novel") else "",
        })
    if not rows:
        print("no arm has both call sets yet", file=sys.stderr)
        return 1
    h = (f"  {'arm':<22}{'mm1':>8}{'mm25':>8}{'ratio':>7}{'jacc':>7}"
         f"{'mm25_only':>11}{'runs mm1/mm25':>15}{'comparable':>18}")
    print("\n" + h); print("  " + "-" * (len(h) - 2))
    for r in rows:
        print(f"  {r['arm']:<22}{r['mm1']:>8,}{r['mm25']:>8,}{r['ratio']:>7.3f}"
              f"{r['jaccard']:>7.3f}"
              f"{r['mm25_only']:>11,}"
              f"{(str(r['runs_mm1'])+'/'+str(r['runs_mm25'])):>15}"
              f"{r['comparable']:>18}")
    out = NEW / "results" / "xspecies_mm1_vs_mm25_calls.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"\n  -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
