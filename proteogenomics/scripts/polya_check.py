#!/usr/bin/env python3
"""Is a salmon quant poly(A)-selected or Ribo-Zero total RNA? (project rule
`feedback_rnaseq_must_be_polya`: SRA `library_selection=cDNA` does NOT prove poly(A) selection.)

TWO EARLIER VERSIONS OF THIS CHECK WERE WRONG, both silently. Recorded because the failure shape
matters more than the fix.

  v1 looked `quant.sf`'s `Name` column up in `tx2biotype.tsv` and summed the TPM in
     rRNA/Mt/misc/sn(o)RNA. Those key spaces do not match -- salmon carries the full pipe-delimited
     GENCODE header (`ENST00000641515.2|ENSG00000186092.7|...`), the table is keyed on the bare
     versioned id -- so every lookup missed, the junk total stayed 0.0, and it printed "0.00%".
     A total-RNA library would have printed the same thing.

  v2 fixed the key split and added a "did enough TPM get a biotype" guard. It reported 97.8% of TPM
     assigned and still 0.00% junk. Also meaningless: **the decoy-aware v49/vM38 salmon indexes are
     built from pc + lncRNA transcripts ONLY.** There is no rRNA, snRNA, snoRNA, misc_RNA, Mt_rRNA or
     Mt_tRNA in the index, so no library can ever put TPM there. The biotype approach is structurally
     incapable of detecting total RNA through this index and would stamp POLY(A)-SELECTED on a
     Ribo-Zero library.

WHAT ACTUALLY DISCRIMINATES, calibrated on this project's own known-good and known-bad libraries
(GSE243134 mouse liver totalRNA is the known-bad; Wang/Janich/THP-1 are known-good):

  1. salmon mapping rate against the pc+lncRNA index. Total RNA is mostly rRNA and intronic
     pre-mRNA, neither of which is in the index, so the rate collapses.
     known-bad 9.75-23.64% (median 17.73)   known-good 81-93%
  2. universe size at TPM >= 1.
     known-bad 3,321-3,534                  known-good 10,431-41,701
  3. TPM CONCENTRATION in a handful of structural transcripts. This is the sharpest signal and the
     only one that survives the index restriction, because the structural contaminants that DO get
     into a pc+lncRNA index are typed `lncRNA`: GSE243134 put 93.6% of TPM into 10 transcripts
     (`n-TKctt14`, tRNA-derived, alone = 82.8%; `Gm59647` = a 7SL/SRP duplicate). A poly(A) library
     is nowhere near that concentrated.

Usage: polya_check.py --quant <quant.sf> [--meta <meta_info.json>] [--top-n 10]
Exit 0 = poly(A), 1 = reject, 2 = could not evaluate.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Thresholds bracket the observed gap between the known-bad and known-good libraries; they are not
# tuned, and the gap they sit in is roughly 4x wide on every axis.
MAP_FAIL, MAP_WARN = 40.0, 65.0          # salmon % mapped
UNIV_FAIL, UNIV_WARN = 6000, 9000        # transcripts at TPM >= 1
# SINGLE-transcript share is the sharpest axis and the thresholds are set from the observed spread:
# known-bad n-TKctt14 79.8%, Wang liver Alb 11.1%, THP-1 MT-CO3 1.0%. Top-10 share is kept as a
# secondary signal but must stay LOOSE -- a healthy tissue with one dominant transcript reaches 34.6%
# (Wang liver), so the 25% warn this started with flagged a known-GOOD library.
TOP1_FAIL, TOP1_WARN = 50.0, 30.0
CONC_FAIL, CONC_WARN = 75.0, 60.0        # % of TPM in the top-N transcripts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quant", required=True)
    ap.add_argument("--meta", default=None,
                    help="salmon aux_info/meta_info.json (default: alongside --quant)")
    ap.add_argument("--top-n", type=int, default=10)
    a = ap.parse_args()

    q = Path(a.quant)
    rows = []
    with open(q) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ni, ti = hdr.index("Name"), hdr.index("TPM")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            rows.append((f[ni], float(f[ti])))
    tot = sum(t for _, t in rows)
    if tot <= 0:
        print(f"ABORT: total TPM is {tot} in {q}", file=sys.stderr)
        return 2

    n_expr = sum(1 for _, t in rows if t >= 1.0)
    top = sorted(rows, key=lambda r: -r[1])[: a.top_n]
    conc = 100.0 * sum(t for _, t in top) / tot

    meta = Path(a.meta) if a.meta else q.parent / "aux_info" / "meta_info.json"
    mapped = None
    if meta.exists():
        try:
            mapped = float(json.loads(meta.read_text())["percent_mapped"])
        except (KeyError, ValueError, json.JSONDecodeError):
            mapped = None

    print(f"  salmon mapped        : {f'{mapped:.2f}%' if mapped is not None else 'UNKNOWN (no meta_info.json)'}")
    print(f"  transcripts TPM >= 1 : {n_expr:,}")
    top1 = 100.0 * top[0][1] / tot if top else 0.0
    print(f"  top-{a.top_n} TPM share    : {conc:.2f}%")
    print(f"  single-tx max share  : {top1:.2f}%")
    print(f"  top 3: " + ", ".join(f"{n.split('|')[5] if '|' in n else n}={t / tot * 100:.1f}%"
                                   for n, t in top[:3]))

    fails, warns = [], []
    if mapped is not None:
        (fails if mapped < MAP_FAIL else warns if mapped < MAP_WARN else []).append(
            f"mapping rate {mapped:.1f}%")
    if n_expr < UNIV_FAIL:
        fails.append(f"universe {n_expr:,} tx")
    elif n_expr < UNIV_WARN:
        warns.append(f"universe {n_expr:,} tx")
    if top1 > TOP1_FAIL:
        fails.append(f"single transcript carries {top1:.1f}% of TPM ({top[0][0].split('|')[5] if '|' in top[0][0] else top[0][0]})")
    elif top1 > TOP1_WARN:
        warns.append(f"single transcript carries {top1:.1f}% of TPM")
    if conc > CONC_FAIL:
        fails.append(f"top-{a.top_n} concentration {conc:.1f}%")
    elif conc > CONC_WARN:
        warns.append(f"top-{a.top_n} concentration {conc:.1f}%")

    if fails:
        print(f"  VERDICT: REJECT -- looks like total RNA ({'; '.join(fails)})")
        return 1
    if warns:
        print(f"  VERDICT: AMBIGUOUS -- inspect before use ({'; '.join(warns)})")
        return 1
    if mapped is None:
        print("  VERDICT: POLY(A)-SELECTED on universe + concentration; mapping rate UNVERIFIED")
        return 0
    print("  VERDICT: POLY(A)-SELECTED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
