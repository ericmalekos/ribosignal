#!/usr/bin/env python3
"""EXPERIMENT (a), CORRECTED: what fraction of each TRAINING transcript's Ribo-seq did mm1 discard?

Uses bamlib.multimap_cost on TOTAL records per transcript (mm1 vs mm25), so isoform expansion
cancels. The earlier NH-split version measured isoform multiplicity and returned a nonsensical
99.8%; see CLAUDE.md "BAM / ALIGNMENT TRAPS" trap 1.

Reports a MANDATORY cross-check against STAR's own per-library multimap rates. If the aggregate
cost is not in the same neighbourhood, the measurement is wrong and the script says so.
"""
from __future__ import annotations
import csv, gzip, json, re, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

NEW = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(NEW / "scripts"))
from bamlib import load_counts_tsv, multimap_cost      # noqa: E402

D = NEW / "data/mm25_diagnostic"


def main():
    mm1, mm25 = defaultdict(int), defaultdict(int)
    for f in sorted((D / "chothani_mm1_counts").glob("*.mm1.tsv.gz")):
        for t, v in load_counts_tsv(f, (1,)).items(): mm1[t] += v
    # nh TSV is tx \t n_unique \t n_multi -- the TOTAL is what is sound (cols 1+2)
    for f in sorted((D / "chothani_nh_counts").glob("*.nh.tsv.gz")):
        for t, v in load_counts_tsv(f, (1, 2)).items(): mm25[t] += v
    print(f"  mm1 tx {len(mm1):,}   mm25 tx {len(mm25):,}")
    tot1, tot25 = sum(mm1.values()), sum(mm25.values())
    agg = 1 - tot1 / tot25
    print(f"  total records: mm1 {tot1:,}   mm25 {tot25:,}   aggregate cost {100*agg:.1f}%")

    # ---- MANDATORY CROSS-CHECK against STAR's own logs (CLAUDE.md rule 5) --------------------
    rates = []
    for p in (NEW / "logs/mm25").glob("cho_*.out"):
        for line in p.read_text().splitlines():
            m = re.search(r"unique=([0-9.]+)% multi=([0-9.]+)%", line)
            if m:
                u, mm = float(m.group(1)), float(m.group(2))
                rates.append(mm / (u + mm))
    star = float(np.mean(rates)) if rates else float("nan")
    print(f"  STAR logs imply multi/(unique+multi) = {100*star:.1f}%  (n={len(rates)} libraries)")
    if not (0.5 * star <= agg <= 2.0 * star):
        print(f"\n  *** CROSS-CHECK FAILED: {100*agg:.1f}% vs STAR's {100*star:.1f}%. "
              f"Measurement is WRONG -- do not report. ***")
        return 1
    print(f"  cross-check PASSES (within 2x of STAR's own rate)\n")

    cost = multimap_cost(mm1, mm25)
    order = (NEW / "data/packed_union/tx_order.txt").read_text().split()
    trained = [t for t in order if mm1.get(t, 0) >= 50]
    a = np.array([cost[t] for t in trained if t in cost])
    print(f"  universe {len(order):,}   trainable (>=50 mm1 records) {len(trained):,}")
    print(f"\n  FRACTION OF SIGNAL mm1 DISCARDED, per trainable transcript")
    for q in (10, 25, 50, 75, 90, 99):
        print(f"    p{q:<3} {np.percentile(a,q):.3f}")
    for thr in (0.2, 0.3, 0.5, 0.8):
        n = int((a > thr).sum())
        print(f"    >{int(thr*100)}% discarded: {n:,} tx ({100*n/len(a):.1f}%)")

    bt = {}
    with open(NEW / "data/tx2biotype.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"): bt[r["tx_id"]] = r.get("transcript_type", "?")
    print(f"\n  {'biotype':<18} {'n':>6} {'median':>8} {'>30% lost':>11} {'>50% lost':>11}")
    byb = defaultdict(list)
    for t in trained:
        if t in cost: byb[bt.get(t, "?")].append(cost[t])
    for b, v in sorted(byb.items(), key=lambda kv: -len(kv[1]))[:5]:
        v = np.array(v)
        print(f"  {b:<18} {len(v):>6} {np.median(v):>8.3f} {100*(v>0.3).mean():>10.1f}% "
              f"{100*(v>0.5).mean():>10.1f}%")

    json.dump({"model": None,
               "source": "mm1 vs mm25 total records/transcript over 74 Chothani training libraries",
               "provenance_note": "bamlib.multimap_cost; cross-checked against STAR per-library "
                                  "multimap rate. Supersedes the NH-split version (trap 1).",
               "aggregate_cost": round(agg, 4), "star_log_rate": round(star, 4),
               "n_trainable": len(a), "median_cost": round(float(np.median(a)), 4),
               "frac_over_30pct": round(float((a > 0.3).mean()), 4),
               "frac_over_50pct": round(float((a > 0.5).mean()), 4)},
              open(D / "training_multimap_cost.json", "w"), indent=2)
    with gzip.open(D / "training_multimap_cost.tsv.gz", "wt") as fh:
        fh.write("tx_id\tmm1_records\tmm25_records\tfrac_discarded\n")
        for t in sorted(trained, key=lambda x: -cost.get(x, 0)):
            if t in cost: fh.write(f"{t}\t{mm1.get(t,0)}\t{mm25.get(t,0)}\t{cost[t]:.4f}\n")
    print(f"\n  wrote training_multimap_cost.{{json,tsv.gz}}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
