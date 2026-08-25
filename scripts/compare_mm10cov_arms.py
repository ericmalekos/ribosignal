#!/usr/bin/env python3
"""Headline test-metric comparison: baseline (3 seeds) vs arm A vs arm B.

Every difference is judged against the seed-noise floor, not against zero.

Test split is identical across all three arms (n=70,883), which is what makes them comparable.
Val is NOT comparable between A and B -- arm B's val has the multimap-corrupted tx removed, so it is
a strictly easier val set and B's val number is inflated by construction. This script deliberately
reports TEST ONLY.

The floor is measured from the three baseline seeds: pearson_median range 0.0098,
frame0_pred_median range 0.0013 (the sensitive discriminator, an order of magnitude
tighter), period_pred_median range 0.0345. A difference smaller than the floor is
not a result.
"""
import argparse
import json
import statistics as st
import sys
from pathlib import Path

R = Path(__file__).resolve().parent.parent
LOTO = R / "results/loto"
BASE = "orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes"
MM10 = "orf_v2_mamba4_onehot_union_noBrain_nokozak_mm10cov"
METRICS = ["pearson_median", "frame0_pred_median", "period_pred_median", "spearman_median"]


def load(run):
    p = LOTO / run / "test_metrics.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    flat = {}

    def fl(o, pre=""):
        for k, v in o.items():
            if isinstance(v, dict):
                fl(v, pre + k + ".")
            elif isinstance(v, int | float):
                flat[pre + k] = v
    fl(d)
    return {m: flat.get(f"all.{m}") for m in METRICS} | {"n": flat.get("all.n")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default="")
    a = ap.parse_args()
    seeds = [load(BASE), load(f"{BASE}_seed1"), load(f"{BASE}_seed2")]
    seeds = [s for s in seeds if s]
    armA = load(f"{MM10}_holdout_Hepatocytes")
    armB = load(f"{MM10}_exclmm50_holdout_Hepatocytes")
    if not seeds:
        sys.exit("no baseline seeds found")

    floor = {m: (max(v) - min(v)) for m in METRICS
             if len(v := [s[m] for s in seeds if s.get(m) is not None]) > 1}
    bmean = {m: st.mean([s[m] for s in seeds if s.get(m) is not None]) for m in METRICS}

    allns = sorted({s["n"] for s in seeds} | {x["n"] for x in (armA, armB) if x})
    print(f"  baseline seeds: {len(seeds)}   test n: {allns}")
    print("  (a single n across all rows is required; differing n means")
    print("   the arms are not comparable)\n")

    print(f"  {'metric':22s} {'base mean':>10s} {'armA':>9s} {'armB':>9s} "
          f"{'A-base':>9s} {'B-base':>9s} {'B-A':>9s} {'floor':>8s} verdict")
    out = {}
    for m in METRICS:
        b = bmean.get(m)
        A = armA.get(m) if armA else None
        B = armB.get(m) if armB else None
        f = floor.get(m, float("nan"))
        def fmt(x):
            return f"{x:9.4f}" if isinstance(x, float) else f"{'--':>9s}"
        dA = (A - b) if (A is not None and b is not None) else None
        dB = (B - b) if (B is not None and b is not None) else None
        dBA = (B - A) if (A is not None and B is not None) else None
        verdict = "pending"
        if dBA is not None and f == f:
            sig = [n for n, d in (("A", dA), ("B", dB), ("B-A", dBA))
                   if d is not None and abs(d) > f]
            verdict = ("exceeds floor: " + ",".join(sig)) if sig else "all within seed noise"
        print(f"  {m:22s} {fmt(b)[1:]:>10s} {fmt(A)} {fmt(B)} {fmt(dA)} {fmt(dB)} {fmt(dBA)} "
              f"{f:8.4f} {verdict}")
        out[m] = {"baseline_mean": b, "armA": A, "armB": B, "d_A": dA, "d_B": dB,
                  "d_BA": dBA, "seed_floor": f if f == f else None, "verdict": verdict}
    if not armA or not armB:
        print("\n  NOTE: one or both arms have no test_metrics.json yet (training unfinished).")
    if a.json:
        Path(a.json).write_text(json.dumps(
            {"model": ["mamba4"], "source": "scripts/compare_mm10cov_arms.py over "
             "results/loto/*/test_metrics.json",
             "note": "TEST ONLY. Val is not comparable A-vs-B (B's val excludes the corrupted tx). "
                     "Differences judged against the 3-seed baseline noise floor.",
             "metrics": out}, indent=2) + "\n")
        print(f"  wrote {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
