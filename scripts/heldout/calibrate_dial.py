#!/usr/bin/env python3
"""O1 / Check 6 -- expose the CDS-anchored stringency dial.

Calls de novo ORFs from a predicted P-site profile at a user-chosen canonical-sensitivity level and
attaches a per-ORF confidence. Unlike the evaluation-side anchoring in orf_call_metrics.py, this
needs **no observed Ribo-seq**, so it is usable on a transcript set that has never been
ribosome-profiled (which is the point of the model).

  --cds-recall 0.90  ->  "call at whatever stringency retains 90% of the canonical ORFs I could
                          recover at all on this data"

Two recall denominators, and why the dial uses the second
--------------------------------------------------------
  absolute : called annotated-CDS tx / ALL CDS-bearing tx in the dump.
             Honest but NOT dial-able: it saturates well below 1 (0.35 on mouse liver) because most
             annotated transcripts are simply not expressed in the tissue, so no stringency setting
             can ever recover them. Reported for transparency.
  relative : called annotated-CDS tx / annotated-CDS tx called at the LOOSEST theta in the grid.
             This is what --cds-recall targets. It spans 0..1 by construction, is monotone in
             theta, and reads plainly: "how much canonical sensitivity am I giving up, relative
             to the un-thresholded caller, in exchange for precision?" Expression-limited
             transcripts drop out of numerator and denominator together, so the dial measures
             stringency alone.

How it works
  1. Sweep the Poisson count-scale theta. At each theta, run the RiboCode drop-in on the predicted
     profile (variant pred_preddepth, --pred_poisson).
  2. Score each theta by relative CDS recall (above). Interpolate theta* onto --cds-recall.
  3. Re-run at theta* with N independent Poisson draws.
  4. Emit the union of calls with **confidence = fraction of the N draws that called the ORF**.
     Poisson noise is the calibration lever (Check 4), so stability under re-draw is the natural
     per-ORF confidence: an ORF surviving every draw is robust, one appearing in 2/8 rides the
     noise.

Outputs (in --out):
  calls_calibrated.tsv   the theta* call set + `confidence` and `n_draws` columns
  calibration.json       both recall curves, theta*, per-class counts, confidence histogram

Env: the RiboCode env python (it shells out to ribocode_dropin.py).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
DROPIN = NEW / "scripts/ribocode_dropin.py"


def coord_key(orf_id):
    parts = orf_id.rsplit("_", 3)
    return "_".join(parts[-3:]) if len(parts) >= 4 else orf_id


def read_calls(path):
    """-> {coord_key: row dict}; empty dict if the file is missing."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(hdr):
                continue
            out[coord_key(p[hdr.index("ORF_ID")])] = dict(zip(hdr, p, strict=False))
    return out


def cds_transcripts(tx2cds, dumped):
    """Dumped transcripts that HAVE an annotated CDS -> the absolute-recall denominator."""
    keep = set()
    with open(tx2cds) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {k: h.index(k) for k in ("tx_id", "cds_len")}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if int(f[ix["cds_len"]]) > 0 and f[ix["tx_id"]] in dumped:
                keep.add(f[ix["tx_id"]])
    return keep


def cds_hit_tx(calls):
    """Transcripts whose ANNOTATED CDS was called."""
    return {r["transcript_id"] for r in calls.values() if r.get("ORF_type") == "annotated"}


def run_dropin(py, profiles, annot, out_dir, theta, seed, min_aa, pval):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [py, str(DROPIN), "--profiles", str(profiles), "--annot", str(annot),
           "--variant", "pred_preddepth", "--out", str(out_dir),
           "--min_aa", str(min_aa), "--pval", str(pval),
           "--pred_scale", str(theta), "--pred_poisson", "--seed", str(seed)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stderr[-2500:])
        raise SystemExit(f"ribocode_dropin failed at theta={theta} seed={seed}")
    return out_dir / "pred_preddepth_collapsed.txt"


def interp_theta(curve, target):
    """Linear interpolation of theta at a target recall on a (theta, recall) curve."""
    pts = sorted(curve, key=lambda t: t[0])
    for (t0, r0), (t1, r1) in zip(pts, pts[1:], strict=False):
        lo, hi = min(r0, r1), max(r0, r1)
        if lo <= target <= hi and r1 != r0:
            return round(t0 + (t1 - t0) * (target - r0) / (r1 - r0), 5), "interpolated"
    best = min(pts, key=lambda tr: abs(tr[1] - target))
    return best[0], "clamped (target outside the swept range)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--annot", required=True)
    ap.add_argument("--tx2cds", required=True,
                    help="tx_id/cds_len table for the absolute-recall denominator (reporting)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cds-recall", type=float, default=0.90, dest="cds_recall",
                    help="target RELATIVE canonical recall (vs the loosest theta). Default 0.90.")
    ap.add_argument("--thetas", default="0.02,0.05,0.1,0.2,0.5,1.0")
    ap.add_argument("--replicates", type=int, default=8)
    ap.add_argument("--min_aa", type=int, default=20)
    ap.add_argument("--pval", type=float, default=0.05)
    ap.add_argument("--python", default="/private/groups/carpenterlab/emalekos/conda_envs/"
                                        "ribocode/bin/python3")
    ap.add_argument("--replay", action="store_true",
                    help="reuse existing sweep/draw dirs instead of re-running the caller")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    import numpy as np
    dumped = {str(t) for t in np.load(a.profiles, allow_pickle=False)["tx_ids"]}
    denom_abs = cds_transcripts(a.tx2cds, dumped)
    print(f"dumped tx {len(dumped):,}; CDS-bearing {len(denom_abs):,}", file=sys.stderr)

    # ---- 1. sweep theta ----
    thetas = sorted(float(x) for x in a.thetas.split(","))
    hits, ncalls = {}, {}
    for th in thetas:
        d = out / f"sweep/theta_{th}"
        f = d / "pred_preddepth_collapsed.txt"
        if not (a.replay and f.exists()):
            f = run_dropin(a.python, a.profiles, a.annot, d, th, 0, a.min_aa, a.pval)
        c = read_calls(f)
        hits[th], ncalls[th] = cds_hit_tx(c), len(c)

    # ---- 2. score: relative recall (dial-able) + absolute recall (context) ----
    loosest = thetas[-1]
    denom_rel = hits[loosest]
    if not denom_rel:
        raise SystemExit("no annotated-CDS calls at the loosest theta; cannot calibrate")
    curve_rel, curve_abs = [], []
    for th in thetas:
        rel = len(hits[th] & denom_rel) / len(denom_rel)
        abs_ = (len(hits[th] & denom_abs) / len(denom_abs)) if denom_abs else float("nan")
        curve_rel.append((th, round(rel, 5)))
        curve_abs.append((th, round(abs_, 5)))
        print(f"  theta={th:<6} rel-CDS-recall={rel:.4f}  abs={abs_:.4f}  calls={ncalls[th]:,}",
              file=sys.stderr)

    theta_star, how = interp_theta(curve_rel, a.cds_recall)
    print(f"\ntheta* = {theta_star}  ({how})  for target relative CDS recall {a.cds_recall}",
          file=sys.stderr)

    # ---- 3/4. replicate draws at theta*; confidence = call stability ----
    seen, rows = Counter(), {}
    for s in range(a.replicates):
        d = out / f"draws/seed_{s}"
        f = d / "pred_preddepth_collapsed.txt"
        if not (a.replay and f.exists()):
            f = run_dropin(a.python, a.profiles, a.annot, d, theta_star, s, a.min_aa, a.pval)
        c = read_calls(f)
        seen.update(c.keys())
        rows.update(c)
        print(f"  draw {s}: {len(c):,} calls", file=sys.stderr)

    hdr = list(next(iter(rows.values())).keys()) if rows else []
    with open(out / "calls_calibrated.tsv", "w") as w:
        w.write("\t".join([*hdr, "confidence", "n_draws"]) + "\n")
        for k, r in sorted(rows.items(), key=lambda kv: -seen[kv[0]]):
            w.write("\t".join([r.get(c, "") for c in hdr]
                              + [f"{seen[k] / a.replicates:.4f}", str(a.replicates)]) + "\n")

    hit_star = cds_hit_tx(rows)
    summary = {
        "target_relative_cds_recall": a.cds_recall,
        "theta_star": theta_star,
        "theta_selection": how,
        "replicates": a.replicates,
        "relative_cds_recall_curve": [{"theta": t, "recall": r} for t, r in curve_rel],
        "absolute_cds_recall_curve": [{"theta": t, "recall": r} for t, r in curve_abs],
        "relative_denominator_tx": len(denom_rel),
        "absolute_denominator_tx": len(denom_abs),
        "achieved_relative_cds_recall": round(len(hit_star & denom_rel) / len(denom_rel), 5),
        "achieved_absolute_cds_recall": (round(len(hit_star & denom_abs) / len(denom_abs), 5)
                                         if denom_abs else None),
        "calls_union": len(rows),
        "calls_by_class": dict(Counter(r.get("ORF_type", "?") for r in rows.values())),
        "calls_at_confidence_1.0": sum(1 for k in rows if seen[k] == a.replicates),
        "calls_at_confidence_ge_0.5": sum(1 for k in rows if seen[k] >= a.replicates / 2),
        "confidence_histogram": {str(round(seen[k] / a.replicates, 3)): v for k, v in
                                 sorted(Counter(round(seen[k] / a.replicates, 3)
                                                for k in rows).items())},
        "inputs": {"profiles": a.profiles, "annot": a.annot, "tx2cds": a.tx2cds,
                   "min_aa": a.min_aa, "pval": a.pval},
    }
    json.dump(summary, open(out / "calibration.json", "w"), indent=2)
    print(f"\nwrote {out}/calls_calibrated.tsv ({len(rows):,} ORFs; "
          f"{summary['calls_at_confidence_1.0']:,} at confidence 1.0)", file=sys.stderr)
    print(f"wrote {out}/calibration.json", file=sys.stderr)


if __name__ == "__main__":
    main()
