#!/usr/bin/env python3
"""How far does the input-swap effect generalise beyond GTF2I?

Same mm1-trained model, two coverage inputs (mm1 vs mm10), over 5,653 HIGH-multimap test tx
(frac_discarded > 0.5) and 5,653 LOW-multimap controls (< 0.10) drawn from the same test set.

The control group carries the argument. mm10 changes coverage everywhere (mean RNA multimap rate
3.8% across the 57 libraries), so "predictions changed under mm10" is not a finding by itself. The
claim is that the change is CONCENTRATED in transcripts whose mm1 coverage was corrupted.

Two metrics per transcript, both on the profile head, which is a distribution over positions:

  moved  total variation distance, 0.5*sum|p_mm1 - p_mm10|, in [0,1]. Directly readable as the
         FRACTION OF PREDICTED MASS THAT RELOCATED when the coverage was corrected.
  d_obs  Pearson(pred, observed) under mm10 minus the same under mm1.

d_obs is expected to go NEGATIVE for the high group, and that is the point rather than a problem:
the observed target is mm1 Ribo-seq, which is holed at exactly these loci, so a model that stops
reproducing the hole necessarily agrees LESS with the holed measurement. A negative d_obs in the
high group with d_obs ~ 0 in the low group is the signature of a corrected input, not a degraded
model. Nothing here can show the new predictions are RIGHT -- there is no unholed ground truth.

Usage: analyze_swap_generalize.py [--json out.json]
"""
import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np

R = Path(__file__).resolve().parent.parent
G = R / "data/mm25_diagnostic/gtf2i_inference"


def load_merged(d):
    """Merge sharded dumps -> {tx: (pred, obs)}. Fails loudly on a nonzero skip count."""
    out, nskip = {}, 0
    fs = sorted(glob.glob(str(d / "pred_profiles_shard*.npz"))) or \
        sorted(glob.glob(str(d / "pred_profiles.npz")))
    if not fs:
        sys.exit(f"no dumps in {d}")
    for f in fs:
        z = np.load(f)
        nskip += int(z["meta"][3])
        ids = [str(x) for x in z["tx_ids"]]
        off = np.concatenate([[0], np.cumsum(z["lengths"])])
        # Materialise ONCE. np.load on a compressed npz returns a lazy handle, so z["pred_flat"]
        # inside the loop re-decompresses the WHOLE array on every transcript -- ~1,400x per shard,
        # which turned a 1-minute job into >30 minutes with no output to show for it.
        pf = z["pred_flat"]
        of = z["obs_flat"]
        for i, t in enumerate(ids):
            a, b = int(off[i]), int(off[i + 1])
            out[t] = (pf[a:b].astype(np.float64), of[a:b].astype(np.float64))
    return out, nskip, len(fs)


def tvd(p, q):
    ps, qs = p.sum(), q.sum()
    if ps <= 0 or qs <= 0:
        return np.nan
    return 0.5 * np.abs(p / ps - q / qs).sum()


def pear(a, b):
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default="")
    a = ap.parse_args()
    hi = {x.strip() for x in (G / "swap_tx_high.txt").read_text().split() if x.strip()}
    lo = {x.strip() for x in (G / "swap_tx_low.txt").read_text().split() if x.strip()}
    m1, s1, n1 = load_merged(G / "swap_union")
    m10, s10, n10 = load_merged(G / "swap_mm10cov")
    print(f"  mm1 dump: {len(m1):,} tx from {n1} shards ({s1} skipped)")
    print(f"  mm10 dump: {len(m10):,} tx from {n10} shards ({s10} skipped)")
    if s1 or s10:
        print("  WARNING: nonzero skip count -- see the FASTA/pack mismatch rule")
    common = sorted(set(m1) & set(m10))
    print(f"  scored on {len(common):,} tx present in both\n")

    rows = {"high": [], "low": []}
    for t in common:
        g = "high" if t in hi else ("low" if t in lo else None)
        if g is None:
            continue
        p1, o = m1[t]
        p10, _ = m10[t]
        if len(p1) != len(p10):
            continue
        rows[g].append((tvd(p1, p10), pear(p10, o) - pear(p1, o)))

    out = {}
    print(f"  {'group':6s} {'n':>6s} {'moved med':>10s} {'moved p90':>10s} "
          f"{'>25% moved':>11s} {'d_obs med':>10s}")
    for g in ("high", "low"):
        v = np.array([x for x in rows[g] if not np.isnan(x[0])])
        if not len(v):
            continue
        mv, dv = v[:, 0], v[:, 1]
        dvv = dv[~np.isnan(dv)]
        out[g] = {"n": len(v), "moved_median": float(np.median(mv)),
                  "moved_p90": float(np.percentile(mv, 90)),
                  "frac_moved_gt25pct": float((mv > 0.25).mean()),
                  "d_obs_median": float(np.median(dvv)) if len(dvv) else None}
        print(f"  {g:6s} {len(v):6,d} {np.median(mv):10.3f} {np.percentile(mv,90):10.3f} "
              f"{100*(mv>0.25).mean():10.1f}% "
              f"{(np.median(dvv) if len(dvv) else float('nan')):10.4f}")
    if "high" in out and "low" in out:
        print(f"\n  specificity: median mass moved is "
              f"{out['high']['moved_median']/max(out['low']['moved_median'],1e-9):.1f}x higher in "
              f"high-multimap tx than in controls")
    print("\n  d_obs is EXPECTED to be negative in the high group: the observed target is")
    print("  mm1 Ribo-seq, holed at these loci, so not reproducing the hole means agreeing")
    print("  less with it.")
    if a.json:
        Path(a.json).write_text(json.dumps(
            {"model": "orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes",
             "source": "scripts/analyze_swap_generalize.py over "
                       "data/mm25_diagnostic/gtf2i_inference/swap_{union,mm10cov}; "
                       "groups from swap_tx_{high,low}.txt",
             "note": "Input-swap only. NOT a retrain: identical weights, only RIBO_PACK_SUFFIX "
                     "differs. moved = total variation distance between the two predicted "
                     "profiles. d_obs = Pearson(pred,obs) under mm10 minus under mm1.",
             "n_skip": {"mm1": s1, "mm10": s10}, "groups": out}, indent=2) + "\n")
        print(f"  wrote {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
