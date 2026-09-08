#!/usr/bin/env python3
"""Score the demo: how close is the predicted profile to the real one, and do the ORF
calls made on it match the calls made on real Ribo-seq?

Two questions, kept separate because they fail differently:

  profile agreement   per-transcript Pearson between predicted and observed P-sites.
                      This is the model's own metric and is directly comparable to
                      release/orf_v2_*/test_metrics.json, which reports 0.6595 (attn) and
                      0.6799 (mamba4) as 3-seed medians over 70,883 held-out transcripts.
  ORF-call overlap    precision / recall / F1 of the calls RiboCode makes on the PREDICTED
                      density against the calls it makes on the OBSERVED density, same
                      caller and same parameters, so the only difference is the input.

Spearman here averages tied ranks. That is not a detail: the project's stored
`spearman_median` was computed with ordinal ranks and no tie averaging, which on profiles
that are 74-95% zeros ranks position rather than signal and drove the reported value
negative. Those stored numbers are renamed `spearman_median_INVALID_ordinal_ties`
throughout. This computes it correctly, so it is comparable to nothing in the release --
only to itself.

    python scripts/demo/score_demo.py --pred-dir pred --calls-dir calls --out RESULTS.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_dropin_calls import load_calls  # noqa: E402  the project's own parser

ARCHS = ("attn", "mamba4")


def rank_avg(x):
    """Ranks with TIES AVERAGED. The thing the stored metric got wrong."""
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x), dtype=np.float64)
    r[order] = np.arange(len(x), dtype=np.float64)
    xs = x[order]
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[i]:
            j += 1
        if j > i:
            r[order[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return r


def pearson(a, b):
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def periodicity(v):
    """Fraction of signal in the dominant 3-nt frame; 1/3 means none."""
    t = v.sum()
    if t <= 0:
        return np.nan
    f = [v[i::3].sum() for i in range(3)]
    return float(max(f) / t)


def profile_scores(npz_path, min_signal=50):
    z = np.load(npz_path, allow_pickle=False)
    off = np.concatenate([[0], np.cumsum(z["lengths"])])
    pr, sp, pp, po, n = [], [], [], [], 0
    for j in range(len(z["lengths"])):
        p = z["pred_flat"][off[j]:off[j + 1]].astype(np.float64)
        o = z["obs_flat"][off[j]:off[j + 1]].astype(np.float64)
        if o.sum() < min_signal:
            continue
        n += 1
        pr.append(pearson(p, o))
        sp.append(pearson(rank_avg(p), rank_avg(o)))
        pp.append(periodicity(p))
        po.append(periodicity(o))
    med = lambda v: float(np.nanmedian(v)) if v else float("nan")  # noqa: E731
    return {"n_scored": n, "min_signal": min_signal,
            "pearson_median": med(pr), "spearman_median_tie_corrected": med(sp),
            "period_pred_median": med(pp), "period_obs_median": med(po)}


def prf(pred, real):
    tp = len(set(pred) & set(real))
    p = tp / len(pred) if pred else 0.0
    r = tp / len(real) if real else 0.0
    return {"n_pred": len(pred), "n_real": len(real), "tp": tp,
            "precision": p, "recall": r,
            "f1": (2 * p * r / (p + r)) if (p + r) else 0.0}


def call_scores(real_path, pred_path, max_pval=0.05, min_len=90):
    real = load_calls(str(real_path), key_mode="genomic", max_pval=max_pval, min_len=min_len)
    pred = load_calls(str(pred_path), key_mode="genomic", max_pval=max_pval, min_len=min_len)
    out = {"overall": prf(pred, real), "by_type": {}}
    types = {v["type"] for v in real.values()} | {v["type"] for v in pred.values()}
    for t in sorted(types):
        rt = {k for k, v in real.items() if v["type"] == t}
        pt = {k for k, v in pred.items() if v["type"] == t}
        if rt or pt:
            out["by_type"][t] = prf(pt, rt)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pred-dir", default="pred")
    ap.add_argument("--calls-dir", default="calls")
    ap.add_argument("--out", default="RESULTS.json")
    ap.add_argument("--min-signal", type=int, default=50)
    a = ap.parse_args()

    res = {"profile": {}, "calls": {}}
    for arch in ARCHS:
        npz = Path(a.pred_dir) / arch / "pred_profiles.npz"
        if npz.is_file():
            res["profile"][arch] = profile_scores(npz, a.min_signal)
    for arch in ARCHS:
        real = Path(a.calls_dir) / arch / "real_collapsed.txt"
        if not real.is_file():
            continue
        for label, pred in ((f"{arch}:pred_preddepth",
                             Path(a.calls_dir) / arch / "pred_preddepth_collapsed.txt"),
                            (f"{arch}:pred_preddepth+poisson0.05",
                             Path(a.calls_dir) / f"{arch}_poisson" / "pred_preddepth_collapsed.txt")):
            if pred.is_file():
                res["calls"][label] = call_scores(real, pred)

    print("\n" + "=" * 78)
    print("PROFILE AGREEMENT  (predicted vs observed P-sites, per transcript)")
    print("=" * 78)
    if not res["profile"]:
        print("  no pred_profiles.npz found")
    for arch, s in res["profile"].items():
        print(f"  {arch:<8} n={s['n_scored']:<6} "
              f"pearson {s['pearson_median']:+.4f}   "
              f"spearman(tie-corrected) {s['spearman_median_tie_corrected']:+.4f}")
        print(f"           periodicity  pred {s['period_pred_median']:.4f}  "
              f"obs {s['period_obs_median']:.4f}   (1/3 = none)")
    print("\n  release/orf_v2_*/test_metrics.json, for reference, reports pearson_median")
    print("  0.6585 (attn) / 0.6851 (mamba4) on 70,883 held-out transcripts, seed 0.")
    print("  These are a different, much smaller transcript set: not a reproduction of it.")

    print("\n" + "=" * 78)
    print("ORF CALLS  (RiboCode on predicted density vs on OBSERVED density)")
    print("=" * 78)
    if not res["calls"]:
        print("  no call sets found")
    for label, s in res["calls"].items():
        o = s["overall"]
        print(f"\n  {label}")
        print(f"    overall   P {o['precision']:.3f}  R {o['recall']:.3f}  F1 {o['f1']:.3f}"
              f"   ({o['tp']} of {o['n_pred']} predicted, {o['n_real']} real)")
        for t, v in sorted(s["by_type"].items(), key=lambda kv: -kv[1]["n_real"])[:7]:
            print(f"      {t:<16} P {v['precision']:.3f}  R {v['recall']:.3f}  "
                  f"F1 {v['f1']:.3f}   (real {v['n_real']}, pred {v['n_pred']})")

    Path(a.out).write_text(json.dumps(res, indent=2) + "\n")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
