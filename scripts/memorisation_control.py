#!/usr/bin/env python3
"""Phase 2: does the LOTO transcript overlap inflate held-out-tissue performance?

LOTO holds out a TISSUE, not transcripts. Measured on the Hepatocytes holdout: of 33,918 test
transcripts, 27,007 (79.6%) appeared in training gradients via some held-in tissue and 6,813 (20.1%)
in validation, leaving only 98 (0.3%) never seen in any held-in tissue. For a given transcript the
SEQUENCE and ORF TRACK are byte-identical between its train and test appearances -- only coverage and
target differ. So the model has seen the exact input sequence of 99.7% of its "held-out" test set.

Cross-species mouse (zero transcript overlap) scores within ~0.02-0.04 F1 of Hepatocytes, which
argues the overlap is not inflating much. But mouse also differs in species, tissue and library, so
it does not isolate the variable. These 98 transcripts are the ONLY clean within-human,
within-tissue, same-library isolation available.

Test: per-transcript agreement between predicted and observed profile, for the never-seen transcripts
versus a MATCHED sample of seen ones. Matching is on log10 observed CDS counts and transcript length,
because both drive profile agreement on their own and an unmatched comparison would mostly measure
depth rather than novelty.

n is small. Bootstrap CIs are reported and no per-class claim is made from this alone.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np


def per_tx_agreement(npz, min_counts=50):
    """Pearson r between predicted and observed per-nt profile, per transcript."""
    z = np.load(npz, allow_pickle=True)
    tx = [t.decode() if isinstance(t, bytes) else str(t) for t in z["tx_ids"]]
    lens = np.asarray(z["lengths"], dtype=np.int64)
    off = np.concatenate([[0], np.cumsum(lens)])
    obs, pred = np.asarray(z["obs_flat"]), np.asarray(z["pred_flat"])
    out = {}
    for i, t in enumerate(tx):
        a, b = off[i], off[i + 1]
        o = np.asarray(obs[a:b], dtype=np.float64)
        p = np.asarray(pred[a:b], dtype=np.float64)
        if o.sum() < min_counts or o.std() == 0 or p.std() == 0:
            continue
        out[t] = dict(r=float(np.corrcoef(o, p)[0, 1]), counts=float(o.sum()), length=int(b - a))
    return out


def boot_ci(vals, n=5000, seed=0):
    rng = np.random.default_rng(seed)
    v = np.asarray(vals, dtype=float)
    if len(v) < 2:
        return (float("nan"), float("nan"))
    m = rng.choice(v, size=(n, len(v)), replace=True).mean(axis=1)
    return (float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True)
    ap.add_argument("--holdout", default="Hepatocytes")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-match", type=int, default=10, help="matched seen tx per unseen tx")
    a = ap.parse_args()

    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    import dataset as D
    train, val, _, test_tx = D.loto_split(a.holdout)
    seen = set()
    for _, txs in train:
        seen.update(txs)
    for _, txs in val:
        seen.update(txs)
    unseen_set = {t for t in test_tx if t not in seen}
    print(f"  holdout={a.holdout}: {len(test_tx):,} test tx, "
          f"{len(unseen_set):,} never seen in any held-in tissue")

    agree = per_tx_agreement(pathlib.Path(a.npz))
    u = {t: v for t, v in agree.items() if t in unseen_set}
    s = {t: v for t, v in agree.items() if t not in unseen_set}
    print(f"  scorable (>=50 obs counts): unseen {len(u):,}, seen {len(s):,}")
    if len(u) < 5:
        print("  ABORT: too few scorable unseen transcripts to compare", file=sys.stderr)
        return 1

    # Match on log10 counts and length -- both drive profile agreement independently of novelty, so
    # an unmatched comparison would mostly report a depth difference.
    sk = list(s)
    sf = np.array([[np.log10(s[t]["counts"] + 1), np.log10(s[t]["length"] + 1)] for t in sk])
    used, matched = set(), []
    for t, v in u.items():
        q = np.array([np.log10(v["counts"] + 1), np.log10(v["length"] + 1)])
        d = np.linalg.norm(sf - q, axis=1)
        taken = 0
        for j in np.argsort(d):
            if sk[j] in used:
                continue
            used.add(sk[j]); matched.append(sk[j]); taken += 1
            if taken >= a.n_match:
                break

    ur = [v["r"] for v in u.values()]
    mr = [s[t]["r"] for t in matched]
    allr = [v["r"] for v in s.values()]
    res = dict(
        holdout=a.holdout, npz=str(a.npz),
        n_unseen=len(ur), n_matched_seen=len(mr), n_all_seen=len(allr),
        mean_r_unseen=float(np.mean(ur)), ci_unseen=boot_ci(ur),
        mean_r_matched_seen=float(np.mean(mr)), ci_matched_seen=boot_ci(mr),
        mean_r_all_seen=float(np.mean(allr)),
        delta_unseen_minus_matched=float(np.mean(ur) - np.mean(mr)),
        median_counts_unseen=float(np.median([v["counts"] for v in u.values()])),
        median_counts_matched=float(np.median([s[t]["counts"] for t in matched])),
    )
    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2))
    print(f"  unseen        n={res['n_unseen']:<5} mean r={res['mean_r_unseen']:.4f} "
          f"CI [{res['ci_unseen'][0]:.4f}, {res['ci_unseen'][1]:.4f}]")
    print(f"  matched seen  n={res['n_matched_seen']:<5} mean r={res['mean_r_matched_seen']:.4f} "
          f"CI [{res['ci_matched_seen'][0]:.4f}, {res['ci_matched_seen'][1]:.4f}]")
    print(f"  all seen      n={res['n_all_seen']:<5} mean r={res['mean_r_all_seen']:.4f}")
    print(f"  delta (unseen - matched seen) = {res['delta_unseen_minus_matched']:+.4f}")
    print(f"  median obs counts: unseen {res['median_counts_unseen']:,.0f} vs "
          f"matched {res['median_counts_matched']:,.0f}")
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
