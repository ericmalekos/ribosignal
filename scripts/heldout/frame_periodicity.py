#!/usr/bin/env python3
"""ORF-aware model-quality metric: predicted 3-nt FRAME-0 fraction within annotated CDS.

val_pearson (per-nt Pearson) is dominated by coverage magnitude and saturates, so it does NOT
track ORF-call quality (attn4 beats attn2 on ORF calls at equal val_pearson). RiboCode's detectORF
is a periodicity test, so the quantity that matters is how sharply the PREDICTED per-nt profile
concentrates on frame 0 (first nt of each codon) inside real coding sequence -- a candidate
checkpoint-selection / model-ranking signal that val_pearson misses.

Per annotated-CDS ORF (from the run's own dropin/real_collapsed.txt: transcript_id + 1-based
ORF_tstart/tstop, tstart defines frame 0), the model's predicted profile for that tx
(dropin/pred_profiles.npz) is summed by frame within the CDS: frame0_frac = sum(pred at
(p-tstart0)%3==0) / sum(pred over CDS). The OBSERVED profile's frame0_frac is the per-dataset
ceiling (same across models). Reported: signal-weighted pooled frame0_frac (headline), unweighted
per-tx mean, and the frame 0/1/2 split.

Usage: frame_periodicity.py --run <run_dir> [--min_cds_signal 50] [--label name]
"""
import argparse
import json

import numpy as np


def load_cds(path):
    """-> {transcript_id: (tstart0, tstop)} for annotated ORFs (tstart0 = 0-based CDS start)."""
    out = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[ci["ORF_type"]] != "annotated":
                continue
            tx = p[ci["transcript_id"]]
            ts, te = int(p[ci["ORF_tstart"]]) - 1, int(p[ci["ORF_tstop"]])   # 0-based [ts, te)
            out[tx] = (ts, te)   # one canonical CDS per tx (annotated)
    return out


def frame_sums(profile, ts, te):
    """(f0,f1,f2) sums of profile[ts:te] by (pos-ts)%3."""
    seg = profile[ts:te]
    if seg.size == 0:
        return 0.0, 0.0, 0.0
    ph = (np.arange(seg.size)) % 3          # pos-ts == arange, so frame = arange%3
    return float(seg[ph == 0].sum()), float(seg[ph == 1].sum()), float(seg[ph == 2].sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--min_cds_signal", type=float, default=50.0,
                    help="skip CDS with fewer OBSERVED real P-sites than this (too noisy); "
                         "pred is a normalized shape, so it is not thresholded")
    ap.add_argument("--label", default=None)
    ap.add_argument("--out-json", default=None, dest="out_json")
    a = ap.parse_args()

    cds = load_cds(f"{a.run}/dropin/real_collapsed.txt")
    d = np.load(f"{a.run}/dropin/pred_profiles.npz", allow_pickle=True)
    tx_ids = d["tx_ids"]
    lengths = d["lengths"].astype(np.int64)
    off = np.zeros(len(lengths) + 1, dtype=np.int64)
    off[1:] = np.cumsum(lengths)
    pred, obs = d["pred_flat"], d["obs_flat"]
    idx = {t: i for i, t in enumerate(tx_ids)}

    # pooled (signal-weighted) + per-tx frame-0 fractions
    pf = np.zeros(3)   # pooled predicted frame sums (f0,f1,f2)
    of = np.zeros(3)   # pooled observed frame sums
    per_pred, per_obs = [], []
    n = 0
    for tx, (ts, te) in cds.items():
        i = idx.get(tx)
        if i is None:
            continue
        s, e = off[i], off[i + 1]
        L = int(lengths[i])
        te2 = min(te, L)
        if ts >= te2:
            continue
        pfr = np.array(frame_sums(pred[s:e], ts, te2))
        ofr = np.array(frame_sums(obs[s:e], ts, te2))
        # threshold on OBSERVED real P-sites (pred is a normalized shape); need nonzero pred too
        if ofr.sum() < a.min_cds_signal or pfr.sum() <= 0:
            continue
        pf += pfr
        of += ofr
        per_pred.append(pfr[0] / pfr.sum())
        per_obs.append(ofr[0] / ofr.sum())
        n += 1

    pt = sum(pf) or 1.0
    ot = sum(of) or 1.0
    rec = {
        "label": a.label or a.run.split("/")[-1],
        "run": a.run, "n_cds": n, "min_cds_signal": a.min_cds_signal,
        "pred_frame0_frac_pooled": round(pf[0] / pt, 4),
        "pred_frame_split": [round(float(x) / pt, 3) for x in pf],
        "pred_frame0_frac_mean": round(float(np.mean(per_pred)) if per_pred else 0.0, 4),
        "obs_frame0_frac_pooled": round(of[0] / ot, 4),   # ceiling (same real data all models)
        "obs_frame0_frac_mean": round(float(np.mean(per_obs)) if per_obs else 0.0, 4),
        "pred_over_obs": round((pf[0] / pt) / (of[0] / ot), 4) if of[0] else None,
    }
    if a.out_json:
        json.dump(rec, open(a.out_json, "w"), indent=2)
    print(f"{rec['label']:28s}  n_CDS={n:6d}  pred_frame0={rec['pred_frame0_frac_pooled']:.3f} "
          f"(mean {rec['pred_frame0_frac_mean']:.3f}, split {rec['pred_frame_split']})  "
          f"obs_frame0={rec['obs_frame0_frac_pooled']:.3f}  pred/obs={rec['pred_over_obs']}")
    return rec


if __name__ == "__main__":
    main()
