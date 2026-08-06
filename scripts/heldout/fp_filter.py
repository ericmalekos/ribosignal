#!/usr/bin/env python3
"""O3 FP-filter -- re-score RiboCode candidate ORFs on model-predicted-density SHAPE features to
drop the smooth-spurious non-canonical calls that RiboCode's periodicity test passes.

RiboCode's single Wilcoxon frame test over-calls short non-canonical ORFs on the model's smooth
predicted density (a spurious weak ORF still looks cleanly periodic). This filter re-scores each
candidate with a logistic regression over the O3 shape features (codon uniformity PME, coverage
Gini/CV, 5' ramp, 3'-of-stop drop, in-frame fraction, length) and drops the ones that look
artificially uniform. It is a POST-RiboCode step, BAM-free: it runs on the same predicted
{tx: np.array} density the drop-in already uses.

Trained/scored PER ORF class (annotated / uORF / dORF / novel / ...) because base rates and feature
distributions differ; a pooled "non-canonical" model is also fit as the fallback for rare classes.

Ceiling caveat (O3 Phase 1): the model smooths REAL ORFs too, so shape features are only weakly
discriminative (novel CV AUC ~0.67). This is the modest post-hoc win; the higher-ceiling fix is the
anti-smoothing model (train_loto_noBrain_peaky.sbatch), re-featurized with orf_features.py.

Subcommands
  train  --features o3_features.tsv --out fp_filter.joblib
  apply  --model fp_filter.joblib --profiles pred_profiles.npz --candidates collapsed.txt
         --out filtered.txt [--target-recall 0.9 | --threshold 0.5] [--only novel,uORF,dORF]
"""
import argparse
import json
import os
import sys

import joblib
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orf_features import features, read_collapsed  # noqa: E402

FEAT = ["frac_f0", "f1max", "pme", "gini", "cv", "max_to_mean", "frac_zero",
        "ramp5", "drop3", "len_codons"]
# non-canonical classes the filter cleans (annotated/CDS is trusted, left alone by default)
NONCANON = {"uORF", "dORF", "novel", "Overlap_uORF", "Overlap_dORF", "internal"}
RECALLS = [0.95, 0.9, 0.8, 0.7, 0.5]


def load_features(path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append(f)
    otype = np.array([r[ci["ORF_type"]] for r in rows])
    X = np.array([[float(r[ci[k]]) for k in FEAT] for r in rows], dtype=np.float64)
    y = np.array([int(r[ci["label"]]) for r in rows], dtype=np.int64)
    return otype, X, y


def fit_one(X, y, seed=0):
    """5-fold CV out-of-fold probabilities + a final full-data model. Returns (scaler, clf, oof)."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler
    oof = np.zeros(len(y))
    n_splits = min(5, int(y.sum()), int((y == 0).sum()))
    if n_splits >= 2:
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        for tr, te in skf.split(X, y):
            sc = StandardScaler().fit(X[tr])
            clf = LogisticRegression(max_iter=2000, class_weight="balanced")
            clf.fit(sc.transform(X[tr]), y[tr])
            oof[te] = clf.predict_proba(sc.transform(X[te]))[:, 1]
    scaler = StandardScaler().fit(X)
    clf = LogisticRegression(max_iter=2000, class_weight="balanced")
    clf.fit(scaler.transform(X), y)
    return scaler, clf, oof


def pr_at_recalls(y, score):
    """precision achieved when the threshold is set to reach each target recall of the positives."""
    order = np.argsort(-score)
    ys = y[order]
    tp = np.cumsum(ys)
    kept = np.arange(1, len(ys) + 1)
    prec = tp / kept
    rec = tp / max(1, ys.sum())
    out = {}
    for R in RECALLS:
        idx = np.searchsorted(rec, R)          # rec is monotone non-decreasing
        idx = min(idx, len(rec) - 1)
        out[R] = {"precision": float(prec[idx]),
                  "threshold": float(score[order][idx]),
                  "kept_frac": float(kept[idx] / len(ys))}
    return out


def auc_auprc(y, score):
    from sklearn.metrics import average_precision_score, roc_auc_score
    if y.sum() == 0 or y.sum() == len(y):
        return float("nan"), float("nan")
    return float(roc_auc_score(y, score)), float(average_precision_score(y, score))


def cmd_train(a):
    otype, X, y = load_features(a.features)
    classes = sorted(set(otype))
    model = {"feat": FEAT, "classes": {}, "pooled": None}
    print(f"loaded {len(y):,} ORFs, {len(classes)} classes from {a.features}\n")
    print(f"{'class':16s} {'n':>7s} {'real':>6s} {'base':>6s} {'AUC':>6s} {'AUPRC':>6s} "
          f"{'P@.9':>6s} {'P@.7':>6s} {'P@.5':>6s} {'lift@.5':>7s}")
    report = {}
    for cls in classes:
        m = otype == cls
        if m.sum() < 20 or y[m].sum() < 5 or (y[m] == 0).sum() < 5:
            continue
        sc, clf, oof = fit_one(X[m], y[m], seed=a.seed)
        auc, ap = auc_auprc(y[m], oof)
        base = float(y[m].mean())
        pr = pr_at_recalls(y[m], oof)
        p9, p7, p5 = pr[0.9]["precision"], pr[0.7]["precision"], pr[0.5]["precision"]
        lift = p5 / base if base > 0 else float("nan")
        print(f"{cls:16s} {m.sum():7d} {int(y[m].sum()):6d} {base:6.3f} {auc:6.3f} {ap:6.3f} "
              f"{p9:6.3f} {p7:6.3f} {p5:6.3f} {lift:7.2f}")
        model["classes"][cls] = {"scaler": sc, "clf": clf, "pr": pr, "base": base,
                                 "auc": auc, "auprc": ap}
        report[cls] = {"n": int(m.sum()), "real": int(y[m].sum()), "base": base,
                       "auc": auc, "auprc": ap,
                       "pr": {str(k): v for k, v in pr.items()}}
    # pooled non-canonical model (fallback for rare classes + a data-richer alternative)
    mnc = np.array([o in NONCANON for o in otype])
    if mnc.sum() >= 20 and y[mnc].sum() >= 5:
        sc, clf, oof = fit_one(X[mnc], y[mnc], seed=a.seed)
        auc, ap = auc_auprc(y[mnc], oof)
        base = float(y[mnc].mean())
        pr = pr_at_recalls(y[mnc], oof)
        lift = pr[0.5]["precision"] / base if base > 0 else float("nan")
        print(f"{'POOLED_noncanon':16s} {mnc.sum():7d} {int(y[mnc].sum()):6d} {base:6.3f} "
              f"{auc:6.3f} {ap:6.3f} {pr[0.9]['precision']:6.3f} {pr[0.7]['precision']:6.3f} "
              f"{pr[0.5]['precision']:6.3f} {lift:7.2f}")
        model["pooled"] = {"scaler": sc, "clf": clf, "pr": pr, "base": base,
                           "auc": auc, "auprc": ap}
        report["POOLED_noncanon"] = {"n": int(mnc.sum()), "real": int(y[mnc].sum()),
                                     "base": base, "auc": auc, "auprc": ap,
                                     "pr": {str(k): v for k, v in pr.items()}}
    joblib.dump(model, a.out)
    json.dump(report, open(os.path.splitext(a.out)[0] + "_report.json", "w"), indent=2)
    print(f"\nsaved {a.out} (+ _report.json)")
    print("lift@.5 > 1 => the filter concentrates real ORFs vs the RiboCode base rate; "
          "P@.9 vs base = precision if we keep 90% of real non-canonical ORFs.")


def cmd_apply(a):
    model = joblib.load(a.model)
    only = set(a.only.split(",")) if a.only else None
    npz = np.load(a.profiles, allow_pickle=True)
    tx = np.array([str(x) for x in npz["tx_ids"]])
    lengths = npz["lengths"].astype(np.int64)
    off = np.concatenate([[0], np.cumsum(lengths)])
    pred = npz["pred_flat"]
    row_of = {t: i for i, t in enumerate(tx)}

    cand, ci = read_collapsed(a.candidates)
    kept, dropped = [], {}
    n_scored = n_nofeat = 0
    for r in cand:
        cls = r[ci["ORF_type"]]
        # classes we do not filter (CDS/annotated, or anything not in --only) pass through
        if (only is not None and cls not in only) or \
           (only is None and cls not in model["classes"] and model["pooled"] is None):
            kept.append(r)
            continue
        entry = model["classes"].get(cls) or model["pooled"]
        if entry is None:
            kept.append(r)
            continue
        txid = r[ci["transcript_id"]]
        i = row_of.get(txid)
        if i is None:
            kept.append(r)          # cannot score -> keep (do not silently drop)
            continue
        ts = int(r[ci["ORF_tstart"]])
        olen = int(r[ci["ORF_length"]])
        s = off[i] + (ts - 1)
        prof = np.asarray(pred[s:s + olen], dtype=np.float64)
        fv = features(prof) if prof.size >= 9 else None
        if fv is None:
            n_nofeat += 1
            kept.append(r)
            continue
        x = np.array([[fv[k] for k in model["feat"]]], dtype=np.float64)
        score = float(entry["clf"].predict_proba(entry["scaler"].transform(x))[0, 1])
        if a.threshold is not None:
            thr = a.threshold
        else:
            thr = entry["pr"][a.target_recall]["threshold"]
        n_scored += 1
        if score >= thr:
            kept.append(r)
        else:
            dropped[cls] = dropped.get(cls, 0) + 1

    hdr = None
    with open(a.candidates) as fh:
        hdr = fh.readline()
    with open(a.out, "w") as o:
        o.write(hdr)
        for r in kept:
            o.write("\t".join(r) + "\n")
    tot_drop = sum(dropped.values())
    at = (f"threshold={a.threshold}" if a.threshold is not None
          else f"target_recall={a.target_recall}")
    print(f"scored {n_scored:,} candidates; kept {len(kept):,}, dropped {tot_drop:,} ({at})")
    for cls in sorted(dropped):
        print(f"  dropped {dropped[cls]:5d}  {cls}")
    if n_nofeat:
        print(f"  ({n_nofeat} too-short-to-featurize kept)")
    print(f"wrote {a.out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("train")
    t.add_argument("--features", required=True)
    t.add_argument("--out", required=True)
    t.add_argument("--seed", type=int, default=0)
    t.set_defaults(func=cmd_train)
    p = sub.add_parser("apply")
    p.add_argument("--model", required=True)
    p.add_argument("--profiles", required=True)
    p.add_argument("--candidates", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--target-recall", type=float, default=0.9, choices=RECALLS,
                   dest="target_recall")
    p.add_argument("--threshold", type=float, default=None,
                   help="override the target-recall threshold with a fixed score cutoff")
    p.add_argument("--only", default="uORF,dORF,novel,Overlap_uORF,Overlap_dORF,internal",
                   help="comma ORF types to filter; others pass through (default: non-canonical)")
    p.set_defaults(func=cmd_apply)
    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
