#!/usr/bin/env python3
"""Comprehensive per-model RiboCode ORF-call metrics on a holdout, per the standing evaluation rule
(feedback_orf_call_eval_standing + feedback_two_arm_orf_calling). Manuscript-grade: emits a JSON + a
markdown block, and records every input path + count so results are fully reproducible.

Scores a model's RiboCode drop-in calls against that run's OWN observed (real) calls on the holdout,
split by ORF class (annotated=CDS, uORF, novel, dORF, non-canonical aggregate, all):
  - pred_obsdepth vs real : SHAPE -- predicted profile at the REAL per-tx depth (did the model
                            get the ORF distribution right, independent of the count head).
  - pred_preddepth vs real: STANDALONE theta=1 deterministic -- the fully de novo call set (shape x
                            count head). precision/recall/F1 + composition + over-call ratio.
  - Poisson arm (optional): two-arm partner -- from a theta sweep, the CDS-anchored operating
                            point (theta whose CDS recall ~ --cds-recall) + per-class precision.

Coord-key matched on gstart_gstop_len (annotation-version robust). Reference = the model's own
dropin/real_collapsed.txt (the holdout's observed RiboCode calls on that model's universe).

Usage: orf_call_metrics.py --run <dir> --label "<model>" [--universe fibroblast|union]
       [--poisson-sweep <dir>] [--cds-recall 0.90] [--out-json f] [--append-md f]
"""
import argparse
import glob
import json
import os

CLASSES = ["annotated", "uORF", "novel", "dORF"]


def coord_key(orf_id):
    parts = orf_id.rsplit("_", 3)
    return "_".join(parts[-3:]) if len(parts) >= 4 else orf_id


def load_calls(path):
    """-> {coord_key: ORF_type}, or None if the file is absent."""
    if not path or not os.path.exists(path):
        return None
    d = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            d[coord_key(p[ci["ORF_ID"]])] = p[ci["ORF_type"]]
    return d


def sel(d, cls):
    if cls == "all":
        return set(d)
    if cls == "non-canonical":
        return {k for k, t in d.items() if t != "annotated"}
    return {k for k, t in d.items() if t == cls}


def load_lengths(path):
    """-> {coord_key: aa_length}. ORF_length is nt incl. stop; aa = nt//3 - 1."""
    if not path or not os.path.exists(path):
        return {}
    out = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            out[coord_key(p[ci["ORF_ID"]])] = int(p[ci["ORF_length"]]) // 3 - 1
    return out


def min_aa_sweep(pred, pred_len, real, real_len, thresholds):
    """For each min-aa floor, restrict BOTH predicted and observed calls to ORFs >= that aa
    length and recompute per-class precision/recall/F1. Shows how the metrics depend on the
    ORF-length threshold (the tiny-ORF tail is uORF-heavy and the most over-call-prone)."""
    sweep = []
    for thr in thresholds:
        pk = {k for k in pred if pred_len.get(k, 0) >= thr}
        rk = {k for k in real if real_len.get(k, 0) >= thr}
        pf = {k: pred[k] for k in pk}
        rf = {k: real[k] for k in rk}
        row = {"min_aa": thr}
        for cls in ["all", "annotated", "uORF", "novel", "dORF", "non-canonical"]:
            row[cls] = prf(sel(pf, cls), sel(rf, cls))
        sweep.append(row)
    return sweep


def prf(pred, ref):
    inter = len(pred & ref)
    p = inter / len(pred) if pred else 0.0
    r = inter / len(ref) if ref else 0.0
    f = 2 * p * r / (p + r) if (p + r) else 0.0
    return {"precision": round(p, 4), "recall": round(r, 4), "f1": round(f, 4),
            "n_pred": len(pred), "n_ref": len(ref), "n_match": inter}


def score_variant(pred, real):
    """per-class precision/recall/F1 of a predicted call set vs the observed (real) calls."""
    out = {}
    for cls in ["all", "annotated", "uORF", "novel", "dORF", "non-canonical"]:
        m = prf(sel(pred, cls), sel(real, cls))
        if cls in CLASSES or cls in ("all", "non-canonical"):
            m["overcall_ratio"] = round(m["n_pred"] / m["n_ref"], 3) if m["n_ref"] else None
        out[cls] = m
    return out


def best_val(run):
    h = f"{run}/history.json"
    if not os.path.exists(h):
        return None, None
    hist = json.load(open(h))
    vps = [(e.get("val", {}).get("all", {}).get("pearson_median"), e["epoch"]) for e in hist]
    vps = [(v, e) for v, e in vps if v is not None]
    return (max(vps) if vps else (None, None))


def at_floor(d, len_map, floor):
    """restrict a {coord_key: type} call set to ORFs >= floor aa."""
    if d is None:
        return None
    return {k: v for k, v in d.items() if len_map.get(k, 0) >= floor}


def cds_anchored_point(sweep_dir, real, target_recall, floor=0):
    """Poisson theta whose CDS recall (vs real CDS) is nearest target; its per-class score.
    Both the theta calls and the CDS-recall anchor are taken at the primary >= floor aa cut."""
    if not sweep_dir or not os.path.isdir(sweep_dir):
        return None
    cand = []
    for d in sorted(glob.glob(f"{sweep_dir}/theta_*"),
                    key=lambda p: float(p.rsplit("_", 1)[1])):
        f = f"{d}/pred_preddepth_collapsed.txt"
        calls = at_floor(load_calls(f), load_lengths(f), floor)
        if not calls:
            continue
        th = float(d.rsplit("_", 1)[1])
        cds_rec = prf(sel(calls, "annotated"), sel(real, "annotated"))["recall"]
        cand.append((th, cds_rec, calls))
    if not cand:
        return None
    th, cds_rec, calls = min(cand, key=lambda x: abs(x[1] - target_recall))
    return {"theta": th, "cds_recall": round(cds_rec, 4),
            "target_cds_recall": target_recall, "score": score_variant(calls, real)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--universe", default="fibroblast")
    ap.add_argument("--poisson-sweep", default=None, dest="poisson_sweep")
    ap.add_argument("--cds-recall", type=float, default=0.90, dest="cds_recall")
    ap.add_argument("--primary-min-aa", type=int, default=20, dest="primary_min_aa",
                    help="primary ORF-length floor (aa) for the headline P/R/F1 tables; calls are "
                         "restricted to ORFs >= this length (default 20, RiboCode's own default)")
    ap.add_argument("--min-aa-sweep", default="5,10,15,20,30,40", dest="min_aa_sweep",
                    help="comma list of min ORF aa lengths; predicted+observed calls both "
                         "restricted to ORFs >= each floor, P/R/F1 recomputed (standalone theta=1)")
    ap.add_argument("--out-json", default=None, dest="out_json")
    ap.add_argument("--append-md", default=None, dest="append_md")
    a = ap.parse_args()

    rpath = f"{a.run}/dropin/real_collapsed.txt"
    opath = f"{a.run}/dropin/pred_obsdepth_collapsed.txt"
    ppath = f"{a.run}/dropin/pred_preddepth_collapsed.txt"
    real0, obs0, pred0 = load_calls(rpath), load_calls(opath), load_calls(ppath)
    if real0 is None or pred0 is None:
        raise SystemExit(f"missing drop-in calls under {a.run}/dropin/ (real/pred_preddepth)")
    real_len, obs_len, pred_len = load_lengths(rpath), load_lengths(opath), load_lengths(ppath)
    thresholds = [int(x) for x in a.min_aa_sweep.split(",") if x.strip()]
    # sweep uses the FULL call sets; primary tables use the >= primary-floor cut
    sweep = min_aa_sweep(pred0, pred_len, real0, real_len, thresholds)
    pf = a.primary_min_aa
    real = at_floor(real0, real_len, pf)
    obs = at_floor(obs0, obs_len, pf)
    pred = at_floor(pred0, pred_len, pf)
    bv, be = best_val(a.run)

    rec = {
        "label": a.label, "run": a.run, "universe": a.universe, "primary_min_aa": pf,
        "val_pearson": bv, "val_epoch": be,
        "inputs": {"real": rpath, "pred_obsdepth": opath, "pred_preddepth": ppath,
                   "poisson_sweep": a.poisson_sweep},
        "n_real": {c: len(sel(real, c)) for c in ["all", "annotated", "uORF", "novel", "dORF"]},
        "obsdepth_vs_real": score_variant(obs, real) if obs else None,
        "preddepth_theta1_vs_real": score_variant(pred, real),
        "poisson_cds_anchored": cds_anchored_point(a.poisson_sweep, real, a.cds_recall, pf),
        "min_aa_sweep": sweep,
    }
    if a.out_json:
        json.dump(rec, open(a.out_json, "w"), indent=2)

    # markdown block
    def row(name, s):
        c, u, n, dd = s["annotated"], s["uORF"], s["novel"], s["dORF"]
        return (f"| {name} | {c['precision']:.3f}/{c['recall']:.3f}/{c['f1']:.3f} | "
                f"{u['precision']:.3f}/{u['recall']:.3f} | "
                f"{n['precision']:.3f}/{n['recall']:.3f} | "
                f"{dd['precision']:.3f}/{dd['recall']:.3f} | {s['non-canonical']['f1']:.3f} |")
    md = [f"\n### {a.label}  ({a.universe} universe)",
          f"- run: `{a.run}`  |  val_pearson {bv:.4f} (e{be})" if bv else f"- run: `{a.run}`",
          f"- **primary floor = ORFs >= {pf} aa** (headline tables); sweep below spans full range.",
          f"- reference = observed Hepatocytes calls (>= {pf} aa): {rec['n_real']['all']:,} "
          f"(CDS {rec['n_real']['annotated']:,}, uORF {rec['n_real']['uORF']:,}, "
          f"novel {rec['n_real']['novel']:,}, dORF {rec['n_real']['dORF']:,})",
          "",
          "| arm | CDS P/R/F1 | uORF P/R | novel P/R | dORF P/R | non-canon F1 |",
          "|---|---|---|---|---|---|"]
    if obs:
        md.append(row("pred_obsdepth (shape @ real depth)", rec["obsdepth_vs_real"]))
    md.append(row("pred_preddepth (standalone theta=1)", rec["preddepth_theta1_vs_real"]))
    p1 = rec["preddepth_theta1_vs_real"]
    md.append(f"- standalone theta=1 over-call: total {p1['all']['n_pred']:,} calls vs "
              f"{p1['all']['n_ref']:,} real (all {p1['all']['overcall_ratio']}x; "
              f"novel {p1['novel']['overcall_ratio']}x, dORF {p1['dORF']['overcall_ratio']}x)")
    pz = rec["poisson_cds_anchored"]
    if pz:
        md.append(row(f"Poisson CDS-anchored (th={pz['theta']}, CDSrec {pz['cds_recall']:.2f})",
                      pz["score"]))

    # min-aa length sweep (standalone theta=1): how P/R/F1 shift as tiny ORFs are excluded
    def prf3(m):
        return f"{m['precision']:.3f}/{m['recall']:.3f}/{m['f1']:.3f}"
    md += ["",
           "_Length sweep_ (standalone theta=1; predicted+observed calls both restricted to "
           "ORFs >= min_aa; P/R/F1). CDS length-invariant (ref); tiny-ORF tail is uORF-heavy:",
           "| min_aa | n_pred/n_real | CDS F1 | uORF P/R/F1 | novel P/R/F1 | dORF P/R/F1 "
           "| non-canon P/R/F1 |",
           "|---|---|---|---|---|---|---|"]
    for s in sweep:
        md.append(f"| {s['min_aa']} | {s['all']['n_pred']:,}/{s['all']['n_ref']:,} | "
                  f"{s['annotated']['f1']:.3f} | {prf3(s['uORF'])} | {prf3(s['novel'])} | "
                  f"{prf3(s['dORF'])} | {prf3(s['non-canonical'])} |")
    md_txt = "\n".join(md) + "\n"
    if a.append_md:
        with open(a.append_md, "a") as f:
            f.write(md_txt)
    print(md_txt)


if __name__ == "__main__":
    main()
