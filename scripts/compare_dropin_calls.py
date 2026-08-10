#!/usr/bin/env python3
"""RiboCode drop-in test (Task 14) comparison: do the ORF calls RiboCode makes on the model's
PREDICTED per-nt density track the calls it makes on the real experimental density?

Reads the three drop-in RiboCode outputs (real / pred_obsdepth / pred_preddepth _collapsed.txt from
ribocode_dropin.py) plus the OFFICIAL <Tissue>_collapsed.txt, restricts every call set to the
model's held-out test transcripts (the tx in pred_profiles.npz), and reports, keyed by (tx, tstart):
  * real vs official       -- harness validation: does running detectORF on the real packed counts
                              reproduce the official per-tissue calls? (overlap on shared tx)
  * pred_obsdepth vs real  -- PRIMARY: is the predicted profile SHAPE (at real depth) enough to
                              recover the calls? precision / recall / F1, plus recall by ORF_type
  * pred_preddepth vs real -- standalone drop-in (predicted shape AND predicted depth, no real data)
Plus, on the matched calls, the agreement of RiboCode's own confidence (-log10 pval_combined) and
frame0 P-site magnitude, and the ORF_type confusion. cas12a env (numpy only). Writes
<out>/dropin_metrics.json + a printed table; the figure is a separate script.

Usage: compare_dropin_calls.py --dropin_dir <dir> --official <Tissue>_collapsed.txt \
    --profiles <pred_profiles.npz> [--out <dir>] [--adj_pval 0.05]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

NONCANON = {"uORF", "Overlap_uORF", "dORF", "Overlap_dORF", "novel", "internal"}
TX2GENE = ("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/"
           "riboseq_signal_model/data/tx2biotype.tsv")


def load_tx2gene(path=TX2GENE):
    t2g = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(header)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            t2g[f[ci["tx_id"]]] = f[ci["gene_id"]]
    return t2g


def _num(f, ci, col):
    try:
        return float(f[ci[col]])
    except (ValueError, KeyError):
        return float("nan")


def load_calls(path, key_mode="genomic", keep_tx=None, keep_genes=None, max_pval=None, min_len=0,
               enrich=None, min_enrichment=0.0):
    """Parse a RiboCode *_collapsed.txt. Returns {key: {type, length, adj, logp, f0}} filtered to
    raw pval_combined <= max_pval and ORF_length >= min_len.

    If enrich is given (a callable enrich(tx, a, e) -> predicted enrichment over uniform = mean
    predicted probability over the ORF x transcript length), calls with enrichment < min_enrichment
    are dropped. This removes the diffuse low-magnitude softmax leak into the 3'UTR (the spurious
    dORF over-calls) while keeping the concentrated CDS/uORF peaks; pass it ONLY for the predicted
    variants (the real/official reference calls are never enrichment-filtered).

    key_mode="genomic" keys by (gene_id, ORF_gstop) -- the genomic ORF locus, INVARIANT to which
    isoform RiboCode picks as the collapse representative. The collapse chooses the longest ORF per
    genomic stop (detectORF.py:440-442, density only breaks exact-start ties), so restricting the
    transcript universe (34k test tx here vs 509k official) reshuffles which isoform carries a call;
    genomic keying neutralizes that artifact. Restricted to keep_genes. key_mode="transcript" keys
    by (transcript_id, ORF_tstart), restricted to keep_tx (isoform-level, universe-sensitive).

    pval_combined is per-ORF and tx-set-independent (adjusted_pval's BH burden depends on how many
    tx were run -- so it must NOT be the threshold). min_len drops the short-ORF regime where the
    caller's periodicity test is fragile and the model over-calls (dORFs); default 90 nt."""
    calls = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(header)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            tx = f[ci["transcript_id"]]
            gene = f[ci["gene_id"]] if "gene_id" in ci else "NA"
            if key_mode == "transcript":
                if keep_tx is not None and tx not in keep_tx:
                    continue
            elif keep_genes is not None and gene not in keep_genes:
                continue
            try:
                tstart = int(f[ci["ORF_tstart"]])
                length = int(f[ci["ORF_length"]])
            except (ValueError, KeyError):
                continue
            if length < min_len:
                continue
            pc = _num(f, ci, "pval_combined")
            if max_pval is not None and not (pc <= max_pval):
                continue
            if enrich is not None and min_enrichment > 0:
                try:
                    enr = enrich(tx, tstart - 1, int(f[ci["ORF_tstop"]]))
                except (ValueError, KeyError):
                    enr = float("nan")
                if not (enr >= min_enrichment):
                    continue
            key = (tx, tstart) if key_mode == "transcript" else (gene, f[ci["ORF_gstop"]])
            logp = -np.log10(pc) if pc == pc and pc > 0 else np.nan
            calls[key] = {
                "type": f[ci["ORF_type"]] if "ORF_type" in ci else "NA", "length": length,
                "adj": _num(f, ci, "adjusted_pval"), "logp": logp,
                "f0": _num(f, ci, "Psites_sum_frame0"),
            }
    return calls


def prf(pred_keys, real_keys):
    inter = len(pred_keys & real_keys)
    p = inter / len(pred_keys) if pred_keys else float("nan")
    r = inter / len(real_keys) if real_keys else float("nan")
    f1 = 2 * p * r / (p + r) if (p == p and r == r and p + r > 0) else float("nan")
    return {"n_pred": len(pred_keys), "n_real": len(real_keys), "n_match": inter,
            "precision": p, "recall": r, "f1": f1}


def compare(pred, real, label):
    pk, rk = set(pred), set(real)
    out = {"label": label, **prf(pk, rk)}
    # recall by real ORF_type (of the real calls of each type, how many pred recovered)
    by_type = {}
    for grp_name, grp in (("annotated", {"annotated"}), ("noncanonical", NONCANON)):
        rk_g = {k for k in rk if real[k]["type"] in grp}
        by_type[grp_name] = {"n_real": len(rk_g),
                             "recall": len(pk & rk_g) / len(rk_g) if rk_g else float("nan")}
    per_type = {}
    for k in rk:
        per_type.setdefault(real[k]["type"], [0, 0])
        per_type[real[k]["type"]][1] += 1
        if k in pk:
            per_type[real[k]["type"]][0] += 1
    out["recall_grouped"] = by_type
    out["recall_per_type"] = {t: {"recovered": a, "n_real": b, "recall": a / b}
                              for t, (a, b) in sorted(per_type.items(), key=lambda kv: -kv[1][1])}
    # agreement on matched calls: RiboCode confidence + frame0 magnitude + type concordance
    matched = pk & rk
    if matched:
        pl = np.array([pred[k]["logp"] for k in matched])
        rl = np.array([real[k]["logp"] for k in matched])
        good = np.isfinite(pl) & np.isfinite(rl)
        out["matched_logp_pearson"] = (float(np.corrcoef(pl[good], rl[good])[0, 1])
                                       if good.sum() > 2 else float("nan"))
        pf = np.log1p([pred[k]["f0"] for k in matched])
        rf = np.log1p([real[k]["f0"] for k in matched])
        gf = np.isfinite(pf) & np.isfinite(rf)
        out["matched_logf0_pearson"] = (float(np.corrcoef(pf[gf], rf[gf])[0, 1])
                                        if gf.sum() > 2 else float("nan"))
        out["matched_type_concordance"] = float(
            np.mean([pred[k]["type"] == real[k]["type"] for k in matched]))
    return out


def build_loader(profiles, key="genomic", tx2gene=TX2GENE, max_pval=0.05, min_len=90,
                 min_enrichment=0.5):
    """Return (lc, keep_tx, keep_genes) for one dump's pred_profiles.npz.

    `lc(path, is_pred=False)` loads a *_collapsed.txt under THIS dump's filtering: restricted to the
    model's test transcripts (or their genes), pval <= max_pval, ORF >= min_len nt, and -- on the
    predicted side only -- mean predicted density over the ORF >= min_enrichment x uniform.

    Factored out of main() so other scorers can reuse the exact same filtering rather than reimplement
    it. Two tables that disagree because one of them quietly skipped the enrichment filter is a
    manuscript-level hazard, and copying the constants across files does not prevent it.
    """
    npz = np.load(profiles, allow_pickle=False)
    keep_tx = set(npz["tx_ids"].tolist())
    _lengths = npz["lengths"]
    _pflat = npz["pred_flat"]
    _off = np.concatenate([[0], np.cumsum(_lengths)])
    _idx = {str(t): i for i, t in enumerate(npz["tx_ids"])}

    def enrich(tx, a, e):
        i = _idx.get(str(tx))
        if i is None:
            return float("nan")
        p = _pflat[_off[i]:_off[i + 1]]
        if a < 0 or e > len(p) or e <= a:
            return float("nan")
        return float(p[a:e].mean() * len(p))

    keep_genes = None
    if key == "genomic":
        t2g = load_tx2gene(tx2gene)
        keep_genes = {t2g[t] for t in keep_tx if t in t2g}
        if keep_tx and not keep_genes:
            raise SystemExit(
                f"ERROR: 0 of {len(keep_tx):,} test tx matched the tx2gene map {tx2gene}. It is "
                "likely for the wrong assembly (e.g. the human map against mouse tx). Pass the "
                "matching tx2gene, or use key='transcript'.")

    def lc(path, is_pred=False):
        return load_calls(path, key, keep_tx, keep_genes, max_pval, min_len,
                          enrich if is_pred else None, min_enrichment)

    return lc, keep_tx, keep_genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dropin_dir", required=True)
    # Optional: not every held-out dataset has a call set from RiboCode run over the FULL observed
    # data outside the pack. Wang and GSE243134 liver do; Janich does not. The two *_vs_official rows
    # are harness validation, not the primary metric -- pred_* vs real is self-contained -- so omitting
    # official drops those two rows rather than blocking the comparison.
    ap.add_argument("--official", default=None)
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--max_pval", type=float, default=0.05,
                    help="keep calls with raw pval_combined <= this in every set (tx-set-indep; "
                         "collapsed.txt is already filtered here, so 0.05 = use membership)")
    ap.add_argument("--min_len", type=int, default=90,
                    help="drop ORFs shorter than this (nt) in every set; the short-ORF regime is "
                         "where the predicted density over-calls (spurious dORFs). Default 90.")
    ap.add_argument("--key", choices=["genomic", "transcript"], default="genomic",
                    help="match calls by (gene_id, ORF_gstop) [genomic: invariant to the collapse "
                         "representative choice] or (transcript_id, ORF_tstart) [transcript: "
                         "isoform-level, universe-sensitive]. Default genomic.")
    ap.add_argument("--min_enrichment", type=float, default=0.5,
                    help="drop PREDICTED calls whose mean predicted density over the ORF is below "
                         "this multiple of uniform (mean prob x L). Removes the diffuse 3'UTR leak "
                         "(spurious dORFs). Pred variants only; 0 = off. Default 0.5.")
    ap.add_argument("--tx2gene", default=TX2GENE,
                    help="transcript_id -> gene_id map (tx2biotype.tsv schema) used by --key genomic. "
                         "Default is the human v49 table; pass the matching assembly's table for a "
                         "non-human held-out (e.g. data/mouse_tx2biotype.tsv for the mouse vM38 "
                         "held-out), else every tx misses the map and all calls are filtered out.")
    args = ap.parse_args()
    d = Path(args.dropin_dir)
    out = Path(args.out) if args.out else d
    out.mkdir(parents=True, exist_ok=True)

    lc, keep_tx, keep_genes = build_loader(
        args.profiles, key=args.key, tx2gene=args.tx2gene, max_pval=args.max_pval,
        min_len=args.min_len, min_enrichment=args.min_enrichment)
    if args.key == "genomic":
        print(f"test tx: {len(keep_tx):,}; test genes: {len(keep_genes):,} "
              f"(tx2gene={args.tx2gene})", file=sys.stderr)
    else:
        print(f"test transcripts (model inference set): {len(keep_tx):,}", file=sys.stderr)
    real = lc(d / "real_collapsed.txt")
    pred_obs = lc(d / "pred_obsdepth_collapsed.txt", is_pred=True)
    pred_pred = lc(d / "pred_preddepth_collapsed.txt", is_pred=True)
    official = lc(args.official) if args.official else None
    off_str = (f"official={len(official):,}" if official is not None
               else "official=NONE (--official not given; the two *_vs_official rows are skipped)")
    print(f"calls (key={args.key}, pval<={args.max_pval}, len>={args.min_len}nt, "
          f"pred enrichment>={args.min_enrichment}): real={len(real):,} "
          f"pred_obsdepth={len(pred_obs):,} pred_preddepth={len(pred_pred):,} "
          f"{off_str}", file=sys.stderr)

    res = {
        "test_tx": len(keep_tx), "max_pval": args.max_pval, "min_len": args.min_len,
        "key": args.key, "min_enrichment": args.min_enrichment,
        "pred_obsdepth_vs_real": compare(pred_obs, real, "pred_obsdepth_vs_real"),
        "pred_preddepth_vs_real": compare(pred_pred, real, "pred_preddepth_vs_real"),
    }
    if official is not None:
        res["real_vs_official"] = compare(real, official, "real_vs_official")
        res["pred_obsdepth_vs_official"] = compare(pred_obs, official, "pred_obsdepth_vs_official")
    (out / "dropin_metrics.json").write_text(json.dumps(res, indent=2))

    print("\n=== RiboCode drop-in: predicted-density ORF calls vs real-density calls "
          "(held-out Hepatocytes) ===")
    print(f"{'comparison':<28} {'nPred':>7} {'nReal':>7} {'match':>7} "
          f"{'prec':>6} {'recall':>7} {'F1':>6} {'logpR':>6} {'typeConc':>8}")
    for key in ("real_vs_official", "pred_obsdepth_vs_real", "pred_preddepth_vs_real",
                "pred_obsdepth_vs_official"):
        c = res.get(key)
        if c is None:      # no --official: those two rows were not computed
            continue
        print(f"{key:<28} {c['n_pred']:>7,} {c['n_real']:>7,} {c['n_match']:>7,} "
              f"{c['precision']:>6.3f} {c['recall']:>7.3f} {c['f1']:>6.3f} "
              f"{c.get('matched_logp_pearson', float('nan')):>6.2f} "
              f"{c.get('matched_type_concordance', float('nan')):>8.2f}")
    print("\n-- pred_obsdepth vs real: recall by ORF type --")
    for t, v in res["pred_obsdepth_vs_real"]["recall_per_type"].items():
        print(f"  {t:<16} {v['recovered']:>6,}/{v['n_real']:<6,} recall {v['recall']:.3f}")
    print(f"\nwrote {out / 'dropin_metrics.json'}", file=sys.stderr)


if __name__ == "__main__":
    main()
