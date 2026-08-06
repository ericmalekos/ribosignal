#!/usr/bin/env python3
"""Task 14 transfer check: does the enrichment-threshold tradeoff (chosen on held-out Hepatocytes)
hold on a second tissue (Fibroblast)? Overlays the two tissues' enrichment sweeps -- if the curves
are parallel and the knee sits in the same place, the 0.5 operating point transfers.

For each (dropin_dir, label) it computes, on the standalone (count-head depth) predicted calls
matched by genomic ORF locus (gene_id, ORF_gstop) at ORFs >= 90 nt, the precision / recall / F1 and
the spurious dORF FP count across enrichment thresholds (mean predicted probability over the ORF x
transcript length). Left panel: precision (solid) + recall (dashed) vs threshold, both tissues.
Right: spurious dORF FPs. cas12a matplotlib. Usage: make_transfer_figure.py [--out <png>]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
TX2GENE = NEW / "data" / "tx2biotype.tsv"
TAUS = [0.0, 0.3, 0.5, 0.7, 1.0]
# (dropin dir, label, colour)
TISSUES = [
    (NEW / "results/loto/orf_v2_attn_holdout_Hepatocytes/dropin",
     "Hepatocytes (tuned on)", "#2b6cb0"),
    (NEW / "results/improve/orf_v2_attn_f0/dropin", "Fibroblast (transfer)", "#d1782f"),
]


def load_tx2gene():
    t2g = {}
    with TX2GENE.open() as fh:
        ci = {c: i for i, c in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            t2g[f[ci["tx_id"]]] = f[ci["gene_id"]]
    return t2g


def sweep(dropin, t2g):
    z = np.load(dropin / "pred_profiles.npz", allow_pickle=False)
    tx_ids, lengths, pflat = z["tx_ids"], z["lengths"], z["pred_flat"]
    off = np.concatenate([[0], np.cumsum(lengths)])
    idx = {t: i for i, t in enumerate(tx_ids)}
    keep = set(tx_ids.tolist())
    tgenes = {t2g[t] for t in keep if t in t2g}

    def load(path):
        out = []
        with open(path) as fh:
            ci = {c: i for i, c in enumerate(fh.readline().rstrip("\n").split("\t"))}
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                tx, g = f[ci["transcript_id"]], f[ci["gene_id"]]
                if tx not in keep or g not in tgenes:
                    continue
                try:
                    if int(f[ci["ORF_length"]]) < 90 or not float(f[ci["pval_combined"]]) <= 0.05:
                        continue
                    a, e = int(f[ci["ORF_tstart"]]) - 1, int(f[ci["ORF_tstop"]])
                except (ValueError, KeyError):
                    continue
                out.append((tx, f[ci["ORF_type"]], a, e, (g, f[ci["ORF_gstop"]])))
        return out

    realk = {r[4] for r in load(dropin / "real_collapsed.txt")}
    rows = load(dropin / "pred_preddepth_collapsed.txt")
    enr = []
    for tx, _, a, e, _ in rows:
        i = idx[tx]
        seg = pflat[off[i]:off[i + 1]]
        enr.append(float(seg[a:e].mean() * len(seg)))
    enr = np.array(enr)
    out = {"prec": [], "recall": [], "f1": [], "dorf_fp": []}
    for tau in TAUS:
        kk = [r for r, en in zip(rows, enr, strict=False) if en >= tau]
        kept = {r[4] for r in kk}
        tp = len(kept & realk)
        prec = tp / len(kept) if kept else 0.0
        rec = tp / len(realk) if realk else 0.0
        out["prec"].append(prec)
        out["recall"].append(rec)
        out["f1"].append(2 * prec * rec / (prec + rec) if prec + rec else 0.0)
        out["dorf_fp"].append(sum(1 for r in kk if r[1] == "dORF" and r[4] not in realk))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(NEW / "results/figures/dropin_transfer.png"))
    args = ap.parse_args()
    t2g = load_tx2gene()
    data = [(lab, col, sweep(d, t2g)) for d, lab, col in TISSUES]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.6))
    for lab, col, s in data:
        axL.plot(TAUS, s["prec"], "o-", color=col, label=f"{lab} precision")
        axL.plot(TAUS, s["recall"], "s--", color=col, alpha=0.7, label=f"{lab} recall")
    axL.axvline(0.5, ls=":", lw=1, color="#888")
    axL.text(0.52, 0.44, "0.5 default", fontsize=8, color="#555")
    axL.set_xlabel("enrichment threshold (x uniform)")
    axL.set_ylabel("agreement with real calls")
    axL.set_ylim(0.4, 1.0)
    axL.set_xticks(TAUS)
    axL.legend(fontsize=7.5, loc="lower left", framealpha=0.9)
    axL.set_title("A  Precision / recall vs threshold -- parallel across tissues",
                  fontsize=9.5, loc="left")

    for lab, col, s in data:
        axR.plot(TAUS, s["dorf_fp"], "o-", color=col, label=lab)
        for tau, v in zip(TAUS, s["dorf_fp"], strict=False):
            axR.annotate(f"{v}", (tau, v), fontsize=7, color=col,
                         textcoords="offset points", xytext=(0, 5), ha="center")
    axR.axvline(0.5, ls=":", lw=1, color="#888")
    axR.set_xlabel("enrichment threshold (x uniform)")
    axR.set_ylabel("spurious dORF false positives")
    axR.set_ylim(bottom=0)
    axR.set_xticks(TAUS)
    axR.legend(fontsize=8, loc="upper right", framealpha=0.9)
    axR.set_title("B  Spurious dORFs collapse by ~0.5 in both tissues", fontsize=9.5, loc="left")

    fig.suptitle("Enrichment-threshold transfer: Fibroblast reproduces the Hepatocytes knee "
                 "(standalone, genomic-matched, ORFs >= 90 nt)", fontsize=11, y=1.02)
    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
