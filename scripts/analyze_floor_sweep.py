#!/usr/bin/env python3
"""Task 14 FP-reduction sweep: how a per-position density floor (zero predicted positions below
floor_mult x uniform, ribocode_dropin.py --floor_mult) trades precision for recall, and how far it
knocks down the spurious 3'UTR dORF over-calls.

For each pred variant x floor level it scores the drop-in ORF calls against the real-density calls
on the genomic ORF locus (gene_id, ORF_gstop), restricted to the test-tx genes, ORFs >= 90 nt,
pval_combined <= 0.05 -- the same contract as compare_dropin_calls.py --key genomic. Reports
precision / recall / F1, the spurious dORF FP count, and per-class recall (annotated / uORF / dORF),
and writes a tradeoff figure. cas12a matplotlib.

Usage: analyze_floor_sweep.py [--dropin_dir <dir>] [--profiles <npz>] [--out <png>]
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
VARIANTS = ["pred_obsdepth", "pred_preddepth"]
FLOORS = [0.0, 0.5, 1.0]


def load_tx2gene():
    t2g = {}
    with TX2GENE.open() as fh:
        ci = {c: i for i, c in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            t2g[f[ci["tx_id"]]] = f[ci["gene_id"]]
    return t2g


def load(path, keep_genes, min_len=90, max_pval=0.05):
    """{(gene_id, ORF_gstop): ORF_type} over test-gene calls passing the filters."""
    calls = {}
    with open(path) as fh:
        ci = {c: i for i, c in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            g = f[ci["gene_id"]]
            if g not in keep_genes:
                continue
            try:
                if int(f[ci["ORF_length"]]) < min_len:
                    continue
                if not float(f[ci["pval_combined"]]) <= max_pval:
                    continue
            except (ValueError, KeyError):
                continue
            calls[(g, f[ci["ORF_gstop"]])] = f[ci["ORF_type"]]
    return calls


def fname(variant, floor):
    return f"{variant}_collapsed.txt" if floor <= 0 else f"{variant}_floor{floor:g}_collapsed.txt"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dropin_dir",
                    default=str(NEW / "results/loto/orf_v2_attn_holdout_Hepatocytes/dropin"))
    ap.add_argument("--profiles", default=None)
    ap.add_argument("--out", default=str(NEW / "results/figures/dropin_floor_sweep.png"))
    args = ap.parse_args()
    d = Path(args.dropin_dir)
    profiles = args.profiles or str(d / "pred_profiles.npz")

    keep_tx = set(np.load(profiles, allow_pickle=False)["tx_ids"].tolist())
    t2g = load_tx2gene()
    keep_genes = {t2g[t] for t in keep_tx if t in t2g}
    real = load(d / "real_collapsed.txt", keep_genes)
    realk = set(real)
    real_by_type = {}
    for k, t in real.items():
        real_by_type.setdefault(t, set()).add(k)

    print(f"test genes {len(keep_genes):,}; real calls {len(realk):,} "
          f"(dORF {len(real_by_type.get('dORF', ())):,})\n")
    hdr = (f"{'variant':<15}{'floor':>6}{'nPred':>8}{'prec':>7}{'recall':>8}{'F1':>7}"
           f"{'dORF_FP':>9}{'rec_ann':>8}{'rec_uORF':>9}{'rec_dORF':>9}")
    print(hdr)
    results = {}
    for v in VARIANTS:
        results[v] = []
        for fm in FLOORS:
            fp = d / fname(v, fm)
            if not fp.exists():
                print(f"{v:<15}{fm:>6}  MISSING {fp.name}")
                results[v].append(None)
                continue
            pred = load(fp, keep_genes)
            pk = set(pred)
            m = len(pk & realk)
            prec = m / len(pk) if pk else float("nan")
            rec = m / len(realk) if realk else float("nan")
            f1 = 2 * prec * rec / (prec + rec) if prec + rec else float("nan")
            dorf_fp = sum(1 for k, t in pred.items() if t == "dORF" and k not in realk)

            def rct(t, pk=pk):
                rk = real_by_type.get(t, set())
                return len(pk & rk) / len(rk) if rk else float("nan")
            rec_ann, rec_uorf, rec_dorf = rct("annotated"), rct("uORF"), rct("dORF")
            results[v].append({"floor": fm, "prec": prec, "recall": rec, "f1": f1,
                               "dorf_fp": dorf_fp, "n": len(pk)})
            print(f"{v:<15}{fm:>6}{len(pk):>8,}{prec:>7.3f}{rec:>8.3f}{f1:>7.3f}"
                  f"{dorf_fp:>9,}{rec_ann:>8.3f}{rec_uorf:>9.3f}{rec_dorf:>9.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    for ax, v in zip(axes, VARIANTS, strict=False):
        pts = [r for r in results[v] if r]
        if not pts:
            continue
        x = [r["floor"] for r in pts]
        for key, col, lab in [("prec", "#2b6cb0", "precision"), ("recall", "#d1782f", "recall"),
                              ("f1", "#2ca02c", "F1")]:
            ax.plot(x, [r[key] for r in pts], "o-", color=col, label=lab)
            for r in pts:
                ax.annotate(f"{r[key]:.2f}", (r["floor"], r[key]), fontsize=6.5,
                            textcoords="offset points", xytext=(0, 4), ha="center", color=col)
        ax2 = ax.twinx()
        ax2.plot(x, [r["dorf_fp"] for r in pts], "s--", color="#888", label="dORF FPs")
        ax2.set_ylabel("spurious dORF FPs", color="#888")
        ax2.set_ylim(bottom=0)
        for r in pts:
            ax2.annotate(f"{r['dorf_fp']}", (r["floor"], r["dorf_fp"]), fontsize=6.5, color="#666",
                         textcoords="offset points", xytext=(0, -10), ha="center")
        ax.set_xlabel("density floor (x uniform)")
        ax.set_ylabel("agreement with real calls")
        ax.set_ylim(0, 1.0)
        ax.set_xticks(FLOORS)
        ax.set_title(f"{v}", fontsize=10, loc="left")
        ax.legend(fontsize=8, loc="lower left")
    fig.suptitle("Per-position density floor vs the spurious-dORF over-call "
                 "(Hepatocytes, genomic-locus matched, ORFs >= 90 nt)", fontsize=11, y=1.02)
    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
