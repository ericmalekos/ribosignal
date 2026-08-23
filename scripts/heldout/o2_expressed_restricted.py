#!/usr/bin/env python3
"""O2 re-scored on EXPRESSED transcripts only.

Why: both O2 dumps used `--tx_list`, which bypasses the `min_signal` gate that every other eval
applies, so the model predicted across the entire 221,835-tx Wang universe (essentially all of
vM38). The reference Ribo-seq call sets only ever span ~16.5k transcripts, so ~85% of the
transcripts carrying model calls are ones the reference could never match. Those count as false
positives even where the reference experiment had no coverage at all, which deflates precision and
hits the non-canonical class hardest.

This restricts ALL FOUR call sets (both observed, both model) to a shared expressed-transcript
space and recomputes the ceiling and both arms, so the comparison is like-for-like.

Expression is derived from the packs' own per-nt coverage (summed depth / length, normalised to
1e6), because Wang has no salmon quant: its RNA fastqs were processed into coverage and then
deleted. The proxy is validated against Janich's real salmon quant, recovering 95.4% of the true
TPM>=1 set, but it is more permissive overall, so two thresholds are reported: proxy>=1, and a
salmon-calibrated cut reproducing Janich's real 11,510-tx count.

cas12a env (numpy). Run: python o2_expressed_restricted.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
ECH = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/biotype_probe/expression_context_human")

WANG_OBS = NEW / "data/heldout_psites/mouse_wang_liver/mouse_wang_liver_collapsed.txt"
JANICH_OBS = NEW / "results/merged_liver_ribocode/merged_liver_collapsed.txt"
WANG_MODEL = NEW / "results/o2_liver/wang_nokozak/dropin/pred_preddepth_collapsed.txt"
JANICH_MODEL = NEW / "results/o2_liver/janich_wanguni/dropin/pred_preddepth_collapsed.txt"
WANG_PACK = NEW / "data/packed_heldout_mouse_wang_liver"
JANICH_PACK = NEW / "data/packed_heldout_mouse_janich_liver_wanguni"
SALMON = NEW / "data/tpm/salmon_quant_janich_liver/janich_merged_quant.sf"
SALMON_WANG = NEW / "data/tpm/salmon_quant_wang_liver/wang_merged_quant.sf"
OUT = NEW / "results/o2_liver/o2_expressed_restricted.md"


def salmon_tpm(path):
    """Real salmon TPM, keyed by bare versioned tx id (quant.sf carries full GENCODE headers)."""
    d = {}
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.split("\t")
            d[p[0].split("|")[0]] = float(p[3])
    return d


def coord_key(orf_id):
    p = orf_id.rsplit("_", 3)
    return "_".join(p[-3:]) if len(p) >= 4 else orf_id


def load(path):
    """-> list of (coord_key, ORF_type, transcript_id)."""
    out = []
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        oi, ti, ty = h.index("ORF_ID"), h.index("transcript_id"), h.index("ORF_type")
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            out.append((coord_key(p[oi]), p[ty], p[ti]))
    return out


def cov_tpm(pack):
    off = np.load(pack / "offsets.npy")
    ln = np.load(pack / "lengths.npy").astype(np.float64)
    cov = np.load(pack / "coverage.npy", mmap_mode="r")
    tot = np.add.reduceat(np.asarray(cov, dtype=np.float64), off[:-1])
    rate = np.divide(tot, ln, out=np.zeros_like(tot), where=ln > 0)
    tpm = rate / rate.sum() * 1e6
    tx = [x.strip() for x in open(pack / "tx_order.txt")]
    return dict(zip(tx, tpm, strict=False))


def f1(pred, ref):
    i = len(pred & ref)
    p = i / len(pred) if pred else 0.0
    r = i / len(ref) if ref else 0.0
    return p, r, (2 * p * r / (p + r) if (p + r) else 0.0)


def sel(rows, cls, keep_tx):
    """coord-key set for a class, restricted to keep_tx (None = no restriction)."""
    out = set()
    for k, t, tx in rows:
        if keep_tx is not None and tx not in keep_tx:
            continue
        if cls == "all" or (cls == "cds" and t == "annotated") or \
           (cls == "noncan" and t != "annotated"):
            out.add(k)
    return out


def main():
    wo, jo = load(WANG_OBS), load(JANICH_OBS)
    wm, jm = load(WANG_MODEL), load(JANICH_MODEL)
    wt, jt = cov_tpm(WANG_PACK), cov_tpm(JANICH_PACK)

    real = {}
    with open(SALMON) as f:
        f.readline()
        for line in f:
            p = line.split("\t")
            real[p[0].split("|")[0]] = float(p[3])       # strip GENCODE pipe header
    n_real = sum(1 for v in real.values() if v >= 1)
    cal = sorted((v for v in jt.values() if v > 0), reverse=True)
    cal_thr = cal[n_real - 1] if n_real <= len(cal) else 1.0

    L = ["# O2 re-scored on EXPRESSED transcripts only\n",
         "Both O2 dumps used `--tx_list`, bypassing the `min_signal` gate every other eval applies, so the "
         "model predicted over the whole 221,835-tx annotation while the reference Ribo-seq sets span only "
         "~16.5k transcripts. Restricting all four call sets to a shared expressed space makes the "
         "comparison like-for-like.\n",
         f"Expression is a coverage-derived proxy (Wang has no salmon quant). Validated on Janich: recovers "
         f"95.4% of its real TPM>=1 set. `calibrated` = the proxy cut reproducing Janich's real "
         f"{n_real:,}-transcript count (proxy TPM >= {cal_thr:.2f}).\n"]

    # REAL salmon TPM for both datasets (Wang quantified 2026-07-31; Janich rebuilt 2026-07-30)
    real_w = salmon_tpm(SALMON_WANG) if SALMON_WANG.exists() else None
    keep_real = None
    if real_w is not None:
        keep_real = ({t for t, v in real_w.items() if v >= 1}
                     & {t for t, v in real.items() if v >= 1})
        L.append(f"\n**Real salmon TPM>=1**: Wang {sum(1 for v in real_w.values() if v>=1):,}, "
                 f"Janich {n_real:,}, shared **{len(keep_real):,}**. This is the primary result; the "
                 f"coverage-proxy rows below are kept for comparison.\n")

    rows = []
    variants = [("unrestricted", None)]
    if keep_real is not None:
        variants.append(("REAL salmon TPM >= 1 (both)", "REAL"))
    variants += [("proxy TPM >= 1", 1.0), (f"proxy calibrated (>= {cal_thr:.2f})", cal_thr)]
    for lbl, thr in variants:
        if thr is None:
            keep = None
            ntx = "all"
        elif thr == "REAL":
            keep = keep_real
            ntx = f"{len(keep):,}"
        else:
            keep = {t for t, v in wt.items() if v >= thr} & {t for t, v in jt.items() if v >= thr}
            ntx = f"{len(keep):,}"
        L.append(f"\n## {lbl}  (shared expressed transcripts: {ntx})\n")
        L.append("| class | ceiling F1(Wang,Janich) | Wang_model vs Janich_obs | % of ceiling | "
                 "Janich_model vs Wang_obs | % of ceiling |")
        L.append("|---|---:|---:|---:|---:|---:|")
        for cls, name in (("all", "ALL"), ("cds", "CDS"), ("noncan", "non-canonical")):
            w, j = sel(wo, cls, keep), sel(jo, cls, keep)
            _, _, ceil = f1(w, j)
            _, _, a = f1(sel(wm, cls, keep), j)      # Wang model vs the OTHER real set
            _, _, b = f1(sel(jm, cls, keep), w)      # Janich model vs the OTHER real set
            L.append(f"| {name} | {ceil:.3f} | {a:.3f} | {a/ceil*100:.0f}% | {b:.3f} | {b/ceil*100:.0f}% |"
                     if ceil else f"| {name} | 0 | {a:.3f} | -- | {b:.3f} | -- |")
            rows.append((lbl, name, ceil, a, b))
        nm = len(sel(wm, "all", keep))
        nr = len(sel(jo, "all", keep))
        L.append(f"\nWang_model calls {nm:,} vs Janich_obs {nr:,} "
                 f"(over-call {nm/nr:.1f}x)" if nr else "")
    md = "\n".join(L) + "\n"
    OUT.write_text(md)
    print(md)
    json.dump([{"restriction": r[0], "class": r[1], "ceiling": r[2],
                "wang_model_vs_janich": r[3], "janich_model_vs_wang": r[4]} for r in rows],
              open(str(OUT).replace(".md", ".json"), "w"), indent=2)
    print(f"wrote {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
