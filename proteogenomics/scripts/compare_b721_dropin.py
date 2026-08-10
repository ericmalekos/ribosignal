#!/usr/bin/env python3
"""B721.221 ground-truth check: do the model's PREDICTED ORF calls recover the ones RiboCode makes
on REAL Ribo-seq from the same cell line?

This is the first model-vs-measured-translation comparison in the project. Every prior
proteogenomics result scored the model against enumerated or coding-potential NULLS, which tests
"is the model better than naive enumeration", not "does the model find what is actually
translated". B721.221 has both arms in the same parental line: Sarkizova RNA-seq (drives the
prediction) and Ouspenskaia Ribo-seq (327 M unique footprints, gives the measured calls).

TRANSCRIPT SPACE (standing rule): the two arms are NOT called over the same transcript set. The
measured call set comes from RiboCode over the full prepared v49 annotation; the predicted arms run
only over the model's universe restricted to B721-expressed transcripts. Scoring is therefore
GENOMIC-keyed (gene_id, ORF_gstop) and restricted to the genes present in the predicted universe --
otherwise the model is charged for every gene it was never given. n_ref is reported for both the
raw and the restricted reference so the restriction is auditable.

Usage:
  compare_b721_dropin.py --measured <B721_collapsed.txt> --profiles <pred_profiles.npz> \
      --arms mamba4/theta1=<dir> ... --out <dir> [--min_len 90] [--pval 0.05]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from compare_dropin_calls import NONCANON, load_calls, prf  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--measured", required=True, help="RiboCode *_collapsed.txt from real Ribo-seq")
    ap.add_argument("--profiles", required=True, help="pred_profiles.npz (defines the model universe)")
    ap.add_argument("--arms", nargs="+", required=True,
                    help="label=<dropin dir> (uses pred_preddepth_collapsed.txt) or "
                         "label=<path to a *_collapsed.txt>. The file form lets the split-half "
                         "CEILING be scored by this same code path -- keying, thresholds and the "
                         "gene restriction must be identical to the drop-in or the ceiling is not "
                         "comparable to what it is meant to bound.")
    ap.add_argument("--tx2biotype", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min_len", type=int, default=90)
    ap.add_argument("--pval", type=float, default=0.05)
    a = ap.parse_args()

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    # The model universe: transcripts the model actually predicted a profile for.
    with np.load(a.profiles, allow_pickle=True) as z:
        keep_tx = {str(t).split()[0] for t in z["tx_ids"]} if "tx_ids" in z else set(z.files)
    t2g = {}
    src = a.tx2biotype or (Path(__file__).resolve().parents[2] / "data" / "tx2biotype.tsv")
    with open(src) as fh:
        ci = {c: i for i, c in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            t2g[f[ci["tx_id"]]] = f[ci["gene_id"]]
    keep_genes = {t2g[t] for t in keep_tx if t in t2g}
    print(f"model universe: {len(keep_tx)} tx -> {len(keep_genes)} genes")

    kw = dict(key_mode="genomic", keep_genes=keep_genes, max_pval=a.pval, min_len=a.min_len)
    meas_all = load_calls(a.measured, key_mode="genomic", max_pval=a.pval, min_len=a.min_len)
    meas = load_calls(a.measured, **kw)
    print(f"measured calls: {len(meas_all)} genome-wide -> {len(meas)} in the model's gene space")

    res = {"n_ref_genomewide": len(meas_all), "n_ref_restricted": len(meas),
           "n_tx_universe": len(keep_tx), "n_gene_universe": len(keep_genes),
           "min_len": a.min_len, "pval": a.pval, "arms": {}}
    rk = set(meas)
    rows = []
    for spec in a.arms:
        label, _, d = spec.partition("=")
        p = Path(d)
        f = p if p.is_file() else p / "pred_preddepth_collapsed.txt"
        if not f.exists():
            print(f"  MISSING {label}: {f}", file=sys.stderr)
            continue
        pred = load_calls(f, **kw)
        pk = set(pred)
        m = prf(pk, rk)
        # split the measured side: canonical CDS vs everything non-canonical
        for grp, want in (("annotated", {"annotated"}), ("noncanonical", NONCANON)):
            rg = {k for k in rk if meas[k]["type"] in want}
            pg = {k for k in pk if pred[k]["type"] in want}
            m[f"recall_{grp}"] = len(pk & rg) / len(rg) if rg else float("nan")
            m[f"n_real_{grp}"] = len(rg)
            m[f"n_pred_{grp}"] = len(pg)
            m[f"precision_{grp}"] = len(pg & rk) / len(pg) if pg else float("nan")
        res["arms"][label] = m
        rows.append((label, m))

    hdr = (f"{'arm':<16}{'n_pred':>8}{'n_meas':>8}{'match':>8}{'prec':>7}{'rec':>7}{'F1':>7}"
           f"{'rec_ann':>9}{'rec_nc':>8}{'prec_nc':>9}")
    print("\n" + hdr)
    print("-" * len(hdr))
    for label, m in rows:
        print(f"{label:<16}{m['n_pred']:>8}{m['n_real']:>8}{m['n_match']:>8}"
              f"{m['precision']:>7.3f}{m['recall']:>7.3f}{m['f1']:>7.3f}"
              f"{m['recall_annotated']:>9.3f}{m['recall_noncanonical']:>8.3f}"
              f"{m['precision_noncanonical']:>9.3f}")
    (out / "b721_dropin_metrics.json").write_text(json.dumps(res, indent=2))
    print(f"\nwrote {out / 'b721_dropin_metrics.json'}")


if __name__ == "__main__":
    main()
