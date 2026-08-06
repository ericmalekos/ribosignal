#!/usr/bin/env python3
"""Check 2/3: build the CDS-anchored operating curve from a depth sweep, and read off the
non-canonical ORF yield at each CDS operating point (the user-facing stringency dial).

Reads:
  <run_dir>/dropin/real_collapsed.txt                         observed calls (the reference/truth)
  <run_dir>/dropin_sweep/theta_*/pred_preddepth_collapsed.txt predicted calls at each effective depth
Writes:
  <run_dir>/dropin_sweep/cds_anchored_curve.tsv   per-theta metrics
  <run_dir>/dropin_sweep/dial_summary.txt         non-canonical yield at CDS recall {0.8,0.9,1.0}

CDS (ORF_type == annotated) is the trusted anchor: at each theta we measure CDS recall/precision,
then report how many non-canonical ORFs (uORF/novel/dORF) come along and their precision vs the
observed calls. Lower CDS recall = stricter = fewer, higher-confidence novel ORFs.
"""
import glob
import os
import sys
from collections import Counter

run_dir = sys.argv[1]
sweep_name = sys.argv[2] if len(sys.argv) > 2 else "dropin_sweep"   # allow poisson/minaa sweep dirs
sweep = f"{run_dir}/{sweep_name}"
NONCAN = ["uORF", "novel", "dORF", "Overlap_uORF", "Overlap_dORF", "internal"]


def load(path):
    d = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            d[p[ci["ORF_ID"]]] = p[ci["ORF_type"]]
    return d


real = load(f"{run_dir}/dropin/real_collapsed.txt")
real_cds = {i for i, t in real.items() if t == "annotated"}
real_nc = {t: {i for i, tt in real.items() if tt == t} for t in NONCAN}


def metrics(pred):
    pred_cds = {i for i, t in pred.items() if t == "annotated"}
    inter = len(pred_cds & real_cds)
    cds_rec = inter / len(real_cds) if real_cds else 0.0
    cds_pre = inter / len(pred_cds) if pred_cds else 0.0
    row = dict(n_pred=len(pred), cds_recall=cds_rec, cds_precision=cds_pre,
               calls_over_cds=len(pred) / len(real_cds) if real_cds else 0.0)
    for t in NONCAN:
        pt = {i for i, tt in pred.items() if tt == t}
        row[f"{t}_n"] = len(pt)
        row[f"{t}_prec"] = (len(pt & real_nc[t]) / len(pt)) if pt else 0.0
    return row


rows = []
for d in sorted(glob.glob(f"{sweep}/theta_*"), key=lambda p: float(p.split("theta_")[1])):
    f = f"{d}/pred_preddepth_collapsed.txt"
    if not os.path.exists(f):
        continue
    theta = float(d.split("theta_")[1])
    r = metrics(load(f))
    r["theta"] = theta
    rows.append(r)

cols = ["theta", "n_pred", "cds_recall", "cds_precision", "calls_over_cds",
        "uORF_n", "uORF_prec", "novel_n", "novel_prec", "dORF_n", "dORF_prec"]
with open(f"{sweep}/cds_anchored_curve.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

# dial: nearest theta to each CDS-recall target
L = [f"CDS-anchored dial ({run_dir.split('/')[-1]}); real CDS={len(real_cds)}\n",
     f"{'CDS_recall':>10s}{'theta':>7s}{'CDS_prec':>9s}{'calls/CDS':>10s}"
     f"{'uORF(prec)':>14s}{'novel(prec)':>14s}{'dORF(prec)':>13s}"]
for r in rows:
    u = f"{r['uORF_n']}({r['uORF_prec']:.2f})"
    nv = f"{r['novel_n']}({r['novel_prec']:.2f})"
    dd = f"{r['dORF_n']}({r['dORF_prec']:.2f})"
    L.append(f"{r['cds_recall']:>10.3f}{r['theta']:>7.2f}{r['cds_precision']:>9.3f}"
             f"{r['calls_over_cds']:>10.2f}{u:>14s}{nv:>14s}{dd:>13s}")
report = "\n".join(L) + "\n\nRead: lower CDS_recall (stricter) -> fewer novel ORFs at higher precision.\n"
open(f"{sweep}/dial_summary.txt", "w").write(report)
print(report)
print(f"wrote {sweep}/cds_anchored_curve.tsv + dial_summary.txt")
