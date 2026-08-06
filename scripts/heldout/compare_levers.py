#!/usr/bin/env python3
"""Check 4: compare calibration levers at MATCHED CDS recall. Each lever has a CDS-anchored curve
(cds_anchored_curve.tsv). At CDS-recall targets we interpolate the non-canonical precision; the
lever giving the highest novel/uORF precision at the same CDS recall wins.

Levers:
  depth    dropin_sweep         (a) scale predicted depth (--pred_scale)
  poisson  dropin_sweep_poisson (b) same, Poisson-sampled (noise injection)
  pval     dropin_sweep_pval    (d) fixed depth, tighten RiboCode p-value
  obsdepth single point         (c) isotonic-count ceiling ~= pred_obsdepth (shape at REAL per-tx depth)

Usage: compare_levers.py <run_dir>
"""
import sys

run = sys.argv[1]
TARGETS = [0.80, 0.85, 0.90]
LEVERS = [("depth", "dropin_sweep"), ("poisson", "dropin_sweep_poisson"), ("pval", "dropin_sweep_pval")]


def read_curve(path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            rows.append({k: float(p[ci[k]]) for k in ("cds_recall", "uORF_prec", "novel_prec", "novel_n")})
    return sorted(rows, key=lambda r: r["cds_recall"])


def interp(rows, target, key):
    """linear-interpolate key at cds_recall=target."""
    below = [r for r in rows if r["cds_recall"] <= target]
    above = [r for r in rows if r["cds_recall"] >= target]
    if not below or not above:
        return None
    lo, hi = below[-1], above[0]
    if hi["cds_recall"] == lo["cds_recall"]:
        return lo[key]
    w = (target - lo["cds_recall"]) / (hi["cds_recall"] - lo["cds_recall"])
    return lo[key] + w * (hi[key] - lo[key])


def load_calls(path):
    d = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            d[p[ci["ORF_ID"]]] = p[ci["ORF_type"]]
    return d


curves = {}
for name, sub in LEVERS:
    try:
        curves[name] = read_curve(f"{run}/{sub}/cds_anchored_curve.tsv")
    except FileNotFoundError:
        curves[name] = None

# isotonic ceiling: pred_obsdepth single operating point (shape at REAL per-tx depth)
real = load_calls(f"{run}/dropin/real_collapsed.txt")
obsd = load_calls(f"{run}/dropin/pred_obsdepth_collapsed.txt")
rc = {i for i, t in real.items() if t == "annotated"}
rn = {i for i, t in real.items() if t == "novel"}
oc = {i for i, t in obsd.items() if t == "annotated"}
on = {i for i, t in obsd.items() if t == "novel"}
obs_cds_recall = len(oc & rc) / len(rc) if rc else 0.0
obs_novel_prec = len(on & rn) / len(on) if on else 0.0

L = ["Check 4: non-canonical precision at matched CDS recall (higher = better lever)\n"]
L.append(f"{'CDS_recall':>10s}" + "".join(f"{n+' novel/uORF':>20s}" for n, _ in LEVERS))
for t in TARGETS:
    cells = ""
    for name, _ in LEVERS:
        c = curves[name]
        if c is None:
            cells += f"{'(pending)':>20s}"
            continue
        nv, uu = interp(c, t, "novel_prec"), interp(c, t, "uORF_prec")
        cells += f"{(f'{nv:.2f}/{uu:.2f}' if nv is not None else 'n/a'):>20s}"
    L.append(f"{t:>10.2f}{cells}")
L.append("")
L.append(f"isotonic CEILING (pred_obsdepth, shape at REAL depth): "
         f"CDS_recall={obs_cds_recall:.2f}  novel_prec={obs_novel_prec:.2f}")
L.append("  = the best a per-tx count calibration could reach; a fitted isotonic map (tissue A->B) "
         "is the runnable approximation (Check 5 transfer).")
report = "\n".join(L) + "\n"
open(f"{run}/dropin_sweep/lever_comparison.txt", "w").write(report)
print(report)
