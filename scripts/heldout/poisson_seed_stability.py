#!/usr/bin/env python3
"""Aggregate the Poisson CDS-anchored dial across seeds (0,1,2) to confirm it is stable and put
error bars on it for the supplemental figure. Reads each seed's cds_anchored_curve.tsv and reports,
per theta, the mean and min..max across seeds of CDS recall, novel/uORF precision, and novel count.

Usage: poisson_seed_stability.py <run_dir>
"""
import os
import sys

run = sys.argv[1]
SEEDS = {"0": f"{run}/dropin_sweep_poisson",
         "1": f"{run}/dropin_sweep_poisson_seed1",
         "2": f"{run}/dropin_sweep_poisson_seed2"}
KEYS = ("cds_recall", "novel_prec", "uORF_prec", "novel_n")


def read(path):
    d = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[round(float(f[ci["theta"]]), 4)] = {k: float(f[ci[k]]) for k in KEYS}
    return d


data = {}
for s, sub in SEEDS.items():
    p = f"{sub}/cds_anchored_curve.tsv"
    if os.path.exists(p):
        data[s] = read(p)
if not data:
    raise SystemExit("no seed curves found")

thetas = sorted(set().union(*[set(d) for d in data.values()]))
L = [f"Poisson dial seed-stability across seeds {sorted(data)} ({run.split('/')[-1]})\n",
     f"{'theta':>6s}{'CDS_recall (min..max)':>26s}{'novel_prec (min..max)':>26s}"
     f"{'uORF_prec':>18s}{'novel_n':>16s}"]
for th in thetas:
    row = {}
    for k in KEYS:
        vals = [data[s][th][k] for s in data if th in data[s]]
        row[k] = (sum(vals) / len(vals), min(vals), max(vals))
    cr, npc, up, nn = row["cds_recall"], row["novel_prec"], row["uORF_prec"], row["novel_n"]
    L.append(f"{th:>6.2f}"
             f"{f'{cr[0]:.3f} ({cr[1]:.3f}..{cr[2]:.3f})':>26s}"
             f"{f'{npc[0]:.2f} ({npc[1]:.2f}..{npc[2]:.2f})':>26s}"
             f"{f'{up[0]:.2f} ({up[1]:.2f}..{up[2]:.2f})':>18s}"
             f"{f'{nn[0]:.0f} ({nn[1]:.0f}..{nn[2]:.0f})':>16s}")
report = "\n".join(L) + "\n\nStable if min..max ranges are tight (Poisson draw averages out over many ORFs).\n"
open(f"{run}/dropin_sweep_poisson/seed_stability.txt", "w").write(report)
print(report)
