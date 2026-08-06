#!/usr/bin/env python3
"""Compare RiboCode drop-in ORF calls across two or more models on the SAME holdout, to see
how the ORF-call behaviour changes between models (standing eval rule). For each run dir it
reads dropin/{real,pred_obsdepth,pred_preddepth}_collapsed.txt and reports, vs that run's own
`real` reference: detection precision/recall/F1, and the ORF-type composition + over-call ratio
of the standalone `pred_preddepth` calls (the deployment scenario, where over-calling shows up).

Usage: compare_dropin_models.py <label1>=<run_dir1> <label2>=<run_dir2> ...
"""
import sys
from collections import Counter


def load(run_dir, variant):
    f = f"{run_dir}/dropin/{variant}_collapsed.txt"
    d = {}
    with open(f) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            d[p[ci["ORF_ID"]]] = p[ci["ORF_type"]]
    return d


def prf(pred, real):
    inter = len(pred & real)
    p = inter / len(pred) if pred else 0.0
    r = inter / len(real) if real else 0.0
    return p, r, (2 * p * r / (p + r) if (p + r) else 0.0)


models = []
for arg in sys.argv[1:]:
    label, run = arg.split("=", 1)
    models.append((label, run))

TYPES = ["annotated", "uORF", "Overlap_uORF", "novel", "internal", "dORF", "Overlap_dORF"]
print(f"{'model':22s}{'real':>8s}{'obsF1':>8s}{'predF1':>8s}{'predP':>8s}{'predR':>8s}")
comp = {}
for label, run in models:
    real = load(run, "real")
    po = load(run, "pred_obsdepth")
    pp = load(run, "pred_preddepth")
    _, _, of1 = prf(set(po), set(real))
    pp_p, pp_r, pp_f1 = prf(set(pp), set(real))
    comp[label] = (Counter(real.values()), Counter(pp.values()))
    print(f"{label:22s}{len(real):>8d}{of1:>8.3f}{pp_f1:>8.3f}{pp_p:>8.3f}{pp_r:>8.3f}")

print("\nStandalone (pred_preddepth) over-call ratio vs real, by ORF_type:")
print(f"{'ORF_type':16s}" + "".join(f"{lab[:14]:>16s}" for lab, _ in models))
for t in TYPES:
    cells = ""
    for lab, _ in models:
        real_c, pp_c = comp[lab]
        rc, pc = real_c.get(t, 0), pp_c.get(t, 0)
        cells += f"{f'{pc}/{rc}={pc/max(rc,1):.1f}x':>16s}"
    print(f"{t:16s}{cells}")
print("\nlower over-call ratio on non-canonical (uORF/novel/dORF) = fewer spurious calls")
