#!/usr/bin/env python3
"""O2 two-arm: score each Poisson-calibrated theta of the Wang de-novo calls vs the Janich ceiling.

For every theta in the sweep it reports, split CDS vs non-canonical:
  - call counts (total / CDS / non-canonical)
  - CDS RECALL vs Wang_obs (the CDS-anchored dial: pick the operating point by this)
  - F1 / precision / recall vs Janich_obs (the O2 de-novo bar; ceiling = Wang_obs vs Janich_obs)
so we can see whether Poisson calibration closes the standard-arm gap (CDS 0.621, non-canon 0.067)
toward the ceiling (CDS 0.735, non-canon 0.505).
"""
import glob
import os

NEW = "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model"
ECH = ("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/"
       "biotype_probe/expression_context_human")
WANG_OBS = f"{NEW}/data/heldout_psites/mouse_wang_liver/mouse_wang_liver_collapsed.txt"
JANICH_OBS = f"{ECH}/data/ribocode_mouse_liver/Liver_5samp/Liver_5samp_collapsed.txt"
SWEEP = f"{NEW}/results/o2_liver/wang_nokozak/dropin_sweep_poisson"


def load(path):
    d = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            oid = p[ci["ORF_ID"]]
            parts = oid.rsplit("_", 3)
            key = "_".join(parts[-3:]) if len(parts) >= 4 else oid
            d[key] = p[ci["ORF_type"]]
    return d


def cds(d):
    return {k for k, t in d.items() if t == "annotated"}


def noncan(d):
    return {k for k, t in d.items() if t != "annotated"}


def prf(pred, ref):
    i = len(pred & ref)
    p = i / len(pred) if pred else 0.0
    r = i / len(ref) if ref else 0.0
    f = 2 * p * r / (p + r) if (p + r) else 0.0
    return p, r, f


def main():
    wang, jan = load(WANG_OBS), load(JANICH_OBS)
    wc, jc, wn, jn = cds(wang), cds(jan), noncan(wang), noncan(jan)
    # references / bars
    _, _, cds_ceil = prf(wc, jc)
    _, _, nc_ceil = prf(wn, jn)
    print(f"refs: Wang_obs CDS={len(wc)} nc={len(wn)}  Janich_obs CDS={len(jc)} nc={len(jn)}")
    print(f"CEILING (Wang_obs vs Janich_obs):   CDS F1={cds_ceil:.3f}  non-canon F1={nc_ceil:.3f}")
    print("STANDARD ARM (theta=1, no Poisson):  CDS F1=0.621  non-canon F1=0.067  (o2_summary)")
    print()
    thetas = sorted(glob.glob(f"{SWEEP}/theta_*"),
                    key=lambda p: float(p.rsplit("_", 1)[1]))
    hdr = (f"{'theta':>6s} {'calls':>7s} {'cdsN':>6s} {'ncN':>7s} | "
           f"{'CDSrecWang':>10s} | {'CDS_F1j':>8s} {'CDS_Pj':>7s} {'CDS_Rj':>7s} | "
           f"{'nc_F1j':>7s} {'nc_Pj':>6s}")
    print(hdr)
    print("-" * len(hdr))
    for d in thetas:
        th = d.rsplit("_", 1)[1]
        f = f"{d}/pred_preddepth_collapsed.txt"
        if not os.path.exists(f):
            print(f"{th:>6s}  (pending)")
            continue
        m = load(f)
        mc, mn = cds(m), noncan(m)
        cds_r_wang = len(mc & wc) / len(wc) if wc else 0.0    # CDS recall vs Wang (the anchor)
        cp, cr, cf = prf(mc, jc)                              # CDS vs Janich (O2)
        _, _, nf = prf(mn, jn)
        npj = len(mn & jn) / len(mn) if mn else 0.0
        print(f"{th:>6s} {len(m):>7d} {len(mc):>6d} {len(mn):>7d} | "
              f"{cds_r_wang:>10.3f} | {cf:>8.3f} {cp:>7.3f} {cr:>7.3f} | {nf:>7.3f} {npj:>6.3f}")
    print("\nCDS-anchored operating point = theta whose CDSrecWang ~ the ceiling Wang->Janich")
    print("CDS recall (0.756); read its CDS_F1j / nc_F1j vs ceiling (CDS 0.735, nc 0.505).")


if __name__ == "__main__":
    main()
