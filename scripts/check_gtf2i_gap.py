#!/usr/bin/env python3
"""Pre-registered locus check for the mm10-coverage retrain arms (methods.md 2026-08-23).

GTF2I ENST00000901263.1 has a 1,617 nt window (transcript positions 1682:3299) where mm1 RNA
coverage reads 8.9% of the transcript's CDS mean and mm10 reads 138.3%. The deployed mm1 model
paints that window near-empty. The question this answers: does a model trained on mm10 coverage
stop doing that?

This is the check that the aggregate metrics cannot answer. The mean RNA multimap rate is 3.8%
across the 57 libraries, so mm10 is a small global change concentrated at loci like this one --
val Pearson can be flat while the locus behaviour changes completely, or vice versa.

Reads the npz written by dump_pred_profiles.py. Reports, per run, the predicted density in the gap
as a percentage of predicted CDS density, alongside the observed target and the RNA input for
reference. Both reference rows are model-free.

Usage:
  check_gtf2i_gap.py --dump <run>/dropin/pred_profiles.npz [--dump ...] [--json out.json]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

R = Path(__file__).resolve().parent.parent
TX = "ENST00000901263.1"
POSTURES = R / "figures/X_gtf2i_isoforms/X_gtf2i_rna_postures_values.json"


def tx_from_npz(path, tx):
    """(pred, obs) per-nt for one tx, or (None, None). Also returns n_skip from the dump."""
    z = np.load(path, allow_pickle=False)
    ids = [str(x) for x in z["tx_ids"]]
    n_skip = int(z["meta"][3]) if "meta" in z else -1
    if tx not in ids:
        return None, None, n_skip
    i = ids.index(tx)
    off = np.concatenate([[0], np.cumsum(z["lengths"])])
    a, b = int(off[i]), int(off[i + 1])
    return z["pred_flat"][a:b], z["obs_flat"][a:b], n_skip


def label_for(p):
    """<run>/dropin/pred_profiles.npz -> the run name; any other layout -> the containing dir."""
    return p.parent.parent.name if p.parent.name.startswith("dropin") else p.parent.name


def pct(arr, gap, cds):
    g = float(np.asarray(arr[gap[0]:gap[1]], dtype=np.float64).mean())
    c = float(np.asarray(arr[cds[0]:cds[1]], dtype=np.float64).mean())
    return (100 * g / c) if c else float("nan"), g, c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", action="append", default=[], required=True,
                    help="path to a run's pred_profiles.npz (repeatable)")
    ap.add_argument("--json", default="")
    a = ap.parse_args()

    pv = json.loads(POSTURES.read_text())
    gap, cds = tuple(pv["gap"]), tuple(pv["cds"])
    print(f"  {TX}  gap {gap[0]}:{gap[1]} ({gap[1]-gap[0]} nt)  CDS {cds[0]}:{cds[1]}\n")
    print("  MODEL-FREE reference (RNA input, from the BAM-level posture analysis):")
    for k in ("mm1", "mm10"):
        print(f"    RNA {k:5s} gap = {pv['postures'][k]['gap_pct_of_cds']:6.1f}% of CDS mean")
    print(f"    observed Ribo in gap = {pv['observed_ribo_in_gap']} of "
          f"{pv['observed_ribo_total']} P-sites  <- the LABEL is holed here too\n")

    rows = []
    print(f"  {'run':52s} {'gap%CDS':>9s} {'n_skip':>7s}")
    for d in a.dump:
        p = Path(d)
        if not p.exists():
            print(f"  {label_for(p)[:52]:52s} {'MISSING':>9s} {'-':>7s}")
            continue
        pred, obs, n_skip = tx_from_npz(p, TX)
        if pred is None:
            print(f"  {label_for(p)[:52]:52s} {'TX ABSENT':>9s} {n_skip:>7d}")
            rows.append({"run": label_for(p), "gap_pct_of_cds": None, "n_skip": n_skip})
            continue
        gp, g, c = pct(pred, gap, cds)
        print(f"  {label_for(p)[:52]:52s} {gp:8.1f}% {n_skip:>7d}")
        rows.append({"run": label_for(p), "gap_pct_of_cds": round(gp, 1),
                     "gap_mean": g, "cds_mean": c, "n_skip": n_skip})
    print("\n  n_skip > 0 means the dump silently dropped tx on a onehot/FASTA length mismatch;")
    print("  a nonzero value invalidates any cross-run comparison built on this dump.")

    if a.json:
        Path(a.json).write_text(json.dumps(
            {"model": [r["run"] for r in rows],
             "source": "scripts/check_gtf2i_gap.py over dump_pred_profiles.py npz; "
                       "reference from figures/X_gtf2i_isoforms/X_gtf2i_rna_postures_values.json",
             "transcript": TX, "gap": list(gap), "cds": list(cds),
             "rna_input_gap_pct": {k: pv["postures"][k]["gap_pct_of_cds"] for k in ("mm1", "mm10")},
             "observed_ribo_in_gap": pv["observed_ribo_in_gap"],
             "observed_ribo_total": pv["observed_ribo_total"],
             "runs": rows}, indent=2) + "\n")
        print(f"  wrote {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
