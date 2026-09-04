#!/usr/bin/env python3
"""Emit the three publication tables from results/orfcaller_comparison_chr1.json.

T1  cross-tool agreement on the OBSERVED profile (no model involved)
T2  predicted vs observed, per tool, per class, at theta = 1
T3  the theta = 1 -> Poisson-calibrated delta, per tool, per class

Absent classes are written as the empty string, never 0: RiboTaper has no `internal` category and
called no `Overlap_uORF`, which is not the same statement as calling zero of something it can call.
"""
from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

NEW = Path(__file__).resolve().parents[2]
SRC = NEW / "results/orfcaller_comparison_chr1.json"
OUT = NEW / "results/orfcaller_tables"
TOOLS = ["RiboCode", "Ribo-TISH", "RiboTaper"]
CLASSES = ["canonical", "uORF", "Overlap_uORF", "dORF", "Overlap_dORF", "internal", "novel"]
# classes a tool cannot emit at all, as opposed to emitted-zero
NO_CATEGORY = {("RiboTaper", "internal"), ("RiboTaper", "Overlap_uORF")}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    d = json.loads(SRC.read_text())
    counts, pvr = d["counts"], d["pred_vs_real"]

    # ---- T1: cross-tool on observed ------------------------------------------------------
    ov = d.get("observed_overlap")
    rows = []
    for cl in CLASSES:
        r = {"class": cl}
        for t in TOOLS:
            n = counts[t]["real"].get(cl)
            r[f"n_{t}"] = "" if (t, cl) in NO_CATEGORY else (n if n else 0)
        if ov:
            for a, b in itertools.combinations(TOOLS, 2):
                k = f"{a}|{b}|{cl}"
                r[f"shared_{a}_{b}"] = ov.get(k, {}).get("shared", "")
                r[f"jaccard_{a}_{b}"] = ov.get(k, {}).get("jaccard", "")
            r["shared_all_three"] = ov.get(f"ALL3|{cl}", "")
        rows.append(r)
    write(OUT / "T1_cross_tool_observed.tsv", rows)

    # ---- T2: predicted vs observed at theta = 1 -------------------------------------------
    rows = []
    for cl in CLASSES + ["ALL"]:
        r = {"class": cl}
        for t in TOOLS:
            m = pvr.get(f"{t}|pred_obsdepth", {}).get(cl)
            blank = (t, cl) in NO_CATEGORY
            for f in ("P", "R", "F1"):
                r[f"{f}_{t}"] = "" if (blank or m is None) else round(m[f], 4)
            for f in ("TP", "FP", "FN", "n_real"):
                r[f"{f}_{t}"] = "" if (blank or m is None) else m[f]
        rows.append(r)
    write(OUT / "T2_pred_vs_observed_theta1.tsv", rows)

    # ---- T3: Poisson-calibration delta -----------------------------------------------------
    rows = []
    for cl in CLASSES + ["ALL"]:
        for t in TOOLS:
            a = pvr.get(f"{t}|pred_obsdepth", {}).get(cl)
            b = pvr.get(f"{t}|pred_preddepth", {}).get(cl)
            if a is None or b is None or (t, cl) in NO_CATEGORY:
                continue
            rows.append({
                "class": cl, "tool": t,
                "P_theta1": round(a["P"], 4), "P_poisson": round(b["P"], 4),
                "dP": round(b["P"] - a["P"], 4),
                "R_theta1": round(a["R"], 4), "R_poisson": round(b["R"], 4),
                "dR": round(b["R"] - a["R"], 4),
                "F1_theta1": round(a["F1"], 4), "F1_poisson": round(b["F1"], 4),
                "dF1": round(b["F1"] - a["F1"], 4),
                "FP_theta1": a["FP"], "FP_poisson": b["FP"],
                "FP_removed": a["FP"] - b["FP"],
                "TP_theta1": a["TP"], "TP_poisson": b["TP"],
                "TP_lost": a["TP"] - b["TP"], "n_real": a["n_real"]})
    write(OUT / "T3_poisson_calibration_delta.tsv", rows)
    return 0


def write(path: Path, rows):
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"wrote {path}  ({len(rows)} rows)")


if __name__ == "__main__":
    raise SystemExit(main())
