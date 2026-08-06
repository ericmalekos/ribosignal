#!/usr/bin/env python3
"""CDS-anchor the macrophage MODEL-DB selection threshold instead of hard-coding pred_frame0 >= 0.5.

WHY THIS EXISTS
---------------
The macrophage MODEL database is built by `build_a549_dbs.py --thresh 0.5`, i.e. every candidate ORF
whose model-predicted in-frame fraction clears a FIXED 0.5 is admitted. That 0.5 was never calibrated
against anything, and it is not comparable across checkpoints: two models can put the same biological
ORF on opposite sides of 0.5 purely because their score distributions differ. Comparing an attn DB to a
mamba4 DB at a shared fixed cut therefore confounds "which model ranks ORFs better" with "which model
happens to output larger numbers".

This applies the project's CDS-anchoring principle (docs/count_head_calibration_checks.md) to the caller
the macrophage path actually uses. Annotated CDS is the trusted anchor: choose the threshold so that
ANNOTATED-CANONICAL recall hits a target, then read off what that stringency does to the non-canonical
classes. Anchoring on recall rather than on the raw score makes the operating point commensurate across
checkpoints.

Like `scripts/heldout/calibrate_dial.py`, this needs NO observed Ribo-seq -- the anchor is the GENCODE
CDS annotation carried in `candidates.tsv` as `orf_class == "canonical"`. That matters here because the
macrophage populations have no ribosome profiling at all, which is the entire point of the model.

DIFFERENCE FROM calibrate_dial.py: that script anchors the Poisson/theta depth dial routed through
RiboCode's frame test. This one anchors the model's own per-ORF pred_frame0 score. They are different
callers on purpose -- Task 19b showed RiboCode's caller collapses to CTG on predicted profiles, which is
why `enumerate_score_orfs.py` scores ORFs directly from the predicted profile instead.

OUTPUT
------
  <out>/f0_curve_<pop>.tsv     per-threshold canonical recall + per-class survivor counts
  <out>/f0_calibration.json    chosen threshold per population + the pooled threshold
Both arms are reported per the standing two-arm rule: the legacy fixed 0.5 AND the CDS-anchored cut.

cas12a env (stdlib only).
"""
import argparse
import csv
import json
import statistics
from pathlib import Path

NOVEL_CLASSES = ("uORF", "dORF", "lncRNA_orf", "nterm_ext")
GRID = [round(x / 100, 2) for x in range(0, 100, 2)] + [0.99, 1.0]


def load(pop_tsv):
    """-> list of (orf_class, pred_frame0). One row per candidate ORF."""
    rows = []
    with open(pop_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                rows.append((r["orf_class"], float(r["pred_frame0"])))
            except (KeyError, ValueError):
                continue
    return rows


def curve(rows):
    """Per-threshold canonical recall and per-class survivor counts."""
    canon = [f for (k, f) in rows if k == "canonical"]
    n_canon = len(canon)
    out = []
    for th in GRID:
        rec = sum(1 for f in canon if f >= th) / n_canon if n_canon else 0.0
        per = {k: sum(1 for (kk, f) in rows if kk == k and f >= th) for k in NOVEL_CLASSES}
        out.append({"thresh": th, "cds_recall": rec, "n_canonical": sum(1 for f in canon if f >= th),
                    **{f"n_{k}": per[k] for k in NOVEL_CLASSES},
                    "n_novel_total": sum(per.values())})
    return out, n_canon


def pick(cur, target):
    """Largest threshold whose canonical recall is still >= target (most stringent cut meeting the anchor).

    Recall is monotone non-increasing in the threshold, so scanning for the last qualifying grid point
    gives the tightest cut that still retains the target fraction of annotated CDS.
    """
    ok = [c for c in cur if c["cds_recall"] >= target]
    return max(ok, key=lambda c: c["thresh"]) if ok else cur[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--orfs-dir", required=True,
                    help="dir holding <pop>/candidates.tsv (e.g. .../orfs_mamba4_union)")
    ap.add_argument("--out", required=True, help="output dir for the curve + calibration json")
    ap.add_argument("--cds-recall", type=float, default=0.90, help="anchor target (default 0.90)")
    ap.add_argument("--legacy-thresh", type=float, default=0.5,
                    help="the uncalibrated threshold to report alongside (default 0.5)")
    ap.add_argument("--label", default="", help="label recorded in the json (e.g. mamba4_union)")
    a = ap.parse_args()

    od = Path(a.out); od.mkdir(parents=True, exist_ok=True)
    orfs = Path(a.orfs_dir)
    pops = sorted(p.name for p in orfs.iterdir() if (p / "candidates.tsv").exists())
    if not pops:
        raise SystemExit(f"no <pop>/candidates.tsv under {orfs}")

    res, chosen = {}, []
    print(f"{'population':<18}{'n_canon':>9}{'thresh*':>9}{'CDSrec':>8}{'novel@th*':>11}"
          f"{'novel@' + str(a.legacy_thresh):>12}{'ratio':>8}")
    print("-" * 75)
    for pop in pops:
        rows = load(orfs / pop / "candidates.tsv")
        cur, n_canon = curve(rows)
        with open(od / f"f0_curve_{pop}.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(cur[0].keys()), delimiter="\t")
            w.writeheader(); w.writerows(cur)
        best = pick(cur, a.cds_recall)
        legacy = min(cur, key=lambda c: abs(c["thresh"] - a.legacy_thresh))
        ratio = best["n_novel_total"] / legacy["n_novel_total"] if legacy["n_novel_total"] else float("nan")
        res[pop] = {"n_canonical_total": n_canon, "anchored": best, "legacy": legacy}
        chosen.append(best["thresh"])
        print(f"{pop:<18}{n_canon:>9,}{best['thresh']:>9.2f}{best['cds_recall']:>8.3f}"
              f"{best['n_novel_total']:>11,}{legacy['n_novel_total']:>12,}{ratio:>8.2f}")

    pooled = statistics.median(chosen)
    print("-" * 75)
    print(f"  median anchored threshold across {len(pops)} populations: {pooled:.2f}  "
          f"(target CDS recall {a.cds_recall:.2f}; legacy fixed cut {a.legacy_thresh})")
    print(f"  spread: {min(chosen):.2f} - {max(chosen):.2f}")
    (od / "f0_calibration.json").write_text(json.dumps(
        {"label": a.label, "orfs_dir": str(orfs), "cds_recall_target": a.cds_recall,
         "legacy_thresh": a.legacy_thresh, "pooled_threshold": pooled,
         "per_population_threshold": {p: res[p]["anchored"]["thresh"] for p in pops},
         "per_population": res}, indent=2))
    print(f"  wrote {od}/f0_calibration.json")


if __name__ == "__main__":
    main()
