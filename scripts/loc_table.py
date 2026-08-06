#!/usr/bin/env python3
"""Render the localization eval (eval_localization.py) across runs as compact comparison tables.

Reads every <run>/localization_metrics.json under the given roots (default: the within-tissue
fold-0 improve runs + any LOTO runs) and prints two tables:

  1. Discrimination -- AUROC separating RiboCode-translated ORFs from untranslated candidate AUG
     ORFs, for the predicted frame0 (periodicity) and predicted density scores, with the
     length-controlled frame0 AUROC, the ORF-length floor, and the observed-profile ceiling.
     Shown for ALL translated ORFs and for the NON-CANONICAL subset (uORF/dORF/novel/...).
  2. Fidelity -- on translated ORFs only, median predicted vs observed frame0 and the
     predicted-vs-observed frame0 Pearson (does the model reproduce the periodicity it was
     never explicitly trained on).

Usage: loc_table.py [run_dir ...]   (each dir must contain localization_metrics.json)
       loc_table.py                 (auto-discovers results/improve/*_f0 + results/loto/*)
Writes results/localization_summary.tsv.
"""
import json
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")


def discover():
    runs = []
    for d in sorted((NEW / "results" / "improve").glob("*_f0")):
        if (d / "localization_metrics.json").exists():
            runs.append(d)
    for d in sorted((NEW / "results" / "loto").glob("*_holdout_*")):
        if (d / "localization_metrics.json").exists():
            runs.append(d)
    return runs


def fmt(x):
    return f"{x:.3f}" if isinstance(x, (int, float)) and x == x else "  .  "


def main():
    argv = sys.argv[1:]
    runs = [Path(a) for a in argv] if argv else discover()
    runs = [r if r.is_absolute() else NEW / r for r in runs]
    metas = []
    for r in runs:
        f = r / "localization_metrics.json"
        if f.exists():
            metas.append(json.loads(f.read_text()))
    if not metas:
        print("no localization_metrics.json found", file=sys.stderr)
        return 1

    out = ["\t".join(["run", "tissue", "stratum", "n_pos", "n_neg",
                      "pred_frame0", "pred_frame0_lenctrl", "pred_density",
                      "length_floor", "obs_frame0_ceiling"])]
    print("=== Discrimination: AUROC translated-vs-untranslated candidate AUG ORFs ===")
    print(f"{'run':26} {'tissue':11} {'stratum':13} {'nP':>5} {'nN':>6} "
          f"{'pF0':>6} {'pF0lc':>6} {'pDen':>6} {'lenFl':>6} {'obsCeil':>7}")
    for m in metas:
        for stratum in ("all", "noncanonical"):
            d = m["discrimination"].get(stratum)
            if not d:
                continue
            print(f"{m['run']:26.26} {m['tissue']:11.11} {stratum:13} "
                  f"{d['n_pos']:>5} {d['n_neg']:>6} "
                  f"{fmt(d['auroc_pred_frame0']):>6} "
                  f"{fmt(d.get('auroc_pred_frame0_lengthctrl')):>6} "
                  f"{fmt(d['auroc_pred_density']):>6} "
                  f"{fmt(d['auroc_orf_length_floor']):>6} "
                  f"{fmt(d['auroc_obs_frame0_ceiling']):>7}")
            out.append("\t".join([m["run"], m["tissue"], stratum,
                       str(d["n_pos"]), str(d["n_neg"]),
                       fmt(d["auroc_pred_frame0"]),
                       fmt(d.get("auroc_pred_frame0_lengthctrl")),
                       fmt(d["auroc_pred_density"]),
                       fmt(d["auroc_orf_length_floor"]),
                       fmt(d["auroc_obs_frame0_ceiling"])]))

    print("\n=== Fidelity: predicted vs observed periodicity on translated ORFs ===")
    print(f"{'run':26} {'stratum':13} {'n':>5} {'medPredF0':>9} {'medObsF0':>9} {'F0pearson':>9}")
    for m in metas:
        for stratum in ("all_positives", "noncanonical", "annotated"):
            fb = m["fidelity"].get(stratum)
            if not fb or "median_pred_frame0" not in fb:
                continue
            print(f"{m['run']:26.26} {stratum:13} {fb['n']:>5} "
                  f"{fmt(fb['median_pred_frame0']):>9} {fmt(fb['median_obs_frame0']):>9} "
                  f"{fmt(fb['frame0_pearson']):>9}")

    (NEW / "results" / "localization_summary.tsv").write_text("\n".join(out) + "\n")
    print(f"\nwrote {NEW / 'results' / 'localization_summary.tsv'}  ({len(metas)} run(s))")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
