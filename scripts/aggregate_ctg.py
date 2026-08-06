#!/usr/bin/env python3
"""Consolidate the start-codon-stratified drop-in metrics (dropin_ctg_metrics.json) across runs into a
single ATG-vs-CTG table, to see whether the Hepatocytes CTG F1 ~0.41 holds up across backends, an
independent human dataset (Ruiz-Orera), cross-species mouse (Wang), and stronger mixers (attn4, mamba).
Reads each run's results/<run>/dropin_ctg/dropin_ctg_metrics.json. cas12a env (stdlib only)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model")


def pct_ctg(counts):
    tot = sum(counts.values()) or 1
    return 100.0 * counts.get("CTG", 0) / tot


def uorf_frac(breakdown):
    tot = sum(breakdown.values()) or 1
    return 100.0 * (breakdown.get("uORF", 0) + breakdown.get("Overlap_uORF", 0)) / tot


def main():
    runs = sys.argv[1:]
    rows = []
    for r in runs:
        p = NEW / "results" / r / "dropin_ctg" / "dropin_ctg_metrics.json"
        if not p.exists():
            print(f"MISSING {p}", file=sys.stderr)
            continue
        d = json.loads(p.read_text())
        atg = d["pred_obsdepth_vs_real_ATG"]["strict"]
        ctg = d["pred_obsdepth_vs_real_CTG"]["strict"]
        ctg_sa = d["pred_preddepth_vs_real_CTG"]["strict"]
        rows.append({
            "run": r,
            "atg_f1": atg["f1"],
            "ctg_nreal": ctg["n_real"], "ctg_npred": ctg["n_pred"], "ctg_match": ctg["n_match"],
            "ctg_prec": ctg["precision"], "ctg_rec": ctg["recall"], "ctg_f1": ctg["f1"],
            "ctg_sa_f1": ctg_sa["f1"],
            "pct_ctg_real": pct_ctg(d["codon_counts"]["real"]),
            "pct_ctg_pred": pct_ctg(d["codon_counts"]["pred_obsdepth"]),
            "uorf_of_ctg": uorf_frac(d["real_CTG_type_breakdown"]),
        })

    def short(r):
        return (r.replace("loto/orf_v2_attn_", "").replace("heldout/", "").replace("ablation/", "")
                .replace("_holdout_Hepatocytes", "").replace("_onehot", "/onehot"))

    print("\n================ CTG drop-in across runs (pred_obsdepth vs real, genomic key) ================")
    print(f"{'run':<26} {'ATG_F1':>6} | {'nRealC':>6} {'nPredC':>6} {'match':>5} "
          f"{'precC':>6} {'recC':>6} {'F1C':>6} {'F1saC':>6} | {'%CTGr':>6} {'%uORF':>6}")
    for r in rows:
        print(f"{short(r['run']):<26} {r['atg_f1']:>6.3f} | {r['ctg_nreal']:>6} {r['ctg_npred']:>6} "
              f"{r['ctg_match']:>5} {r['ctg_prec']:>6.3f} {r['ctg_rec']:>6.3f} {r['ctg_f1']:>6.3f} "
              f"{r['ctg_sa_f1']:>6.3f} | {r['pct_ctg_real']:>5.1f}% {r['uorf_of_ctg']:>5.0f}%")
    if rows:
        cf = [r["ctg_f1"] for r in rows]
        af = [r["atg_f1"] for r in rows]
        print(f"\nATG F1 range {min(af):.3f}-{max(af):.3f}; CTG F1 range {min(cf):.3f}-{max(cf):.3f} "
              f"(mean {sum(cf)/len(cf):.3f})")
        out = NEW / "results/ctg_across_runs.json"
        out.write_text(json.dumps(rows, indent=2))
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
