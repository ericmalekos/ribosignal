#!/usr/bin/env python3
"""Aggregate the mokapot stochasticity replicates: for each immunopeptidome and each seed's rescore_seed<K>
dir, compute the class-FDR novel-peptide count per DB (reusing compare_rescored's own load_rank1 +
class_fdr_novel), then report mean/SD/min/max over seeds. Writes results/mokapot_stochastic.json.
Usage: aggregate_mokapot_stochastic.py [N]"""
from __future__ import annotations

import glob
import json
import os
import statistics as st
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
sys.path.insert(0, str(NEW / "proteogenomics/scripts"))
from compare_rescored import class_fdr_novel, load_rank1  # noqa: E402

DATASETS = {"HBL1": "HBL1_pilot", "DoHH2": "DoHH2_pilot", "SUDHL4": "SUDHL4_pilot"}
DBS = ["model", "null"]


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    out = {"n_requested": n, "datasets": {}}
    for name, pilot in DATASETS.items():
        # discover which seed dirs actually completed all DBs
        seed_dirs = sorted(glob.glob(str(NEW / f"proteogenomics/data/{pilot}/rescore_seed*")))
        per_db = {db: [] for db in DBS}
        used_seeds = 0
        for sd in seed_dirs:
            tag = os.path.basename(sd)                       # rescore_seed3
            ok = all(glob.glob(str(NEW / f"proteogenomics/data/{pilot}/{tag}/{db}/*.psms.tsv")) for db in DBS)
            if not ok:
                continue
            used_seeds += 1
            for db in DBS:
                psms = load_rank1(db, pilot, tag)
                keep, t, d = class_fdr_novel(psms, 0.01)
                per_db[db].append(len(keep))
        ds = {"n_seeds": used_seeds, "per_db": {}}
        for db in DBS:
            v = per_db[db]
            # novel-ORF DB size (deterministic; for the discovery-rate = peptides per 1k ORFs)
            cm = NEW / f"proteogenomics/data/{pilot}/db/{db}_class_map.tsv"
            n_orfs = (sum(1 for _ in open(cm)) - 1) if cm.exists() else 0
            if v:
                ds["per_db"][db] = {
                    "counts": v, "mean": st.mean(v),
                    "sd": st.stdev(v) if len(v) > 1 else 0.0,
                    "min": min(v), "max": max(v), "n_orfs": n_orfs}
        out["datasets"][name] = ds
        row = "  ".join(f"{db} {st.mean(per_db[db]):.1f}+/-{(st.stdev(per_db[db]) if len(per_db[db])>1 else 0):.1f} "
                        f"(n={len(per_db[db])}, {min(per_db[db])}-{max(per_db[db])})"
                        for db in DBS if per_db[db])
        print(f"{name}: {row}")
    dst = NEW / "results/mokapot_stochastic.json"
    dst.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()
