#!/usr/bin/env python3
"""Deterministic DB-tradeoff from the RAW MSFragger search output (search_<db>/*.tsv) -- no ms2rescore, no
mokapot, so every number is reproducible. For each dataset x DB (canonical baseline, model, null):

- rank1_canonical : # spectra whose MSFragger rank-1 (hit_rank==1) is a canonical target. Drop vs baseline =
                    canonical PSMs mechanically DISPLACED when novel ORFs enter the DB (COST 1).
- global 1% target-decoy FDR on hyperscore over ALL rank-1 PSMs -> confident canonical + novel + the score
  THRESHOLD + the decoy count. A bigger decoy DB (null ~3x model) inflates the decoy competition and raises
  the threshold = FDR RECALIBRATION (COST 2), shown in the fully controlled search space.

Writes results/raw_search_fdr.json. cas12a env (stdlib only)."""
from __future__ import annotations

import csv
import glob
import json
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
DATASETS = {"HBL1": "HBL1_pilot", "DoHH2": "DoHH2_pilot", "SUDHL4": "SUDHL4_pilot"}


def classify(prot_field):
    prots = [p for p in prot_field.split(";") if p]
    tgt = [p for p in prots if not (p.startswith("REV_") or p.startswith("rev_"))]
    if not tgt:
        return "decoy"
    return "novel" if all("nuORF|" in p for p in tgt) else "canonical"


def load_rank1(pilot, db):
    """deterministic rank-1 PSMs: list of (hyperscore, class)."""
    out = []
    for f in glob.glob(str(NEW / f"proteogenomics/data/{pilot}/search_{db}/*.tsv")):
        with open(f) as fh:
            rd = csv.reader(fh, delimiter="\t")
            hdr = next(rd)
            ci = {h: i for i, h in enumerate(hdr)}
            ir, ip, ih = ci["hit_rank"], ci["proteins"], ci["hyperscore"]
            for row in rd:
                if row[ir] != "1":
                    continue
                try:
                    out.append((float(row[ih]), classify(row[ip])))
                except (ValueError, IndexError):
                    pass
    return out


def global_fdr(psms, fdr=0.01):
    """1% target-decoy FDR on hyperscore. Returns confident canonical, novel, threshold, total decoys."""
    ps = sorted(psms, key=lambda z: -z[0])
    t = d = 0
    thr = None
    for h, k in ps:
        if k == "decoy":
            d += 1
        else:
            t += 1
        if t > 0 and d / t <= fdr:
            thr = h
    ndec = sum(1 for _, k in ps if k == "decoy")
    if thr is None:
        return 0, 0, None, ndec
    canon = sum(1 for h, k in ps if h >= thr and k == "canonical")
    novel = sum(1 for h, k in ps if h >= thr and k == "novel")
    return canon, novel, thr, ndec


def main():
    out = {"datasets": {}}
    for name, pilot in DATASETS.items():
        ds = {"per_db": {}}
        for db in ("canonical", "model", "null"):
            psms = load_rank1(pilot, db)
            r1_canon = sum(1 for _, k in psms if k == "canonical")
            c, nv, thr, ndec = global_fdr(psms)
            ds["per_db"][db] = {"rank1_canonical": r1_canon, "conf_canonical": c, "conf_novel": nv,
                                "fdr_threshold": thr, "n_decoy_rank1": ndec, "n_rank1": len(psms)}
        b = ds["per_db"]["canonical"]
        for db in ("model", "null"):
            ds["per_db"][db]["rank1_canonical_displaced"] = ds["per_db"][db]["rank1_canonical"] - b["rank1_canonical"]
            ds["per_db"][db]["conf_canonical_delta"] = ds["per_db"][db]["conf_canonical"] - b["conf_canonical"]
        ds["baseline_rank1_canonical"] = b["rank1_canonical"]
        ds["baseline_conf_canonical"] = b["conf_canonical"]
        out["datasets"][name] = ds
        print(f"{name}: base rank1 {b['rank1_canonical']} conf {b['conf_canonical']} (thr {b['fdr_threshold']})")
        for db in ("model", "null"):
            e = ds["per_db"][db]
            print(f"  {db}: rank1-canon {e['rank1_canonical']} ({e['rank1_canonical_displaced']:+}) | "
                  f"conf-canon {e['conf_canonical']} ({e['conf_canonical_delta']:+}) conf-novel {e['conf_novel']} | "
                  f"1%thr {e['fdr_threshold']:.2f} decoys {e['n_decoy_rank1']}")
    (NEW / "results/raw_search_fdr.json").write_text(json.dumps(out, indent=2) + "\n")
    print("wrote results/raw_search_fdr.json")


if __name__ == "__main__":
    main()
