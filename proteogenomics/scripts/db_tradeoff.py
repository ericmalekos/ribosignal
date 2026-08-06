#!/usr/bin/env python3
"""DB-tradeoff accounting: GAIN (confident novel peptides) vs PAY (canonical/CDS PSMs lost), decomposed into
the TWO mechanisms by which a bloated search DB hurts, so the model's mitigation is explicit.

Metrics per dataset x DB, over the 8 replicate rescoring runs where relevant:
- n_orfs                : novel ORFs added to the DB = the OFF-TARGET / DECOY LOAD (deterministic). The null
                          is ~3x the model. Each added ORF also adds a REV decoy, so this is what inflates the
                          FDR denominator.
- novel_*               : confident novel peptides at 1% class-specific FDR (from mokapot_stochastic.json;
                          median + IQR + raw counts; sparse + stochastic).
- search_canonical      : DETERMINISTIC MSFragger search rank-1 canonical PSMs (before_rescoring_rank==1,
                          canon_t). Displacement vs baseline = spectra whose top match is STOLEN by a novel
                          ORF. Mechanism 1 (mechanical displacement).
- conf_canonical_*       : canonical PSMs at 1% GLOBAL mokapot q (canon_t, q<=0.01), median + IQR over 8 runs.
                          Mechanism 2 (FDR recalibration): a larger decoy DB raises the q-threshold, so fewer
                          confident canonical survive -- in the ACTUAL analysis space (mokapot partly
                          compensates for the large canonical class, so this is smaller than the raw
                          displacement; the recalibration bites hardest on the sparse novel class = the lower
                          null discovery RATE in the gain panel).

Writes results/db_tradeoff.json. cas12a env.
"""
from __future__ import annotations

import ast
import csv
import glob
import json
import statistics as st
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
sys.path.insert(0, str(NEW / "proteogenomics/scripts"))
from compare_rescored import klass  # noqa: E402

DATASETS = {"HBL1": "HBL1_pilot", "DoHH2": "DoHH2_pilot", "SUDHL4": "SUDHL4_pilot"}
QCUT = 0.01


def parse_accs(s):
    try:
        return tuple(ast.literal_eval(s))
    except Exception:
        return tuple(a.strip(" '\"[]") for a in s.split(",") if a.strip(" '\"[]"))


def scan(pilot, db, seed_tag):
    """One pass over a seed's psms.tsv: (search_rank1_canon, conf_canon@q<=0.01)."""
    search_canon = 0
    conf_canon = set()
    for f in glob.glob(str(NEW / f"proteogenomics/data/{pilot}/{seed_tag}/{db}/*.psms.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            a = parse_accs(r["protein_list"]); k = klass(a)
            if r.get("provenance:before_rescoring_rank") in ("1", "1.0") and k == "canon_t":
                search_canon += 1
            if r.get("rank") in ("1", "1.0") and k == "canon_t":
                try:
                    if float(r["qvalue"]) <= QCUT:
                        conf_canon.add(r["peptidoform"] + "|" + r["spectrum_id"])
                except (ValueError, KeyError):
                    pass
    return search_canon, len(conf_canon)


def stats(v):
    if not v:
        return None
    v = sorted(v)
    return {"median": v[len(v)//2], "iqr": [v[len(v)//4], v[(3*len(v))//4]], "min": min(v), "max": max(v)}


def main():
    runs = json.load(open(NEW / "results/mokapot_stochastic.json"))["datasets"]
    out = {"datasets": {}}
    for name, pilot in DATASETS.items():
        seeds = [Path(p).name for p in sorted(glob.glob(str(NEW / f"proteogenomics/data/{pilot}/rescore_seed*")))]
        ds = {"per_db": {}}
        for db in ("canonical", "model", "null"):
            search_c, conf_c = [], []
            for tag in seeds:
                sc, cc = scan(pilot, db, tag)
                search_c.append(sc); conf_c.append(cc)
            entry = {"search_canonical": st.median(search_c),
                     "conf_canonical": stats(conf_c)}
            if db != "canonical":
                rd = runs[name]["per_db"][db]
                c = sorted(rd["counts"])
                entry.update(n_orfs=rd["n_orfs"], novel_counts=rd["counts"],
                             novel_median=c[len(c)//2], novel_iqr=[c[len(c)//4], c[(3*len(c))//4]],
                             novel_min=min(c), novel_max=max(c))
            ds["per_db"][db] = entry
        base_search = ds["per_db"]["canonical"]["search_canonical"]
        base_conf = ds["per_db"]["canonical"]["conf_canonical"]["median"]
        for db in ("model", "null"):
            ds["per_db"][db]["canonical_displaced"] = ds["per_db"][db]["search_canonical"] - base_search
            ds["per_db"][db]["conf_canonical_delta"] = ds["per_db"][db]["conf_canonical"]["median"] - base_conf
        ds["baseline_search_canonical"] = base_search
        ds["baseline_conf_canonical"] = base_conf
        out["datasets"][name] = ds
        m, n = ds["per_db"]["model"], ds["per_db"]["null"]
        print(f"{name}: base search {base_search} conf {base_conf}")
        for db, e in (("model", m), ("null", n)):
            print(f"  {db}: novel {e['novel_median']} (IQR {e['novel_iqr']}) | ORFs {e['n_orfs']:,} | "
                  f"CDS displaced {e['canonical_displaced']:+} | conf-canon FDR-recal {e['conf_canonical_delta']:+} "
                  f"(median {e['conf_canonical']['median']})")
    dst = NEW / "results/db_tradeoff.json"
    dst.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()
