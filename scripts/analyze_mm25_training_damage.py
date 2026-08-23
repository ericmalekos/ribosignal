#!/usr/bin/env python3
"""EXPERIMENT (a): how much of the TRAINING SET's Ribo-seq signal did --outFilterMultimapNmax 1
discard, and where?

Per transcript, from the mm25 re-alignment of all 74 Chothani training libraries:

    multimap_frac = n_multi / (n_unique + n_multi)

n_unique is what mm1 kept; n_multi is what it threw away. A transcript at 0.60 had 60% of its
Ribo-seq removed before the model ever saw it -- and the profile head is trained on the SHAPE of
what remained, so it is being taught that those positions are empty. That is worse than noise.

Scope: transcripts the model is actually TRAINED on, not the whole universe. A transcript with no
training signal cannot be damaged in a way that matters here.

  analyze_mm25_training_damage.py
"""
from __future__ import annotations
import csv, gzip, json
from collections import defaultdict
from pathlib import Path
import numpy as np

NEW = Path(__file__).resolve().parents[1]
CNT = NEW / "data/mm25_diagnostic/chothani_nh_counts"
OUT = NEW / "data/mm25_diagnostic"


def main():
    man = {}
    with open(NEW / "data/external/chothani_ribo_refetch/manifest.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = list(r.values())
            man[k[0]] = k[1].split("_")[0] if len(k) > 1 else "?"

    uni = defaultdict(int); mul = defaultdict(int)
    per_tis = defaultdict(lambda: [0, 0])
    files = sorted(CNT.glob("*.nh.tsv.gz"))
    for f in files:
        srr = f.name.split(".")[0]; tis = man.get(srr, "?")
        with gzip.open(f, "rt") as fh:
            for line in fh:
                t, u, m = line.rstrip("\n").split("\t")
                u, m = int(u), int(m)
                uni[t] += u; mul[t] += m
                per_tis[tis][0] += u; per_tis[tis][1] += m
    print(f"  {len(files)} libraries, {len(uni):,} transcripts with any signal")

    print(f"\n  PER TISSUE -- share of Ribo-seq reads mm1 discarded")
    print(f"  {'tissue':<14} {'unique':>14} {'multi':>14} {'% discarded':>12}")
    for t, (u, m) in sorted(per_tis.items(), key=lambda kv: -kv[1][1] / max(kv[1][0] + kv[1][1], 1)):
        print(f"  {t:<14} {u:>14,} {m:>14,} {100*m/max(u+m,1):>11.1f}%")

    # restrict to transcripts the model actually TRAINS on: >=50 pooled mm1 P-sites (min_train_signal)
    order = (NEW / "data/packed_union/tx_order.txt").read_text().split()
    keep = set(order)
    trained = {t for t in keep if uni.get(t, 0) >= 50}
    print(f"\n  universe {len(keep):,}   with >=50 mm1 reads (trainable) {len(trained):,}")

    bt = {}
    with open(NEW / "data/tx2biotype.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            bt[r["tx_id"]] = r.get("transcript_type", "?")

    frac = {t: mul[t] / (uni[t] + mul[t]) for t in trained if (uni[t] + mul[t]) > 0}
    a = np.array(list(frac.values()))
    print(f"\n  MULTIMAP FRACTION over {len(a):,} trainable transcripts")
    for q in (10, 25, 50, 75, 90, 95, 99):
        print(f"    p{q:<3} {np.percentile(a,q):.3f}")
    for thr in (0.2, 0.3, 0.5, 0.8):
        n = int((a > thr).sum())
        print(f"    >{int(thr*100)}% of signal discarded: {n:,} tx ({100*n/len(a):.1f}%)")

    print(f"\n  BY BIOTYPE (trainable transcripts)")
    print(f"  {'biotype':<18} {'n':>7} {'median frac':>12} {'>50% lost':>11}")
    byb = defaultdict(list)
    for t, v in frac.items():
        byb[bt.get(t, "?")].append(v)
    for b, v in sorted(byb.items(), key=lambda kv: -len(kv[1]))[:6]:
        v = np.array(v)
        print(f"  {b:<18} {len(v):>7} {np.median(v):>12.3f} {100*(v>0.5).mean():>10.1f}%")

    json.dump({"model": None,
               "source": "mm25 re-alignment of 74 Chothani training libraries; "
                         "data/mm25_diagnostic/chothani_nh_counts/*.nh.tsv.gz",
               "provenance_note": "DIAGNOSTIC, mm25 violates the project Ribo-seq rule on purpose. "
                                  "multimap_frac = discarded/(kept+discarded) per transcript.",
               "n_libraries": len(files), "n_trainable_tx": len(a),
               "median_multimap_frac": round(float(np.median(a)), 4),
               "frac_over_50pct_lost": round(float((a > 0.5).mean()), 4),
               "per_tissue_pct_discarded": {t: round(100*m/max(u+m,1), 2)
                                            for t, (u, m) in per_tis.items()}},
              open(OUT / "training_damage_summary.json", "w"), indent=2)
    with gzip.open(OUT / "training_multimap_frac.tsv.gz", "wt") as fh:
        fh.write("tx_id\tn_unique\tn_multi\tmultimap_frac\n")
        for t in sorted(frac, key=lambda x: -frac[x]):
            fh.write(f"{t}\t{uni[t]}\t{mul[t]}\t{frac[t]:.4f}\n")
    print(f"\n  wrote training_damage_summary.json + training_multimap_frac.tsv.gz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
