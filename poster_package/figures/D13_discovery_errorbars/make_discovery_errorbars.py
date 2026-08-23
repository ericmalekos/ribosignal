#!/usr/bin/env python3
"""Figure D13: the model-selected search DB vs the 3-frame null DB as a DB-design tradeoff, computed
entirely from the DETERMINISTIC MSFragger search (1% target-decoy FDR on hyperscore -- no mokapot, no
rescoring stochasticity), so every number is reproducible. Three panels tell one story: predicting which
ORFs to include keeps the DB compact, which mitigates the FDR-recalibration penalty a naive 3-frame null
pays.

(a) GAIN -- discovery efficiency: confident novel peptides per 1,000 novel ORFs (1% FDR). The null finds
    more RAW novel peptides but needs ~3x the DB; per ORF the model is 2-3.4x more efficient, because the
    null's off-target ORFs add decoys that dilute the FDR without proportional discovery.
(b) COST -- canonical/CDS displacement: canonical PSMs whose MSFragger rank-1 is STOLEN by a novel ORF. The
    null displaces 1.5-3.1x more canonical rank-1 than the model.
(c) DRIVER -- decoy load: novel ORFs added to the DB (= added decoy sequences). The null is ~3x the model;
    this decoy inflation is the FDR RECALIBRATION that drives both the lower efficiency (a) and the larger
    displacement (b). The model mitigates by predicting a compact, enriched ORF set.

Reads results/raw_search_fdr.json (+ results/mokapot_stochastic.json for the novel-ORF counts). cas12a env.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/D13_discovery_errorbars"
RAW = NEW / "results/raw_search_fdr.json"
RUNS = NEW / "results/mokapot_stochastic.json"
MODEL_C, NULL_C = "#2C6FBB", "#B0B0B0"
ORDER = ["HBL1", "DoHH2", "SUDHL4"]
LABEL = {"HBL1": "HBL-1", "DoHH2": "DoHH2", "SUDHL4": "SU-DHL-4"}


def main():
    raw = json.load(open(RAW))["datasets"]
    norf = {n: {db: json.load(open(RUNS))["datasets"][n]["per_db"][db]["n_orfs"]
                for db in ("model", "null")} for n in ORDER}
    x = np.arange(len(ORDER)); w = 0.38
    fig, (axa, axb, axc) = plt.subplots(1, 3, figsize=(14.4, 4.4))

    def grouped(ax, valfn, ylabel, title, fmt, pct_base=None):
        for db, c, off in [("model", MODEL_C, -w/2), ("null", NULL_C, +w/2)]:
            vals = [valfn(name, db) for name in ORDER]
            ax.bar(x + off, vals, w, color=c, edgecolor="white", lw=0.5,
                   label=("Model-selected DB" if db == "model" else "3-frame null DB"))
            for xi, name, v in zip(x + off, ORDER, vals):
                extra = f"\n({100*v/pct_base(name):+.0f}%)" if pct_base else ""
                ax.text(xi, v + (np.sign(v) or 1) * max(abs(v)*0.02, 0.01*max(abs(v) for v in vals)),
                        fmt.format(v) + extra, ha="center",
                        va=("bottom" if v >= 0 else "top"), fontsize=6.8,
                        color=(MODEL_C if db == "model" else "#555"))
        ax.set_xticks(x); ax.set_xticklabels([LABEL[n] for n in ORDER], fontsize=8.6)
        ax.set_ylabel(ylabel, fontsize=8.6)
        ax.legend(frameon=False, fontsize=7.2, loc="best")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(title, fontsize=9, loc="left")
        ax.margins(y=0.16)

    # (a) discovery efficiency = confident novel / 1k ORFs
    grouped(axa, lambda n, db: raw[n]["per_db"][db]["conf_novel"] * 1000.0 / norf[n][db],
            "Novel peptides per 1,000 ORFs\n(1% target-decoy FDR)",
            "a  GAIN: discovery efficiency\n(model 2-3.4x more novel per ORF)", "{:.2f}")
    axa.set_ylim(bottom=0)

    # (b) canonical rank-1 displaced (negative)
    grouped(axb, lambda n, db: raw[n]["per_db"][db]["rank1_canonical_displaced"],
            "Canonical/CDS PSMs displaced vs baseline\n(MSFragger search rank-1)",
            "b  COST: CDS displacement\n(null steals 1.5-3.1x more)", "{:+,.0f}",
            pct_base=lambda n: raw[n]["baseline_rank1_canonical"])
    axb.axhline(0, color="#333", lw=0.8)

    # (c) decoy load = novel ORFs added (thousands)
    grouped(axc, lambda n, db: norf[n][db] / 1000.0,
            "Novel ORFs added to DB (thousands)\n= added decoy sequences",
            "c  DRIVER: decoy load\n(null ~3x -> FDR recalibration penalty)", "{:.0f}k")
    axc.set_ylim(bottom=0)

    fig.suptitle("Predicting which ORFs to include keeps the DB compact: 2-3.4x higher discovery efficiency, "
                 "far less CDS displacement, ~1/3 the decoy load vs the naive 3-frame null",
                 fontsize=9.8, y=1.02)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"D13_discovery_errorbars.{ext}", dpi=300, bbox_inches="tight")

    summary = {}
    for name in ORDER:
        summary[name] = {"baseline_rank1_canonical": raw[name]["baseline_rank1_canonical"]}
        for db in ("model", "null"):
            e = raw[name]["per_db"][db]
            summary[name][db] = {"conf_novel": e["conf_novel"], "n_orfs": norf[name][db],
                                 "efficiency_per_1k": e["conf_novel"] * 1000.0 / norf[name][db],
                                 "rank1_canonical_displaced": e["rank1_canonical_displaced"],
                                 "fdr_threshold": e["fdr_threshold"]}
    (HERE / "D13_values.json").write_text(json.dumps(summary, indent=2) + "\n")
    for name in ORDER:
        m, n = raw[name]["per_db"]["model"], raw[name]["per_db"]["null"]
        em = m["conf_novel"]*1000.0/norf[name]["model"]; en = n["conf_novel"]*1000.0/norf[name]["null"]
        print(f"  {LABEL[name]}: eff model {em:.2f} vs null {en:.2f} ({em/en:.1f}x) | CDS displaced "
              f"model {m['rank1_canonical_displaced']:+} vs null {n['rank1_canonical_displaced']:+} | "
              f"ORFs {norf[name]['model']/1000:.0f}k/{norf[name]['null']/1000:.0f}k")


if __name__ == "__main__":
    main()
