#!/usr/bin/env python3
"""Poster panel P11: the proteogenomic lift across every mouse macrophage population.

Three databases per population, searched against identical spectra in all 12 macrophage populations:
the model's ORFs, a naive AUG null, and a near-cognate null.

CURRENT DEFAULTS: the ATTENTION model at a 30-AA ORF floor (2026-08-15). Rule 5 sets the tryptic
floor at 30 aa -- macrophages are trypsin/termini=2 -- and the transformer is the model the
manuscript emphasises. `--source` / `--arm` reproduce the superseded 7-aa mamba panel for comparison.

TWO HEADLINE NUMBERS FELL WHEN THE FLOOR WAS CORRECTED, AND THAT IS THE HONEST DIRECTION.
At 7 aa the null carried 235,686 sequences, ~65% of them sub-30-aa ORFs the ORF-call track excludes
BY RULE. Applying the correct floor shrinks the null to ~84,000 and the ratios fall with it:

    DB size vs null      156x smaller  ->   78x smaller
    discovery density    111x higher   ->   43x higher

The earlier figures were inflated by counting null sequences that should never have been in the
comparison. The model's own net benefit went the other way: above the GENCODE-only baseline in
6 of 12 populations (median +3 peptides) versus 5 of 12 (median -8) before.

THE CLAIM IS COST-NEUTRALITY, NOT A GAIN. Median +3 total peptides out of ~55,000 is indistinguishable
from zero, and well inside the +-50 frozen-vs-unfrozen search-parameter swing measured on BMDM. What
is NOT within noise is the null: it loses ~1,100 peptides in EVERY population, 0 of 12 above baseline.

ABSOLUTE COUNTS ARE SHOWN even where they disfavour the model (panel A): the nulls find MORE novel
peptides in raw terms because they search ~78x more sequences. Density is the per-sequence axis on
which the model wins, and panel B is the axis an experimentalist actually cares about.

Novel columns are at 1% CLASS-SPECIFIC FDR (novel targets vs their own decoys); canonical columns and
dPSM at 1% global FDR (feedback_db_tradeoff_table_columns).

cas12a env. Regenerate: `python3 make_bmdm_proteomics.py`
                        `python3 make_bmdm_proteomics.py --source <7aa frozen_reports> --arm model_predicted`
"""
from __future__ import annotations

import json
import statistics as st
import pathlib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/P11_bmdm_proteomics"
# DEFAULTS ARE THE 30-AA ATTENTION RESULT (2026-08-15). Rule 5 sets the tryptic floor at 30 aa, and
# the transformer is the model the manuscript emphasises. The 7-aa mamba tree is retained and
# selectable via --source/--arm so the two can be compared, but it is SUPERSEDED for macrophages:
# it violates the standing floor.
SRCDIR = NEW / "proteogenomics/data/macrophage_tissue/pgx_xsubtype_aa30/reports"
MODEL_ARM = "model_predicted_attn"
MODELCOL = "#2b6cb0"
NULLCOL = "#a0aec0"
NULLNCCOL = "#cbd5e0"


def load():
    rows = []
    for f in sorted(SRCDIR.glob("report_*.json")):
        j = json.loads(f.read_text())
        a = j["arms"]
        b, m, n = a.get("gencode"), a.get(MODEL_ARM), a.get("null_atg")
        if not all([b, m, n]):
            print(f"  SKIP {f.stem}: missing an arm")
            continue
        base = b["gencode_peptides"] + b["novel_peptides"]
        rec = dict(pop=f.stem.replace("report_", ""), baseline_total=base)
        for tag, d in [("model", m), ("null_atg", n), ("null_nc", a.get("null_nc"))]:
            if d is None:
                continue
            tot = d["gencode_peptides"] + d["novel_peptides"]
            rec[tag] = dict(novel=d["novel_peptides"], canon=d["gencode_peptides"], total=tot,
                            dtot=tot - base, db=d["db_novel_seqs"], dpsm=d["d_gencode_psms"],
                            density=1000 * d["novel_peptides"] / d["db_novel_seqs"]
                            if d["db_novel_seqs"] else 0.0)
        rows.append(rec)
    return rows


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=str(SRCDIR),
                    help="report dir. Default = the 30-aa tree; pass the 7-aa frozen_reports to "
                         "reproduce the superseded mamba panel.")
    ap.add_argument("--arm", default=MODEL_ARM, choices=["model_predicted_attn", "model_predicted"],
                    help="which model arm to plot (attn = transformer, default)")
    a_ = ap.parse_args()
    globals()["SRCDIR"] = pathlib.Path(a_.source)
    globals()["MODEL_ARM"] = a_.arm
    if not SRCDIR.exists():
        raise SystemExit(f"missing {SRCDIR}")
    rows = load()
    if not rows:
        raise SystemExit("no population reports found")
    # Sorted by the honest axis, so panel B reads as a gradient rather than noise.
    rows.sort(key=lambda r: -r["model"]["dtot"])
    pops = [r["pop"] for r in rows]
    x = np.arange(len(rows))
    n_pop = len(rows)

    mdt = [r["model"]["dtot"] for r in rows]
    ndt = [r["null_atg"]["dtot"] for r in rows]
    dbr = [r["null_atg"]["db"] / r["model"]["db"] for r in rows]
    densr = [r["model"]["density"] / r["null_atg"]["density"] for r in rows]
    mwin = sum(1 for v in mdt if v > 0)
    nwin = sum(1 for v in ndt if v > 0)

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(14.0, 4.1), constrained_layout=True)
    w = 0.27

    # A: absolute novel peptides -- the model finds FEWER, shown plainly
    for i, (tag, col, lab) in enumerate([("model", MODELCOL, "model (Poisson-calibrated)"),
                                         ("null_atg", NULLCOL, "null: every AUG"),
                                         ("null_nc", NULLNCCOL, "null: + near-cognate")]):
        v = [r.get(tag, {}).get("novel", np.nan) for r in rows]
        axA.bar(x + (i - 1) * w, v, w, color=col, label=lab, zorder=3)
    axA.set_ylabel("novel peptides identified", fontsize=8.5)
    axA.set_title("A  Novel peptides (absolute) -- the model finds FEWER", fontsize=9, loc="left")
    axA.legend(fontsize=6.8, frameon=False, loc="upper right")
    mn = [r["model"]["novel"] for r in rows]
    axA.text(0.02, 0.90, f"model {min(mn)}-{max(mn)} per population",
             transform=axA.transAxes, fontsize=7.5, color="#2b6cb0", fontweight="bold")

    # B: the honest axis -- what an experimentalist actually takes home
    for i, (tag, col, lab) in enumerate([("model", MODELCOL, "model"),
                                         ("null_atg", NULLCOL, "null: every AUG")]):
        v = [r[tag]["dtot"] for r in rows]
        axB.bar(x + (i - 0.5) * (w + 0.06), v, w + 0.06, color=col, label=lab, zorder=3)
    axB.axhline(0, color="black", lw=1.1, zorder=4)
    axB.set_ylabel("$\\Delta$ TOTAL unique peptides\nvs searching GENCODE alone", fontsize=8.5)
    axB.set_title("B  Cost to the proteome -- the honest axis", fontsize=9, loc="left")
    axB.set_yscale("symlog", linthresh=50)
    axB.text(0.5, 0.10,
             f"model above zero in {mwin}/{n_pop}  (median {st.median(mdt):+.0f} of ~{rows[0]['baseline_total']/1000:.0f},000)\n"
             f"null above zero in {nwin}/{n_pop}  (median {st.median(ndt):+,.0f})",
             transform=axB.transAxes, ha="center", fontsize=7.6, fontweight="bold")
    axB.legend(fontsize=6.8, frameon=False, loc="upper right")

    # C: yield per sequence searched -- why a small curated DB is the right trade
    axC.bar(x, densr, color=MODELCOL, zorder=3)
    axC.set_ylabel("discovery density, model / null\n(novel peptides per 1,000 DB sequences)",
                   fontsize=8.5)
    axC.set_title("C  Yield per sequence searched", fontsize=9, loc="left")
    axC.axhline(1, color="black", lw=1.0, ls="--", zorder=4)
    axC.set_ylim(0, max(densr) * 1.22)   # headroom: at 1.0 the caption sat on the tallest bar
    axC.text(0.5, 0.92, f"median {st.median(densr):.0f}x higher "
                        f"({min(densr):.0f}-{max(densr):.0f}x), on a DB "
                        f"{st.median(dbr):.0f}x smaller",
             transform=axC.transAxes, ha="center", fontsize=8, fontweight="bold", color="#c05621")

    for ax in (axA, axB, axC):
        ax.set_xticks(x)
        ax.set_xticklabels(pops, fontsize=6.8, rotation=40, ha="right")
        ax.grid(axis="y", alpha=0.25, zorder=0)
        ax.tick_params(labelsize=7.5)

    fig.suptitle(f"Mouse macrophages, {n_pop} populations: the model's ORF database is COST-NEUTRAL "
                 f"(median {st.median(mdt):+.0f} peptides) while the naive null loses "
                 f"~{abs(st.median(ndt)):,.0f} in every one", fontsize=9.5)
    HERE.mkdir(parents=True, exist_ok=True)
    fig.savefig(HERE / "P11_bmdm_proteomics.pdf")
    fig.savefig(HERE / "P11_bmdm_proteomics.png", dpi=200)
    plt.close(fig)

    out = {"source_dir": str(SRCDIR), "search_parameters": "frozen",
           "n_populations": n_pop, "populations": rows,
           "model_beats_baseline_n": mwin, "null_atg_beats_baseline_n": nwin,
           "median_model_dtotal": st.median(mdt), "median_null_atg_dtotal": st.median(ndt),
           "db_size_ratio_null_over_model": {"median": st.median(dbr),
                                             "range": [min(dbr), max(dbr)]},
           "density_ratio_model_over_null": {"median": st.median(densr),
                                             "range": [min(densr), max(densr)]},
           "supersedes": "BMDM-only version of this panel, whose +31 total-peptide gain is the best "
                         "of 12 populations and was presented as if general"}
    (HERE / "P11_values.json").write_text(json.dumps(out, indent=2))

    print(f"  wrote {HERE/'P11_bmdm_proteomics.pdf'}  ({n_pop} populations)")
    print(f"  TOTAL unique peptides vs GENCODE-only: model beats it {mwin}/{n_pop} "
          f"(median {st.median(mdt):+.0f}); null_atg {nwin}/{n_pop} (median {st.median(ndt):+,.0f})")
    print(f"  DB {st.median(dbr):.0f}x smaller (median), density {st.median(densr):.0f}x higher")
    print(f"  {'population':<18}{'novel_m':>8}{'novel_n':>9}{'dtot_m':>8}{'dtot_n':>9}"
          f"{'dbratio':>9}{'densratio':>11}")
    for r, dr, dn in zip(rows, dbr, densr):
        print(f"  {r['pop']:<18}{r['model']['novel']:>8}{r['null_atg']['novel']:>9}"
              f"{r['model']['dtot']:>+8}{r['null_atg']['dtot']:>+9,}{dr:>8.0f}x{dn:>10.0f}x")


if __name__ == "__main__":
    main()
