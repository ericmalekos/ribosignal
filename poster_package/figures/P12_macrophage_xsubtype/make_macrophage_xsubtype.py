#!/usr/bin/env python3
"""P12 -- the full macrophage sweep: a PREDICTED ORF database vs a REAL Ribo-seq one, 12 populations.

P11 asks whether the model beats naive enumeration. This asks the harder and more useful question:
does a predicted database substitute for one built from actual Ribo-seq? The data to answer it has
existed as a table (`pgx_xsubtype/crosssubtype_table.md`) with no figure.

CURRENT DEFAULTS: the ATTENTION model at a 30-AA ORF floor (2026-08-15), matching P11 and P13. Rule 5
sets the tryptic floor at 30 aa -- macrophages are trypsin/termini=2 -- and the transformer is the
model the manuscript emphasises. `--source` / `--arm` reproduce the superseded 7-aa mamba panel.

Three databases, searched against identical spectra in each of 12 mouse macrophage populations
(18 mzML fractions each, so the sweep is depth-balanced and cross-population comparison is fair):

  model_predicted_attn     the model's ORFs, called PER POPULATION from that population's own data
  model_ribocode_bmdm_nt   REAL Ribo-seq, RiboCode on BMDM untreated (506 seqs at 30 aa)
  model_ribocode_bmdm      REAL Ribo-seq, RiboCode on BMDM NT + LPS (1,265 seqs at 30 aa)

DATABASE SIZES MOVE WITH THE FLOOR, so the panel measures them rather than printing them from this
docstring: the same two Ribo-seq databases are 699 and 2,173 sequences at 7 aa and 506 and 1,265 at
30 aa. The suptitle likewise derives its comparison and announces a reversal if one ever occurs.

THE ASYMMETRY THAT DEFINES THE COMPARISON. The two Ribo-seq databases are BMDM-derived and SHARED by
every row -- one fixed database applied to all 12 populations. The model's database is rebuilt per
population. So BMDM is the MATCHED case for the Ribo-seq arms and the other 11 are a TRANSFER test,
and it is drawn separately for that reason. Every summary statistic here is computed on the 11
TRANSFER populations, because including the matched row would flatter the Ribo-seq arms on their own
home ground.

WHAT IT SHOWS, STATED PLAINLY BECAUSE IT IS NOT FLATTERING. On the 11 transfer populations a FIXED
506-sequence database derived from BMDM Ribo-seq matches or beats the per-population predicted
database on every axis: comparable novel peptides (median 26 vs 24), better on total unique peptides
(median +18 vs -3; above baseline in 8/11 vs 5/11), and 2.4x the discovery density (51.4 vs 21.3 per
1,000 sequences) on a database less than half the size (506 vs a median 1,091).

CORRECTING THE FLOOR DID NOT RESCUE THE MODEL HERE, and that is worth stating because the 30-AA
rebuild helped it elsewhere. The 7-aa mamba version of this panel read 27 vs 26 novel, +4 vs -8
total, 3.0x density. Moving to 30 aa and the attention model narrowed the density gap (3.0x -> 2.4x)
but WIDENED the total-peptide gap (+4/-8 -> +18/-3). The conclusion is unchanged and the panel is
still the honest counterweight to P11.

The model's case in this application is therefore NOT that it beats Ribo-seq. It is that it does not
require Ribo-seq: the predicted arm needs only RNA-seq plus sequence, while every Ribo-seq arm here
needed a ribosome profiling experiment in some macrophage. A reader should leave with that
distinction, not with "predicted beats measured".

CAVEATS THAT TRAVEL WITH THE PANEL:
  * The real-Ribo-seq arms are ATG-only BY CONSTRUCTION -- RiboCode reports N-terminal extensions
    without testing whether the upstream start is used. Their ncStart is structurally zero, a
    property of the caller and not evidence about non-AUG initiation.
  * Differences in delta-total here are of the same order as the frozen-vs-unfrozen search-parameter
    swing measured on BMDM (+31 vs -21 for the same databases and spectra). Read panel B as
    "indistinguishable from zero for all three arms", not as a ranking.
  * Frozen search parameters, matching P11 and F2b so the three panels sit on one footing.

Novel columns at 1% CLASS-SPECIFIC FDR; canonical columns at 1% global FDR.

cas12a env. Regenerate: `python3 make_macrophage_xsubtype.py`.
"""
from __future__ import annotations

import json
import statistics as st
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/P12_macrophage_xsubtype"
# DEFAULTS ARE THE 30-AA ATTENTION TREE (2026-08-15). Rule 5 puts the tryptic floor at 30 aa --
# macrophages are trypsin/termini=2 -- and the transformer is the model the manuscript emphasises.
# --source / --arm reproduce the superseded 7-aa mamba panel for comparison.
SRCDIR = NEW / "proteogenomics/data/macrophage_tissue/pgx_xsubtype_aa30/reports"
MODEL_ARM = "model_predicted_attn"
MATCHED = "BMDM"          # the population the Ribo-seq databases were built from
# DATABASE SIZES ARE MEASURED AT PLOT TIME, NEVER WRITTEN INTO THE LABEL. The 7-aa version of this
# panel hardcoded "699 seqs" and "2,173 seqs"; at the 30-aa floor the SAME two databases are 506 and
# 1,265, so a fixed label would have silently misreported them when the source tree changed.
ARM_META = [("MODEL", "model (predicted, per population)", "#2b6cb0"),
            ("model_ribocode_bmdm_nt", "real Ribo-seq, BMDM NT", "#c05621"),
            ("model_ribocode_bmdm", "real Ribo-seq, BMDM NT+LPS", "#dd9a63")]


def arms():
    """ARM_META with the model slot resolved to whichever model arm was selected."""
    return [(MODEL_ARM if k == "MODEL" else k, lab, col) for k, lab, col in ARM_META]


def dblab(lab, key, rows, tr):
    """Append the measured DB size. One value if the DB is fixed across populations, else a median.

    This is also the panel's evidence for the asymmetry it is built on: the Ribo-seq arms collapse to
    a single number because one database is reused everywhere, while the model's does not.
    """
    sizes = {r[key]["db"] for r in rows}
    if len(sizes) == 1:
        return f"{lab} ({sizes.pop():,} seqs, fixed)"
    return f"{lab} (median {int(st.median(r[key]['db'] for r in tr)):,} seqs)"


def load():
    rows = []
    for f in sorted(SRCDIR.glob("report_*.json")):
        a = json.loads(f.read_text())["arms"]
        if not all(k in a for k, _, _ in arms()) or "gencode" not in a:
            missing = [k for k, _, _ in arms() if k not in a]
            print(f"  SKIP {f.stem}: missing {missing or ['gencode']}")
            continue
        b = a["gencode"]
        base = b["gencode_peptides"] + b["novel_peptides"]
        r = {"pop": f.stem.replace("report_", ""), "baseline_total": base}
        for k, _, _ in arms():
            d = a[k]
            tot = d["gencode_peptides"] + d["novel_peptides"]
            r[k] = dict(novel=d["novel_peptides"], total=tot, dtot=tot - base,
                        db=d["db_novel_seqs"],
                        density=1000 * d["novel_peptides"] / d["db_novel_seqs"]
                        if d["db_novel_seqs"] else 0.0)
        rows.append(r)
    return rows


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=str(SRCDIR),
                    help="report dir. Default = the 30-aa tree; pass the 7-aa frozen_reports to "
                         "reproduce the superseded mamba panel.")
    ap.add_argument("--arm", default=MODEL_ARM,
                    choices=["model_predicted_attn", "model_predicted"],
                    help="which model arm to plot (attn = transformer, default)")
    a_ = ap.parse_args()
    globals()["SRCDIR"] = Path(a_.source)
    globals()["MODEL_ARM"] = a_.arm
    ARMS = arms()
    if not SRCDIR.exists():
        raise SystemExit(f"missing {SRCDIR}")
    rows = load()
    if not rows:
        raise SystemExit("no population reports found")
    # Matched population LAST, after a visual break, so it is never read as one of the transfer rows.
    rows.sort(key=lambda r: (r["pop"] == MATCHED, r["pop"]))
    pops = [r["pop"] for r in rows]
    x = np.arange(len(rows))
    mi = pops.index(MATCHED) if MATCHED in pops else None
    tr = [r for r in rows if r["pop"] != MATCHED]

    def med(key, field):
        return st.median(r[key][field] for r in tr)

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15.4, 4.5), constrained_layout=True)
    w = 0.27
    for i, (key, lab, col) in enumerate(ARMS):
        off = (i - 1) * w
        axA.bar(x + off, [r[key]["novel"] for r in rows], w, color=col,
                label=dblab(lab, key, rows, tr), zorder=3)
        axB.bar(x + off, [r[key]["dtot"] for r in rows], w, color=col, zorder=3)
        axC.bar(x + off, [r[key]["density"] for r in rows], w, color=col, zorder=3)

    axA.set_ylabel("novel peptides identified", fontsize=8.5)
    axA.set_title("A  Novel peptides (absolute)", fontsize=9, loc="left")
    axA.legend(fontsize=6.6, frameon=False, loc="upper left")
    # Headroom for the legend AND the MATCHED flag; at the default limit both sat on the bars.
    axA.set_ylim(0, max(r[k]["novel"] for r in rows for k, _, _ in ARMS) * 1.42)
    axA.text(0.50, 0.02, f"transfer medians (n={len(tr)}): "
                         + " / ".join(f"{med(k,'novel'):.0f}" for k, _, _ in ARMS),
             transform=axA.transAxes, ha="center", fontsize=7.2, style="italic", color="0.3")

    axB.axhline(0, color="black", lw=1.1, zorder=4)
    axB.set_ylabel("$\\Delta$ TOTAL unique peptides\nvs searching GENCODE alone", fontsize=8.5)
    axB.set_title("B  Cost to the proteome -- all three straddle zero", fontsize=9, loc="left")
    nz = {k: sum(1 for r in tr if r[k]["dtot"] > 0) for k, _, _ in ARMS}
    axB.text(0.5, 0.04,
             "above zero, transfer only:  " + "   ".join(
                 f"{lab.split('(')[0].strip()} {nz[k]}/{len(tr)}" for k, lab, _ in ARMS),
             transform=axB.transAxes, ha="center", fontsize=7.0, fontweight="bold")
    _bv = [r[k]["dtot"] for r in rows for k, _, _ in ARMS]
    axB.set_ylim(min(_bv) * 1.30, max(_bv) * 2.6)   # room for the caveat above, the tally below
    axB.text(0.5, 0.97, "differences here are the size of the search-parameter swing\n"
                        "(BMDM: +31 frozen vs -21 unfrozen, same DBs) -- do not rank on this",
             transform=axB.transAxes, ha="center", va="top", fontsize=6.8, style="italic",
             color="#c05621")

    axC.set_ylabel("novel peptides per 1,000 DB sequences", fontsize=8.5)
    axC.set_title("C  Discovery density -- the REAL Ribo-seq DB is denser", fontsize=9, loc="left")
    axC.set_ylim(0, max(r[k]["density"] for r in rows for k, _, _ in ARMS) * 1.22)
    d_model = med(MODEL_ARM, "density")
    d_ribo = med("model_ribocode_bmdm_nt", "density")
    dens_ratio = d_ribo / d_model if d_model else float("nan")
    axC.text(0.5, 0.97, f"transfer medians: model {d_model:.1f}  vs  "
                        f"real NT {d_ribo:.1f}  ({dens_ratio:.1f}x)",
             transform=axC.transAxes, ha="center", va="top", fontsize=7.6, fontweight="bold",
             color="#c05621")

    for ax in (axA, axB, axC):
        ax.set_xticks(x)
        ax.set_xticklabels(pops, fontsize=6.6, rotation=42, ha="right")
        ax.grid(axis="y", alpha=0.25, zorder=0)
        ax.tick_params(labelsize=7.5)
        if mi is not None:
            # The Ribo-seq DBs came from BMDM, so that column is their home ground, not a transfer
            # test. Separated by a rule and labelled rather than silently averaged in.
            ax.axvline(mi - 0.5, color="0.35", lw=1.1, ls=":", zorder=2)
            ax.get_xticklabels()[mi].set_color("#c05621")
            ax.get_xticklabels()[mi].set_fontweight("bold")
    if mi is not None:
        axA.text(mi, axA.get_ylim()[1] * 0.98, "MATCHED\n(Ribo-seq DB\nis from here)",
                 ha="center", va="top", fontsize=6.2, color="#c05621", fontweight="bold")

    # THE HEADLINE IS DERIVED FROM THE NUMBERS, NOT REPRINTED. At the 7-aa floor a fixed 699-seq
    # Ribo-seq DB beat the model ~3x on density; both the DB size and the direction can change with
    # the floor or the model arm, so a flip is detected and announced rather than silently contradicted
    # by the panel underneath it.
    nt_db = rows[0]["model_ribocode_bmdm_nt"]["db"]
    m_novel, r_novel = med(MODEL_ARM, "novel"), med("model_ribocode_bmdm_nt", "novel")
    if dens_ratio >= 1:
        verdict = (f"a fixed {nt_db:,}-seq Ribo-seq DB reaches {r_novel:.0f} novel peptides vs the "
                   f"model's {m_novel:.0f} and is {dens_ratio:.1f}x denser")
        tail = "the model's case is NOT needing Ribo-seq, not beating it"
    else:
        verdict = (f"the model reaches {m_novel:.0f} novel peptides vs {r_novel:.0f} for a fixed "
                   f"{nt_db:,}-seq Ribo-seq DB and is {1/dens_ratio:.1f}x denser")
        tail = "NOTE: this REVERSES the 7-aa result -- do not reuse the older caption"
    fig.suptitle(f"Mouse macrophages, {len(rows)} populations: a PREDICTED ORF database vs a REAL "
                 f"Ribo-seq one\nOn the {len(tr)} transfer populations {verdict} -- {tail}",
                 fontsize=9.2)
    HERE.mkdir(parents=True, exist_ok=True)
    fig.savefig(HERE / "P12_macrophage_xsubtype.pdf")
    fig.savefig(HERE / "P12_macrophage_xsubtype.png", dpi=200)
    plt.close(fig)

    out = {"source_dir": str(SRCDIR), "search_parameters": "frozen",
           "model_arm": MODEL_ARM, "min_aa": 30 if "aa30" in str(SRCDIR) else 7,
           "density_ratio_ribo_over_model": dens_ratio,
           "n_populations": len(rows), "matched_population": MATCHED,
           "n_transfer": len(tr), "mzml_per_population": 18,
           "arms": {k: lab for k, lab, _ in ARMS}, "populations": rows,
           "transfer_summary": {k: dict(novel_median=med(k, "novel"),
                                        dtot_median=med(k, "dtot"),
                                        dtot_above_zero=nz[k],
                                        density_median=med(k, "density"),
                                        db_seqs=rows[0][k]["db"] if "ribocode" in k else None)
                                for k, _, _ in ARMS},
           "caveats": [
               "Ribo-seq DBs are BMDM-derived and SHARED across all rows; the model's DB is "
               "per-population. BMDM is the matched case, the other 11 are a transfer test.",
               "Real-Ribo-seq arms are ATG-only by construction (RiboCode reports N-terminal "
               "extensions without testing the upstream start), so ncStart is structurally zero.",
               "delta-total differences are the size of the frozen-vs-unfrozen search-parameter "
               "swing measured on BMDM (+31 vs -21); panel B is not a ranking."]}
    (HERE / "P12_values.json").write_text(json.dumps(out, indent=2))

    print(f"  wrote {HERE/'P12_macrophage_xsubtype.pdf'}  ({len(rows)} populations, "
          f"{len(tr)} transfer + {MATCHED} matched)")
    print(f"\n  TRANSFER-ONLY medians (n={len(tr)}):")
    print(f"    {'arm':<40}{'novel':>7}{'dTOTAL':>8}{'>0':>7}{'density':>9}{'DB':>8}")
    for k, lab, _ in ARMS:
        db = rows[0][k]["db"] if "ribocode" in k else int(st.median(r[k]["db"] for r in tr))
        print(f"    {lab:<40}{med(k,'novel'):>7.0f}{med(k,'dtot'):>+8.0f}"
              f"{nz[k]:>4}/{len(tr)}{med(k,'density'):>9.1f}{db:>8,}")
    print(f"\n  {MATCHED} (matched, NOT in the medians): "
          + ", ".join(f"{lab.split('(')[0].strip()} novel="
                      f"{next(r for r in rows if r['pop']==MATCHED)[k]['novel']} "
                      f"dtot={next(r for r in rows if r['pop']==MATCHED)[k]['dtot']:+}"
                      for k, lab, _ in ARMS))


if __name__ == "__main__":
    main()
