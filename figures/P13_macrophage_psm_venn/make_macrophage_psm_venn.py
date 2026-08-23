"""P13 -- per-population Venn of NOVEL PSMs: model-predicted DB vs a real BMDM-NT Ribo-seq DB.

P12 shows the two databases deliver comparable novel-peptide YIELD across 12 macrophage populations.
That leaves the question counts cannot answer: are they finding the SAME spectra? Two arms can each
report 27 novel PSMs and either agree completely or not at all.

One Venn per macrophage population, on NOVEL PSMs at each arm's own 1% class-specific FDR.

THE UNIT IS (spec_id, peptide), NOT spec_id. A scan counts as shared only when BOTH arms assigned it
the SAME peptide. Scans passing in both arms but assigned DIFFERENT peptides are reported under each
circle as `diff-pep` and are deliberately NOT in the intersection: they are rank-1 competition
changing the winner as the search space changes, and folding them in would inflate agreement. A
two-circle Venn cannot draw that category, which is exactly why it is annotated rather than dropped.

EACH ARM CARRIES ITS OWN FDR CUT, and that is correct -- each controls its own novel-class error
rate. It does mean a spectrum can be inside one circle and outside the other purely by threshold.
That is a real difference in what a database delivers at a fixed error rate, not an artifact.

BMDM IS THE MATCHED CASE and is drawn with a coloured frame: the Ribo-seq database was built from
BMDM Ribo-seq, so it is that arm's home ground. The other 11 are the transfer test, and the summary
line reports transfer-only totals.

Counts come from `proteogenomics/scripts/macro_novel_psm_venn.py`, which imports the pipeline's own
`pgx.report.cut` rather than reimplementing FDR, verifies its single-pass loader against
`pgx.score_compare.passing_psms`, and cross-checks every arm against the frozen report's
`novel_psms`.

cas12a env. Regenerate:
  python3 proteogenomics/scripts/macro_novel_psm_venn.py     # ~5 min, writes novel_psm_venn.json
  cd figures/P13_macrophage_psm_venn && python3 make_macrophage_psm_venn.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
HERE = NEW / "figures/P13_macrophage_psm_venn"
# DEFAULT IS THE 30-AA ATTENTION TREE (2026-08-15), matching P11 and P12. --source selects the
# superseded 7-aa mamba venn.
SRC = NEW / "proteogenomics/data/macrophage_tissue/pgx_xsubtype_aa30/novel_psm_venn.json"
MODELCOL = "#2b6cb0"
RIBOCOL = "#c05621"


def ribo_db_seqs(J):
    """Measure the Ribo-seq arm's database size from the reports beside the venn JSON.

    The legend used to assert "699 seqs". That is the 7-aa number; the same database is 506 at the
    30-aa floor, so asserting it would have mislabelled the panel the moment the source tree moved.
    Returns None if the reports are not locatable, and the legend then omits the size rather than
    guessing.
    """
    root = Path(J["search_root"]).parent          # .../pgx_xsubtype*/frozen -> .../pgx_xsubtype*
    for sub in ("reports", "frozen_reports"):
        for f in sorted((root / sub).glob("report_*.json")) if (root / sub).is_dir() else []:
            arm = json.loads(f.read_text())["arms"].get(J["arms"]["B"])
            if arm and arm.get("db_novel_seqs"):
                return arm["db_novel_seqs"]
    return None


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=str(SRC),
                    help="novel_psm_venn.json. Default = the 30-aa attention venn.")
    a_ = ap.parse_args()
    globals()["SRC"] = Path(a_.source)
    if not SRC.exists():
        raise SystemExit(f"missing {SRC}\nRun proteogenomics/scripts/macro_novel_psm_venn.py first.")
    J = json.loads(SRC.read_text())
    rows = J["populations"]
    # Transfer populations first, matched last, so the home-ground case reads as the exception.
    rows.sort(key=lambda r: (r["matched"], r["population"]))
    n = len(rows)
    ncol = 4
    nrow = -(-n // ncol)

    fig, axes = plt.subplots(nrow, ncol, figsize=(3.5 * ncol, 3.7 * nrow))
    axes = axes.ravel()
    for ax, r in zip(axes, rows):
        a_only, b_only, shared = r["a_only"], r["b_only"], r["shared"]
        if a_only + b_only + shared == 0:
            ax.text(0.5, 0.5, "no novel PSMs\nin either arm", ha="center", va="center",
                    fontsize=8, color="0.45", transform=ax.transAxes)
        else:
            v = venn2(subsets=(a_only, b_only, shared), ax=ax,
                      set_colors=(MODELCOL, RIBOCOL), alpha=0.55)
            venn2_circles(subsets=(a_only, b_only, shared), ax=ax, lw=0.8, color="0.35")
            for lab in ("A", "B"):
                t = v.get_label_by_id(lab)
                if t:
                    t.set_fontsize(6.8)
                    t.set_text("")          # set labels drawn once in the figure legend instead
            for rid in ("10", "01", "11"):
                t = v.get_label_by_id(rid)
                if t:
                    t.set_fontsize(9)
                    t.set_fontweight("bold")
            # matplotlib_venn draws NO label for an empty region, so a zero intersection renders as
            # two disjoint circles with nothing between them -- visually identical to "not measured".
            # Zero overlap is the headline finding here, so it is stated explicitly.
            if shared == 0:
                ax.text(0.5, 0.5, "0 shared", transform=ax.transAxes, ha="center", va="center",
                        fontsize=8.5, fontweight="bold", color="#a02020")
        ttl = r["population"] + ("  (MATCHED)" if r["matched"] else "")
        ax.set_title(ttl, fontsize=9.5, fontweight="bold",
                     color=RIBOCOL if r["matched"] else "0.15")
        # The category a 2-circle Venn cannot draw, kept visible.
        ax.text(0.5, -0.01, f"J={r['jaccard']:.2f}   diff-pep={r['same_scan_diff_peptide']}",
                transform=ax.transAxes, ha="center", va="top", fontsize=7.5, color="0.35")
        if r["matched"]:
            # venn2 calls set_axis_off(), so making the spines visible is a no-op -- the earlier
            # version of this silently drew nothing. An explicit unclipped rectangle does render.
            ax.add_patch(plt.Rectangle((-0.02, -0.10), 1.04, 1.16, transform=ax.transAxes,
                                       fill=False, edgecolor=RIBOCOL, lw=1.8, clip_on=False,
                                       zorder=5))
    for ax in axes[n:]:
        ax.axis("off")

    t = J["transfer_totals"]
    tr = [r for r in rows if not r["matched"]]
    shared_pop = sum(1 for r in tr if r["shared"] > 0)
    nseq = ribo_db_seqs(J)
    mlab = "model-predicted database (per population)"
    if "attn" in J["arms"]["A"]:
        mlab += ", attention model"
    fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=10, alpha=0.55, color=MODELCOL,
                                   label=mlab),
                        plt.Line2D([], [], marker="o", ls="", ms=10, alpha=0.55, color=RIBOCOL,
                                   label="real Ribo-seq database, BMDM NT"
                                         + (f" ({nseq:,} seqs, fixed)" if nseq else " (fixed)"))],
               loc="lower center", ncol=2, frameon=False, fontsize=9.5, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle(
        f"Novel PSMs, {J['n_populations']} mouse macrophage populations: a predicted ORF database "
        f"and a real Ribo-seq one find LARGELY DIFFERENT spectra\n"
        f"transfer totals (n={len(tr)}, BMDM excluded): model {t['A_n']:,} PSMs, "
        f"Ribo-seq {t['B_n']:,}, shared {t['shared']:,}  "
        f"({100*t['shared']/max(1, t['A_n']+t['B_n']-t['shared']):.1f}% of the union; "
        f"any overlap at all in {shared_pop}/{len(tr)} populations)",
        fontsize=10.5)
    fig.tight_layout(rect=(0, 0.045, 1, 0.93), h_pad=2.6)
    HERE.mkdir(parents=True, exist_ok=True)
    fig.savefig(HERE / "P13_macrophage_psm_venn.pdf")
    fig.savefig(HERE / "P13_macrophage_psm_venn.png", dpi=200)
    plt.close(fig)

    out = {"source": str(SRC), "arms": J["arms"], "fdr": J["fdr"], "unit": J["unit"],
           "matched_population": J["matched_population"], "n_populations": n,
           "transfer_totals": t, "transfer_union": t["A_n"] + t["B_n"] - t["shared"],
           "transfer_shared_frac_of_union":
               t["shared"] / max(1, t["A_n"] + t["B_n"] - t["shared"]),
           "populations_with_any_overlap": shared_pop, "n_transfer": len(tr),
           "populations": rows}
    (HERE / "P13_values.json").write_text(json.dumps(out, indent=2))

    with open(HERE / "P13_macrophage_psm_venn.tsv", "w") as fh:
        cols = ["population", "matched", "A_n", "B_n", "shared", "a_only", "b_only",
                "same_scan_diff_peptide", "jaccard", "A_cut", "B_cut",
                "shared_peptides", "a_only_peptides", "b_only_peptides"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    print(f"  wrote {HERE/'P13_macrophage_psm_venn.pdf'}  ({n} populations)")
    print(f"  {'population':<18}{'model':>7}{'ribo':>7}{'shared':>8}{'mdl-only':>10}"
          f"{'ribo-only':>11}{'diff-pep':>10}{'J':>7}")
    for r in rows:
        print(f"  {r['population']:<18}{r['A_n']:>7}{r['B_n']:>7}{r['shared']:>8}{r['a_only']:>10}"
              f"{r['b_only']:>11}{r['same_scan_diff_peptide']:>10}{r['jaccard']:>7.3f}"
              + ("  *MATCHED" if r["matched"] else ""))
    print(f"\n  TRANSFER (n={len(tr)}): shared {t['shared']:,} of a "
          f"{t['A_n']+t['B_n']-t['shared']:,} union = "
          f"{100*t['shared']/max(1, t['A_n']+t['B_n']-t['shared']):.1f}%; "
          f"any overlap in {shared_pop}/{len(tr)} populations")


if __name__ == "__main__":
    main()
