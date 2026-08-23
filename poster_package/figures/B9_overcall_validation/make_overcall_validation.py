#!/usr/bin/env python3
"""B9 -- are the model's experiment-unsupported ORF calls false positives, or deep discoveries?

The mouse-liver UpSet (P9) shows 1,942 ORFs called by the model and by NO single Ribo-seq experiment.
Its bar is annotated with the answer, but the MEASUREMENT that produced it -- four strata and a
gene-coverage split -- had no figure. This is that figure, and it is a rigor panel: its subject is a
failure mode, not a capability.

THE ARBITER. `results/merged_liver_ribocode` is a JOINT RiboCode call over all 28 canonically
aligned mouse-liver libraries (Janich 5 + GSE243134 21 + Wang 2), 26,388 raw calls, +77% to +102%
over any single dataset. Not a union of three call sets: 28 per-library P-site offsets feed ONE
periodicity test on the pooled signal, so an ORF with weak but consistent signal everywhere can clear
significance here while failing in each dataset alone.

PANEL A IS THE WHOLE ARGUMENT, AND IT ONLY WORKS BECAUSE OF THE OTHER TWO BARS. A recovery rate for
the model's extras means nothing alone -- a deeper reference recovers more of everything. The two
reference strata are what make 7.7% interpretable:

  ANCHOR   model calls that ARE experiment-supported. 99.3%. Proves the merged reference is not the
           limitation: it finds essentially every model call some experiment saw.
  CONTROL  observed calls made by exactly ONE experiment and missed by the other two. 74.5%. These
           are known-real and equally unsupported by the rest of the 3x3, so they measure how often
           "unsupported by the other arms" merely means "too deep for them".

Against those, the model's extras at 7.7% are a ~10x shortfall. They are false positives.

PANEL B: TWO ARCHITECTURES AGREEING IS NOT EVIDENCE. The 5-set UpSet family was built expecting the
shared bar to identify real-but-shallow ORFs. It does not: calls made by BOTH mamba4 and attn and no
experiment are corroborated at 8.7%, indistinguishable from either alone. Cross-architecture
consensus is a shared systematic error and must never be used as a confidence filter.

PANEL C CLOSES THE "NO DATA THERE" OBJECTION. A missed call is only FP evidence if the merged
reference HAD signal at that locus. Splitting on whether merged called any ORF on the same GENE
separates "the specific ORF was rejected on its merits" (881 calls, still only 16.9%) from "the gene
shows no translation at all" (1,061 calls). Both halves say false positive, and since the universe is
already salmon TPM >= 1, the second half is RNA-expressed transcripts with no ribosome signal --
exactly the lncRNA over-call mode. The 0% recovery in that stratum is DEFINITIONAL (recovery requires
a merged call on the gene), so the panel plots its SIZE, never its rate.

Per-class rates are shown because lncRNA and uORF fail differently, and because `internal` is the one
class where the CONTROL is itself weak (33-54%) -- the one class where a low model number is
genuinely ambiguous.

cas12a env. Regenerate:
  python3 scripts/score_merged_liver_overcall.py
  cd figures/B9_overcall_validation && python3 make_overcall_validation.py
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
HERE = NEW / "figures/B9_overcall_validation"
OC = NEW / "results/merged_liver_ribocode/overcall"
ANCHOR, CONTROL, MODEL, SHARED = "#2f855a", "#2b6cb0", "#c05621", "#9b2c2c"
CLASSES = [("annotated", "annotated\nCDS"), ("novel", "lncRNA\n(novel)"), ("uORF", "uORF"),
           ("dORF", "dORF"), ("internal", "internal")]


def main():
    mf = OC / "overcall_meta.json"
    if not mf.exists():
        raise SystemExit(f"missing {mf}\nRun scripts/score_merged_liver_overcall.py first.")
    M = json.loads(mf.read_text())
    rows = {(r["kind"], r["source"]): r for r in M["rows"]}
    ctrl = [r for r in M["rows"] if r["kind"].startswith("experiment-unique")]
    c_rate = sum(r["recovered"] for r in ctrl) / max(1, sum(r["n"] for r in ctrl))
    c_n = sum(r["n"] for r in ctrl)

    by_class = {}
    for line in (OC / "overcall_recovery_by_class.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        by_class.setdefault(f[0], {})[f[1]] = dict(n=int(f[2]), rec=int(f[3]), rate=float(f[4]))

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(14.6, 4.5), constrained_layout=True)

    # A: the four strata
    bars = [("ANCHOR\nmodel calls that ARE\nexperiment-supported",
             rows[("model, experiment-supported (anchor)", "mamba4")], ANCHOR),
            ("CONTROL\nobserved calls unique\nto ONE experiment",
             dict(n=c_n, rate=c_rate), CONTROL),
            ("THE QUESTION\nmodel calls NO\nexperiment supports",
             rows[("model-unique (question)", "mamba4")], MODEL)]
    x = np.arange(len(bars))
    v = [100 * b[1]["rate"] for b in bars]
    axA.bar(x, v, 0.62, color=[b[2] for b in bars], zorder=3)
    for xi, (lab, d, _c) in zip(x, bars):
        axA.annotate(f"{100*d['rate']:.1f}%", (xi, 100 * d["rate"]), xytext=(0, 3),
                     textcoords="offset points", ha="center", fontsize=11, fontweight="bold")
        axA.annotate(f"n={d['n']:,}", (xi, 2), ha="center", fontsize=7.5, color="white",
                     fontweight="bold")
    axA.axhline(100 * c_rate, color=CONTROL, lw=1.2, ls="--", zorder=4)
    axA.set_xticks(x)
    axA.set_xticklabels([b[0] for b in bars], fontsize=7.8)
    axA.set_ylabel("corroborated by the 28-library merged reference (%)", fontsize=8.5)
    axA.set_title("A  A ~10x shortfall against the control", fontsize=9.5, loc="left")
    axA.set_ylim(0, 118)
    axA.grid(axis="y", alpha=0.25, zorder=0)
    axA.text(0.5, 0.90, "the extras are FALSE POSITIVES", transform=axA.transAxes, ha="center",
             fontsize=9.5, fontweight="bold", color=MODEL)

    # B: cross-architecture agreement does not rescue them
    b2 = [("mamba4\nalone", rows[("model-unique (question)", "mamba4")], MODEL),
          ("attn\nalone", rows[("model-unique (question)", "attn")], MODEL),
          ("BOTH models,\nno experiment", rows[("model-unique, BOTH models", "mamba4+attn")], SHARED)]
    x2 = np.arange(len(b2))
    axB.bar(x2, [100 * b[1]["rate"] for b in b2], 0.6, color=[b[2] for b in b2], zorder=3)
    for xi, (lab, d, _c) in zip(x2, b2):
        axB.annotate(f"{100*d['rate']:.1f}%\nn={d['n']:,}", (xi, 100 * d["rate"]), xytext=(0, 3),
                     textcoords="offset points", ha="center", fontsize=9, fontweight="bold")
    axB.axhline(100 * c_rate, color=CONTROL, lw=1.4, ls="--", zorder=4)
    axB.text(len(b2) - 0.45, 100 * c_rate + 1.5, f"CONTROL {100*c_rate:.1f}%", ha="right",
             fontsize=8, color=CONTROL, fontweight="bold")
    axB.set_xticks(x2)
    axB.set_xticklabels([b[0] for b in b2], fontsize=8)
    axB.set_ylabel("corroborated (%)", fontsize=8.5)
    axB.set_title("B  Two architectures agreeing is NOT evidence", fontsize=9.5, loc="left")
    axB.set_ylim(0, 100 * c_rate * 1.45)
    axB.grid(axis="y", alpha=0.25, zorder=0)
    axB.text(0.5, 0.62, "a shared systematic error --\nnever use as a confidence filter",
             transform=axB.transAxes, ha="center", fontsize=8.2, style="italic", color=SHARED)

    # C: per class, control vs model-unique
    xs = np.arange(len(CLASSES))
    w = 0.38
    cg = {c: [by_class[f"obs_unique_{d}"][c]["rate"] for d in ("janich", "gse243134", "wang")
              if c in by_class.get(f"obs_unique_{d}", {})] for c, _ in CLASSES}
    mg = by_class.get("model_unique_mamba4", {})
    axC.bar(xs - w / 2, [100 * (sum(cg[c]) / len(cg[c])) if cg[c] else np.nan for c, _ in CLASSES],
            w, color=CONTROL, label="experiment-unique (CONTROL)", zorder=3)
    axC.bar(xs + w / 2, [100 * mg[c]["rate"] if c in mg else np.nan for c, _ in CLASSES],
            w, color=MODEL, label="model-unique, mamba4", zorder=3)
    for i, (c, _) in enumerate(CLASSES):
        if c in mg:
            axC.annotate(f"{100*mg[c]['rate']:.0f}%", (i + w / 2, 100 * mg[c]["rate"]),
                         xytext=(0, 2), textcoords="offset points", ha="center", fontsize=7)
    axC.set_xticks(xs)
    axC.set_xticklabels([lab for _, lab in CLASSES], fontsize=7.8)
    axC.set_ylabel("corroborated (%)", fontsize=8.5)
    axC.set_title("C  The gap holds in EVERY class, including canonical CDS", fontsize=9.5,
                  loc="left")
    # Dropped below the caveat text, which occupies the top strip of this panel.
    axC.legend(fontsize=7.5, frameon=False, loc="upper right",
               bbox_to_anchor=(1.0, 0.84))
    axC.grid(axis="y", alpha=0.25, zorder=0)
    axC.set_ylim(0, 138)   # headroom for the caveat, which has no clear space among the bars
    # No region among the bars is free at this many classes, so the axis gets explicit headroom
    # and the caveat sits above everything. A bottom strip ran straight across the bars.
    axC.text(0.02, 0.985, "`internal` is the one class where the CONTROL is also weak\n"
                          "-- the one place a low model number is ambiguous",
             transform=axC.transAxes, va="top", fontsize=7.2, style="italic", color="0.35")

    on = M["detail"]["model_unique_mamba4_gene_covered"]
    off = M["detail"]["model_unique_mamba4_gene_uncovered"]
    fig.suptitle(
        f"The model's {rows[('model-unique (question)','mamba4')]['n']:,} experiment-unsupported "
        f"mouse-liver ORF calls are ~{100-100*rows[('model-unique (question)','mamba4')]['rate']:.0f}% "
        f"FALSE POSITIVES, not depth-limited discoveries\n"
        f"and it is not a coverage artifact: {off['n']:,} sit on genes with NO detectable translation "
        f"in 28 pooled libraries; of the {on['n']:,} on demonstrably translated genes, only "
        f"{100*on['rate']:.1f}% are corroborated", fontsize=9.4)
    HERE.mkdir(parents=True, exist_ok=True)
    fig.savefig(HERE / "B9_overcall_validation.pdf")
    fig.savefig(HERE / "B9_overcall_validation.png", dpi=200)
    plt.close(fig)

    out = {# CONVENTION D1: name the architecture, never leave it to be inferred. This panel is NOT
           # single-model: (a) is the mamba4 anchor/control/model-unique triplet, (b) is the
           # cross-model comparison naming mamba4, attn and both-models. A bare "model" key would
           # misrepresent that, so it says so explicitly.
           "model": "mamba4 (panel a) + attn and both-models (panel b)",
           "models_present": ["mamba4", "attn", "mamba4+attn"],
           "panel_a_model": "mamba4",
           "source": str(OC), "merged_libraries": M["n_libraries_merged"],
           "merged_calls_in_space": M["merged_calls_in_space"],
           "transcript_space": M["transcript_space"], "arm": M["arm"],
           "control_rate": c_rate, "control_n": c_n,
           "strata": {k[0] + " / " + k[1]: v for k, v in rows.items()},
           "gene_covered": on, "gene_uncovered": off,
           "by_class": by_class}
    (HERE / "B9_values.json").write_text(json.dumps(out, indent=2))
    print(f"  wrote {HERE/'B9_overcall_validation.pdf'}")
    print(f"  anchor {100*rows[('model, experiment-supported (anchor)','mamba4')]['rate']:.1f}%  "
          f"control {100*c_rate:.1f}%  model-unique "
          f"{100*rows[('model-unique (question)','mamba4')]['rate']:.1f}%  "
          f"both-models {100*rows[('model-unique, BOTH models','mamba4+attn')]['rate']:.1f}%")
    print(f"  gene-covered {on['n']:,} at {100*on['rate']:.1f}%; "
          f"gene-uncovered {off['n']:,} (rate is 0 by construction -- read n)")


if __name__ == "__main__":
    main()
