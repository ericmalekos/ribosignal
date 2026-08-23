#!/usr/bin/env python3
"""Assemble docs/POSTER_LAYOUT.md -- the 4-column poster layout, with LIVE headline numbers.

The panels existed before this script did; what did not exist was a single document saying which
panel goes where, what each one claims, and which claims are NOT safe to make. That gap is how a
superseded number reaches a printed poster: P11's BMDM-only "+31 total peptides" survived in the
figures README for a week after 12-population data contradicted it.

WHY A GENERATOR RATHER THAN A HAND-WRITTEN DOC. Every headline number here is read from the panel's
own `*_values.json` at build time, the same contract every figure in this project follows. A number
on the poster layout therefore cannot drift from the figure it describes -- if a panel is
regenerated with different data, re-running this script changes the layout doc too. The prose
(captions, caveats, column rationale) is editorial and lives in this file; only the numbers are
pulled.

COLUMNS FOLLOW THE REQUESTED NARRATIVE: methods and training data, then held-out human validation,
then the mouse cross-study, then the proteomics lift for mouse macrophages and human HLA.

Each panel entry carries a MUST-STATE line where the honest version of the claim is narrower than
the eye-catching one. Those are not optional caption garnish; several exist because the broad
version was measured and found false.

cas12a env. Regenerate: `python3 make_poster_layout.py`.
"""
from __future__ import annotations

import json
# NOT `as st`: line ~72 binds a local `st` for the P5 standalone-arm dict, which shadows the module
# inside build() and makes any later st.median() an AttributeError.
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
NEW = HERE.parent
OUT = NEW / "docs/POSTER_LAYOUT.md"


# MODEL-AWARE. The package is built --model attn, but this generator loaded the UNSUFFIXED values
# files, which are mamba4 for every panel that renders both. It therefore printed mamba4 numbers on
# an all-attn poster -- F1 0.831-0.913 and non-canonical 0.54-0.60 (mamba4) where attn is 0.831-0.911
# and 0.52-0.59. A doc edit could not fix that; it regenerated wrong on every refresh.
#
# j() now prefers <stem>_<MODEL>.json when it exists and falls back to the unsuffixed file, so a
# panel that renders one architecture is unaffected while every dual-render panel follows --model.
MODEL = "attn"


def j(rel):
    p = HERE / rel
    alt = p.with_name(f"{p.stem}_{MODEL}{p.suffix}")
    chosen = alt if alt.exists() else p
    if chosen.exists():
        d = json.loads(chosen.read_text())
        if isinstance(d, dict):
            d.setdefault("_values_file", chosen.name)
        return d
    return None


def fmt(v, spec=""):
    return "??" if v is None else format(v, spec)


def build():
    print(f"  model = {MODEL} (values files: *_{MODEL}.json where present)")
    P3 = j("P3_training_data/P3_values.json")
    P5 = j("P5_heldout_human/P5_values.json")
    P9 = j("P9_upset_mouse3x3/P9_values.json")
    P11 = j("P11_bmdm_proteomics/P11_values.json")
    P12 = j("P12_macrophage_xsubtype/P12_values.json")
    P13 = j("P13_macrophage_psm_venn/P13_values.json")
    F2b = j("F2b_discovery_forest/F2b_values.json")
    D12 = j("D12_cpat_cpc2/D12_values.json")
    E2 = j("E2_class_ceiling/E2_values.json")
    QC = j("S_riboseq_qc/riboseq_qc_values.json")
    OC = j("../results/merged_liver_ribocode/overcall/overcall_meta.json")

    # --- P3 ---
    p3 = "??"
    if P3:
        tr = [t for t in P3["tissues"] if not t["holdout"]]
        p3 = (f"{len(tr)} human cell types, {sum(t['n_ribo'] for t in tr)} Ribo-seq libraries; "
              f"Hepatocytes ({sum(t['n_ribo'] for t in P3['tissues'] if t['holdout'])}) held out "
              f"entirely. {P3['universe_tx']:,}-transcript universe. "
              f"{P3['depth_spread_train']:.0f}x depth spread across training tissues.")

    # --- P5 ---
    p5, p5_canon = "??", None
    if P5:
        p5_canon = P5.get("cart_source_is_canonical")
        arms = P5["data"]
        st = {k: v["pred_preddepth"] for k, v in arms.items()}
        f1s = [v["f1"] for v in st.values()]
        an = [v["annot_recall"] for v in st.values()]
        nc = [v["noncanon_recall"] for v in st.values()]
        p5 = (f"{len(arms)} independent human held-outs. Standalone F1 {min(f1s):.3f}-{max(f1s):.3f}. "
              f"Annotated-CDS recall {min(an):.2f}-{max(an):.2f}; "
              f"NON-CANONICAL recall {min(nc):.2f}-{max(nc):.2f}.")

    # --- P9 + the merged-reference measurement ---
    p9 = "??"
    if P9 and OC:
        rows = {(r["kind"], r["source"]): r for r in OC["rows"]}
        ctrl = [r for r in OC["rows"] if r["kind"].startswith("experiment-unique")]
        c = 100 * sum(r["recovered"] for r in ctrl) / max(1, sum(r["n"] for r in ctrl))
        mu = rows.get(("model-unique (question)", "mamba4"))
        an_ = rows.get(("model, experiment-supported (anchor)", "mamba4"))
        both = next((r for r in OC["rows"] if r["kind"] == "model-unique, BOTH models"), None)
        p9 = (f"{mu['n']:,} ORFs called by the model and NO experiment. A "
              f"{OC['n_libraries_merged']}-library merged reference corroborates "
              f"{100*mu['rate']:.1f}% of them, against a {c:.1f}% control and a "
              f"{100*an_['rate']:.1f}% anchor -> ~{100-100*mu['rate']:.0f}% are FALSE POSITIVES. "
              f"Both models agreeing does not rescue them ({100*both['rate']:.1f}%).")

    # --- P11 ---
    p11 = "??"
    if P11:
        p11 = (f"{P11['n_populations']} mouse macrophage populations. Database "
               f"{P11['db_size_ratio_null_over_model']['median']:.0f}x smaller and "
               f"{P11['density_ratio_model_over_null']['median']:.0f}x denser than a naive AUG null. "
               f"On TOTAL unique peptides the model is above the GENCODE-only baseline in "
               f"{P11['model_beats_baseline_n']}/{P11['n_populations']} "
               f"(median {P11['median_model_dtotal']:+.0f}); the null in "
               f"{P11['null_atg_beats_baseline_n']}/{P11['n_populations']} "
               f"(median {P11['median_null_atg_dtotal']:+,.0f}).")

    # --- P12 ---
    p12 = "??"
    if P12:
        t = P12["transfer_summary"]
        # THE MODEL ARM IS READ, NOT ASSUMED. This was `t["model_predicted"]` and became a KeyError
        # the moment P12 was rebuilt on the attention arm -- the generator crashed rather than
        # printing a stale number, but a crash is still a drift failure. Older JSONs lack
        # `model_arm`, hence the fallback.
        marm = P12.get("model_arm", "model_predicted")
        m, n = t[marm], t["model_ribocode_bmdm_nt"]
        p12 = (f"{P12['n_populations']} populations x 3 databases, {P12['mzml_per_population']} mzML "
               f"each. On the {P12['n_transfer']} TRANSFER populations (BMDM excluded -- the Ribo-seq "
               f"databases came from it) a FIXED {n['db_seqs']}-sequence real-Ribo-seq database gives "
               f"{n['novel_median']:.0f} novel peptides vs the model's {m['novel_median']:.0f}, "
               f"{n['dtot_median']:+.0f} vs {m['dtot_median']:+.0f} on TOTAL unique peptides "
               f"({n['dtot_above_zero']}/{P12['n_transfer']} vs "
               f"{m['dtot_above_zero']}/{P12['n_transfer']} above baseline), at "
               f"{n['density_median']/m['density_median']:.1f}x the discovery density.")

    # --- P13 ---
    p13 = "??"
    if P13:
        t = P13["transfer_totals"]
        p13 = (f"Novel PSMs, {P13['n_transfer']} transfer populations: model {t['A_n']:,}, "
               f"real Ribo-seq {t['B_n']:,}, shared only {t['shared']:,} "
               f"({100*P13['transfer_shared_frac_of_union']:.0f}% of the union). "
               f"{t['a_only']:,} model-only vs {t['b_only']:,} Ribo-seq-only PSMs -- the two "
               f"databases are COMPLEMENTARY, not redundant.")

    # --- F2b ---
    f2b = "??"
    if F2b:
        rr = F2b.get("records") or F2b.get("points") or []
        if rr:
            ratios = [r["ratio"] for r in rr]
            po = [r for r in rr if r["arm"] == "model_poisson"]
            npos = sum(1 for r in po if r["d_total_peptides"] > 0)
            # statistics.median, NOT sorted(x)[n//2]. With an even n the latter returns the UPPER
            # middle value, which printed 96x here where the true median (and F2b_values.json and
            # F2b's own doc) is 93.6x. n is 20, so it is always even in practice.
            f2b = (f"{len(rr)} points (5 datasets x 2 models x 2 arms), ALL above 1.0; median "
                   f"{statistics.median(ratios):.1f}x, range {min(ratios):.0f}-{max(ratios):.0f}x. "
                   f"On TOTAL unique peptides the Poisson arm is positive in {npos}/{len(po)}.")

    # --- D12 ---
    d12 = "??"
    if D12:
        t = D12["datasets"]
        pos = sum(1 for d in t if t[d]["arms"].get("model_poisson", {}).get("d_total_peptides", -1) > 0)
        # The assay split is GONE with A549 (2026-08-15). This used to read "CPAT/CPC2 win the
        # tryptic proteome, the model wins all three HLA-I sets"; the tryptic half was A549 and is
        # the half the model LOST. Only the winning substrate remains, so the claim is scoped to it.
        d12 = (f"Against CPAT / CPC2, the standard sequence-only selectors, on "
               f"{len(t)} HLA-I immunopeptidomes. The model wins on this substrate. On TOTAL "
               f"unique peptides the mamba4 Poisson arm is positive in {pos}/{len(t)}. "
               f"HLA-I ONLY -- the tryptic comparison, which the model lost, was removed with A549.")

    # --- E2 ---
    # E2_values.json became an OBJECT (model / n_cells_per_arm / aggregation / classes) so a pooled
    # two-architecture mean could no longer read as a single-model result. It was a bare list before;
    # both shapes are accepted so an older package still renders.
    e2 = "??"
    if E2:
        e2_rows = E2["classes"] if isinstance(E2, dict) else E2
        e2_model = E2.get("model", "both") if isinstance(E2, dict) else "unlabelled"
        a = next((r for r in e2_rows if r["orf_class"] == "annotated"), None)
        n = next((r for r in e2_rows if r["orf_class"] == "novel"), None)
        if a and n:
            # The architecture is NAMED. "both" is a mean over 2 models x 9 RNA/Ribo cells, and
            # printing it unqualified beside single-model panels is what made it misleading.
            who = ("pooled over BOTH architectures" if e2_model == "both"
                   else f"{e2_model} only")
            e2 = (f"Model vs the experiment-vs-experiment CEILING, per ORF class ({who}). "
                  f"Annotated CDS {a['model_pred_preddepth_mean']:.3f} against a ceiling of "
                  f"{a['ceiling_min']:.3f}-{a['ceiling_max']:.3f}; lncRNA "
                  f"{n['model_pred_preddepth_mean']:.3f} against "
                  f"{n['ceiling_min']:.3f}-{n['ceiling_max']:.3f}."
                  + (f" attn alone: annotated {a['model_pred_preddepth_mean_attn']:.3f}, "
                     f"lncRNA {n['model_pred_preddepth_mean_attn']:.3f}."
                     if a.get("model_pred_preddepth_mean_attn") else ""))

    qc = "??"
    if QC:
        rows = QC if isinstance(QC, list) else QC.get("datasets", [])
        # f0_frac is a FRACTION in this JSON, not a percent -- reading it as a percent would
        # print "0.7-0.9% frame 0", i.e. below the 33% floor the panel exists to clear.
        f0 = [100 * r["f0_frac"] for r in rows if isinstance(r, dict) and r.get("f0_frac")]
        if f0:
            qc = (f"Ribo-seq library quality for all {len(rows)} dataset arms: frame-0 fraction "
                  f"{min(f0):.1f}-{max(f0):.1f}% against a 33% no-periodicity floor.")

    COLS = [
        ("1. What it is, and what it was trained on",
         "A poster reader's first two questions. Column 1 answers them before any result is shown.",
         [("arch_attn", "figures/arch_attn/", "Model architecture",
           "One-hot sequence + a 5-channel ORF track + RNA-seq coverage -> per-nucleotide P-site "
           "density. ~5M parameters, CPU-runnable at inference.",
           "The STANDALONE arm uses NO Ribo-seq at inference. Say so explicitly -- 'predicts "
           "ribosome density' is ambiguous about it and that ambiguity is the whole application."),
          ("P3", "figures/P3_training_data/", "Training data", p3,
           "Hepatocytes is drawn in a different colour because it is the LOTO holdout behind every "
           "human number on the poster, not because it was trained on. Brain was dropped for low "
           "periodicity before any pack was built."),
          ("S-QC", "figures/S_riboseq_qc/", "Library quality", qc,
           "Establishes that the input data clears a periodicity floor, so a downstream failure is "
           "the model's and not the library's.")]),

        ("2. Held-out human validation",
         "The central claim. Four held-outs that fail in DIFFERENT ways, which is the point: "
         "same-study-withheld, cross-study, cross-lab cell line, and primary engineered cells.",
         [("P5", "figures/P5_heldout_human/", "Four-arm held-out human -- THE core panel", p5,
           "Panel B is the honest half and must not be cropped. Overall F1 is dominated by "
           "annotated CDS, so a single number reads ~0.9 everywhere and hides the real behaviour. "
           "Both calling arms are always reported."
           + ("" if p5_canon else "  *** CAR-T IS PROVISIONAL: cart_source_is_canonical=false. "
              "Do NOT print until the final-recipe rescore lands (the old pack carries ~29% PCR "
              "duplicates). ***")),
          ("A1", "figures/A1_localization_ceiling/", "Translation localization vs the ceiling",
           "Predicted localization BEATS the observed Ribo-seq periodicity ceiling on all-ORFs.",
           "On NON-CANONICAL ORFs it reaches 94-96% of the ceiling but does NOT beat it. The "
           "caption must say both -- the all-ORF number alone overstates it."),
          ("A2", "figures/A2_ribocode_dropin/", "Drop-in ORF calling",
           "The predicted profile dropped into RiboCode reproduces that dataset's own ORF calls, "
           "no Ribo-seq for the query.",
           "State n_ref and the transcript space. A precision computed against an unstated "
           "transcript space is not interpretable.")]),

        ("3. Mouse cross-study: does it transfer, and where does it fail?",
         "Three independent mouse-liver Ribo-seq experiments on one pipeline and one universe, so a "
         "difference is between EXPERIMENTS and not between code paths. This column is where the "
         "poster is honest about failure modes.",
         [("E1", "figures/E1_rna_provenance/", "RNA-seq provenance",
           "Matched vs mismatched RNA input across the 3x3 factorial.",
           "Reads from results/mouse_liver_3x3_canon/scored -- the CANONICAL alignment. Do not "
           "compare these to any previously reported non-canonical number."),
          ("E2", "figures/E2_class_ceiling/", "Per-class ceiling", e2,
           "The ceiling is experiment-vs-experiment agreement, not 1.0. Scoring the model against "
           "1.0 on non-canonical classes charges it for irreproducibility in the reference."),
          ("E3", "figures/E3_shallow_experiment_win/", "Predicted beats measured, below a depth",
           "The predicted profile calls ORFs better than the measured Ribo-seq below ~29M P-sites.",
           "This is the 'why not just do Ribo-seq' answer, and it is depth-conditional. State the "
           "crossover, not just the win."),
          ("P9", "figures/P9_upset_mouse3x3/", "Where the disagreement lives -- and how much is real",
           p9,
           "The 'model only' bar is a MEASUREMENT now, not a guess, and the control must travel "
           "with it: a genuinely real singly-observed ORF corroborates 74.5% of the time, so '8% "
           "corroborated' is only interpretable beside it. NEVER present cross-architecture "
           "agreement as a confidence filter -- it was tested and it is not one.")]),

        ("4. Proteogenomics: what a better ORF database is worth",
         "The application. Mouse macrophages and human HLA-I, searched with model-selected ORF "
         "databases against naive enumeration.",
         [("P11", "figures/P11_bmdm_proteomics/", "Mouse macrophages, 12 populations", p11,
           "The claim is COST-NEUTRALITY, not a gain. The model is above the GENCODE-only baseline "
           f"in {fmt(P11 and P11['model_beats_baseline_n'])} of "
           f"{fmt(P11 and P11['n_populations'])} populations. An earlier BMDM-only version of this "
           "panel claimed the model was 'the only one that increases total unique peptides'; BMDM "
           "is the best of the 12 and that claim does NOT generalise. Do not print it."),
          ("P12", "figures/P12_macrophage_xsubtype/", "Predicted database vs a REAL Ribo-seq one",
           p12,
           "This is the panel that says what the model is NOT. On the "
           f"{fmt(P12 and P12['n_transfer'])} transfer populations a FIXED "
           f"{fmt(P12 and P12['transfer_summary']['model_ribocode_bmdm_nt']['db_seqs'], ',')}"
           "-sequence database built from BMDM Ribo-seq matches the per-population predicted "
           f"database on yield and is ~{fmt(P12 and P12.get('density_ratio_ribo_over_model'), '.1f')}"
           "x denser. Never caption this as the model beating Ribo-seq. The claim is that it does "
           "not REQUIRE Ribo-seq -- cost and applicability, not accuracy. BMDM is the Ribo-seq "
           "arms' home ground and is excluded from every median."),
          ("P13", "figures/P13_macrophage_psm_venn/", "Do the two databases find the SAME spectra?",
           p13,
           "Pair with P12, never alone. P12 says the predicted database does not out-yield "
           "Ribo-seq; P13 says it finds DIFFERENT spectra, which is what makes it worth running at "
           "all. A scan counts as shared only when both arms assigned the SAME peptide; "
           "same-scan-different-peptide totals "
           f"{fmt(P13 and P13['transfer_totals']['same_scan_diff_peptide'])} across all "
           f"{fmt(P13 and P13['n_populations'])} populations, so the disagreement is about which "
           "spectra clear FDR, not about peptide assignment. Each arm has its own FDR cut, so part "
           "of the non-overlap is a threshold effect -- say so."),
          ("F2b", "figures/F2b_discovery_forest/", "Human HLA-I, 4 datasets x 2 models x 2 arms",
           f2b,
           "Panel (a) is a RATE RATIO; panel (b) is the absolute total and must be shown beside it. "
           "HLA-I IMMUNOPEPTIDOMES ONLY since A549 was removed 2026-08-15 -- do not call this a "
           "whole-proteome result. The Poisson arm is now positive in 8/8 because the only two "
           "failures were the A549 points (attn -35, mamba4 +11). A clean 8/8 is LESS informative "
           "than the old 9/10, not more: the dataset where the arms disagreed is gone."),
          ("D12", "figures/D12_cpat_cpc2/", "Against the honest competitor, HLA-I only", d12,
           "CPAT/CPC2 are the comparison a reviewer asks for; 'better than enumerating every ORF' "
           "is not a strong claim on its own. THE ASSAY SPLIT IS GONE: this used to show CPAT/CPC2 "
           "winning the tryptic proteome and the model winning HLA-I, and the tryptic half (A549) "
           "was removed. Caption it as an HLA-I result, never as 'the model beats CPAT/CPC2'. On "
           "HLA-I the model's edge depends on sub-30-aa ORFs -- 27-33% of its discoveries there come "
           "only from short ORFs, versus 0% for CPAT/CPC2.")]),
    ]

    L = ["# Poster layout -- 48 x 36 in landscape, 4 columns, 15 panels",
         "",
         "GENERATED by `figures/make_poster_layout.py`. Every headline number is read from the "
         "panel's own `*_values.json` at build time, so this document cannot drift from the "
         "figures. Prose is editorial and lives in the generator. Re-run it after regenerating any "
         "panel.",
         "",
         "**Checkpoints: the DEPLOYED models** (`results/loto/*_mm1_holdout_Hepatocytes`), not the "
         "final-recipe retrains, which are a separate open question.",
         "",
         "## How to read the MUST STATE lines",
         "",
         "Every panel below carries one. They are not caption garnish: most exist because the "
         "broader, more eye-catching version of the claim was measured and found false. A poster "
         "that drops them is claiming something this project has evidence against.",
         ""]
    for title, rationale, panels in COLS:
        L += [f"## Column {title}", "", rationale, "",
              "| panel | folder | what it shows |", "|---|---|---|"]
        for tag, folder, short, _head, _must in panels:
            L.append(f"| **{tag}** | `{folder}` | {short} |")
        L.append("")
        for tag, folder, short, head, must in panels:
            L += [f"### {tag} -- {short}", "", head, "",
                  f"**MUST STATE:** {must}", ""]

    L += ["## Numbers deliberately NOT on the poster", "",
          "- Any non-canonical F1 measured before the final-recipe alignment redo. The earlier "
          "reference contained ~22% more novel calls than the corrected pipeline produces, so those "
          "values were inflated and are not comparable to anything here.",
          "- P11's BMDM-only proteomics numbers (154x / 227x / 95x with '+31 total peptides'). "
          "Superseded by the 12-population panel.",
          "- Cross-architecture agreement as evidence for an ORF call. Measured at 8.7% "
          "corroboration, indistinguishable from a single model.",
          "- Any ORF-call precision without its transcript space and n_ref stated.",
          "- \"The predicted database beats Ribo-seq.\" It does not, on the 11 macrophage transfer "
          "populations (P12). The claim is that it does not REQUIRE Ribo-seq.",
          "",
          "## Regenerate", "",
          "```bash",
          "PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3",
          "cd figures && $PY make_poster_layout.py",
          "```", ""]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(L))
    print(f"  wrote {OUT}")
    missing = [t for t, v in [("P3", p3), ("P5", p5), ("P9", p9), ("P11", p11), ("P12", p12), ("P13", p13),
                              ("F2b", f2b), ("D12", d12), ("E2", e2), ("S-QC", qc)] if v == "??"]
    if missing:
        print(f"  UNRESOLVED headline(s) (panel values JSON absent or shape changed): {missing}")
    if p5_canon is False:
        print("  WARNING: P5 CAR-T is still the ORIGINAL pack (cart_source_is_canonical=false)")


if __name__ == "__main__":
    import argparse
    _ap = argparse.ArgumentParser()
    _ap.add_argument("--model", default="attn", choices=["attn", "mamba4"],
                     help="which architecture's values files to prefer (default attn = the poster)")
    MODEL = _ap.parse_args().model
    build()
