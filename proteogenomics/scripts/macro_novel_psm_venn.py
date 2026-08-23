#!/usr/bin/env python3
"""Per-population overlap of NOVEL PSMs: model-predicted database vs a real BMDM-NT Ribo-seq one.

P12 compares the two databases on summary counts. This asks whether they are finding the SAME
spectra: a model arm and a Ribo-seq arm can each report 27 novel PSMs and either agree completely or
not at all, and the counts cannot tell those apart.

THE UNIT IS A SPECTRUM, AND MATCHING IS ON (spec_id, peptide). Two arms "share" a PSM only when the
same scan was assigned the SAME peptide in both. Scans that pass in both arms but were assigned
DIFFERENT peptides are counted separately and are NOT overlap -- they are rank-1 competition
changing the winner as the search space changes, and folding them into the intersection would
inflate agreement. This is the same pairing `pgx.score_compare` uses for its matched-spectra test.

FDR IS THE PIPELINE'S OWN. Each arm gets its own 1% CLASS-SPECIFIC cut (novel targets vs
`REV_nuORF|` decoys only) via `pgx.report.cut`, imported rather than reimplemented. Two arms have
DIFFERENT cuts, which is correct -- each controls its own novel-class FDR -- but it means a spectrum
can be "in" one arm and out of the other purely by threshold. That is a real difference in what the
database delivers at a fixed error rate, not an artifact, and it is what the Venn is measuring.

SINGLE PASS, VERIFIED. `load_rank1_spec` is read once per arm and reshaped for `cut()`, which reads
only (score, class) and so is unaffected. The first population is ALSO scored through the pipeline's
own `passing_psms` and the two are asserted equal, so the fast path cannot silently diverge. Every
arm's PSM count is cross-checked against the frozen report's `novel_psms`.

  python3 proteogenomics/scripts/macro_novel_psm_venn.py
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from pgx.report import cut, load_rank1, load_rank1_spec  # noqa: E402
from pgx.score_compare import passing_psms  # noqa: E402

NEW = HERE.parents[1]
# DEFAULTS ARE THE 30-AA ATTENTION TREE (2026-08-15): Rule 5 puts the tryptic floor at 30 aa and the
# transformer is the emphasised model. --xs / --a-arm select the superseded 7-aa mamba tree.
XS = NEW / "proteogenomics/data/macrophage_tissue/pgx_xsubtype_aa30"
SEARCH = XS / "frozen"
REPORTS = XS / "reports"
A_ARM = "model_predicted_attn"
B_ARM = "model_ribocode_bmdm_nt"
MATCHED = "BMDM"        # the population the Ribo-seq database was built from
FDR = 0.01


def novel_psms(d, fdr=FDR):
    """-> ({spec_id: peptide}, cut). Novel PSMs above this arm's own class-specific FDR."""
    rows = load_rank1_spec(d)
    # cut() reads only positions 1 and 2 of a load_rank1 4-tuple (score, class); the peptide and
    # novel-accession slots are ignored, so this reshape is exact rather than an approximation.
    c = cut([(pep, h, k, ()) for (_s, pep, h, k) in rows], {"novel_t"}, {"novel_d"}, fdr)
    if c is None:
        return {}, None
    return {s: pep for (s, pep, h, k) in rows if k == "novel_t" and h >= c}, c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(XS / "novel_psm_venn.json"))
    ap.add_argument("--fdr", type=float, default=FDR)
    ap.add_argument("--xs", default=None, help="pgx_xsubtype tree root (default: the 30-aa one)")
    ap.add_argument("--a-arm", default=None, help="model arm for circle A (default: attn)")
    a = ap.parse_args()
    if a.xs:
        globals()["XS"] = pathlib.Path(a.xs)
        globals()["SEARCH"] = XS / "frozen"
        globals()["REPORTS"] = XS / ("reports" if (XS / "reports").exists() else "frozen_reports")
    if a.a_arm:
        globals()["A_ARM"] = a.a_arm

    pops = sorted(p.name for p in SEARCH.iterdir() if p.is_dir())
    print(f"{len(pops)} populations: {A_ARM} vs {B_ARM}, 1% class-specific FDR\n")
    out, verified = [], False
    for pop in pops:
        rec = {"population": pop, "matched": pop == MATCHED}
        sets = {}
        for tag, arm in (("A", A_ARM), ("B", B_ARM)):
            d = SEARCH / pop / arm
            if not (d / "rank1.tsv.gz").exists() and not list(d.glob("*.tsv")):
                print(f"  SKIP {pop}/{arm}: no search output")
                sets = None
                break
            s, c = novel_psms(d, a.fdr)
            sets[tag] = s
            rec[f"{tag}_arm"], rec[f"{tag}_n"], rec[f"{tag}_cut"] = arm, len(s), c
            # Guard the single-pass reshape against the pipeline's own path, once.
            if not verified:
                P, C = passing_psms(load_rank1(d), a.fdr)
                assert len(P) == len(s) and (C is None) == (c is None) and (C is None or abs(C - c) < 1e-9), \
                    f"fast path diverged on {pop}/{arm}: {len(P)} vs {len(s)}, cut {C} vs {c}"
                print(f"  [verified] {pop}/{arm}: fast path == pgx.passing_psms "
                      f"({len(P)} PSMs, cut {C:.2f})")
                verified = True
            # Cross-check against the frozen report the figures read.
            rp = REPORTS / f"report_{pop}.json"
            if rp.exists():
                want = json.loads(rp.read_text())["arms"].get(arm, {}).get("novel_psms")
                if want is not None and want != len(s):
                    print(f"  WARNING {pop}/{arm}: {len(s)} PSMs here vs {want} in the frozen report")
                rec[f"{tag}_report_novel_psms"] = want
        if sets is None:
            continue
        A, B = sets["A"], sets["B"]
        both_scans = set(A) & set(B)
        same = {s for s in both_scans if A[s] == B[s]}
        diff = both_scans - same
        rec.update(shared=len(same), a_only=len(A) - len(both_scans), b_only=len(B) - len(both_scans),
                   same_scan_diff_peptide=len(diff),
                   union=len(A) + len(B) - len(same),
                   jaccard=len(same) / (len(A) + len(B) - len(same)) if (A or B) else 0.0,
                   shared_peptides=len({A[s] for s in same}),
                   a_only_peptides=len({A[s] for s in A if s not in same}),
                   b_only_peptides=len({B[s] for s in B if s not in same}))
        out.append(rec)
        flag = "  *MATCHED" if rec["matched"] else ""
        print(f"  {pop:<17} model={rec['A_n']:>4}  ribo={rec['B_n']:>4}  "
              f"shared={rec['shared']:>4}  model-only={rec['a_only']:>4}  ribo-only={rec['b_only']:>4}"
              f"  diff-pep={rec['same_scan_diff_peptide']:>3}  J={rec['jaccard']:.3f}{flag}")

    tr = [r for r in out if not r["matched"]]
    payload = {"arms": {"A": A_ARM, "B": B_ARM}, "fdr": a.fdr, "unit": "(spec_id, peptide)",
               "matched_population": MATCHED, "n_populations": len(out),
               "search_root": str(SEARCH), "populations": out,
               "transfer_totals": {k: sum(r[k] for r in tr)
                                   for k in ("A_n", "B_n", "shared", "a_only", "b_only",
                                             "same_scan_diff_peptide")}}
    pathlib.Path(a.out).write_text(json.dumps(payload, indent=2))
    t = payload["transfer_totals"]
    print(f"\n  TRANSFER TOTALS (n={len(tr)}, {MATCHED} excluded): model {t['A_n']:,} PSMs, "
          f"ribo {t['B_n']:,}, shared {t['shared']:,}, model-only {t['a_only']:,}, "
          f"ribo-only {t['b_only']:,}, same-scan-different-peptide {t['same_scan_diff_peptide']:,}")
    print(f"  wrote {a.out}")


if __name__ == "__main__":
    main()
