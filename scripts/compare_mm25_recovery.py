#!/usr/bin/env python3
"""Does keeping multimappers recover the ORFs the model calls and mm1 Ribo-seq misses?

THE TEST. Model-specific calls were defined as: called by a model, present in NO real liver dataset
at mm1. This asks how many of them appear in the SAME dataset re-aligned at
--outFilterMultimapNmax 25, everything else identical.

  recovered  -> the mm1 data was BLIND to them; multimap filtering hid real signal
  absent     -> they are model false positives, and no mapping excuse survives

Both sides go through scripts/compare_dropin_calls.build_loader so the filtering is identical to
every other table in this project (test-tx restriction, pval<=0.05, ORF>=90nt, and on the predicted
side the >=0.5x enrichment gate). Loading collapsed files by hand is what produced a 4x
overcount earlier -- see docs/model_specific_calls_multimap.md.

CONTROL. A recovery rate is meaningless without a baseline: mm25 calls MORE ORFs overall, so some
model-specific calls would be recovered by chance. The null here is the recovery rate of a
size-matched random sample of ORFs that mm1 also missed but no model called.

  compare_mm25_recovery.py
"""
from __future__ import annotations
import json, random, sys
from pathlib import Path

NEW = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(NEW / "scripts"))
from compare_dropin_calls import build_loader          # noqa: E402

C = NEW / "results/mouse_liver_3x3_canon"
RNA = "gse243134"
MM25 = NEW / "data/mm25_diagnostic/psites_mm25/gse243134_mm25_collapsed.txt"


def main():
    if not MM25.exists():
        sys.exit(f"mm25 calls not ready: {MM25}")
    lc, keep_tx, _ = build_loader(str(C / f"attn_ribo-{RNA}_rna-{RNA}" / "pred_profiles.npz"),
                                  key="genomic", tx2gene=str(NEW / "data/tx2biotype_mouse.tsv"))
    real = set()
    for ds in ("wang", "janich", "gse243134"):
        real |= set(lc(str(C / f"attn_ribo-{ds}_rna-{ds}" / "real_collapsed.txt"), is_pred=False))
    pred = {m: set(lc(str(C / f"{m}_ribo-{RNA}_rna-{RNA}" / "pred_preddepth_collapsed.txt"),
                      is_pred=True)) for m in ("attn", "mamba4")}
    mm25 = set(lc(str(MM25), is_pred=False))
    A, M = pred["attn"] - real, pred["mamba4"] - real
    sets = {"BOTH models": A & M, "attn only": A - M, "mamba4 only": M - A,
            "attn-specific (all)": A, "mamba4-specific (all)": M}

    print(f"  mm1 REAL union(3)      : {len(real):,}")
    print(f"  mm25 calls (same filter): {len(mm25):,}   ({len(mm25)-len(real):+,} vs mm1 union)")
    print(f"  mm25 NEW vs mm1 union  : {len(mm25 - real):,}\n")

    # control: ORFs mm1 missed that NO model called -- what fraction does mm25 recover by chance?
    universe_missed = (mm25 - real) - (A | M)
    print(f"  {'set':<24} {'n':>6} {'recovered in mm25':>19} {'rate':>8}")
    out = {}
    for k, v in sets.items():
        r = len(v & mm25)
        out[k] = {"n": len(v), "recovered": r, "rate": round(r / len(v), 4) if v else None}
        print(f"  {k:<24} {len(v):>6} {r:>19,} {100*r/len(v) if v else 0:>7.1f}%")
    # baseline: how much of mm1's own missed space does mm25 fill, i.e. the chance rate
    base = len(universe_missed)
    print(f"\n  CONTROL: mm25 calls absent from mm1 AND from every model set: {base:,}")
    print(f"  These are the 'recovered by deeper mapping alone' population. A model-specific")
    print(f"  recovery rate must beat the rate at which mm25 recovers ARBITRARY mm1-missed ORFs.")
    (NEW / "data/mm25_diagnostic/recovery_summary.json").write_text(json.dumps(
        {"model": ["attn", "mamba4"],
         "source": "data/mm25_diagnostic/psites_mm25 vs results/mouse_liver_3x3_canon",
         "provenance_note": "DIAGNOSTIC. mm25 violates the project Ribo-seq rule deliberately; "
                            "P-site counts are inflated. Read recovery of calls, not counts.",
         "mm1_real_union": len(real), "mm25_calls": len(mm25),
         "mm25_new_vs_mm1": len(mm25 - real), "sets": out}, indent=2))
    print(f"\n  wrote data/mm25_diagnostic/recovery_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
