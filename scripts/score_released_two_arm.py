#!/usr/bin/env python3
"""Two-arm ORF-call scoring for every RELEASED-model held-out dump.

The standing project rule (feedback_two_arm_orf_calling) is that every ORF-calling test reports BOTH
arms, never just one:

  arm 1  theta=1 deterministic   -- the raw standalone prediction
  arm 2  Poisson at the CDS-anchored theta* -- the calibrated prediction

and the anchor for theta* must be RIBO-SEQ-FREE (annotated-CDS relative recall, `pgx.calibrate
--recall-mode relative`), not a match to an observed between-experiment ceiling. Matching a ceiling
requires the very Ribo-seq the standalone claim says you do not need, so a theta picked that way is
not reachable on a tissue without Ribo-seq.

This exists because the two arms had drifted apart in the record: the theta=1 numbers were re-run on
the released checkpoints while the Poisson numbers were still from a stale run at a theta chosen the
ceiling-matched way. One script, one table, both arms, released models only.

Scored against each dump's own `real_collapsed.txt` -- RiboCode on that pack's real counts -- so the
transcript space is identical on both sides by construction and no expression restriction is needed
(feedback_orf_call_transcript_space; n_ref is printed regardless).

FILTERING IS IMPORTED, NOT REIMPLEMENTED. The call-loading and filtering come from
`compare_dropin_calls.build_loader`, the same code path that produces every other drop-in number in
the project: restricted to the dump's test transcripts, pval <= 0.05, ORF >= 90 nt, and predicted
calls additionally required to reach 0.5x uniform mean density over the ORF. A first version of this
script loaded the collapsed files directly and reported attn/Wang theta=1 at F1 0.725 where the
established path reports 0.867 -- same dump, same arm, different convention. Two tables that disagree
by 0.14 on the same quantity is worse than one table.

Usage:
  score_released_two_arm.py [--key genomic|transcript] [--out <json>]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_dropin_calls import build_loader  # noqa: E402

NONCANON = {"uORF", "Overlap_uORF", "dORF", "Overlap_dORF", "novel", "internal"}

# (label, dump dir, species). Every one is on *_union_noBrain_nokozak_mm1_*.
DUMPS = [
    ("attn/Wang",        "results/liver_released/attn_mouse_wang_liver", "mouse"),
    ("mamba4/Wang",      "results/liver_released/mamba4_mouse_wang_liver", "mouse"),
    ("attn/Janich",      "results/liver_released/attn_mouse_janich_liver_decon", "mouse"),
    ("mamba4/Janich",    "results/liver_released/mamba4_mouse_janich_liver_decon", "mouse"),
    ("attn/GSE243134",   "results/liver_released/attn_mouse_gse243134_liver", "mouse"),
    ("mamba4/GSE243134", "results/liver_released/mamba4_mouse_gse243134_liver", "mouse"),
    ("attn/Ruiz-Orera",  "results/heldout/human_ruizorera/released_attn/dropin", "human"),
    ("mamba4/Ruiz-Orera", "results/heldout/human_ruizorera/released_mamba4/dropin", "human"),
]
TX2GENE = {"mouse": ROOT / "data/mouse_tx2biotype.tsv",
           "human": ROOT / "data/tx2biotype.tsv"}
def prf(pred, ref, subset=None):
    """P/R/F1 over call keys, optionally restricted to a set of ORF_type values.

    build_loader's records are dicts ({"type", "length", "adj", "logp", "f0"}), so the class filter
    reads r["type"] -- not the record itself.
    """
    if subset is not None:
        pred = {k: v for k, v in pred.items() if v["type"] in subset}
        ref = {k: v for k, v in ref.items() if v["type"] in subset}
    m = len(set(pred) & set(ref))
    p = m / len(pred) if pred else 0.0
    r = m / len(ref) if ref else 0.0
    return {"n_pred": len(pred), "n_ref": len(ref), "n_match": m, "precision": p, "recall": r,
            "f1": 2 * p * r / (p + r) if p + r else 0.0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--key", default="genomic", choices=["genomic", "transcript"])
    ap.add_argument("--out", default=str(ROOT / "results/released_two_arm_orf_calls.json"))
    a = ap.parse_args()

    rows, missing = [], []
    for label, d, sp in DUMPS:
        dd = ROOT / d
        real = dd / "real_collapsed.txt"
        theta1 = dd / "pred_preddepth_collapsed.txt"
        prof = dd / "pred_profiles.npz"
        # calib/ sits INSIDE the dump dir for the liver runs (results/liver_released/<name>/calib) but
        # BESIDE it for Ruiz-Orera, whose dump dir is .../<name>/dropin while calibrate wrote
        # .../<name>/calib. Check both rather than encode the inconsistency as two DUMPS columns --
        # the earlier version only checked the first and silently reported "Poisson arm not run" for
        # two rows whose calibration had in fact completed.
        tj = next((p for p in (dd / "calib/theta.json", dd.parent / "calib/theta.json")
                   if p.exists()), dd / "calib/theta.json")
        if not (real.exists() and theta1.exists() and prof.exists()):
            missing.append(f"{label}: no real/theta1 collapsed or pred_profiles.npz in {d}")
            continue
        lc, keep_tx, _ = build_loader(prof, key=a.key, tx2gene=str(TX2GENE[sp]))
        R = lc(real)
        arms = {"theta=1": lc(theta1, is_pred=True)}
        theta_star = None
        if tj.exists():
            j = json.loads(tj.read_text())
            theta_star = j["theta"]
            cp = Path(j["collapsed"])
            if cp.exists():
                arms[f"poisson th={theta_star}"] = lc(cp, is_pred=True)
            else:
                missing.append(f"{label}: theta.json points at a missing {cp}")
        else:
            missing.append(f"{label}: no calib/theta.json (Poisson arm not run)")
        for arm, P in arms.items():
            rows.append({"dump": label, "arm": arm, "theta_star": theta_star, "species": sp,
                         "overall": prf(P, R), "non-canonical": prf(P, R, NONCANON),
                         "CDS": prf(P, R, {"annotated"})})

    w = f"{'dump':<20}{'arm':<18}{'nPred':>8}{'nRef':>8}{'P':>7}{'R':>7}{'F1':>7} | " \
        f"{'ncP':>6}{'ncR':>6}{'ncF1':>6}"
    print(w)
    print("-" * len(w))
    prev = None
    for r in rows:
        if prev and prev != r["dump"]:
            print()
        prev = r["dump"]
        o, n = r["overall"], r["non-canonical"]
        print(f"{r['dump']:<20}{r['arm']:<18}{o['n_pred']:>8,}{o['n_ref']:>8,}"
              f"{o['precision']:>7.3f}{o['recall']:>7.3f}{o['f1']:>7.3f} | "
              f"{n['precision']:>6.3f}{n['recall']:>6.3f}{n['f1']:>6.3f}")
    if missing:
        print("\nMISSING:")
        for m in missing:
            print(f"  {m}")
    Path(a.out).write_text(json.dumps({"key": a.key, "rows": rows, "missing": missing}, indent=2))
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
