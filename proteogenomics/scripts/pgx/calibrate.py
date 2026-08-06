#!/usr/bin/env python3
"""Poisson / theta CDS-anchored calibration for the `pgx` pipeline (implements Check 6).

theta scales the PREDICTED effective depth fed to RiboCode (`ribocode_dropin.py --pred_scale`).
Lower theta = shallower simulated experiment = fewer, higher-confidence calls. Check 4 established
Poisson injection as the winning lever (+0.10 novel precision at matched CDS recall), so the
calibrated arm samples Poisson(rate) rather than rounding.

ANCHOR: recall of ANNOTATED CDS among expressed transcripts. This is deliberately self-contained,
because the model predicts translation from sequence + RNA-seq alone; requiring a real Ribo-seq call
set to calibrate would defeat the purpose. Pass --ref-calls to anchor on an observed RiboCode run
instead when one exists.

SELECTION: the SMALLEST theta whose CDS recall still meets --cds-recall. Tightening trades novel
yield for novel precision monotonically (Check 3: novel precision 0.25 -> 0.71 as theta falls), so
the most stringent point that clears the recall floor is the right operating point.

Resumable + array-friendly: each theta writes its own directory and is skipped if already complete.

  # fan out one theta per array task
  python -m pgx.calibrate --profiles p.npz --annot A --species mouse --out S --only-theta 0.05
  # then pick
  python -m pgx.calibrate --profiles p.npz --annot A --species mouse --out S --pick-only
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from . import rc_io
from .refs import species_refs, tool_path

DEFAULT_THETAS = (0.02, 0.05, 0.10, 0.20, 0.35, 0.60, 1.00)
VARIANT = "pred_preddepth"


def theta_dir(out: Path, theta: float) -> Path:
    return Path(out) / f"theta_{theta:g}"


def collapsed_path(out: Path, theta: float) -> Path:
    return theta_dir(out, theta) / f"{VARIANT}_collapsed.txt"


def run_theta(theta, profiles, annot, out, min_aa=5, pval=0.05, seed=0, poisson=True,
              force=False, ribocode_python=None, dropin=None):
    """Run one theta through ribocode_dropin.py. Skips if its collapsed output already exists."""
    dst = collapsed_path(out, theta)
    if dst.exists() and not force:
        print(f"[skip] theta={theta:g} already done -> {dst}", file=sys.stderr)
        return dst
    d = theta_dir(out, theta)
    d.mkdir(parents=True, exist_ok=True)
    cmd = [str(tool_path("ribocode_python", ribocode_python)),
           str(tool_path("ribocode_dropin", dropin)),
           "--profiles", str(profiles), "--annot", str(annot), "--variant", VARIANT,
           "--out", str(d), "--min_aa", str(min_aa), "--pval", str(pval),
           "--pred_scale", f"{theta:g}", "--seed", str(seed)]
    if poisson:
        cmd.append("--pred_poisson")
    print(f"[run ] theta={theta:g}: {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, check=True)
    if not dst.exists():
        raise SystemExit(f"ribocode_dropin produced no output at {dst}")
    return dst


def scored_tx(profiles) -> set:
    """Transcript ids present in the prediction dump (the only ones that can ever be called)."""
    z = np.load(profiles, allow_pickle=False)
    return {str(t) for t in z["tx_ids"]}


def cds_denominator(profiles, tx2cds_path, universe_tx=None):
    """Annotated-CDS transcripts that were actually scored (and in the universe, if given)."""
    denom = set(rc_io.load_tx2cds(tx2cds_path)) & scored_tx(profiles)
    if universe_tx:
        keep = {ln.strip() for ln in open(universe_tx) if ln.strip()}
        denom &= keep
    return denom


def cds_recall(collapsed, denom, qvalue=0.05):
    """Fraction of `denom` transcripts with a significant `annotated` call. -> (recall, hits, n_calls)."""
    hit, n_calls = set(), 0
    for orf in rc_io.read_collapsed(collapsed):
        n_calls += 1
        if orf["orf_type"] == "annotated" and orf["qval"] <= qvalue and orf["tx"] in denom:
            hit.add(orf["tx"])
    return (len(hit) / len(denom) if denom else 0.0), len(hit), n_calls


def class_counts(collapsed, biotype, qvalue=0.05):
    """{pgx_class: n} among significant calls, for the curve's yield columns."""
    out = {}
    for orf in rc_io.read_collapsed(collapsed):
        if orf["qval"] > qvalue:
            continue
        k = rc_io.pgx_class(orf, biotype.get(orf["tx"]))
        if k:
            out[k] = out.get(k, 0) + 1
    return out


def add_relative_recall(curve):
    """Attach cds_recall_rel = hits / (hits at the saturating theta).

    Absolute recall against the full annotation has a hard ceiling well below 1. Measured on BMDM
    the curve plateaus at 0.533 (10,236 of 19,190 CDS-bearing transcripts), so an absolute 0.90
    target is unreachable by construction and the dial would always saturate at theta = 1.

    The cause is ANNOTATION REDUNDANCY, not detection failure -- an earlier version of this comment
    claimed the missing half is "not translated at a detectable level", which is wrong. The
    universe carries 1.90 CDS-bearing isoforms per gene; Ribo-seq signal concentrates on
    essentially one of them, and the caller finds it. At the ceiling, 10,236 transcripts span 9,777
    of 10,080 genes: transcript-level recall 53.3% is **gene-level recall 97.0%**, at 1.05 called
    isoforms per gene. The ceiling is therefore a property of how many isoforms the annotation
    lists, not of the biology or of the model.

    The relative anchor asks the question the dial is actually for: at this stringency, what
    fraction of the CDS this sample CAN support is retained. It normalises away that annotation
    artifact, is self-calibrating across datasets of different depth and isoform density, and on
    BMDM reproduces the operating point Check 5 independently validated as transferable
    (theta ~ 0.05 -> ~90% CDS).
    """
    ceiling = max((r["cds_hits"] for r in curve), default=0)
    for r in curve:
        r["cds_recall_rel"] = round(r["cds_hits"] / ceiling, 4) if ceiling else 0.0
    return ceiling


def saturated(curve, tol=0.02):
    """True if the top two thetas' CDS hits agree within `tol` (so the ceiling is a real plateau)."""
    if len(curve) < 2:
        return False
    top = sorted(curve, key=lambda r: -r["theta"])[:2]
    hi = max(r["cds_hits"] for r in top)
    return hi > 0 and abs(top[0]["cds_hits"] - top[1]["cds_hits"]) / hi <= tol


def pick(curve, target, mode="relative"):
    """Smallest theta meeting the recall floor; else the theta with the highest recall."""
    key = "cds_recall_rel" if mode == "relative" else "cds_recall"
    ok = [r for r in curve if r.get(key, 0.0) >= target]
    if ok:
        return min(ok, key=lambda r: r["theta"]), True
    return max(curve, key=lambda r: r.get(key, 0.0)), False


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profiles", required=True, help="pred_profiles.npz from dump_pred_profiles")
    ap.add_argument("--annot", default=None, help="RiboCode annotation dir (default: species ref)")
    ap.add_argument("--species", required=True, choices=["human", "mouse"])
    ap.add_argument("--out", required=True, help="sweep directory")
    ap.add_argument("--thetas", default=",".join(f"{t:g}" for t in DEFAULT_THETAS))
    ap.add_argument("--only-theta", type=float, default=None,
                    help="run just this theta and exit (sbatch array fan-out)")
    ap.add_argument("--pick-only", action="store_true", help="do not run; read existing sweep output")
    ap.add_argument("--cds-recall", type=float, default=0.90, help="anchor: CDS recall floor")
    ap.add_argument("--recall-mode", default="relative", choices=["relative", "absolute"],
                    help="relative (default): fraction of the CDS the sample can support at all "
                         "(hits / hits at the saturating theta). absolute: fraction of ALL "
                         "annotated CDS, which has a hard sub-1 ceiling and usually cannot reach "
                         "a 0.90 target. See add_relative_recall().")
    ap.add_argument("--orf-qvalue", type=float, default=0.05, help="BH q cutoff on RiboCode calls")
    ap.add_argument("--universe-tx", default=None, help="restrict the denominator to this tx list")
    ap.add_argument("--ref-calls", default=None,
                    help="observed RiboCode *_collapsed.txt; anchor recall on its `annotated` calls "
                         "instead of the annotation (use when real Ribo-seq exists)")
    ap.add_argument("--min-aa", type=int, default=5)
    ap.add_argument("--pval", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--no-poisson", action="store_true", help="deterministic rounding (standard arm)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--tx2cds", default=None)
    ap.add_argument("--tx2biotype", default=None)
    a = ap.parse_args()

    refs = species_refs(a.species, tx2cds=a.tx2cds, tx2biotype=a.tx2biotype,
                        ribocode_annot=a.annot)
    annot = refs["ribocode_annot"]
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    thetas = sorted({float(t) for t in a.thetas.split(",") if t.strip()})

    if a.only_theta is not None:
        run_theta(a.only_theta, a.profiles, annot, out, a.min_aa, a.pval, a.seed,
                  not a.no_poisson, a.force)
        return

    if a.ref_calls:
        denom = {o["tx"] for o in rc_io.read_collapsed(a.ref_calls)
                 if o["orf_type"] == "annotated" and o["qval"] <= a.orf_qvalue}
        denom &= scored_tx(a.profiles)
        anchor = f"observed annotated calls ({Path(a.ref_calls).name})"
    else:
        denom = cds_denominator(a.profiles, refs["tx2cds"], a.universe_tx)
        anchor = "annotated CDS of scored/expressed transcripts"
    if not denom:
        raise SystemExit("empty CDS denominator; check --profiles / --universe-tx / --tx2cds")
    print(f"anchor: {anchor}\nCDS denominator: {len(denom):,} transcripts", file=sys.stderr)

    biotype = rc_io.load_biotype(refs["tx2biotype"])
    curve = []
    for t in thetas:
        cp = collapsed_path(out, t)
        if not cp.exists():
            if a.pick_only:
                print(f"[warn] theta={t:g} missing, skipped", file=sys.stderr)
                continue
            run_theta(t, a.profiles, annot, out, a.min_aa, a.pval, a.seed,
                      not a.no_poisson, a.force)
        rec, hits, n_calls = cds_recall(cp, denom, a.orf_qvalue)
        cc = class_counts(cp, biotype, a.orf_qvalue)
        curve.append({"theta": t, "cds_recall": round(rec, 4), "cds_hits": hits,
                      "n_calls": n_calls, **{k: cc.get(k, 0)
                                             for k in ("uORF", "dORF", "lncRNA_orf")}})
    if not curve:
        raise SystemExit("no theta results found; run without --pick-only first")

    ceiling = add_relative_recall(curve)
    sat = saturated(curve)
    best, met = pick(curve, a.cds_recall, a.recall_mode)
    cols = ["theta", "cds_recall", "cds_recall_rel", "cds_hits", "n_calls",
            "uORF", "dORF", "lncRNA_orf"]
    tsv = out / "theta_curve.tsv"
    with open(tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in curve:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    meta = {"theta": best["theta"], "cds_recall": best["cds_recall"],
            "cds_recall_rel": best["cds_recall_rel"], "recall_mode": a.recall_mode,
            "target_cds_recall": a.cds_recall, "target_met": met, "anchor": anchor,
            "cds_denominator": len(denom), "cds_ceiling": ceiling, "ceiling_saturated": sat,
            "poisson": not a.no_poisson, "orf_qvalue": a.orf_qvalue, "seed": a.seed,
            "variant": VARIANT, "collapsed": str(collapsed_path(out, best["theta"]))}
    (out / "theta.json").write_text(json.dumps(meta, indent=2) + "\n")

    print("\n" + "  ".join(f"{c:>14}" for c in cols))
    for r in curve:
        mark = " <-- chosen" if r["theta"] == best["theta"] else ""
        print("  ".join(f"{r[c]:>14}" for c in cols) + mark)
    print(f"\nCDS ceiling (max hits over the sweep): {ceiling:,} of {len(denom):,} annotated "
          f"({ceiling / len(denom):.3f}); plateau reached: {sat}")
    if not sat:
        print("WARNING: the top two thetas still differ by > 2%, so the ceiling is not a true "
              "plateau and cds_recall_rel understates stringency. Extend --thetas upward.")
    if not met:
        print(f"\nWARNING: no theta reached {a.recall_mode} CDS recall {a.cds_recall}; using the "
              f"best available. Consider a deeper prediction or a lower target.")
    print(f"\ntheta* = {best['theta']:g}  (CDS recall abs {best['cds_recall']}, "
          f"rel {best['cds_recall_rel']})")
    print(f"wrote {tsv}\nwrote {out / 'theta.json'}")


if __name__ == "__main__":
    main()
