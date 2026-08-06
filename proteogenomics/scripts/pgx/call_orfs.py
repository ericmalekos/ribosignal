#!/usr/bin/env python3
"""Two-arm RiboCode ORF calling on the model's PREDICTED signal (`pgx` step 4).

Standing project rule (feedback_two_arm_orf_calling): every ORF-calling test reports BOTH

  standard : pred_preddepth, theta = 1.0, deterministic rounding   (naive standalone prediction)
  poisson  : pred_preddepth, theta = theta*, Poisson sampling      (calibrated de novo recipe)

theta* comes from pgx.calibrate. The Poisson arm at theta* is ALREADY produced by the calibration
sweep, so it is reused in place rather than recomputed -- the sweep directory is the cache.

Calls are filtered at --orf-qvalue on RiboCode's BH-adjusted p (`adjusted_pval`, present only in
`*_collapsed.txt`) and classified with rc_io.pgx_class: uORF / dORF / lncRNA_orf are novel;
`annotated` and `internal` are not (canonical protein and CDS truncation respectively). N-terminal
extensions live in the `annotated` rows and are harvested separately by pgx.extensions.

  python -m pgx.call_orfs --profiles p.npz --species mouse --out <dir> --theta-json <calib>/theta.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import rc_io
from .calibrate import collapsed_path, run_theta
from .refs import species_refs

ARMS = ("standard", "poisson")
COLS = ("orf_id", "tx", "gene_id", "gene_name", "pgx_class", "orf_type", "start0", "end0",
        "start_codon", "aa_len", "psites_f0", "qval", "aaseq")


def resolve_theta(theta_json=None, theta=None):
    if theta is not None:
        return float(theta), None
    if not theta_json:
        raise SystemExit("need --theta or --theta-json")
    meta = json.loads(Path(theta_json).read_text())
    if not meta.get("target_met", True):
        print(f"WARNING: calibration never reached its CDS-recall target "
              f"({meta.get('cds_recall')} < {meta.get('target_cds_recall')})", file=sys.stderr)
    return float(meta["theta"]), meta


def arm_collapsed(arm, theta, profiles, annot, out, calib_dir=None, min_aa=5, pval=0.05,
                  seed=0, force=False):
    """Collapsed RiboCode output for one arm, reusing the calibration sweep where possible.

    The poisson arm at theta* is exactly what the sweep already computed, so the sweep directory
    doubles as its cache. run_theta writes to <dir>/theta_<scale>/ and is itself skip-if-present.
    """
    if arm == "poisson" and calib_dir:
        cached = collapsed_path(Path(calib_dir), theta)
        if cached.exists() and not force:
            print(f"[reuse] poisson arm theta={theta:g} <- {cached}", file=sys.stderr)
            return cached
    scale = 1.0 if arm == "standard" else theta
    d = Path(out) / arm
    d.mkdir(parents=True, exist_ok=True)
    return run_theta(scale, profiles, annot, d, min_aa=min_aa, pval=pval, seed=seed,
                     poisson=(arm == "poisson"), force=force)


def collect(collapsed, biotype, qvalue):
    """-> (novel_rows, stats). novel_rows are significant, non-canonical, non-internal calls."""
    rows, stats = [], {"total": 0, "significant": 0, "annotated": 0, "internal": 0, "nterm_ext": 0}
    for orf in rc_io.read_collapsed(collapsed):
        stats["total"] += 1
        if orf["qval"] > qvalue:
            continue
        stats["significant"] += 1
        if orf["orf_type"] == "annotated":
            stats["annotated"] += 1
            stats["nterm_ext"] += rc_io.is_nterm_extension(orf)
            continue
        if orf["orf_type"] == "internal":
            stats["internal"] += 1
            continue
        k = rc_io.pgx_class(orf, biotype.get(orf["tx"]))
        if not k or not orf["aaseq"]:
            continue
        rows.append({**{c: orf.get(c) for c in COLS if c in orf}, "pgx_class": k})
        stats[k] = stats.get(k, 0) + 1
    return rows, stats


def write_calls(rows, path):
    with open(path, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--species", required=True, choices=["human", "mouse"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--theta-json", default=None, help="theta.json from pgx.calibrate")
    ap.add_argument("--theta", type=float, default=None, help="explicit theta (overrides --theta-json)")
    ap.add_argument("--calib-dir", default=None,
                    help="calibration sweep dir; the poisson arm is reused from it when present "
                         "(defaults to the directory holding --theta-json)")
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--orf-qvalue", type=float, default=0.05)
    ap.add_argument("--min-aa", type=int, default=5)
    ap.add_argument("--pval", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--annot", default=None)
    ap.add_argument("--tx2biotype", default=None)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    refs = species_refs(a.species, tx2biotype=a.tx2biotype, ribocode_annot=a.annot)
    theta, meta = resolve_theta(a.theta_json, a.theta)
    calib_dir = a.calib_dir or (str(Path(a.theta_json).parent) if a.theta_json else None)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    biotype = rc_io.load_biotype(refs["tx2biotype"])
    arms = [x.strip() for x in a.arms.split(",") if x.strip()]

    summary = {}
    for arm in arms:
        if arm not in ARMS:
            raise SystemExit(f"unknown arm {arm!r}; known: {ARMS}")
        cp = arm_collapsed(arm, theta, a.profiles, refs["ribocode_annot"], out, calib_dir,
                           a.min_aa, a.pval, a.seed, a.force)
        rows, stats = collect(cp, biotype, a.orf_qvalue)
        write_calls(rows, out / f"calls_{arm}.tsv")
        stats["collapsed"] = str(cp)
        stats["theta"] = 1.0 if arm == "standard" else theta
        stats["novel_total"] = len(rows)
        summary[arm] = stats
        print(f"[{arm}] theta={stats['theta']:g} calls={stats['total']:,} "
              f"significant={stats['significant']:,} novel={len(rows):,} "
              f"(uORF {stats.get('uORF', 0):,} / dORF {stats.get('dORF', 0):,} / "
              f"lncRNA_orf {stats.get('lncRNA_orf', 0):,})  "
              f"annotated={stats['annotated']:,} of which nterm_ext={stats['nterm_ext']:,}")

    summary["_meta"] = {"theta": theta, "orf_qvalue": a.orf_qvalue, "species": a.species,
                        "profiles": str(a.profiles), "calibration": meta}
    (out / "call_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nwrote {out}/calls_<arm>.tsv + call_summary.json")


if __name__ == "__main__":
    main()
