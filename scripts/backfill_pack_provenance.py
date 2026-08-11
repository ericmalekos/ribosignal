#!/usr/bin/env python3
"""Write a RECONSTRUCTED provenance record for packs that have none.

31 of 60 packs carry no `provenance.json`, so `audit_pack_inputs.py` cannot see them at all. That is
worse than a pack with known-missing inputs: an audit that reports "all clear" while silently skipping
half the packs is misleading in the direction that matters.

This does NOT invent inputs. The source files for these packs are unknown and in several cases
deleted, and guessing paths would be worse than the current gap -- a fabricated provenance would make
an unreproducible pack look reproducible. What it records instead is everything verifiable from the
bytes on disk today:

  * the pack contract (n_tx, sum_L, totals, dtypes) so corruption is detectable later
  * per-file sizes and mtimes
  * `"inputs": null` with `"inputs_known": false` and a reason

so the audit can distinguish three states that are currently conflated:
  1. provenance present, inputs on disk        -> reproducible
  2. provenance present, inputs missing        -> known gap, recoverable if re-fetched
  3. NO provenance at all                      -> invisible; this script converts these to state 2'
                                                  ("inputs unknown"), which is honest and auditable

Existing provenance.json files are never touched.
"""
import argparse
import datetime
import json
import pathlib
import sys

import numpy as np


def summarize(p: pathlib.Path):
    out = {}
    order = p / "tx_order.txt"
    out["n_tx"] = sum(1 for _ in open(order))
    for name in ("target_counts", "coverage", "offsets", "lengths"):
        f = p / f"{name}.npy"
        if not f.exists():
            continue
        a = np.load(f, mmap_mode="r")
        out[name] = dict(shape=list(a.shape), dtype=str(a.dtype),
                         sum=int(np.asarray(a).sum()) if name in ("target_counts", "coverage") else None,
                         bytes=f.stat().st_size)
    cn = p / "coverage_norm.json"
    if cn.exists():
        try:
            out["coverage_norm"] = json.load(open(cn))
        except Exception:
            out["coverage_norm"] = "unparseable"
    tracks = sorted(x.name for x in p.glob("orf_track*.npy"))
    out["orf_tracks"] = tracks
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    root = pathlib.Path(a.data_dir)
    now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")

    wrote, skipped = 0, 0
    for p in sorted(root.glob("packed*")):
        if not (p / "tx_order.txt").exists():
            continue
        prov = p / "provenance.json"
        if prov.exists():
            skipped += 1
            continue
        try:
            s = summarize(p)
        except Exception as e:  # a pack too broken to summarize is itself worth recording
            s = {"error": f"{type(e).__name__}: {e}"}
        rec = {
            "group": p.name.replace("packed_heldout_", "").replace("packed_", "") or "Fibroblast",
            "provenance_reconstructed": True,
            "reconstructed_at": now,
            "reconstructed_by": "scripts/backfill_pack_provenance.py",
            "inputs_known": False,
            "ribo_psites_inputs": None,
            "rna_coverage_inputs": None,
            "why_inputs_unknown": (
                "This pack predates provenance recording. Its source hd5 were never listed and in "
                "several cases have been deleted. Input paths are deliberately NOT guessed: a "
                "fabricated provenance would make an unreproducible pack appear reproducible. "
                "Treat this pack as NOT rebuildable from source until its inputs are re-derived."),
            "verified_contract": s,
            "mtimes": {x.name: datetime.datetime.fromtimestamp(x.stat().st_mtime)
                       .astimezone().isoformat(timespec="seconds")
                       for x in sorted(p.iterdir()) if x.is_file()},
        }
        print(f"  {p.name}: n_tx={s.get('n_tx','?')} "
              f"target={s.get('target_counts',{}).get('sum') if isinstance(s.get('target_counts'), dict) else 'n/a'}")
        if not a.dry_run:
            prov.write_text(json.dumps(rec, indent=2))
        wrote += 1
    print(f"\n  {'would write' if a.dry_run else 'wrote'} {wrote} reconstructed provenance record(s); "
          f"{skipped} pack(s) already had one and were left untouched")
    return 0


if __name__ == "__main__":
    sys.exit(main())
