#!/usr/bin/env python3
"""Check that every pack's recorded source inputs still exist on disk.

WHY THIS EXISTS. In August 2026 the mouse-liver 3x3 factorial was blocked because the five Janich
RIBO-SEQ runs had been deleted -- no BAMs, no `*_psites.hd5`, nothing. A pack cannot be re-projected
onto a different universe without its sources, so an entire experiment was unbuildable, and the loss
was discovered only when the work was attempted.

Two things made it invisible:

  * `data/heldout_bam/mouse_janich_liver/` exists and holds seven BAMs -- but they are the RNA
    samples (`SRR19302xx`), matching the RNA coverage hd5 one for one. The directory name says only
    "janich_liver", so it LOOKS like the Janich BAMs are archived. Pooling those through the Ribo path
    would have produced a "Janich Ribo" arm made of RNA-seq: flat, non-periodic, silently wrong.
  * Nothing ever checked. `provenance.json` records exactly which files a pack was built from, and
    that record was correct -- it was simply never compared against the filesystem.

This closes the second gap. Run it after any cleanup, before planning work that needs re-projection,
and periodically. It is read-only.

Exit status: 0 if every pack has all inputs, 1 if anything is missing.

Usage: audit_pack_inputs.py [--packs-glob 'data/packed*'] [--quiet]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
INPUT_KEYS = ("ribo_psites_inputs", "rna_coverage_inputs", "coverage_inputs", "psites_inputs")


def entry_path(e):
    return Path(e["path"]) if isinstance(e, dict) else Path(str(e))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--packs-glob", default="data/packed*")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()

    packs = sorted(p for p in NEW.glob(a.packs_glob) if p.is_dir())
    bad, noprov, total_missing = [], [], 0
    print(f"auditing {len(packs)} pack(s) under {NEW}/{a.packs_glob}\n")
    for pk in packs:
        prov = pk / "provenance.json"
        if not prov.exists():
            noprov.append(pk.name)
            continue
        try:
            j = json.loads(prov.read_text())
        except Exception as e:
            bad.append((pk.name, f"unreadable provenance.json: {e}", []))
            continue
        missing, present = [], 0
        for k in INPUT_KEYS:
            for e in j.get(k, []) or []:
                p = entry_path(e)
                if p.exists():
                    present += 1
                else:
                    missing.append(f"{k}: {p.name}")
        if missing:
            bad.append((pk.name, f"{present} present, {len(missing)} MISSING", missing))
            total_missing += len(missing)
        elif not a.quiet:
            print(f"  OK    {pk.name:<52} {present} inputs")

    if noprov:
        print(f"\n  {len(noprov)} pack(s) with NO provenance.json -- their sources are unrecorded, "
              "which is its own risk:")
        for n in noprov:
            print(f"    {n}")
    if bad:
        print(f"\n  {len(bad)} pack(s) with MISSING inputs ({total_missing} files):")
        for name, summ, miss in bad:
            print(f"    {name}: {summ}")
            for m in miss[:6]:
                print(f"      - {m}")
            if len(miss) > 6:
                print(f"      ... and {len(miss) - 6} more")
        print("\n  A pack with missing inputs CANNOT be re-projected onto a different universe.")
        print("  Re-fetch from SRA before planning work that needs it.")
        return 1
    print("\n  all packs have every recorded input on disk")
    return 0


if __name__ == "__main__":
    sys.exit(main())
