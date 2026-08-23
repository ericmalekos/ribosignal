#!/usr/bin/env python3
"""Fail if any figure values file does not record which model and which source produced it.

WHY THIS EXISTS. 16 of 47 values files had drifted with no model recorded, and the ones that DID
record it used six different key names (`source_npz`, `source_tsv`, `source_json`, `source_dir`,
`model_arm`, `models`). A convention that is not checked is a convention that drifts, and the cost
is real: a poster panel was captioned with the wrong architecture because an unsuffixed filename
was assumed to mean `attn` when the generator's own default was `mamba4`.

THE CONTRACT. Every `figures/*/<name>_values*.json` must carry, at top level:

  model   -- "attn" | "mamba4" | a list of both | null
             null is legal and means "no model was involved" (a dataset census, a raw-data QC
             panel, an observed-vs-observed comparison). It MUST be explicit: the point is to
             distinguish "no model here" from "nobody wrote it down". Those are different, and
             only one of them is a problem.
  source  -- where the numbers came from: a run dir, results path, or a one-line description.

  usage:  check_figure_provenance.py [--figures DIR]
  exit 0 = all pass, 1 = at least one offender (prints them)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

VALID = {"attn", "mamba4"}


def check(fig_dir: Path) -> int:
    # Retired figures are excluded: they are kept as a record of work that was pulled, so holding
    # them to the live contract would mean either editing a frozen artifact or a permanent red check.
    files = [f for f in sorted(fig_dir.glob("*/*_values*.json"))
             if not f.parent.name.startswith("_retired")]
    if not files:
        print(f"no values files under {fig_dir}", file=sys.stderr)
        return 1
    bad = []
    for p in files:
        try:
            d = json.load(open(p))
        except Exception as e:
            bad.append((p, f"unreadable: {e}"))
            continue
        if not isinstance(d, dict):
            bad.append((p, "top level is not an object"))
            continue
        if "model" not in d:
            bad.append((p, "no 'model' key (use null if genuinely model-free)"))
        else:
            m = d["model"]
            ms = m if isinstance(m, list) else [m]
            unknown = [x for x in ms if x is not None and x not in VALID]
            if unknown:
                bad.append((p, f"unrecognised model {unknown!r}; want {sorted(VALID)} or null"))
        if not d.get("source"):
            bad.append((p, "no 'source' key (where did these numbers come from?)"))
    for p, why in bad:
        print(f"  FAIL  {p.parent.name}/{p.name}: {why}")
    print(f"{len(files) - len({p for p, _ in bad})}/{len(files)} values files carry provenance")
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figures", type=Path, default=Path(__file__).resolve().parents[1] / "figures")
    return check(ap.parse_args().figures)


if __name__ == "__main__":
    raise SystemExit(main())
