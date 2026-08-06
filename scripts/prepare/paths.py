#!/usr/bin/env python3
"""Resolve the riboseq_signal_model project root WITHOUT hardcoding an absolute path.

Resolution order:
  1. $RIBOSEQ_SIGNAL_MODEL_ROOT if set.
  2. Auto-detect: walk up from this file to the first ancestor that contains both a
     `scripts/` and a `data/` directory (the project root).
  3. Fallback: the original on-cluster absolute path (so existing callers keep working).

Every new prepare/ script imports project_root() instead of baking in the absolute path,
which is the single biggest portability fragility in the old prep scripts.
"""
from __future__ import annotations

import os
from pathlib import Path

# Kept only as a last-resort fallback; do NOT rely on it in new code.
_FALLBACK_ROOT = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
                      "RNAZoo/experiments/riboseq_signal_model")


def project_root() -> Path:
    env = os.environ.get("RIBOSEQ_SIGNAL_MODEL_ROOT")
    if env:
        return Path(env).resolve()
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "scripts").is_dir() and (parent / "data").is_dir():
            return parent
    return _FALLBACK_ROOT


if __name__ == "__main__":
    print(project_root())
