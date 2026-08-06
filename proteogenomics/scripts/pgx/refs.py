#!/usr/bin/env python3
"""Species reference + tool path resolution for the `pgx` proteogenomics pipeline.

NO hardcoded absolute paths in callers. Every reference resolves in this order:

  1. explicit CLI argument passed by the caller (highest priority)
  2. environment variable  PGX_<SPECIES>_<KEY>  (species refs) or PGX_<KEY> (tools/shared)
  3. the species registry below, anchored on project_root() from scripts/prepare/paths.py

Supported species: human (GENCODE v49) and mouse (GENCODE vM38).

Usage:
  python -m pgx.refs --species mouse --print          # list + existence-check every path
  from pgx.refs import species_refs, tool_path
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
# scripts/prepare/paths.py owns root resolution; reuse it rather than duplicating the logic.
sys.path.insert(0, str(_HERE.parents[3] / "scripts" / "prepare"))
from paths import project_root  # noqa: E402

ROOT = project_root()
# <...>/RNAZoo_meta/annotations, i.e. two levels above RNAZoo/experiments/riboseq_signal_model
ANN = Path(os.environ.get("PGX_ANNOTATIONS", ROOT.parents[2] / "annotations"))
EXP = ROOT.parent                                    # <...>/RNAZoo/experiments

# Keys that name a single file/dir per species. `tx_fastas` is a list and handled separately.
SPECIES_KEYS = ("gencode_proteome", "tx2cds", "tx2biotype", "ribocode_annot")

REGISTRY = {
    "human": {
        "release": "v49",
        "gencode_proteome": ANN / "gencode_proteins" / "gencode.v49.pc_translations.fa.gz",
        "tx2cds": ROOT / "data" / "tx2cds.tsv",
        "tx2biotype": ROOT / "data" / "tx2biotype.tsv",
        "tx_fastas": [ANN / "gencode.v49.pc_transcripts.fa",
                      ANN / "gencode.v49.lncRNA_transcripts.fa"],
        "ribocode_annot": EXP / "biotype_probe" / "expression_context_human" / "data" / "ribocode_annot",
    },
    "mouse": {
        "release": "vM38",
        "gencode_proteome": ANN / "gencode_proteins" / "gencode.vM38.pc_translations.fa.gz",
        "tx2cds": ROOT / "data" / "mouse_tx2cds.tsv",
        "tx2biotype": ROOT / "data" / "tx2biotype_mouse.tsv",
        "tx_fastas": [ANN / "gencode.vM38.pc_transcripts.fa",
                      ANN / "gencode.vM38.lncRNA_transcripts.fa"],
        "ribocode_annot": ROOT / "data" / "mouse_ribocode_annot",
    },
}

# Tool / env paths: shared across species, same 3-tier resolution via PGX_<KEY>.
TOOLS = {
    "msfragger_jar": Path("/private/groups/carpenterlab/emalekos/conda_envs/msfragger/share/"
                          "msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar"),
    "java": Path("/private/groups/carpenterlab/emalekos/conda_envs/msfragger/bin/java"),
    "ribocode_python": Path("/private/groups/carpenterlab/emalekos/conda_envs/ribocode/bin/python"),
    "python": Path("/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3"),
    # scripts this package shells out to
    "ribocode_dropin": ROOT / "scripts" / "ribocode_dropin.py",
    "build_line_universe": ROOT / "proteogenomics" / "scripts" / "build_line_universe.py",
    "msfragger_params_dir": ROOT / "proteogenomics" / "msfragger",
}


def species_refs(species: str, **overrides) -> dict:
    """Resolved reference paths for `species`. Any keyword overrides win outright.

    Env override name: PGX_<SPECIES>_<KEY> uppercased, e.g. PGX_MOUSE_TX2CDS.
    """
    if species not in REGISTRY:
        raise SystemExit(f"unknown species {species!r}; known: {sorted(REGISTRY)}")
    base = dict(REGISTRY[species])
    for key in SPECIES_KEYS:
        ov = overrides.get(key)
        if ov:
            base[key] = Path(ov)
            continue
        env = os.environ.get(f"PGX_{species.upper()}_{key.upper()}")
        if env:
            base[key] = Path(env)
    ov = overrides.get("tx_fastas")
    if ov:
        base["tx_fastas"] = [Path(p) for p in (ov if isinstance(ov, (list, tuple)) else str(ov).split(","))]
    else:
        env = os.environ.get(f"PGX_{species.upper()}_TX_FASTAS")
        if env:
            base["tx_fastas"] = [Path(p) for p in env.split(",") if p]
    return base


def tool_path(key: str, override=None) -> Path:
    """Resolved path for a tool/helper. Env override: PGX_<KEY> uppercased."""
    if override:
        return Path(override)
    env = os.environ.get(f"PGX_{key.upper()}")
    if env:
        return Path(env)
    if key not in TOOLS:
        raise SystemExit(f"unknown tool {key!r}; known: {sorted(TOOLS)}")
    return TOOLS[key]


def check(species: str, **overrides) -> list:
    """Return [(label, path, exists)] for every resolved path, species refs + tools."""
    refs = species_refs(species, **overrides)
    rows = [(f"{species}.{k}", refs[k], Path(refs[k]).exists()) for k in SPECIES_KEYS]
    rows += [(f"{species}.tx_fastas[{i}]", p, Path(p).exists())
             for i, p in enumerate(refs["tx_fastas"])]
    rows += [(f"tool.{k}", tool_path(k), tool_path(k).exists()) for k in sorted(TOOLS)]
    return rows


def require(species: str, **overrides) -> dict:
    """species_refs() but hard-fails listing every missing path (fail fast, not mid-pipeline)."""
    rows = check(species, **overrides)
    missing = [(lab, p) for lab, p, ok in rows if not ok]
    if missing:
        msg = "\n".join(f"  MISSING {lab}: {p}" for lab, p in missing)
        raise SystemExit(f"unresolved reference paths for species={species}:\n{msg}")
    return species_refs(species, **overrides)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", default="human", choices=sorted(REGISTRY))
    ap.add_argument("--print", action="store_true", help="list every resolved path + existence")
    a = ap.parse_args()
    print(f"project_root : {ROOT}")
    print(f"annotations  : {ANN}")
    print(f"release      : {REGISTRY[a.species]['release']}\n")
    bad = 0
    for lab, p, ok in check(a.species):
        bad += not ok
        print(f"[{'ok ' if ok else 'MISS'}] {lab:<28} {p}")
    print(f"\n{bad} missing" if bad else "\nall paths resolved")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
