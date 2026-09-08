#!/usr/bin/env python3
"""Where the pipeline's data lives.

Every script in this repo used to carry the author's cluster path as a module-level
constant, which made the whole tree unrunnable anywhere else. Paths are resolved here
instead, from, in order of precedence:

  1. an explicit command-line flag,
  2. an environment variable,
  3. a location under the repository itself.

There is no cluster-specific default anywhere. `REPO_ROOT` is derived from this file's
own location, so a clone works from wherever it is checked out.

Environment variables, all optional:

  RIBO_DATA_DIR        root for everything below (default <repo>/data)
  RIBO_PACK_DIR        packed store: tx_order.txt, offsets/lengths/coverage/target_counts
  RIBO_PACK_SUFFIX     select a pack family, e.g. "union" -> packed_union[_<Tissue>]
  RIBO_ORF_TRACK       ORF-candidate track .npy (default <pack>/orf_track.npy)
  RIBO_ONEHOT_FASTA    universe FASTA for the one-hot backend, and the sequence check
  RIBO_SPLIT           chromosome k-fold split JSON
  RIBO_TX2BIOTYPE      tx -> (gene, biotype, chrom, length) TSV
  RIBO_EMB_DIR         root of the per-token FM embedding indexes
  RIBO_TMPDIR          scratch for the BAM projection (default the system temp dir)

`missing()` renders the one error message these all share, so a user who has not set a
variable is told the flag, the variable and the file in one line rather than seeing a
bare FileNotFoundError from three frames down.
"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def env_path(var: str, default: Path) -> Path:
    """`$var` as a Path if it is set and non-empty, else `default`."""
    v = os.environ.get(var, "")
    return Path(v).expanduser() if v else default


_env_path = env_path  # short alias used inside this module


def data_dir() -> Path:
    """Root of the data tree. Defaults inside the repo so a fresh clone resolves."""
    return _env_path("RIBO_DATA_DIR", REPO_ROOT / "data")


def pack_dir() -> Path:
    return _env_path("RIBO_PACK_DIR", data_dir() / "packed")


def onehot_fasta() -> Path:
    return _env_path("RIBO_ONEHOT_FASTA", data_dir() / "universe.fa")


def split_json() -> Path:
    return _env_path("RIBO_SPLIT", data_dir() / "splits" / "fibroblast_chrom_kfold.json")


def tx2biotype() -> Path:
    return _env_path("RIBO_TX2BIOTYPE", data_dir() / "tx2biotype.tsv")


def emb_dir() -> Path:
    return _env_path("RIBO_EMB_DIR", data_dir())


def tmp_dir() -> Path:
    return _env_path("RIBO_TMPDIR", Path(tempfile.gettempdir()))


def missing(path, what: str, env_var: str = "", flag: str = "") -> str:
    """The standard 'you have not pointed me at X' message.

    Names the file, what it is for, and both ways to supply it, because the common case
    is a new user who has the file under a different name in a different directory.
    """
    how = " or ".join(x for x in (f"--{flag.lstrip('-')}" if flag else "", env_var) if x)
    tail = f"; set it with {how}" if how else ""
    return f"{what} not found at {path}{tail}"
