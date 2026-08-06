#!/usr/bin/env python3
"""Index per-token FM embeddings: tx_id -> .npy path, across the 8 chunk dirs.

The extractor writes one {tx_id}_tokens.npy per transcript under
<emb-dir>/chunk_{1..8}/. This flattens them into a single lookup the model data loader
uses to fetch a transcript's (L, d_emb) embedding by id. Also reports coverage of the
universe and any missing transcripts (e.g. dropped over the length cap).

Default indexes RiNALMo (data/rinalmo_token_emb). Pass --emb-dir for another backend, e.g.
  build_embedding_index.py --emb-dir data/orthrus4t_token_emb

Output: <emb-dir>/tx_index.tsv  (tx_id, path)
"""
import argparse
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
UNIV = NEW / "data" / "fibroblast_universe.tsv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb-dir", default=str(NEW / "data" / "rinalmo_token_emb"),
                    help="dir holding chunk_*/<tx>_tokens.npy (default RiNALMo)")
    args = ap.parse_args()
    emb_dir = Path(args.emb_dir)
    if not emb_dir.is_absolute():
        emb_dir = NEW / emb_dir
    out = emb_dir / "tx_index.tsv"

    idx = {}
    for chunk in sorted(emb_dir.glob("chunk_*")):
        for npy in chunk.glob("*_tokens.npy"):
            tx = npy.name[:-len("_tokens.npy")]
            idx[tx] = str(npy)

    with out.open("w") as o:
        o.write("tx_id\tpath\n")
        for tx in sorted(idx):
            o.write(f"{tx}\t{idx[tx]}\n")

    universe = set()
    if UNIV.exists():
        with UNIV.open() as fh:
            fh.readline()
            for line in fh:
                universe.add(line.split("\t", 1)[0])
    have = set(idx)
    missing = universe - have
    print(f"indexed embeddings: {len(idx):,}", file=sys.stderr)
    print(f"universe: {len(universe):,}  with embedding: {len(universe & have):,}  "
          f"missing: {len(missing):,}", file=sys.stderr)
    if missing:
        ex = list(sorted(missing))[:5]
        print(f"  first missing: {ex}", file=sys.stderr)
    print(f"wrote {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
