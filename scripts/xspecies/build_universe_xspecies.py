#!/usr/bin/env python3
"""Build a per-species transcript universe for the cross-species packs.

Follows scripts/heldout/build_mouse_universe.py exactly, generalised off
scripts/xspecies/species.tsv so it works for any of the nine species:

  universe = protein_coding + lncRNA transcripts, mature length <= 10,000 nt,
             that are present on BOTH the P-site and the RNA-coverage transcriptome axes.

No expression prefilter, deliberately, and for the same reason as the mouse build: the
>=50-pooled-P-site scorable floor at eval time defines the scored set, and an untranslated
transcript in the universe is harmless because it simply never becomes a test transcript.
Filtering on TPM here would instead bake a second, different expression criterion into the
axis itself.

The axis intersection is the part that matters. `tx_order` in the pack is this list, and the
ORF track, the one-hot FASTA and the model input all index against it positionally. A
transcript present in the coverage hd5 but absent from the P-site hd5 (or vice versa) would
shift everything downstream of it, and there is no shape check that catches that.

Biotype comes from the NORMALIZED tx2biotype.tsv, so every species is filtered on the same
GENCODE vocabulary rather than on whatever its source GTF happened to call things.

Writes to data/xspecies_refs/:
  <dataset>_universe.txt   one transcript id per line, sorted
  <dataset>_universe.fa    bare-header FASTA (>tx_id)
  <dataset>_universe.tsv   tx_id, length, biotype

Usage:
  build_universe_xspecies.py --species <sp> --ribo-ds <ds> --rna-ds <ds>
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import h5py

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
OUTDIR = NEW / "data" / "xspecies_refs"
MAXLEN = 10000
KEEP_BIOTYPES = {"protein_coding", "lncRNA"}


def iter_fasta(path: Path):
    name, buf = None, []
    with path.open() as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = ln[1:].split("|")[0].split()[0]
                buf = []
            else:
                buf.append(ln.strip())
    if name is not None:
        yield name, "".join(buf)


def axis_ids(path: Path, key: str = "transcript_ids") -> set[str]:
    with h5py.File(path, "r") as h:
        if key not in h:
            sys.exit(f"{path} has no dataset '{key}' (keys: {list(h.keys())[:8]})")
        return set(h[key].asstr()[:])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", required=True)
    ap.add_argument("--ribo-ds", required=True, help="dataset label under data/xspecies_psites/")
    ap.add_argument("--rna-ds", required=True, help="dataset label under data/xspecies_rna_coverage_mm10/")
    ap.add_argument("--out-prefix", default=None)
    args = ap.parse_args()

    prefix = args.out_prefix or args.ribo_ds
    OUTDIR.mkdir(parents=True, exist_ok=True)

    ps_dir = NEW / "data" / "xspecies_psites" / args.ribo_ds
    cov_dir = NEW / "data" / "xspecies_rna_coverage_mm10" / args.rna_ds
    ps = sorted(ps_dir.glob("*_psites.hd5"))
    cv = sorted(cov_dir.glob("*_coverage.hd5"))
    if not ps:
        sys.exit(f"no *_psites.hd5 in {ps_dir}")
    if not cv:
        sys.exit(f"no *_coverage.hd5 in {cov_dir}")

    # Intersect across EVERY file, not just the first: a run whose axis differs would
    # otherwise silently contribute transcripts the others cannot score.
    print(f"  {len(ps)} psite files, {len(cv)} coverage files", file=sys.stderr)
    on_axis = axis_ids(ps[0])
    for p in ps[1:]:
        on_axis &= axis_ids(p)
    n_ps = len(on_axis)
    for c in cv:
        on_axis &= axis_ids(c)
    print(f"  psite axis (intersected) = {n_ps:,};  after coverage = {len(on_axis):,}", file=sys.stderr)
    if not on_axis:
        sys.exit("ABORT: empty axis intersection -- psite and coverage transcriptomes disagree")

    bt_path = NEW / "data" / "annot" / args.species / "tx2biotype.tsv"
    fa_path = NEW / "data" / "annot" / args.species / "transcripts.fa"
    for p in (bt_path, fa_path):
        if not p.exists():
            sys.exit(f"missing {p}")
    biotype = {r["tx_id"]: r["transcript_type"]
               for r in csv.DictReader(bt_path.open(), delimiter="\t")}

    seqs, keep_bt = {}, {}
    n_seen = n_len = n_axis = n_bt = 0
    for tx, seq in iter_fasta(fa_path):
        n_seen += 1
        bt = biotype.get(tx)
        if bt not in KEEP_BIOTYPES:
            n_bt += 1
            continue
        if len(seq) > MAXLEN:
            n_len += 1
            continue
        if tx not in on_axis:
            n_axis += 1
            continue
        seqs[tx] = seq
        keep_bt[tx] = bt
    order = sorted(seqs)
    print(f"  parsed={n_seen:,}  dropped biotype={n_bt:,}  dropped len>{MAXLEN}={n_len:,}  "
          f"dropped off-axis={n_axis:,}  universe={len(order):,}", file=sys.stderr)
    if not order:
        sys.exit("ABORT: empty universe")

    (OUTDIR / f"{prefix}_universe.txt").write_text("\n".join(order) + "\n")
    with (OUTDIR / f"{prefix}_universe.fa").open("w") as fh:
        for tx in order:
            fh.write(f">{tx}\n")
            s = seqs[tx]
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")
    npc = 0
    with (OUTDIR / f"{prefix}_universe.tsv").open("w") as fh:
        fh.write("tx_id\tlength\tbiotype\n")
        for tx in order:
            fh.write(f"{tx}\t{len(seqs[tx])}\t{keep_bt[tx]}\n")
            npc += keep_bt[tx] == "protein_coding"
    print(f"  wrote {OUTDIR}/{prefix}_universe.{{txt,fa,tsv}}  "
          f"pc={npc:,} lncRNA={len(order)-npc:,}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
