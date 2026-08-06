#!/usr/bin/env python3
"""Build the mouse (GENCODE vM38) transcript universe for the cross-species held-out eval.

Universe = all protein_coding + lncRNA transcripts (from the two vM38 transcript FASTAs) with
mature length <= 10000 nt that are also present on the Wang psite + coverage transcriptome axis.
No expression prefilter: the >=50-pooled-P-site scorable floor at eval time defines the scored
set, and untranslated tx in the universe are harmless (they simply never become test tx). One-hot
only, so no FM embeddings are needed -- the point of the cross-species one-hot control.

Writes:
  data/heldout_refs/mouse_wang_universe.txt   one versioned ENSMUST per line (sorted)
  data/heldout_refs/mouse_wang_universe.fa    bare-header FASTA (>tx_id), sequence per tx
  data/heldout_refs/mouse_wang_universe.tsv   tx_id, length, biotype

Usage: build_mouse_universe.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import h5py

BASE = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta")
NEW = BASE / "RNAZoo/experiments/riboseq_signal_model"
ANN = BASE / "annotations"
PC_FA = ANN / "gencode.vM38.pc_transcripts.fa"
LNC_FA = ANN / "gencode.vM38.lncRNA_transcripts.fa"
PSITES = NEW / "data" / "heldout_psites" / "mouse_wang_liver" / "SRR5262890.Aligned.toTranscriptome.out_psites.hd5"
COVER = NEW / "data" / "heldout_rnaseq_coverage" / "mouse_wang_liver" / "SRR5262874_coverage.hd5"
OUTDIR = NEW / "data" / "heldout_refs"
MAXLEN = 10000


def iter_fasta(path):
    """Yield (tx_id, sequence). tx_id = header field 0 (versioned ENSMUST)."""
    name, buf = None, []
    with open(path) as fh:
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


def axis_ids(path, key="transcript_ids"):
    with h5py.File(path, "r") as h:
        return set(h[key].asstr()[:])


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    print("reading Wang transcriptome axis (psites + coverage)...", file=sys.stderr)
    on_axis = axis_ids(PSITES) & axis_ids(COVER)
    print(f"  {len(on_axis):,} tx on both psite + coverage axes", file=sys.stderr)

    seqs, biotype, n_seen, n_len, n_axis = {}, {}, 0, 0, 0
    for fa, bt in ((PC_FA, "protein_coding"), (LNC_FA, "lncRNA")):
        for tx, seq in iter_fasta(fa):
            n_seen += 1
            if len(seq) > MAXLEN:
                n_len += 1
                continue
            if tx not in on_axis:
                n_axis += 1
                continue
            # a tx id appearing in both FASTAs (shouldn't) keeps the first (pc) biotype
            if tx not in seqs:
                seqs[tx] = seq
                biotype[tx] = bt
    order = sorted(seqs)
    print(f"  parsed={n_seen:,}  dropped len>{MAXLEN}={n_len:,}  dropped off-axis={n_axis:,}  "
          f"universe={len(order):,}", file=sys.stderr)

    (OUTDIR / "mouse_wang_universe.txt").write_text("\n".join(order) + "\n")
    with (OUTDIR / "mouse_wang_universe.fa").open("w") as fh:
        for tx in order:
            fh.write(f">{tx}\n")
            s = seqs[tx]
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")
    with (OUTDIR / "mouse_wang_universe.tsv").open("w") as fh:
        fh.write("tx_id\tlength\tbiotype\n")
        npc = 0
        for tx in order:
            fh.write(f"{tx}\t{len(seqs[tx])}\t{biotype[tx]}\n")
            npc += biotype[tx] == "protein_coding"
    print(f"wrote {OUTDIR}/mouse_wang_universe.{{txt,fa,tsv}}  "
          f"pc={npc:,} lncRNA={len(order) - npc:,}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
