#!/usr/bin/env python3
"""Build a tiny synthetic pack so the prediction chain can be run without any real data.

The released checkpoints need four things that are not, and cannot be, in this repository:
a packed store, an ORF-candidate track, a universe FASTA, and RNA-seq coverage. That made
the README quickstart unrunnable end to end and left nothing for CI to execute, so the
loader broke and stayed broken. This writes a miniature version of all of it -- a handful
of invented transcripts with a 5'UTR / CDS / 3'UTR layout and periodic P-site counts --
which is enough to exercise every path from `build_orf_track.py` to `dump_pred_profiles.py`
to `ribocode_dropin.build_density`.

It is a plumbing fixture, NOT data. The sequences are random, so the predictions made on
them mean nothing at all; what it demonstrates is that the code runs and the shapes line
up. Never quote a number obtained from this pack.

    python scripts/make_demo_pack.py --out data/packed_heldout_demo --n-tx 8

Writes, matching what PackedStore and heldout_test_tx read:

    <out>/tx_order.txt          one versioned transcript id per line
    <out>/lengths.npy           (N,) int64 mature transcript lengths
    <out>/offsets.npy           (N+1,) int64 cumulative starts into the ragged arrays
    <out>/coverage.npy          (sum_L,) int32 RNA-seq depth
    <out>/target_counts.npy     (sum_L,) int32 pooled P-sites, the training label
    <out>/coverage_norm.json    the pack's own global mean depth
    <out>/pack_meta.tsv         tx_id, total_psites -- what selects scorable transcripts
    <fasta>                     the matching universe FASTA (default <out>/universe.fa)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

STOPS = ("TAA", "TAG", "TGA")


def make_transcript(rng, utr5, cds_codons, utr3):
    """A sequence with exactly one clean ORF, and P-site counts that sit on it in frame.

    The CDS is built codon by codon from the 61 sense codons so no in-frame stop appears
    before the intended one -- otherwise `build_orf_track.py` would close the ORF early and
    the track would not describe the sequence the counts were drawn against.
    """
    sense = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS]
    seq = ("".join(rng.choice(list("ACGT"), size=utr5))
           + "ATG"
           + "".join(rng.choice(sense, size=cds_codons - 2))
           + rng.choice(list(STOPS))
           + "".join(rng.choice(list("ACGT"), size=utr3)))
    L = len(seq)

    # P-sites: concentrated on the CDS, every third base, with a start peak. Poisson noise
    # so the array is integer counts like a real pooled target, not a smooth curve.
    lam = np.full(L, 0.02)
    cds = np.arange(utr5, utr5 + 3 * cds_codons)
    lam[cds[cds % 3 == utr5 % 3]] = 3.0
    lam[utr5] = 12.0
    counts = rng.poisson(lam).astype(np.int32)

    # Coverage: RNA-seq spans the whole transcript with a mild 5' rise, unlike the P-sites.
    cov = rng.poisson(np.linspace(6.0, 14.0, L)).astype(np.int32)
    return seq, counts, cov


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="pack directory to create")
    ap.add_argument("--fasta", default=None, help="universe FASTA (default <out>/universe.fa)")
    ap.add_argument("--n-tx", type=int, default=8)
    ap.add_argument("--min-len", type=int, default=300)
    ap.add_argument("--max-len", type=int, default=900)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    fasta = Path(args.fasta) if args.fasta else out / "universe.fa"
    fasta.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    ids, seqs, counts, covs = [], [], [], []
    for i in range(args.n_tx):
        L = int(rng.integers(args.min_len, args.max_len + 1))
        utr5 = int(rng.integers(30, 90))
        utr3 = int(rng.integers(30, 150))
        cds_codons = max(10, (L - utr5 - utr3) // 3)
        s, c, v = make_transcript(rng, utr5, cds_codons, utr3)
        ids.append(f"ENSTDEMO{i:011d}.1")        # 21 chars: fits the npz's <U25 tx_ids field
        seqs.append(s)
        counts.append(c)
        covs.append(v)

    lengths = np.array([len(s) for s in seqs], dtype=np.int64)
    offsets = np.concatenate([[0], np.cumsum(lengths)]).astype(np.int64)
    coverage = np.concatenate(covs).astype(np.int32)
    target = np.concatenate(counts).astype(np.int32)

    np.save(out / "lengths.npy", lengths)
    np.save(out / "offsets.npy", offsets)
    np.save(out / "coverage.npy", coverage)
    np.save(out / "target_counts.npy", target)
    (out / "tx_order.txt").write_text("\n".join(ids) + "\n")
    (out / "coverage_norm.json").write_text(json.dumps(
        {"global_mean_coverage": float(coverage.sum()) / coverage.size}, indent=2) + "\n")
    with (out / "pack_meta.tsv").open("w") as fh:
        fh.write("tx_id\ttotal_psites\n")
        for t, c in zip(ids, counts):
            fh.write(f"{t}\t{int(c.sum())}\n")
    with fasta.open("w") as fh:
        for t, s in zip(ids, seqs):
            fh.write(f">{t}\n")
            for j in range(0, len(s), 60):
                fh.write(s[j:j + 60] + "\n")

    print(f"{args.n_tx} tx, {int(lengths.sum()):,} nt -> {out}", file=sys.stderr)
    print(f"  FASTA {fasta}", file=sys.stderr)
    print(f"  P-sites {int(target.sum()):,}  min/tx {int(min(c.sum() for c in counts))}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
