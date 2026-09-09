#!/usr/bin/env python3
"""Pool per-sample coverage into the packed store the model reads.

The chain README section 2 describes stops one step short: `rnaseq_coverage.py` writes a
per-sample hd5, `dump_pred_profiles.py` reads a PACK, and nothing published here turned one
into the other. This does.

A pack is six ragged arrays sharing one transcript order, so a transcript's coverage, its
P-site target and its ORF-track rows are all the same slice `offsets[r]:offsets[r+1]`. See
`scripts/dataset.py:PackedStore`, which is the reader this has to satisfy.

    # inference on a new sample: RNA-seq coverage only, no ribosome profiling
    python scripts/build_pack.py --fasta transcripts.fa --coverage cov.hd5 \\
           --out data/packed_heldout_mysample

    # an evaluation arm, where real P-sites exist to score against
    python scripts/build_pack.py --fasta transcripts.fa --coverage cov1.hd5 cov2.hd5 \\
           --psites psites.hd5 --out data/packed_heldout_mysample

Several --coverage files are SUMMED, which is what pooling replicates of one condition
means; they must be the same transcriptome. The transcript set is the intersection of the
FASTA and every hd5, restricted to transcripts whose lengths agree everywhere -- a
disagreement means the BAM and the FASTA came from different annotation releases, and a
pack built through it would shift every profile silently, so those transcripts are dropped
and counted rather than coerced.

WITHOUT --psites there is no ribosome-profiling target: target_counts is all zeros and
pack_meta.tsv records total_psites = 0. That is the honest state of an inference-only
sample, and it means `heldout_test_tx()`'s signal floor would reject every transcript, so
this also writes `expressed_tx.txt` -- pass it to dump_pred_profiles.py as --tx_list.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def read_fasta_lengths(path):
    """{tx_id: length}, streaming, without holding the sequences."""
    lens, name, n = {}, None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lens[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name is not None:
        lens[name] = n
    return lens


def read_hd5(path):
    """{tx_id: per-nt int64 array} from an rnaseq_coverage.py-schema file."""
    import h5py
    with h5py.File(path, "r") as f:
        ids = [t.decode() if isinstance(t, bytes) else str(t) for t in f["transcript_ids"][:]]
        key = "coverage" if "coverage" in f else "psites"
        if key not in f:
            sys.exit(f"{path}: expected a 'coverage' or 'psites' dataset, found "
                     f"{sorted(f.keys())}")
        arrs = f[key][:]
    return {t: np.asarray(a, dtype=np.int64) for t, a in zip(ids, arrs)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True,
                    help="transcriptome FASTA; defines the transcript set and their lengths")
    ap.add_argument("--coverage", required=True, nargs="+",
                    help="one or more rnaseq_coverage.py hd5 files, summed")
    ap.add_argument("--psites", nargs="*", default=[],
                    help="optional per-nt P-site hd5 (same schema), summed. Omit for inference.")
    ap.add_argument("--out", required=True, help="pack directory to write")
    ap.add_argument("--max-length", type=int, default=10000,
                    help="drop transcripts longer than this (default %(default)s). NOT cosmetic: "
                         "the released checkpoints were trained on a universe capped at 10,000 nt, "
                         "and the transformer's attention is O(L^2). chr22's longest transcript is "
                         "37,852 nt, whose attention matrix is 5.3 GiB per head and 42.7 GiB "
                         "across 8, so a single such transcript exhausts any GPU. Set 0 to "
                         "disable, and expect an out-of-memory abort if you do.")
    ap.add_argument("--min-coverage", type=int, default=1,
                    help="a transcript joins expressed_tx.txt at this total coverage or more "
                         "(default %(default)s)")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    fa = read_fasta_lengths(args.fasta)
    print(f"FASTA {args.fasta}: {len(fa):,} transcripts", file=sys.stderr)

    if args.max_length:
        too_long = {t for t, n in fa.items() if n > args.max_length}
        if too_long:
            longest = max(fa[t] for t in too_long)
            print(f"dropping {len(too_long):,} transcripts longer than {args.max_length:,} nt "
                  f"(longest {longest:,}); the released checkpoints never saw one",
                  file=sys.stderr)
            fa = {t: n for t, n in fa.items() if t not in too_long}

    keep = set(fa)
    sources = [("coverage", p) for p in args.coverage] + [("psites", p) for p in args.psites]
    loaded, n_lenmm = [], 0
    for kind, path in sources:
        d = read_hd5(path)
        mism = {t for t, a in d.items() if t in fa and a.shape[0] != fa[t]}
        n_lenmm = max(n_lenmm, len(mism))
        print(f"{kind} {path}: {len(d):,} tx, {len(mism):,} length mismatches vs the FASTA",
              file=sys.stderr)
        keep &= (set(d) - mism)
        loaded.append((kind, d))
    if not keep:
        sys.exit("no transcript is present, with a consistent length, in the FASTA and every "
                 "input file. The BAM and the FASTA are probably from different annotation "
                 "releases.")
    order = sorted(keep)
    if n_lenmm:
        print(f"dropped transcripts whose length disagreed (max {n_lenmm:,} in one file)",
              file=sys.stderr)

    lengths = np.array([fa[t] for t in order], dtype=np.int64)
    offsets = np.concatenate([[0], np.cumsum(lengths)]).astype(np.int64)
    sum_L = int(offsets[-1])
    coverage = np.zeros(sum_L, dtype=np.int64)
    target = np.zeros(sum_L, dtype=np.int64)
    for kind, d in loaded:
        dst = coverage if kind == "coverage" else target
        for r, t in enumerate(order):
            dst[offsets[r]:offsets[r + 1]] += d[t]

    # int32 is the reader's dtype; a pooled deep sample can overflow it, so say so rather
    # than wrapping negative and poisoning every downstream number.
    for name, arr in (("coverage", coverage), ("target_counts", target)):
        if arr.max(initial=0) > np.iinfo(np.int32).max:
            sys.exit(f"{name} exceeds int32 ({arr.max():,}); pool fewer samples per pack")

    np.save(out / "lengths.npy", lengths)
    np.save(out / "offsets.npy", offsets)
    np.save(out / "coverage.npy", coverage.astype(np.int32))
    np.save(out / "target_counts.npy", target.astype(np.int32))
    (out / "tx_order.txt").write_text("\n".join(order) + "\n")
    (out / "coverage_norm.json").write_text(json.dumps(
        {"global_mean_coverage": float(coverage.sum()) / sum_L}, indent=2) + "\n")

    per_tx_cov = np.add.reduceat(coverage, offsets[:-1])
    per_tx_psi = np.add.reduceat(target, offsets[:-1])
    with (out / "pack_meta.tsv").open("w") as fh:
        fh.write("tx_id\tlength\ttotal_coverage\ttotal_psites\n")
        for r, t in enumerate(order):
            fh.write(f"{t}\t{int(lengths[r])}\t{int(per_tx_cov[r])}\t{int(per_tx_psi[r])}\n")

    expressed = [t for r, t in enumerate(order) if per_tx_cov[r] >= args.min_coverage]
    (out / "expressed_tx.txt").write_text("\n".join(expressed) + "\n")

    print(f"\n{len(order):,} tx, {sum_L:,} nt -> {out}", file=sys.stderr)
    print(f"  coverage total {int(coverage.sum()):,}  mean/nt "
          f"{coverage.sum() / sum_L:.3f}", file=sys.stderr)
    print(f"  P-sites total  {int(target.sum()):,}"
          f"{'  (inference-only pack: no --psites given)' if not args.psites else ''}",
          file=sys.stderr)
    print(f"  expressed_tx.txt: {len(expressed):,} tx at coverage >= {args.min_coverage}",
          file=sys.stderr)
    print("\nNext: build the ORF track for this pack, then predict.\n"
          f"  python scripts/build_orf_track.py --mode ext --kozak none "
          f"--pack {out} --fasta {args.fasta}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
