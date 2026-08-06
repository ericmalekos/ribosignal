#!/usr/bin/env python3
"""DEPRECATED: superseded by scripts/prepare/prepare_pack.py (see scripts/prepare/README.md).
Kept only for the legacy Fibroblast build runner.

Pack per-nt P-site target + pooled RNAseq coverage for the expressed universe into
ragged 1-D numpy arrays, so the training data loader needs only numpy (no h5py) and reads
contiguous per-transcript slices (fast, memmap-friendly) instead of random hd5 vlen rows.

For the 36,668 universe transcripts (sorted tx id = canonical pack order), reads each
transcript's per-nt P-site counts from the target hd5 and per-nt read depth from the
pooled coverage hd5 (joining BY id -- the two hd5 share the same transcript set in
different axis orders), asserts the two per-nt lengths agree, and concatenates into:

  data/packed/tx_order.txt        one tx id per line, the pack order
  data/packed/offsets.npy         int64 (N+1,)  row k = arr[offsets[k]:offsets[k+1]]
  data/packed/lengths.npy         int32 (N,)     per-transcript length (= offsets diff)
  data/packed/target_counts.npy   int32 (sum_L,) concatenated per-nt P-site counts
  data/packed/coverage.npy        int32 (sum_L,) concatenated per-nt pooled coverage
  data/packed/pack_meta.tsv       tx_id, length, chrom, transcript_type, total_psites,
                                  total_coverage

Reads rows in hd5-index order (not tx-sorted order) for gzip chunk locality, then
reassembles in the canonical sorted-tx order.

Usage: pack_target_coverage.py [COVERAGE_HD5] [LIMIT]
  COVERAGE_HD5 defaults to the pooled file; pass a per-sample file to dry-run the code
  before pooling finishes (schema is identical).
  LIMIT (int) packs only the first LIMIT universe tx (sorted) into data/packed_dryrun/
  for a fast chain test; omit for the full universe into data/packed/.
"""
import sys
from pathlib import Path

import h5py
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
UNIV = NEW / "data" / "fibroblast_universe.tsv"
TARGET = NEW / "data" / "target" / "Fibroblast_psites_pooled.hd5"
TARGET_SUM = NEW / "data" / "target" / "Fibroblast_psites_summary.tsv"
POOLED_COV = NEW / "data" / "rnaseq_coverage" / "Fibroblast_rnaseq_coverage_pooled.hd5"
OUT = NEW / "data" / "packed"


def universe_ids():
    ids = []
    with UNIV.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ti = hdr.index("tx_id")
        for line in fh:
            ids.append(line.rstrip("\n").split("\t")[ti])
    return sorted(set(ids))


def meta_from_summary(want):
    """tx_id -> (chrom, transcript_type) for the wanted set, from the target summary."""
    d = {}
    with TARGET_SUM.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            t = f[c["tx_id"]]
            if t in want:
                d[t] = (f[c["chrom"]], f[c["transcript_type"]])
    return d


def read_rows_for(hd5_path, value_ds, want):
    """Return dict tx_id -> per-nt int32 array, reading in hd5-index order for locality."""
    with h5py.File(hd5_path, "r") as h:
        ids = h["transcript_ids"].asstr()[:]
        idx = {t: i for i, t in enumerate(ids)}
        picks = sorted((idx[t], t) for t in want if t in idx)
        missing = [t for t in want if t not in idx]
        if missing:
            raise SystemExit(f"{len(missing)} universe tx missing from {hd5_path.name}, "
                             f"e.g. {missing[:3]}")
        ds = h[value_ds]
        out = {}
        for k, (i, t) in enumerate(picks):
            out[t] = np.asarray(ds[i], dtype=np.int32)
            if (k + 1) % 5000 == 0:
                print(f"  {hd5_path.name}: {k + 1}/{len(picks)}", file=sys.stderr)
    return out


def main():
    cov_path = Path(sys.argv[1]) if len(sys.argv) > 1 else POOLED_COV
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    if not cov_path.exists():
        raise SystemExit(f"coverage hd5 not found: {cov_path}")
    out_dir = (NEW / "data" / "packed_dryrun") if limit else OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    order = sorted(set(universe_ids()))
    if limit:
        order = order[:limit]
    want = set(order)
    print(f"universe: {len(order):,} tx  coverage={cov_path.name}  out={out_dir.name}",
          file=sys.stderr)

    print("reading target P-site rows...", file=sys.stderr)
    tgt = read_rows_for(TARGET, "p_sites", want)
    print("reading coverage rows...", file=sys.stderr)
    cov = read_rows_for(cov_path, "coverage", want)

    meta = meta_from_summary(want)
    lengths = np.empty(len(order), dtype=np.int32)
    for k, t in enumerate(order):
        lt, lc = tgt[t].shape[0], cov[t].shape[0]
        if lt != lc:
            raise SystemExit(f"length mismatch {t}: target {lt} vs coverage {lc}")
        lengths[k] = lt
    offsets = np.zeros(len(order) + 1, dtype=np.int64)
    offsets[1:] = np.cumsum(lengths.astype(np.int64))
    total = int(offsets[-1])
    print(f"sum_L = {total:,} nt  (int32 arrays ~ {2 * total * 4 / 1e9:.2f} GB)",
          file=sys.stderr)

    tgt_cat = np.empty(total, dtype=np.int32)
    cov_cat = np.empty(total, dtype=np.int32)
    for k, t in enumerate(order):
        s, e = offsets[k], offsets[k + 1]
        tgt_cat[s:e] = tgt[t]
        cov_cat[s:e] = cov[t]

    np.save(out_dir / "offsets.npy", offsets)
    np.save(out_dir / "lengths.npy", lengths)
    np.save(out_dir / "target_counts.npy", tgt_cat)
    np.save(out_dir / "coverage.npy", cov_cat)
    (out_dir / "tx_order.txt").write_text("\n".join(order) + "\n")
    with (out_dir / "pack_meta.tsv").open("w") as m:
        m.write("tx_id\tlength\tchrom\ttranscript_type\ttotal_psites\ttotal_coverage\n")
        for k, t in enumerate(order):
            chrom, ttype = meta.get(t, ("NA", "unknown"))
            s, e = offsets[k], offsets[k + 1]
            m.write(f"{t}\t{lengths[k]}\t{chrom}\t{ttype}\t"
                    f"{int(tgt_cat[s:e].sum())}\t{int(cov_cat[s:e].sum())}\n")

    print(f"wrote {out_dir}/ : offsets, lengths, target_counts, coverage, tx_order, "
          f"pack_meta", file=sys.stderr)
    print(f"  transcripts={len(order):,}  sum_L={total:,}  "
          f"target_total={int(tgt_cat.sum()):,}  coverage_total={int(cov_cat.sum()):,}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
