#!/usr/bin/env python3
"""DEPRECATED: superseded by scripts/prepare/prepare_pack.py (see scripts/prepare/README.md).
Kept only for the legacy per-tissue LOTO build runner.

Pack one tissue's per-nt P-site target + pooled RNAseq coverage for the SHARED expressed
universe into ragged 1-D numpy arrays. Generalized from pack_target_coverage.py for LOTO.

Usage: pack_target_coverage_tissue.py <Tissue>

The universe (data/fibroblast_universe.tsv, 36,668 tx) is FIXED across tissues: the FM
embeddings and the ORF track are sequence-based and tissue-independent, and per-tissue
expression enters through the coverage input, not by changing the transcript set. So every
per-tissue pack shares an identical tx_order / offsets / lengths / sum_L, which means the
one shared data/packed/orf_track*.npy and the shared per-token embeddings align 1:1 to
every tissue pack. Only target_counts, coverage, and the per-tissue global-mean depth differ.

Output (data/packed_<Tissue>/):
  tx_order.txt, offsets.npy, lengths.npy   (identical across tissues -- same universe)
  target_counts.npy   int32 (sum_L,)  this tissue's per-nt pooled P-sites
  coverage.npy        int32 (sum_L,)  this tissue's per-nt pooled RNAseq coverage
  coverage_norm.json  {"global_mean_coverage": <this tissue's mean per-nt depth>}
  pack_meta.tsv       tx_id length chrom transcript_type total_psites total_coverage
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
UNIV = NEW / "data" / "fibroblast_universe.tsv"


def universe_ids():
    ids = []
    with UNIV.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ti = hdr.index("tx_id")
        for line in fh:
            ids.append(line.rstrip("\n").split("\t")[ti])
    return sorted(set(ids))


def meta_from_summary(summary_path, want):
    d = {}
    with summary_path.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            t = f[c["tx_id"]]
            if t in want:
                d[t] = (f[c["chrom"]], f[c["transcript_type"]])
    return d


def read_rows_for(hd5_path, value_ds, want):
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
    if len(sys.argv) != 2:
        print("usage: pack_target_coverage_tissue.py <Tissue>", file=sys.stderr)
        return 2
    tissue = sys.argv[1]
    target = NEW / "data" / "target" / f"{tissue}_psites_pooled.hd5"
    target_sum = NEW / "data" / "target" / f"{tissue}_psites_summary.tsv"
    pooled_cov = NEW / "data" / "rnaseq_coverage" / f"{tissue}_rnaseq_coverage_pooled.hd5"
    out_dir = NEW / "data" / f"packed_{tissue}"
    for p in (target, target_sum, pooled_cov):
        if not p.exists():
            raise SystemExit(f"missing input for {tissue}: {p}")
    out_dir.mkdir(parents=True, exist_ok=True)

    order = sorted(set(universe_ids()))
    want = set(order)
    print(f"tissue={tissue}  universe={len(order):,} tx  out={out_dir.name}", file=sys.stderr)

    print("reading target P-site rows...", file=sys.stderr)
    tgt = read_rows_for(target, "p_sites", want)
    print("reading coverage rows...", file=sys.stderr)
    cov = read_rows_for(pooled_cov, "coverage", want)

    meta = meta_from_summary(target_sum, want)
    lengths = np.empty(len(order), dtype=np.int32)
    for k, t in enumerate(order):
        lt, lc = tgt[t].shape[0], cov[t].shape[0]
        if lt != lc:
            raise SystemExit(f"length mismatch {t}: target {lt} vs coverage {lc}")
        lengths[k] = lt
    offsets = np.zeros(len(order) + 1, dtype=np.int64)
    offsets[1:] = np.cumsum(lengths.astype(np.int64))
    total = int(offsets[-1])
    print(f"sum_L = {total:,} nt", file=sys.stderr)

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

    # per-tissue depth normalizer: global mean per-nt coverage over this tissue's pack
    global_mean = float(cov_cat.sum(dtype=np.int64)) / int(cov_cat.size)
    (out_dir / "coverage_norm.json").write_text(json.dumps({
        "global_mean_coverage": global_mean,
        "tissue": tissue,
        "sum_L": total,
        "coverage_total": int(cov_cat.sum(dtype=np.int64)),
    }, indent=2) + "\n")

    with (out_dir / "pack_meta.tsv").open("w") as m:
        m.write("tx_id\tlength\tchrom\ttranscript_type\ttotal_psites\ttotal_coverage\n")
        for k, t in enumerate(order):
            chrom, ttype = meta.get(t, ("NA", "unknown"))
            s, e = offsets[k], offsets[k + 1]
            m.write(f"{t}\t{lengths[k]}\t{chrom}\t{ttype}\t"
                    f"{int(tgt_cat[s:e].sum())}\t{int(cov_cat[s:e].sum())}\n")

    print(f"wrote {out_dir}/  global_mean_coverage={global_mean:.3f}", file=sys.stderr)
    print(f"  transcripts={len(order):,}  sum_L={total:,}  "
          f"target_total={int(tgt_cat.sum()):,}  coverage_total={int(cov_cat.sum()):,}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
