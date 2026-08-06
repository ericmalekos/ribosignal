#!/usr/bin/env python3
"""DEPRECATED: superseded by scripts/prepare/prepare_pack.py (see scripts/prepare/README.md).
Kept only for the legacy per-tissue LOTO build runner.

Pool one tissue's per-sample RNAseq per-nt coverage hd5 into a single tissue file.
Generalized from pool_rnaseq_coverage.py for the LOTO experiment.

Usage: pool_rnaseq_coverage_tissue.py <Tissue>

The tissue's RNAseq SRRs are taken from data/loto_rnaseq_srr_tissue.tsv (SRR<TAB>tissue).
Element-wise sum across the tissue's per-sample coverage hd5 (positional accumulation with
an axis-order assert, int64 accumulate, int32 store-guard). All per-sample files share the
STAR @SQ axis, so accumulation is positional.

Output:
  data/rnaseq_coverage/<Tissue>_rnaseq_coverage_pooled.hd5   (transcript_ids + coverage)
  data/rnaseq_coverage/<Tissue>_rnaseq_coverage_summary.tsv
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
COV_DIR = NEW / "data" / "rnaseq_coverage" / "per_sample"
SRR_TISSUE = NEW / "data" / "loto_rnaseq_srr_tissue.tsv"
TX2B = NEW / "data" / "tx2biotype.tsv"

CHUNK = 10000


def read_tx_ids(h):
    raw = h["transcript_ids"][:]
    return np.array([t.decode() if isinstance(t, (bytes, np.bytes_)) else str(t)
                     for t in raw])


def load_tx2biotype(path):
    d = {}
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[col["tx_id"]]] = (f[col["gene_id"]], f[col["gene_name"]],
                                  f[col["chrom"]], f[col["transcript_type"]])
    return d


def tissue_srrs(tissue):
    srrs = []
    with SRR_TISSUE.open() as fh:
        for line in fh:
            srr, t = line.rstrip("\n").split("\t")
            if t == tissue:
                srrs.append(srr)
    return srrs


def main():
    if len(sys.argv) != 2:
        print("usage: pool_rnaseq_coverage_tissue.py <Tissue>", file=sys.stderr)
        return 2
    tissue = sys.argv[1]
    out_hd5 = NEW / "data" / "rnaseq_coverage" / f"{tissue}_rnaseq_coverage_pooled.hd5"
    out_tsv = NEW / "data" / "rnaseq_coverage" / f"{tissue}_rnaseq_coverage_summary.tsv"

    srrs = tissue_srrs(tissue)
    files = [COV_DIR / f"{s}_coverage.hd5" for s in srrs]
    missing = [str(f) for f in files if not f.exists()]
    assert not missing, f"missing per-sample coverage for {tissue}: {missing}"
    assert files, f"no RNAseq SRRs mapped to tissue {tissue} in {SRR_TISSUE}"
    print(f"tissue={tissue}  samples={len(files)}", file=sys.stderr)

    canon = None
    acc = None
    per_sample = {}

    for k, fp in enumerate(files):
        srr = fp.name.split("_")[0]
        with h5py.File(fp, "r") as h:
            tx = read_tx_ids(h)
            if canon is None:
                canon = tx
                acc = [None] * len(tx)
            else:
                assert np.array_equal(tx, canon), f"axis order mismatch in {fp.name}"
            v = h.attrs.get("mapped_sense_counted", None)
            per_sample[srr] = int(v) if v is not None else -1
            ds = h["coverage"]
            n = ds.shape[0]
            for s in range(0, n, CHUNK):
                block = ds[s:s + CHUNK]
                for j, arr in enumerate(block):
                    i = s + j
                    if acc[i] is None:
                        acc[i] = np.asarray(arr, dtype=np.int64)
                    else:
                        a = np.asarray(arr, dtype=np.int64)
                        assert a.shape[0] == acc[i].shape[0], \
                            f"per-nt length mismatch at tx#{i} ({canon[i]}) in {fp.name}"
                        acc[i] += a
        print(f"[{k + 1:2d}/{len(files)}] {fp.name}  sense_counted={per_sample[srr]:,}",
              file=sys.stderr)

    n_tx = len(canon)
    maxv = max((int(a.max()) for a in acc if a.size), default=0)
    store_dt = np.int32 if maxv < 2 ** 31 else np.int64
    print(f"transcripts={n_tx:,}  max_pooled_cov={maxv:,}  "
          f"store_dtype={np.dtype(store_dt).name}", file=sys.stderr)

    str_dt = h5py.string_dtype("utf-8")
    obj = np.empty(n_tx, dtype=object)
    for i, a in enumerate(acc):
        obj[i] = a.astype(store_dt, copy=False)
    with h5py.File(out_hd5, "w") as o:
        o.create_dataset("transcript_ids", data=np.array(canon, dtype=object), dtype=str_dt)
        o.create_dataset("coverage", data=obj, dtype=h5py.vlen_dtype(store_dt),
                         compression="gzip", compression_opts=4)
        o.attrs["tissue"] = tissue
        o.attrs["n_samples"] = len(files)
        o.attrs["libtype"] = "ISR"
        o.attrs["mapped_sense_total"] = int(sum(v for v in per_sample.values() if v >= 0))
        o.attrs["mapped_sense_per_sample"] = json.dumps(per_sample)
    print(f"wrote {out_hd5}", file=sys.stderr)

    tx2b = load_tx2biotype(TX2B)
    grand = 0
    nz_tx = 0
    with out_tsv.open("w") as t:
        t.write("tx_id\tlength\ttotal_coverage\tn_nonzero_nt\tmax_cov\t"
                "gene_id\tgene_name\tchrom\ttranscript_type\n")
        for i, a in enumerate(acc):
            tx = canon[i]
            tot = int(a.sum(dtype=np.int64))
            grand += tot
            if tot > 0:
                nz_tx += 1
            gid, gname, chrom, ttype = tx2b.get(tx, ("NA", "NA", "NA", "unknown"))
            t.write(f"{tx}\t{a.shape[0]}\t{tot}\t{int((a > 0).sum())}\t"
                    f"{int(a.max()) if a.size else 0}\t{gid}\t{gname}\t{chrom}\t{ttype}\n")
    print(f"wrote {out_tsv}", file=sys.stderr)
    print(f"grand_total_coverage={grand:,}  transcripts_with_coverage={nz_tx:,}/{n_tx:,}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
