#!/usr/bin/env python3
"""SUPERSEDED by scripts/prepare/prepare_pack.py (see scripts/prepare/README.md), but still
called by pool_rnaseq_coverage.sbatch and pool_rnaseq_coverage_tissue.py, so it is maintained
rather than retired.

De-fragilised 2026-08-19 (task #94): the hardcoded project root, the `SRR15513*` glob and the
hard `assert len(files) == 32` are gone. Every path plus the glob is now a CLI argument, the
root resolves via $RIBOSEQ_SIGNAL_MODEL_ROOT or auto-detection, and a wrong file count is a
warning via --expect-n rather than a crash. Defaults reproduce the historical Fibroblast
32-sample behaviour exactly, so existing callers are unaffected.

Pool the per-sample RNAseq per-nt coverage hd5 for one tissue into a single file.

Element-wise sum across samples (positional accumulation with an axis-order assert,
int64 accumulate, int32 store-guard), mirroring build_fibroblast_psite_target.py. All
per-sample files share the STAR @SQ axis (same index), so accumulation is positional.

Output:
  data/rnaseq_coverage/Fibroblast_rnaseq_coverage_pooled.hd5
     transcript_ids  vlen-str    (STAR @SQ order; join to the P-site target by tx-id)
     coverage        vlen-int32  (per-nt pooled read depth, 0-based, length == tx length)
     attrs: tissue, n_samples, libtype, mapped_sense_total,
            mapped_sense_per_sample (JSON srr->int)
  data/rnaseq_coverage/Fibroblast_rnaseq_coverage_summary.tsv
     tx_id length total_coverage n_nonzero_nt max_cov gene_id gene_name chrom transcript_type
"""
import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np

# Project root is resolved, never baked in: $RIBOSEQ_SIGNAL_MODEL_ROOT, else auto-detected by
# walking up to the ancestor holding scripts/ and data/. Same helper the prepare/ toolkit uses.
sys.path.insert(0, str(Path(__file__).resolve().parent / "prepare"))
from paths import project_root  # noqa: E402

CHUNK = 10000


def build_args():
    """Every path and the sample glob are CLI arguments.

    Defaults reproduce the historical Fibroblast/32-sample behaviour EXACTLY, so the existing
    callers (pool_rnaseq_coverage.sbatch, pool_rnaseq_coverage_tissue.py) keep working unchanged.
    """
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=None,
                    help="project root (default: $RIBOSEQ_SIGNAL_MODEL_ROOT or auto-detect)")
    ap.add_argument("--tissue", default="Fibroblast")
    ap.add_argument("--cov-dir", type=Path, default=None,
                    help="dir of per-sample *_coverage.hd5 (default <root>/data/rnaseq_coverage/per_sample)")
    ap.add_argument("--tx2biotype", type=Path, default=None,
                    help="default <root>/data/tx2biotype.tsv")
    ap.add_argument("--out-hd5", type=Path, default=None)
    ap.add_argument("--out-tsv", type=Path, default=None)
    ap.add_argument("--glob", default="SRR15513*_coverage.hd5",
                    help="per-sample coverage glob (default matches the Chothani Fibroblast runs)")
    ap.add_argument("--expect-n", type=int, default=None,
                    help="if given, WARN when the file count differs. This replaces a hard "
                         "`assert len(files)==32`, which made the script unusable on any other "
                         "dataset and is exactly the fragility task #94 exists to remove.")
    a = ap.parse_args()
    root = a.root or project_root()
    a.cov_dir = a.cov_dir or root / "data" / "rnaseq_coverage" / "per_sample"
    a.tx2biotype = a.tx2biotype or root / "data" / "tx2biotype.tsv"
    a.out_hd5 = a.out_hd5 or root / "data" / "rnaseq_coverage" / f"{a.tissue}_rnaseq_coverage_pooled.hd5"
    a.out_tsv = a.out_tsv or root / "data" / "rnaseq_coverage" / f"{a.tissue}_rnaseq_coverage_summary.tsv"
    return a


def read_tx_ids(h):
    raw = h["transcript_ids"][:]
    return np.array([t.decode() if isinstance(t, bytes | np.bytes_) else str(t)
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


def main():
    a = build_args()
    files = sorted(a.cov_dir.glob(a.glob))
    if not files:
        sys.exit(f"no per-sample coverage hd5 matching {a.glob!r} under {a.cov_dir}")
    if a.expect_n is not None and len(files) != a.expect_n:
        print(f"WARNING: expected {a.expect_n} per-sample coverage hd5, found {len(files)}",
              file=sys.stderr)
    print(f"[{a.tissue}] pooling {len(files)} per-sample coverage hd5 from {a.cov_dir}")

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
        print(f"[{k + 1:2d}/32] {fp.name}  sense_counted={per_sample[srr]:,}",
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
    with h5py.File(a.out_hd5, "w") as o:
        o.create_dataset("transcript_ids", data=np.array(canon, dtype=object), dtype=str_dt)
        o.create_dataset("coverage", data=obj, dtype=h5py.vlen_dtype(store_dt),
                         compression="gzip", compression_opts=4)
        o.attrs["tissue"] = a.tissue
        o.attrs["n_samples"] = len(files)
        o.attrs["libtype"] = "ISR"
        o.attrs["mapped_sense_total"] = int(sum(v for v in per_sample.values() if v >= 0))
        o.attrs["mapped_sense_per_sample"] = json.dumps(per_sample)
    print(f"wrote {a.out_hd5}", file=sys.stderr)

    tx2b = load_tx2biotype(a.tx2biotype)
    grand = 0
    nz_tx = 0
    with a.out_tsv.open("w") as t:
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
    print(f"wrote {a.out_tsv}", file=sys.stderr)
    print(f"grand_total_coverage={grand:,}  transcripts_with_coverage={nz_tx:,}/{n_tx:,}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
