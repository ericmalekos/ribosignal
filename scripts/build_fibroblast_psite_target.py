#!/usr/bin/env python3
"""DEPRECATED: superseded by scripts/prepare/prepare_pack.py --ribo-psites (see
scripts/prepare/README.md). Kept only for the legacy Fibroblast build runner.

Step 2: build the pooled per-nt Ribo-seq P-site TARGET for Fibroblast (posture A).

Pools the 32 RiboCode per-nucleotide P-site hd5 (transcriptome coords, GENCODE v49)
element-wise across samples. Each input hd5 already has per-sample per-read-length
P-site offsets + sense-only strand baked in (RiboCode metaplots/pre_config), and
rRNA/tRNA/miRNA transcripts + cross-gene multimappers removed; within-gene isoform
multimappers are retained and each counted (RiboCode default = multimapper posture A,
user-confirmed).

Output (mirrors the RiboCode hd5 schema, so one reader serves inputs and target):
  data/target/Fibroblast_psites_pooled.hd5
     transcript_ids  vlen-str    (versioned ENST, canonical axis order)
     p_sites         vlen-int32  (per-nt pooled counts, 0-based, length == tx length)
     attrs: tissue, n_samples, n_transcripts, psites_number_total,
            psites_number_per_sample (JSON srr->int), multimapper_posture, source_glob
  data/target/Fibroblast_psites_summary.tsv
     tx_id length total_psites n_nonzero_nt max_psite gene_id gene_name chrom transcript_type

The per-nt array for transcript i aligns 1:1 (0-based, position-for-position) with the
RiNALMo per-token embedding {tx}_tokens.npy of shape (L,1280) and with future per-nt
RNAseq coverage on the same transcript-nt axis.
"""
import json
import os
import sys
from pathlib import Path

import h5py
import numpy as np

# ECH removed 2026-08-19: the deleted biotype_probe tree. Consumers repointed to NEW.
# Root resolves via $RIBOSEQ_SIGNAL_MODEL_ROOT or auto-detection (task #94); never baked in.
sys.path.insert(0, str(Path(__file__).resolve().parent / "prepare"))
from paths import project_root  # noqa: E402

NEW = project_root()
RC_DIR = NEW / "data" / "ribocode_per_tissue" / "Fibroblast"
TX2B = NEW / "data" / "tx2biotype.tsv"
OUT_HD5 = NEW / "data" / "target" / "Fibroblast_psites_pooled.hd5"
OUT_TSV = NEW / "data" / "target" / "Fibroblast_psites_summary.tsv"

SRC_GLOB = os.environ.get("PSITE_GLOB", "SRR15513*_psites.hd5")
CHUNK = 10000
TISSUE = "Fibroblast"
POSTURE = "A_ribocode_isoform_multimapped"


def read_tx_ids(h):
    raw = h["transcript_ids"][:]
    return np.array([t.decode() if isinstance(t, bytes | np.bytes_) else str(t)
                     for t in raw])


def load_tx2biotype(path):
    """tx_id -> (gene_id, gene_name, chrom, transcript_type)."""
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
    files = sorted(RC_DIR.glob(SRC_GLOB))
    # De-fragilised 2026-08-19 (task #94): was a hard `assert len(files) == 32`, which made this
    # script unusable on any cohort but the 32-sample Chothani Fibroblast set. Override the glob
    # with $PSITE_GLOB and the expected count with $PSITE_EXPECT_N; defaults reproduce the
    # historical behaviour exactly.
    if not files:
        sys.exit(f"no P-site hd5 matching {SRC_GLOB!r} under {SRC_DIR}")
    _expect = os.environ.get("PSITE_EXPECT_N", "32")
    if _expect and len(files) != int(_expect):
        print(f"WARNING: expected {_expect} P-site hd5, found {len(files)}", file=sys.stderr)
    OUT_HD5.parent.mkdir(parents=True, exist_ok=True)

    canon = None
    acc = None                       # positional list of int64 per-nt arrays
    per_sample = {}

    for k, fp in enumerate(files):
        srr = fp.name.split(".")[0]
        with h5py.File(fp, "r") as h:
            tx = read_tx_ids(h)
            if canon is None:
                canon = tx
                acc = [None] * len(tx)
            else:
                assert np.array_equal(tx, canon), f"axis order mismatch in {fp.name}"
            v = h.attrs.get("psites_number", None)
            per_sample[srr] = int(v) if v is not None else -1
            ds = h["p_sites"]
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
        print(f"[{k + 1:2d}/32] {fp.name}  psites={per_sample[srr]:,}", file=sys.stderr)

    n_tx = len(canon)
    assert all(a is not None for a in acc), "some transcripts never initialised"

    # storage dtype: keep int32 unless a single pooled nt exceeds int32 (it will not)
    maxv = max((int(a.max()) for a in acc if a.size), default=0)
    store_dt = np.int32 if maxv < 2 ** 31 else np.int64
    print(f"transcripts={n_tx:,}  max_pooled_per_nt={maxv:,}  "
          f"store_dtype={np.dtype(store_dt).name}", file=sys.stderr)

    # ---- write pooled hd5 (RiboCode-compatible schema) ----
    str_dt = h5py.string_dtype("utf-8")
    obj = np.empty(n_tx, dtype=object)
    for i, a in enumerate(acc):
        obj[i] = a.astype(store_dt, copy=False)
    with h5py.File(OUT_HD5, "w") as o:
        o.create_dataset("transcript_ids", data=np.array(canon, dtype=object),
                         dtype=str_dt)
        o.create_dataset("p_sites", data=obj, dtype=h5py.vlen_dtype(store_dt),
                         compression="gzip", compression_opts=4)
        o.attrs["tissue"] = TISSUE
        o.attrs["n_samples"] = len(files)
        o.attrs["n_transcripts"] = n_tx
        o.attrs["psites_number_total"] = int(sum(v for v in per_sample.values() if v >= 0))
        o.attrs["psites_number_per_sample"] = json.dumps(per_sample)
        o.attrs["multimapper_posture"] = POSTURE
        o.attrs["source_glob"] = str(RC_DIR / SRC_GLOB)
    print(f"wrote {OUT_HD5}", file=sys.stderr)

    # ---- write per-transcript summary TSV ----
    tx2b = load_tx2biotype(TX2B)
    grand = 0
    nz_tx = 0
    with OUT_TSV.open("w") as t:
        t.write("tx_id\tlength\ttotal_psites\tn_nonzero_nt\tmax_psite\t"
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
    print(f"wrote {OUT_TSV}", file=sys.stderr)
    print(f"grand_total_psites={grand:,}  transcripts_with_signal={nz_tx:,}/{n_tx:,}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
