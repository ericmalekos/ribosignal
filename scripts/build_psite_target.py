#!/usr/bin/env python3
"""DEPRECATED: superseded by scripts/prepare/prepare_pack.py --ribo-psites (see
scripts/prepare/README.md). Kept only for the legacy per-tissue LOTO build runner.

Build the pooled per-nt Ribo-seq P-site TARGET for one tissue (posture A). Generalized
from build_fibroblast_psite_target.py for the leave-one-tissue-out (LOTO) experiment.

Usage: build_psite_target.py <Tissue>

Pools the tissue's RiboCode per-nucleotide P-site hd5 (transcriptome coords, GENCODE v49)
element-wise across its samples. Each input hd5 already has per-sample per-read-length
P-site offsets + sense-only strand baked in (RiboCode metaplots/pre_config), and
rRNA/tRNA/miRNA transcripts + cross-gene multimappers removed; within-gene isoform
multimappers are retained and each counted (RiboCode default = multimapper posture A).
Chrom axis is the RiboCode canonical order (identical across all tissues + samples), so
the pooled per-nt array for transcript i aligns 1:1 with the shared FM per-token embedding
and the per-tissue RNAseq coverage on the same transcript-nt axis.

Output (mirrors the RiboCode hd5 schema, so one reader serves inputs and target):
  data/target/<Tissue>_psites_pooled.hd5   (transcript_ids vlen-str + p_sites vlen-int32)
  data/target/<Tissue>_psites_summary.tsv  (tx_id length total_psites n_nonzero_nt
                                            max_psite gene_id gene_name chrom transcript_type)
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

ECH = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/biotype_probe/expression_context_human")
NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
TX2B = NEW / "data" / "tx2biotype.tsv"

SRC_GLOB = "*_psites.hd5"
CHUNK = 10000
POSTURE = "A_ribocode_isoform_multimapped"


def read_tx_ids(h):
    raw = h["transcript_ids"][:]
    return np.array([t.decode() if isinstance(t, (bytes, np.bytes_)) else str(t)
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
    if len(sys.argv) != 2:
        print("usage: build_psite_target.py <Tissue>", file=sys.stderr)
        return 2
    tissue = sys.argv[1]
    rc_dir = ECH / "data" / "ribocode_per_tissue" / tissue
    out_hd5 = NEW / "data" / "target" / f"{tissue}_psites_pooled.hd5"
    out_tsv = NEW / "data" / "target" / f"{tissue}_psites_summary.tsv"

    files = sorted(rc_dir.glob(SRC_GLOB))
    assert files, f"no {SRC_GLOB} in {rc_dir}"
    out_hd5.parent.mkdir(parents=True, exist_ok=True)
    print(f"tissue={tissue}  samples={len(files)}  dir={rc_dir}", file=sys.stderr)

    canon = None
    acc = None                       # positional list of int64 per-nt arrays
    per_sample = {}

    for k, fp in enumerate(files):
        srr = fp.name.split(".")[0].split("_")[0]
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
        print(f"[{k + 1:2d}/{len(files)}] {fp.name}  psites={per_sample[srr]:,}",
              file=sys.stderr)

    n_tx = len(canon)
    assert all(a is not None for a in acc), "some transcripts never initialised"

    maxv = max((int(a.max()) for a in acc if a.size), default=0)
    store_dt = np.int32 if maxv < 2 ** 31 else np.int64
    print(f"transcripts={n_tx:,}  max_pooled_per_nt={maxv:,}  "
          f"store_dtype={np.dtype(store_dt).name}", file=sys.stderr)

    str_dt = h5py.string_dtype("utf-8")
    obj = np.empty(n_tx, dtype=object)
    for i, a in enumerate(acc):
        obj[i] = a.astype(store_dt, copy=False)
    with h5py.File(out_hd5, "w") as o:
        o.create_dataset("transcript_ids", data=np.array(canon, dtype=object),
                         dtype=str_dt)
        o.create_dataset("p_sites", data=obj, dtype=h5py.vlen_dtype(store_dt),
                         compression="gzip", compression_opts=4)
        o.attrs["tissue"] = tissue
        o.attrs["n_samples"] = len(files)
        o.attrs["n_transcripts"] = n_tx
        o.attrs["psites_number_total"] = int(sum(v for v in per_sample.values() if v >= 0))
        o.attrs["psites_number_per_sample"] = json.dumps(per_sample)
        o.attrs["multimapper_posture"] = POSTURE
        o.attrs["source_glob"] = str(rc_dir / SRC_GLOB)
    print(f"wrote {out_hd5}", file=sys.stderr)

    tx2b = load_tx2biotype(TX2B)
    grand = 0
    nz_tx = 0
    with out_tsv.open("w") as t:
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
    print(f"wrote {out_tsv}", file=sys.stderr)
    print(f"grand_total_psites={grand:,}  transcripts_with_signal={nz_tx:,}/{n_tx:,}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
