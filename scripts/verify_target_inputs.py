#!/usr/bin/env python3
"""Step 1: verify the cleaning posture of the Fibroblast P-site target inputs.

Confirms, WITHOUT rebuilding anything, that the 32 RiboCode per-nt P-site hd5 that
will be pooled into the target are already ncRNA-clean and carry real signal:

  1. rRNA/tRNA/miRNA/etc. drop-list transcripts (ncrna_filter_tx.txt) that are ON the
     transcriptome axis should pool to ~0 P-sites (reads over them were dropped upstream).
  2. Housekeeping positive controls (ACTB, GAPDH) should pool to very large counts
     (proves the reader sees signal, i.e. the ~0 above is real, not an all-zero bug).
  3. Per-tx hd5 array length == GENCODE mature length (tx2biotype) for the controls
     (proves the per-nt axis is the transcript-nt axis embeddings/RNAseq align to).

Light I/O: reads transcript_ids once + fancy-indexes only the drop-list + control rows
across the 32 files. Safe on the head node. Writes a compact report to stdout.
"""
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

ECH = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/biotype_probe/expression_context_human")
NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
RC_DIR = NEW / "data" / "ribocode_per_tissue" / "Fibroblast"
NCRNA = NEW / "data" / "ncrna_tx_human_v49.txt"
TX2B = NEW / "data" / "tx2biotype.tsv"

CONTROLS = ("ACTB", "GAPDH")   # by gene_name; pool all isoforms present on the axis


def read_tx_ids(h):
    raw = h["transcript_ids"][:]
    return np.array([t.decode() if isinstance(t, bytes | np.bytes_) else str(t)
                     for t in raw])


def load_tx2biotype(path):
    """tx_id -> dict(gene_name, gene_type, transcript_type, length)."""
    d = {}
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[col["tx_id"]]] = {
                "gene_name": f[col["gene_name"]],
                "gene_type": f[col["gene_type"]],
                "transcript_type": f[col["transcript_type"]],
                "length": int(f[col["length"]]),
            }
    return d


def pool_rows(files, idx_sorted):
    """P-sites per selected row (idx_sorted), pooled across all files.

    Returns (pooled_total_per_row int64[n], lengths_from_first_file int64[n]).
    """
    idx = np.asarray(idx_sorted, dtype=np.int64)
    n = len(idx)
    pooled = np.zeros(n, dtype=np.int64)
    lengths = np.zeros(n, dtype=np.int64)
    for k, fp in enumerate(files):
        if n == 0:
            break
        with h5py.File(fp, "r") as h:
            rows = h["p_sites"][idx]
        if k == 0:
            lengths = np.array([r.shape[0] for r in rows], dtype=np.int64)
        pooled += np.array([int(r.sum()) for r in rows], dtype=np.int64)
    return pooled, lengths


def main():
    files = sorted(RC_DIR.glob("SRR15513*_psites.hd5"))
    assert len(files) == 32, f"expected 32 Fibroblast hd5, found {len(files)}"

    tx2b = load_tx2biotype(TX2B)

    with h5py.File(files[0], "r") as h:
        tx_arr = read_tx_ids(h)
    pos = {tx: i for i, tx in enumerate(tx_arr)}

    # per-sample sequencing depth (P-site reads), for context
    depths = []
    for fp in files:
        with h5py.File(fp, "r") as h:
            v = h.attrs.get("psites_number", None)
        depths.append(int(v) if v is not None else -1)
    depth_total = sum(d for d in depths if d >= 0)

    # ---- drop-list (ncRNA) partition ----
    drop = [ln.strip() for ln in NCRNA.read_text().splitlines() if ln.strip()]
    present = [tx for tx in drop if tx in pos]
    absent = [tx for tx in drop if tx not in pos]
    present_idx = sorted(pos[tx] for tx in present)
    nc_pooled, _ = pool_rows(files, present_idx)
    nc_types = Counter(tx2b.get(tx, {}).get("gene_type", "unknown") for tx in present)

    # ---- positive controls ----
    ctrl_tx = {g: [] for g in CONTROLS}
    for tx, m in tx2b.items():
        if m["gene_name"] in ctrl_tx and tx in pos:
            ctrl_tx[m["gene_name"]].append(tx)

    # ---- report ----
    out = []
    out.append("=" * 70)
    out.append("STEP 1  Fibroblast P-site target input verification")
    out.append("=" * 70)
    out.append(f"input hd5 files            : {len(files)}  ({files[0].parent})")
    out.append(f"axis transcripts           : {len(tx_arr):,}")
    out.append(f"per-sample P-site depth    : min {min(depths):,}  "
               f"max {max(depths):,}  total {depth_total:,}")
    out.append("")
    out.append("-- ncRNA drop-list (rRNA/tRNA/miRNA/sno/sn...) posture --")
    out.append(f"drop-list transcripts      : {len(drop):,}")
    out.append(f"  present on axis          : {len(present):,}  (expect present-but-zero)")
    out.append(f"  absent from axis         : {len(absent):,}  (filtered at read level, "
               f"never in GENCODE tx set)")
    out.append("present ncRNA by gene_type : " +
               ", ".join(f"{t}={c}" for t, c in nc_types.most_common()))
    if len(present):
        n_nonzero = int((nc_pooled > 0).sum())
        out.append(f"pooled P-sites over present ncRNA : sum={int(nc_pooled.sum()):,}  "
                   f"max_single_tx={int(nc_pooled.max()):,}  nonzero_tx={n_nonzero}/{len(present)}")
        if n_nonzero:
            order = np.argsort(nc_pooled)[::-1][:5]
            worst = [(tx_arr[present_idx[i]], int(nc_pooled[i])) for i in order]
            out.append("  top residual ncRNA (tx, pooled): " +
                       "; ".join(f"{t}={c}" for t, c in worst))
    out.append("")
    out.append("-- positive controls (housekeeping, expect very large) --")
    for g in CONTROLS:
        txs = ctrl_tx[g]
        if not txs:
            out.append(f"{g:<6}: no isoforms found on axis (check gene_name)")
            continue
        idx = sorted(pos[tx] for tx in txs)
        pooled, lengths = pool_rows(files, idx)
        top = int(np.argmax(pooled))
        top_tx = tx_arr[idx[top]]
        # length consistency: hd5 array length vs GENCODE mature length
        len_ok = sum(1 for tx, L in zip([tx_arr[i] for i in idx], lengths, strict=False)
                     if tx2b.get(tx, {}).get("length", -1) == int(L))
        out.append(f"{g:<6}: isoforms={len(txs)}  pooled_total={int(pooled.sum()):,}  "
                   f"top_isoform={top_tx} ({int(pooled[top]):,})  "
                   f"len_match(hd5==GENCODE)={len_ok}/{len(txs)}")
    out.append("=" * 70)

    report = "\n".join(out)
    print(report)


if __name__ == "__main__":
    main()
