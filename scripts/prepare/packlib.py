#!/usr/bin/env python3
"""Path-free pool + pack primitives for riboseq_signal_model model-input packs.

Lifted and generalized from scripts/heldout/build_heldout_pack.py (pool_universe + the
7-file writer) and scripts/pack_target_coverage_tissue.py (the 6-column pack_meta variant).
NO hardcoded project paths, tissue names, or sample-count asserts -- everything is a function
argument, so one library serves the training packs, the held-out packs, and any future
user-supplied dataset.

The pack contract this writes (verified against on-disk bytes; all row-aligned to tx_order):
  tx_order.txt        n_tx lines, versioned tx id per line (defines row order)
  offsets.npy         int64 (n_tx+1,)   [0, cumsum(lengths)], offsets[-1] == sum_L
  lengths.npy         int32 (n_tx,)     per-tx length
  target_counts.npy   int32 (sum_L,)    concatenated per-nt P-sites (zeros if coverage-only)
  coverage.npy        int32 (sum_L,)    concatenated per-nt RNA depth
  coverage_norm.json  JSON              at least {"global_mean_coverage": <float>} (load-bearing)
  pack_meta.tsv       header + per-tx   >=4 cols (tx_id,length,total_psites,total_coverage);
                                        6 cols when a tx->biotype map is supplied
The ORF track is NOT written here (it is a separate sequence-derived artifact built by
build_orf_track.py against this pack's tx_order/offsets, located at load time via
$RIBO_ORF_TRACK).
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np

SCORABLE_MIN_PSITES = 50   # a tx is "scorable" for held-out eval if pooled P-sites >= this


def _log(msg: str) -> None:
    print(msg, file=sys.stderr)


def read_ids(h) -> np.ndarray:
    """Versioned transcript ids from an hd5 (RiboCode psites or rnaseq_coverage output)."""
    return h["transcript_ids"].asstr()[:]


def resolve_hd5(specs, pattern: str):
    """Expand a list of file/dir path specs into a sorted list of hd5 files.

    Each spec may be a single hd5 file (used verbatim) or a directory (globbed by `pattern`).
    Accepts str or Path. Returns [] for an empty/None spec.
    """
    if not specs:
        return []
    files = []
    for spec in specs:
        p = Path(spec)
        if p.is_dir():
            files.extend(sorted(p.glob(pattern)))
        elif p.exists():
            files.append(p)
        else:
            raise FileNotFoundError(f"input path does not exist: {p}")
    return files


def pool_per_nt(files, value_key: str, want: set):
    """tx -> int64 pooled per-nt array, summed across `files`, restricted to `want`.

    Joins strictly by versioned tx id (the coverage hd5 @SQ order may differ from the
    RiboCode order). Reads rows in on-disk order per file (sorted picks) so vlen reads stay
    sequential. Raises ValueError on a per-nt length mismatch for the same tx across files.
    """
    acc = {}
    for k, fp in enumerate(files):
        fp = Path(fp)
        with h5py.File(fp, "r") as h:
            ids = read_ids(h)
            idx = {t: i for i, t in enumerate(ids)}
            picks = sorted((idx[t], t) for t in want if t in idx)
            ds = h[value_key]
            for n, (i, t) in enumerate(picks):
                a = np.asarray(ds[i], dtype=np.int64)
                if t in acc:
                    if acc[t].shape[0] != a.shape[0]:
                        raise ValueError(f"per-nt length mismatch {t} in {fp.name}: "
                                         f"{a.shape[0]} vs {acc[t].shape[0]}")
                    acc[t] += a
                else:
                    acc[t] = a
                if (n + 1) % 8000 == 0:
                    _log(f"  [{k + 1}/{len(files)}] {fp.name}: {n + 1}/{len(picks)}")
        _log(f"pooled {value_key} from {fp.name} ({len(picks)} universe rows)")
    return acc


def load_universe(universe_tx=None, ref_pack=None):
    """Return (order, ref_len_of).

    ref_pack: reuse an existing pack's tx_order + lengths verbatim (so a shared ORF track /
      embeddings / one-hot align 1:1); ref_len_of = {tx: length} enforces annotation identity.
    universe_tx: a fresh tx-id list (one versioned id per line); ref_len_of = None (offsets
      computed fresh from the pooled arrays).
    Exactly one of the two must be given.
    """
    if (universe_tx is None) == (ref_pack is None):
        raise ValueError("provide exactly one of universe_tx / ref_pack")
    if ref_pack is not None:
        ref_pack = Path(ref_pack)
        order = (ref_pack / "tx_order.txt").read_text().split()
        ref_len = np.load(ref_pack / "lengths.npy")
        ref_len_of = {t: int(ref_len[i]) for i, t in enumerate(order)}
    else:
        order = Path(universe_tx).read_text().split()
        ref_len_of = None
    if len(set(order)) != len(order):
        raise ValueError("duplicate tx id in universe")
    return order, ref_len_of


def load_tx_biotype(tx2biotype_path):
    """tx_id -> (chrom, transcript_type) from a build_tx2biotype.py TSV (for 6-col pack_meta)."""
    d = {}
    with Path(tx2biotype_path).open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[c["tx_id"]]] = (f[c["chrom"]], f[c["transcript_type"]])
    return d


def md5(path, blocksize: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(blocksize), b""):
            h.update(block)
    return h.hexdigest()


def write_pack(out_dir, order, tgt, cov, *, group, species,
               biotype=None, ref_len_of=None, n_ribo=0, n_rna=0, extra_norm=None):
    """Write the per-pack contract files into out_dir. Returns the coverage_norm dict.

    order      : list of versioned tx ids (row order)
    tgt, cov   : {tx -> int64 per-nt array} pooled P-sites and RNA coverage
    biotype    : optional {tx -> (chrom, transcript_type)}; when given, pack_meta is 6-col
                 (matching the tissue packs, which are the biotype source for eval)
    ref_len_of : optional {tx -> length}; asserts every pooled length matches (annotation guard)
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    miss_t = [t for t in order if t not in tgt]
    miss_c = [t for t in order if t not in cov]
    if miss_t or miss_c:
        raise ValueError(f"universe tx absent from inputs: {len(miss_t)} in target "
                         f"(e.g. {miss_t[:3]}), {len(miss_c)} in coverage (e.g. {miss_c[:3]})")

    lengths = np.empty(len(order), dtype=np.int32)
    for k, t in enumerate(order):
        lt, lc = tgt[t].shape[0], cov[t].shape[0]
        if lt != lc:
            raise ValueError(f"length mismatch {t}: target {lt} vs coverage {lc}")
        if ref_len_of is not None and lt != ref_len_of[t]:
            raise ValueError(f"pooled length {t}={lt} != ref_pack {ref_len_of[t]} "
                             "(annotation mismatch -> orf_track/embeddings would misalign)")
        lengths[k] = lt
    offsets = np.zeros(len(order) + 1, dtype=np.int64)
    offsets[1:] = np.cumsum(lengths.astype(np.int64))
    total = int(offsets[-1])

    tgt_cat = np.empty(total, dtype=np.int32)
    cov_cat = np.empty(total, dtype=np.int32)
    for k, t in enumerate(order):
        s, e = offsets[k], offsets[k + 1]
        tgt_cat[s:e] = tgt[t].astype(np.int32)
        cov_cat[s:e] = cov[t].astype(np.int32)

    np.save(out_dir / "offsets.npy", offsets)
    np.save(out_dir / "lengths.npy", lengths)
    np.save(out_dir / "target_counts.npy", tgt_cat)
    np.save(out_dir / "coverage.npy", cov_cat)
    (out_dir / "tx_order.txt").write_text("\n".join(order) + "\n")

    global_mean = float(cov_cat.sum(dtype=np.int64)) / int(cov_cat.size) if cov_cat.size else 0.0
    norm = {
        "global_mean_coverage": global_mean,
        "group": group,
        "species": species,
        "n_ribo_samples": int(n_ribo),
        "n_rna_samples": int(n_rna),
        "sum_L": total,
        "coverage_total": int(cov_cat.sum(dtype=np.int64)),
    }
    if extra_norm:
        norm.update(extra_norm)
    (out_dir / "coverage_norm.json").write_text(json.dumps(norm, indent=2) + "\n")

    n_scorable = 0
    with (out_dir / "pack_meta.tsv").open("w") as m:
        if biotype is not None:
            m.write("tx_id\tlength\tchrom\ttranscript_type\ttotal_psites\ttotal_coverage\n")
        else:
            m.write("tx_id\tlength\ttotal_psites\ttotal_coverage\n")
        for k, t in enumerate(order):
            s, e = offsets[k], offsets[k + 1]
            tp = int(tgt_cat[s:e].sum())
            tc = int(cov_cat[s:e].sum())
            if tp >= SCORABLE_MIN_PSITES:
                n_scorable += 1
            if biotype is not None:
                chrom, ttype = biotype.get(t, ("NA", "unknown"))
                m.write(f"{t}\t{lengths[k]}\t{chrom}\t{ttype}\t{tp}\t{tc}\n")
            else:
                m.write(f"{t}\t{lengths[k]}\t{tp}\t{tc}\n")

    _log(f"wrote {out_dir}/  global_mean_coverage={global_mean:.3f}")
    _log(f"  transcripts={len(order):,}  sum_L={total:,}  "
         f"target_total={int(tgt_cat.sum()):,}  coverage_total={int(cov_cat.sum()):,}  "
         f"scorable(>={SCORABLE_MIN_PSITES} psites)={n_scorable:,}")
    return norm
