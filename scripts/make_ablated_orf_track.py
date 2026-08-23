#!/usr/bin/env python3
"""Write a copy of an ORF track with one or more channels zeroed.

`orf_channel_ablation.sbatch` does this inline for the liver-3x3 evaluation universe. That inline
version cannot be reused for a TRAINING run, because the training packs live on a different universe
(212,647,172 rows for the human union vs 56,349,794 for liver-3x3) and a track is only valid for the
universe it was built against. Hence a tool rather than a third copy of the same four lines.

Channels are (0,1,2) frame-0/1/2 occupancy, (3) start_ext propensity, (4) is_stop.

Copies in chunks through a memmap rather than loading the whole array: the union track is ~2.1 GB and
this is often run on the head node. The output is written to a temp path and renamed atomically, so a
concurrent reader can never observe a torn file.

A sibling `<out>_meta.json` is written by copying the source meta and recording `zeroed_channels`, so
a track cannot end up in a training run with provenance that silently claims it is the full one. Any
`kozak` field is carried through unchanged -- `feedback_kozak_never_set` requires kozak=none, and this
tool must not launder that.

cas12a env. Example:
  python3 make_ablated_orf_track.py --in data/packed_union/orf_track_v2.npy \\
      --zero 0,1,2 --out data/packed_union/orf_track_v2_nof012.npy
"""
import argparse
import json
import os
import pathlib

import numpy as np

CHANNEL_NAMES = ["orf_f0", "orf_f1", "orf_f2", "start_ext", "is_stop"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True)
    ap.add_argument("--zero", required=True, help="comma-separated channel indices, e.g. 0,1,2")
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunk", type=int, default=8_000_000, help="rows per copy chunk")
    a = ap.parse_args()

    src = pathlib.Path(a.src)
    out = pathlib.Path(a.out)
    cols = sorted({int(c) for c in a.zero.split(",") if c.strip() != ""})
    if not cols:
        raise SystemExit("--zero selected no channels")

    arr = np.load(src, mmap_mode="r")
    if arr.ndim != 2:
        raise SystemExit(f"expected a 2-D track, got shape {arr.shape}")
    n, c = arr.shape
    bad = [i for i in cols if i < 0 or i >= c]
    if bad:
        raise SystemExit(f"channel(s) {bad} out of range for a {c}-channel track")
    if len(cols) == c:
        raise SystemExit("refusing to zero EVERY channel -- that is an empty track, not an ablation")

    tmp = out.with_suffix(f".tmp{os.getpid()}.npy")
    dst = np.lib.format.open_memmap(tmp, mode="w+", dtype=arr.dtype, shape=arr.shape)
    for i in range(0, n, a.chunk):
        j = min(i + a.chunk, n)
        blk = np.array(arr[i:j])
        blk[:, cols] = 0
        dst[i:j] = blk
    dst.flush()
    del dst
    os.replace(tmp, out)

    meta_src = src.with_name(src.stem + "_meta.json")
    meta = json.loads(meta_src.read_text()) if meta_src.exists() else {}
    zn = [CHANNEL_NAMES[i] if i < len(CHANNEL_NAMES) else str(i) for i in cols]
    # The copied meta describes the SOURCE track, so its `mean_per_channel` still shows ~0.38 for
    # channels this file zeroes -- a record that contradicts `zeroed_channels` in the same file and
    # reads like the ablation silently failed. The zeroed entries are set to 0.0 here (exact, not
    # sampled: the assert below proves it) and the rest are labelled as inherited.
    mpc = dict(meta.get("mean_per_channel") or {})
    if mpc:
        for nm in zn:
            if nm in mpc:
                mpc[nm] = 0.0
        meta["mean_per_channel"] = mpc
        meta["mean_per_channel_note"] = ("zeroed channels set to 0.0; non-zeroed values are "
                                         "INHERITED from derived_from and were not recomputed")
    meta.update({"derived_from": str(src), "zeroed_channels": cols,
                 "zeroed_channel_names": zn,
                 "shape": [int(n), int(c)], "dtype": str(arr.dtype)})
    out.with_name(out.stem + "_meta.json").write_text(json.dumps(meta, indent=2))

    chk = np.array(np.load(out, mmap_mode="r")[: min(n, 100_000)])
    assert chk[:, cols].sum() == 0, "zeroed channels are not zero in the output"
    keep = [i for i in range(c) if i not in cols]
    print(f"  wrote {out}  shape={arr.shape} dtype={arr.dtype}")
    print(f"  zeroed channels {cols} ({', '.join(meta['zeroed_channel_names'])}); kept {keep}")
    print(f"  kozak in meta: {meta.get('kozak', 'n/a')}")


if __name__ == "__main__":
    main()
