#!/usr/bin/env python3
"""Merge dump_pred_profiles.py shards (pred_profiles_shard{s}of{N}.npz) into one pred_profiles.npz,
so the sharded short-partition CPU dump is transparent to ribocode_dropin.py / assemble_heldout.py.

Usage: merge_pred_shards.py <dropin_dir> <nshards>
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main():
    if len(sys.argv) != 3:
        print("usage: merge_pred_shards.py <dropin_dir> <nshards>", file=sys.stderr)
        return 2
    d = Path(sys.argv[1])
    n = int(sys.argv[2])
    tx_all, pred_all, obs_all, ptot_all, len_all = [], [], [], [], []
    for s in range(n):
        p = d / f"pred_profiles_shard{s}of{n}.npz"
        if not p.exists():
            raise SystemExit(f"missing shard {p}")
        z = np.load(p, allow_pickle=True)
        tx = z["tx_ids"]
        lens = z["lengths"].astype(np.int64)
        off = np.zeros(len(lens) + 1, dtype=np.int64)
        off[1:] = np.cumsum(lens)
        pred, obs, ptot = z["pred_flat"], z["obs_flat"], z["pred_total"]
        for k, t in enumerate(tx):
            a, b = off[k], off[k + 1]
            tx_all.append(str(t))
            pred_all.append(pred[a:b])
            obs_all.append(obs[a:b])
            ptot_all.append(float(ptot[k]))
            len_all.append(int(lens[k]))
        print(f"  shard {s}: {len(tx):,} tx", file=sys.stderr)

    order = np.argsort(tx_all)
    tx_s = np.array([tx_all[i] for i in order], dtype="<U25")
    lengths = np.array([len_all[i] for i in order], dtype=np.int64)
    pred_flat = (np.concatenate([pred_all[i] for i in order]).astype(np.float32)
                 if order.size else np.zeros(0, np.float32))
    obs_flat = (np.concatenate([obs_all[i] for i in order]).astype(np.int32)
                if order.size else np.zeros(0, np.int32))
    ptot = np.array([ptot_all[i] for i in order], dtype=np.float64)
    dst = d / "pred_profiles.npz"
    np.savez_compressed(dst, tx_ids=tx_s, lengths=lengths, pred_flat=pred_flat,
                        obs_flat=obs_flat, pred_total=ptot,
                        meta=np.array(["merged", str(n)], dtype="<U64"))
    print(f"wrote {dst}: {len(tx_s):,} tx from {n} shards", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
