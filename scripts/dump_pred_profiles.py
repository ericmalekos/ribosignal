#!/usr/bin/env python3
"""Dump per-nt PREDICTED P-site profiles for a run's held-out transcripts, for the RiboCode drop-in
test (Task 14): can the model's predicted profile stand in for experimental Ribo-seq inside the real
ORF caller?

For every held-out test transcript the model scores, this writes three per-nt arrays keyed by
versioned ENST (0-based, length == transcript length, so they align 1:1 to the RiboCode annotation
DB transcript lengths and to the FM embeddings):
  pred_profile  softmax(profile_logits) over the transcript (sums to ~1)   -- the predicted SHAPE
  obs_counts    the real pooled P-site target counts (integers)            -- the reference density
  pred_total    scalar per tx = expm1(pred_logcount) from the count head   -- the predicted DEPTH

ribocode_dropin.py then builds RiboCode's tpsites_sum dict from these (real / pred x obs-depth /
pred-depth) and calls detectORF.main directly. Mirrors eval_localization.py's model + data load
exactly (same PackedStore / RiboDataset / loto_split / TokenBudgetSampler / RiboSignalModel), so the
profiles are identical to the ones the localization eval scored. GPU (rinalmo SIF). Writes
<out>/pred_profiles.npz (numpy-only; the SIF lacks h5py).

Usage: dump_pred_profiles.py --run <run_dir> [--out <dir>] [--device cuda]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dataset import (  # noqa: E402
    BACKENDS,
    PackedStore,
    RiboDataset,
    TokenBudgetSampler,
    collate_pad,
    heldout_pack_dir,
    heldout_test_tx,
    load_split,
    loto_split,
    tissue_pack_dir,
)
from model import RiboSignalModel  # noqa: E402

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
FASTA = NEW / "data" / "fibroblast_universe.fa"


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = ln[1:].split()[0]
                chunks = []
            else:
                chunks.append(ln.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--out", default=None, help="default <run>/dropin (or <run>/dropin_<heldout>)")
    ap.add_argument("--heldout", default=None,
                    help="external held-out dataset name (packed_heldout_<name>); applies the "
                         "trained model to that dataset's every scorable tx (no train/val split)")
    ap.add_argument("--min_signal", type=int, default=50,
                    help="held-out scorable floor: min pooled P-sites for a test tx")
    ap.add_argument("--tx_list", default=None,
                    help="explicit tx list to score (for predict-only coverage packs with no Ribo-seq "
                         "target, e.g. proteogenomics A549: score the RNA-seq-expressed tx, not the "
                         "Ribo-seq-scorable ones). Only used with --heldout.")
    ap.add_argument("--fasta", default=None,
                    help="override the seq-verification FASTA (cross-species held-out universe)")
    ap.add_argument("--nshards", type=int, default=1,
                    help="split the test tx into N interleaved shards (for short-partition CPU "
                         "runs that cannot dump the whole set in one wall-clock window)")
    ap.add_argument("--shard", type=int, default=0, help="which shard [0, nshards)")
    ap.add_argument("--budget", type=int, default=24000)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    run = Path(args.run)
    if not run.is_absolute():
        run = NEW / run
    default_out = (f"dropin_{args.heldout}" if args.heldout else "dropin")
    out = Path(args.out) if args.out else (run / default_out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((run / "args.json").read_text())
    device = args.device if (args.device != "cuda" or torch.cuda.is_available()) else "cpu"
    backend = cfg.get("emb_backend", "rinalmo")

    if args.heldout:
        tissue = args.heldout
        store = PackedStore(pack=heldout_pack_dir(args.heldout), tx_index=BACKENDS[backend])
        if args.tx_list:
            wl = [t for t in Path(args.tx_list).read_text().split() if t]
            test_ids = store.usable(wl)
            print(f"tx_list: {len(wl)} requested -> {len(test_ids)} usable in pack", file=sys.stderr)
        else:
            test_ids = store.usable(heldout_test_tx(args.heldout, args.min_signal))
    elif cfg.get("holdout"):
        tissue = cfg["holdout"]
        store = PackedStore(pack=tissue_pack_dir(tissue), tx_index=BACKENDS[backend])
        _, _, _, test_ids = loto_split(
            tissue, train_tissues=cfg.get("resolved_train_tissues"),
            min_signal=cfg.get("min_train_signal", 50), val_fold=cfg.get("val_fold", 0))
    else:
        tissue = "Fibroblast"
        store = PackedStore(tx_index=BACKENDS[backend])
        _, _, test_ids = load_split(store, test_fold=cfg.get("test_fold", 0))

    # Interleaved sharding for short-partition CPU runs: shard s takes test_ids[s::nshards], so
    # long transcripts (the O(L^2)-attention cost) spread evenly across shards. Merge later.
    if args.nshards > 1:
        test_ids = sorted(test_ids)[args.shard::args.nshards]

    seqs = read_fasta(Path(args.fasta) if args.fasta else FASTA)
    print(f"run={run.name} tissue={tissue} backend={backend} test_tx={len(test_ids)} "
          f"shard={args.shard}/{args.nshards} device={device}", file=sys.stderr)

    use_orf = cfg.get("use_orf_track", False)
    ds = RiboDataset(store, test_ids, input_mode=cfg.get("input_mode", "both"),
                     use_orf=use_orf, cov_norm=cfg.get("cov_norm", "raw"))
    lens = [store.length_of(t) for t in test_ids]
    eval_budget = min(args.budget, cfg.get("budget", args.budget))
    samp = TokenBudgetSampler(lens, budget=eval_budget, shuffle=False)
    d_emb = store.load_emb(test_ids[0]).shape[1]
    extra_in = store.orf_channels if use_orf else 0
    model = RiboSignalModel(d_emb=d_emb, channels=cfg["channels"], n_blocks=cfg["n_blocks"],
                            dropout=0.0, extra_in=extra_in,
                            n_attn_layers=cfg.get("n_attn_layers", 0),
                            n_heads=cfg.get("n_heads", 8),
                            mixer=cfg.get("mixer", "transformer"),
                            mamba_d_state=cfg.get("mamba_d_state", 16),
                            mamba_d_conv=cfg.get("mamba_d_conv", 4),
                            mamba_expand=cfg.get("mamba_expand", 2),
                            learn_start_context=cfg.get("learn_start_context", False),
                            fm_to_mixer=cfg.get("fm_to_mixer", False)).to(device)
    model.load_state_dict(torch.load(run / "best.pt", map_location=device))
    model.eval()

    tx_ids, pred_profile, obs_counts, pred_total = [], [], [], []
    n_skip = 0
    with torch.no_grad():
        for batch in torch.utils.data.DataLoader(ds, batch_sampler=samp,
                                                 collate_fn=collate_pad, num_workers=4):
            feats = batch["feats"].to(device)
            mask = batch["mask"].to(device)
            logits, pred_lc = model(feats, mask)
            logp = torch.log_softmax(logits.masked_fill(~mask, -1e9), dim=1)
            P = torch.exp(logp).cpu().numpy()
            ptot = torch.expm1(pred_lc).clamp(min=0).cpu().numpy()  # predicted total P-sites per tx
            counts = batch["counts"].numpy()
            for i, tx in enumerate(batch["tx"]):
                L = int(batch["lengths"][i])
                seq = seqs.get(tx)
                if seq is None or len(seq) != L:
                    n_skip += 1
                    continue
                tx_ids.append(tx)
                pred_profile.append(P[i, :L].astype(np.float32))
                obs_counts.append(counts[i, :L].astype(np.int32))
                pred_total.append(float(ptot[i]))

    # numpy-only .npz (the rinalmo SIF lacks h5py): variable-length per-tx arrays are stored
    # flat-concatenated with per-tx lengths (offsets = cumsum), reconstructed in ribocode_dropin.py.
    order = np.argsort(tx_ids)
    tx_ids_s = np.array([tx_ids[i] for i in order], dtype="<U25")
    lengths = np.array([len(pred_profile[i]) for i in order], dtype=np.int64)
    if order.size:
        pred_flat = np.concatenate([pred_profile[i] for i in order]).astype(np.float32)
        obs_flat = np.concatenate([obs_counts[i] for i in order]).astype(np.int32)
    else:
        pred_flat = np.zeros(0, np.float32)
        obs_flat = np.zeros(0, np.int32)
    ptot = np.array([pred_total[i] for i in order], dtype=np.float64)
    dst = (out / f"pred_profiles_shard{args.shard}of{args.nshards}.npz"
           if args.nshards > 1 else out / "pred_profiles.npz")
    np.savez_compressed(dst, tx_ids=tx_ids_s, lengths=lengths, pred_flat=pred_flat,
                        obs_flat=obs_flat, pred_total=ptot,
                        meta=np.array([run.name, tissue, backend, str(n_skip)], dtype="<U64"))
    print(f"wrote {dst}: {len(tx_ids_s)} tx ({n_skip} skipped seq/len mismatch)", file=sys.stderr)


if __name__ == "__main__":
    main()
