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
import os
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
    load_split,
    loto_split,
    pack_test_tx,
    tissue_pack_dir,
)
from model import RiboSignalModel  # noqa: E402

import paths  # noqa: E402


def resolve_run(run, arch=None, checkpoint=None, config=None):
    """Find the checkpoint and its training config inside a --run directory.

    The released weights and this loader disagreed on names: HuggingFace publishes
    `attn_best.pt` + `attn_config.json` and `mamba4_best.pt` + `mamba4_config.json`, while
    this script demanded `best.pt` + `args.json`, so the README quickstart could not run
    for anyone. Both layouts are accepted now, plus explicit --checkpoint/--config for a
    directory that follows neither:

      weights/attn_best.pt   + weights/attn_config.json     (HuggingFace / the container)
      weights/mamba4_best.pt + weights/mamba4_config.json
      <run>/best.pt          + <run>/args.json              (a training run directory)

    A directory holding both released checkpoints is ambiguous, so --arch names one.
    Returns (checkpoint_path, config_path).
    """
    run = Path(run)
    ckpt = Path(checkpoint) if checkpoint else None
    cfg = Path(config) if config else None

    if ckpt is None:
        if not run.is_dir():
            sys.exit(f"--run {run} is not a directory (and no --checkpoint given)")
        if arch:
            ckpt = run / f"{arch}_best.pt"
            if not ckpt.is_file():
                sys.exit(f"--arch {arch} but {ckpt} does not exist; "
                         f"found {sorted(q.name for q in run.glob('*best.pt')) or 'no *best.pt'}")
        else:
            found = sorted(run.glob("*_best.pt")) + ([run / "best.pt"]
                                                     if (run / "best.pt").is_file() else [])
            if not found:
                sys.exit(f"no checkpoint in {run}: expected <arch>_best.pt or best.pt. "
                         f"Fetch the released weights with:\n"
                         f"  pip install huggingface_hub && python -c "
                         f"'from huggingface_hub import snapshot_download; "
                         f"snapshot_download(\"emalek/RiboSignal\", local_dir=\"weights\")'")
            if len(found) > 1:
                sys.exit(f"{run} holds {len(found)} checkpoints "
                         f"({', '.join(q.name for q in found)}); pick one with --arch "
                         f"(e.g. --arch attn) or --checkpoint")
            ckpt = found[0]
    if not ckpt.is_file():
        sys.exit(f"checkpoint {ckpt} does not exist")

    if cfg is None:
        stem = ckpt.name[:-len("_best.pt")] if ckpt.name.endswith("_best.pt") else ""
        # the released name first, then the training-run name, then the same prefix
        cands = ([ckpt.parent / f"{stem}_config.json", ckpt.parent / f"{stem}_args.json"]
                 if stem else []) + [ckpt.parent / "args.json", ckpt.parent / "config.json"]
        cfg = next((c for c in cands if c.is_file()), None)
        if cfg is None:
            sys.exit(f"no config for {ckpt.name}: looked for "
                     f"{', '.join(str(c) for c in cands)}. Pass --config.")
    if not cfg.is_file():
        sys.exit(f"config {cfg} does not exist")
    return ckpt, cfg


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


def _enable_mamba_cpu():
    """Run the Mamba mixer on CPU by swapping mamba_ssm's CUDA kernels for the
    pure-PyTorch reference implementations it already ships.

    The fast path is CUDA-only (causal_conv1d_fwd asserts x.is_cuda), but the package
    also ships selective_scan_ref / mamba_inner_ref / causal_conv1d_ref, which are plain
    PyTorch and run anywhere. The catch: each was resolved by name at import time in
    SEVERAL modules, so patching one namespace is not enough -- mamba_inner_ref calls the
    causal_conv1d_fn bound inside selective_scan_interface, which is why patching only
    causal_conv1d_interface still fails.

    Cost: ~7.9 s per 3,000 nt transcript single-threaded, so shard wide (--nshards).
    Verified finite on the real RiboSignalModel (mamba4, d_state 16) at L=1000 and 3000.

    Enabled by RIBO_MAMBA_CPU=1. Off by default: on a GPU the fast path is far quicker.

    This applies where mamba_ssm is INSTALLED but the GPU is not there to run -- the
    container, a CPU node. Where mamba_ssm cannot be installed at all, which is anywhere
    without a CUDA toolchain, there is nothing here to patch and model.py falls back to
    scripts/mamba_ref.py by itself; this becomes a no-op rather than an ImportError.
    """
    try:
        import causal_conv1d.causal_conv1d_interface as cci
        import mamba_ssm.modules.mamba_simple as ms
        import mamba_ssm.ops.selective_scan_interface as ssi
    except ImportError as e:
        print(f"RIBO_MAMBA_CPU=1 but mamba_ssm is not installed ({e}); "
              f"scripts/mamba_ref.py already runs the mixer on CPU, so this is a no-op",
              file=sys.stderr)
        return
    ssi.selective_scan_fn = ssi.selective_scan_ref
    ssi.mamba_inner_fn = ssi.mamba_inner_ref
    ssi.causal_conv1d_fn = cci.causal_conv1d_ref
    cci.causal_conv1d_fn = cci.causal_conv1d_ref
    ms.selective_scan_fn = ssi.selective_scan_ref
    ms.mamba_inner_fn = ssi.mamba_inner_ref
    ms.causal_conv1d_fn = cci.causal_conv1d_ref
    print("RIBO_MAMBA_CPU=1: mamba_ssm CUDA kernels -> reference path", file=sys.stderr)


def main():
    if os.environ.get("RIBO_MAMBA_CPU", "") == "1":
        _enable_mamba_cpu()
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True,
                    help="directory holding the checkpoint and its config: a HuggingFace "
                         "snapshot (attn_best.pt / attn_config.json) or a training run "
                         "(best.pt / args.json). Used as given -- relative paths resolve "
                         "against the working directory, not any other tree.")
    ap.add_argument("--arch", choices=["attn", "mamba4"], default=None,
                    help="which released checkpoint to load when --run holds both")
    ap.add_argument("--checkpoint", default=None,
                    help="explicit .pt path, overriding discovery inside --run")
    ap.add_argument("--config", default=None,
                    help="explicit training-config json, overriding discovery inside --run")
    ap.add_argument("--out", default=None, help="default <run>/dropin (or <run>/dropin_<heldout>)")
    ap.add_argument("--pack", default=None,
                    help="score this pack directory directly, whatever it is called. The "
                         "alternative, --heldout <name>, only reaches "
                         "$RIBO_DATA_DIR/packed_heldout_<name>; a pack you just built with "
                         "build_pack.py need not be named that way.")
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
                    help="FASTA the predicted transcripts are verified against (default "
                         "$RIBO_ONEHOT_FASTA, the same universe the one-hot backend reads)")
    ap.add_argument("--nshards", type=int, default=1,
                    help="split the test tx into N interleaved shards (for short-partition CPU "
                         "runs that cannot dump the whole set in one wall-clock window)")
    ap.add_argument("--shard", type=int, default=0, help="which shard [0, nshards)")
    ap.add_argument("--budget", type=int, default=24000)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    run = Path(args.run)
    ckpt_path, cfg_path = resolve_run(run, args.arch, args.checkpoint, args.config)
    default_out = (f"dropin_{args.heldout}" if args.heldout else "dropin")
    out = Path(args.out) if args.out else (run / default_out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads(cfg_path.read_text())
    device = args.device if (args.device != "cuda" or torch.cuda.is_available()) else "cpu"
    backend = cfg.get("emb_backend", "rinalmo")

    if args.pack or args.heldout:
        pack = Path(args.pack) if args.pack else heldout_pack_dir(args.heldout)
        if not (pack / "tx_order.txt").is_file():
            sys.exit(f"{pack} is not a pack: no tx_order.txt. Build one with\n"
                     f"  python scripts/build_pack.py --fasta <transcripts.fa> "
                     f"--coverage <coverage.hd5> --out {pack}")
        tissue = args.heldout or pack.name
        store = PackedStore(pack=pack, tx_index=BACKENDS[backend])
        if args.tx_list:
            wl = [t for t in Path(args.tx_list).read_text().split() if t]
            test_ids = store.usable(wl)
            print(f"tx_list: {len(wl)} requested -> {len(test_ids)} usable in pack", file=sys.stderr)
        else:
            test_ids = store.usable(pack_test_tx(pack, args.min_signal))
            if not test_ids:
                sys.exit(f"no transcript in {pack} has >= {args.min_signal} pooled P-sites. "
                         f"An inference-only pack has none by construction -- pass "
                         f"--tx_list {pack / 'expressed_tx.txt'} to score the expressed ones.")
    elif cfg.get("holdout"):
        tissue = cfg["holdout"]
        store = PackedStore(pack=tissue_pack_dir(tissue), tx_index=BACKENDS[backend])
        _, _, _, test_ids = loto_split(
            tissue, train_tissues=cfg.get("resolved_train_tissues"),
            min_signal=cfg.get("min_train_signal", 50), val_fold=cfg.get("val_fold", 0))
        if args.tx_list:
            # --tx_list used to be honoured ONLY on the --heldout path, so passing it to a LOTO run
            # was silently ignored and the full test set got dumped instead. Silent, because the
            # only symptom is a much larger dump than asked for. Restrict here too.
            wl = {t for t in Path(args.tx_list).read_text().split() if t}
            keep = [t for t in test_ids if t in wl]
            print(f"tx_list: {len(wl)} requested -> {len(keep)} of {len(test_ids)} test tx kept",
                  file=sys.stderr)
            if not keep:
                sys.exit(f"--tx_list {args.tx_list} matched none of the {len(test_ids)} test tx")
            test_ids = keep
    else:
        tissue = "Fibroblast"
        store = PackedStore(tx_index=BACKENDS[backend])
        _, _, test_ids = load_split(store, test_fold=cfg.get("test_fold", 0))

    # Interleaved sharding for short-partition CPU runs: shard s takes test_ids[s::nshards], so
    # long transcripts (the O(L^2)-attention cost) spread evenly across shards. Merge later.
    if args.nshards > 1:
        test_ids = sorted(test_ids)[args.shard::args.nshards]

    fasta = Path(args.fasta) if args.fasta else paths.onehot_fasta()
    if not fasta.is_file():
        sys.exit(paths.missing(fasta, "sequence FASTA", "RIBO_ONEHOT_FASTA", "fasta"))
    seqs = read_fasta(fasta)
    print(f"run={run.name} ckpt={ckpt_path.name} cfg={cfg_path.name} tissue={tissue} "
          f"backend={backend} test_tx={len(test_ids)} "
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
    # weights_only=True is not optional here: the README tells users to download a
    # checkpoint off the internet and load it, and a pickled state_dict is arbitrary code
    # execution otherwise. Passed explicitly rather than relying on the torch default,
    # which has changed across releases.
    model.load_state_dict(torch.load(ckpt_path, map_location=device, weights_only=True))
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
                        meta=np.array([ckpt_path.stem, tissue, backend, str(n_skip)],
                                      dtype="<U64"))
    print(f"wrote {dst}: {len(tx_ids_s)} tx ({n_skip} skipped seq/len mismatch)", file=sys.stderr)


if __name__ == "__main__":
    main()
