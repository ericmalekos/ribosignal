#!/usr/bin/env python3
"""Leave-one-tissue-out trainer for the per-nt Ribo-seq signal model.

Trains on the held-in Chothani tissues (concatenated) and tests on a held-out tissue the
model never saw -- the cross-tissue generalization objective (methods.md, LOTO_PLAN.md).
Reuses the model, loss, and evaluation from train.py / model.py unchanged; the only new
piece is multi-tissue data assembly: one PackedStore per tissue (each with its own packed
target + coverage + depth normalizer), combined into a ConcatDataset with a token-budget
sampler over the concatenated lengths. The FM embeddings and the ORF track are
tissue-independent and shared (RIBO_ORF_TRACK selects the ORF variant, as in train.py).

Splits (dataset.loto_split): train / val over the held-in tissues (val = a held-out
chromosome fold, gene-disjoint from train for honest early stopping); test = the held-out
tissue's scorable transcripts. A transcript is "scorable" in a tissue if it has
>= --min_train_signal pooled P-sites there (learn/score shape where shape exists), chrM
excluded (project rule).

Example (inside the CUDA SIF):
  RIBO_ORF_TRACK=.../data/packed/orf_track_v2.npy python3 train_loto.py \
    --holdout Hepatocytes --emb_backend rinalmo --cov_norm global_mean \
    --use_orf_track --n_attn_layers 2 --out results/loto/orf_v2_attn_holdout_Hepatocytes
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import ConcatDataset, DataLoader

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dataset import (  # noqa: E402
    BACKENDS,
    PackedStore,
    RiboDataset,
    TokenBudgetSampler,
    collate_pad,
    loto_split,
    tissue_pack_dir,
)
from model import (  # noqa: E402
    RiboSignalModel,
    count_mse,
    profile_entropy_gap,
    profile_multinomial_nll,
)
from train import evaluate, load_biotype  # noqa: E402


def load_cds_map():
    """tx_id -> (cds_start, cds_stop) in tx coords from data/tx2cds.tsv (or $RIBO_TX2CDS): start =
    utr5_len (0-based), stop = utr5_len + cds_len; only tx with a start codon + cds_len>0. Returns
    None if absent -- frame-0 selection is then disabled (training otherwise unchanged)."""
    default = Path(__file__).resolve().parent.parent / "data" / "tx2cds.tsv"
    p = Path(os.environ.get("RIBO_TX2CDS", str(default)))
    if not p.exists():
        return None
    cds = {}
    with open(p) as f:
        ci = {n: i for i, n in enumerate(f.readline().rstrip("\n").split("\t"))}
        for ln in f:
            r = ln.rstrip("\n").split("\t")
            if r[ci["has_start_codon"]] != "1":
                continue
            u5, cl = int(r[ci["utr5_len"]]), int(r[ci["cds_len"]])
            if cl > 0:
                cds[r[ci["tx_id"]]] = (u5, u5 + cl)
    return cds


def make_multi_loader(stores, per_tissue_tx, budget, shuffle, workers, seed=0,
                      input_mode="both", use_orf=False, cov_norm="global_mean"):
    """ConcatDataset over (tissue, tx_ids) pairs + a token-budget sampler on the combined
    lengths. Each tissue's items are drawn from that tissue's own store (its own target,
    coverage, depth mean); embeddings + ORF track are shared."""
    datasets, lens = [], []
    for tissue, tx_ids in per_tissue_tx:
        if not tx_ids:
            continue
        st = stores[tissue]
        datasets.append(RiboDataset(st, tx_ids, input_mode=input_mode, use_orf=use_orf,
                                    cov_norm=cov_norm))
        lens.extend(st.length_of(t) for t in tx_ids)
    if not datasets:
        raise SystemExit("no transcripts in this split (check --min_train_signal / packs)")
    ds = ConcatDataset(datasets)
    samp = TokenBudgetSampler(lens, budget=budget, shuffle=shuffle, seed=seed)
    loader = DataLoader(ds, batch_sampler=samp, collate_fn=collate_pad,
                        num_workers=workers, pin_memory=torch.cuda.is_available(),
                        persistent_workers=workers > 0)
    return loader, samp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--holdout", required=True, help="tissue to hold out for testing")
    ap.add_argument("--train_tissues", default="",
                    help="comma list of held-in tissues (default: all 9 minus --holdout)")
    ap.add_argument("--min_train_signal", type=int, default=50,
                    help="per-tissue: keep tx with >= this many pooled P-sites (train + test)")
    ap.add_argument("--val_fold", type=int, default=0,
                    help="chromosome fold used as the (gene-disjoint) validation set")
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--weight_decay", type=float, default=1e-2)
    ap.add_argument("--channels", type=int, default=256)
    ap.add_argument("--n_blocks", type=int, default=10)
    ap.add_argument("--dropout", type=float, default=0.1)
    ap.add_argument("--budget", type=int, default=16000)
    ap.add_argument("--count_weight", type=float, default=0.1)
    ap.add_argument("--peakiness_weight", type=float, default=0.0,
                    help="O3 anti-smoothing: weight on the entropy-gap penalty "
                         "H(pred)-H(obs) (0 = off, the deployed default; try 0.05-0.2 to "
                         "push the model toward structured/peaked profiles).")
    ap.add_argument("--warmup", type=int, default=2)
    ap.add_argument("--patience", type=int, default=6)
    ap.add_argument("--eval_cap", type=int, default=2500)
    ap.add_argument("--num_workers", type=int, default=8)
    ap.add_argument("--max_tx_per_tissue", type=int, default=0,
                    help="cap (random, seeded) scorable tx per tissue for tissue balance + "
                         "cost; 0 = all scorable. Deep tissues would otherwise dominate.")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--resume", action="store_true",
                    help="resume from <out>/last.pt (model + optimizer + epoch) if present -- "
                         "recovers a run interrupted mid-training (e.g. a filesystem-quota kill)")
    ap.add_argument("--input_mode", default="both", choices=["both", "emb", "cov"])
    ap.add_argument("--cov_norm", default="global_mean", choices=["raw", "global_mean"])
    ap.add_argument("--emb_backend", default="rinalmo", choices=list(BACKENDS))
    ap.add_argument("--use_orf_track", action="store_true")
    ap.add_argument("--n_attn_layers", type=int, default=0,
                    help="number of global-context mixer layers after the CNN body (0 = pure CNN)")
    ap.add_argument("--n_heads", type=int, default=8)
    ap.add_argument("--fm_to_mixer", action="store_true",
                    help="feed the FM embedding DIRECTLY to the mixer via a dedicated projection and "
                         "a zero-init gate, in addition to the normal in_proj->CNN path. Tests whether "
                         "the shared 1x1 in_proj bottleneck (1280->256 for RiNALMo) is what suppresses "
                         "the FM lift. Starts numerically identical to the baseline; the learned gate "
                         "magnitude reports how much the model actually wants the un-convolved emb.")
    ap.add_argument("--mixer", default="transformer", choices=["transformer", "mamba"],
                    help="global-context mixer: transformer (default) or a bidirectional Mamba "
                         "stack (needs mamba_ssm; run inside the Orthrus/HydraRNA CUDA SIF)")
    ap.add_argument("--mamba_d_state", type=int, default=16)
    ap.add_argument("--mamba_d_conv", type=int, default=4)
    ap.add_argument("--mamba_expand", type=int, default=2)
    ap.add_argument("--learn_start_context", action="store_true",
                    help="Q3 Kozak ablation: replace the ORF-track's precomputed start-context weight "
                         "with a learnable conv over the one-hot window [-6,+4] (onehot backend only; "
                         "pair with the no-Kozak ORF track). The 4x10 kernel is a learned Kozak matrix.")
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = args.device if (args.device != "cuda" or torch.cuda.is_available()) else "cpu"
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    train_tissues = ([t for t in args.train_tissues.split(",") if t]
                     if args.train_tissues else None)
    train_split, val_split, holdout, test_tx = loto_split(
        args.holdout, train_tissues=train_tissues,
        min_signal=args.min_train_signal, val_fold=args.val_fold)
    resolved_train = [t for t, _ in train_split]
    if args.max_tx_per_tissue:
        # balance tissues + cap cost: seeded-random subsample each tissue to at most N tx
        rng = np.random.default_rng(args.seed)

        def cap(ids, n):
            if not ids or len(ids) <= n:
                return sorted(ids)
            return sorted(rng.choice(ids, size=n, replace=False).tolist())

        train_split = [(t, cap(ids, args.max_tx_per_tissue)) for t, ids in train_split]
        vcap = max(300, args.max_tx_per_tissue // 5)
        val_split = [(t, cap(ids, vcap)) for t, ids in val_split]

    # one store per tissue (held-in + holdout); embeddings + ORF track shared
    tissues = resolved_train + [holdout]
    stores = {t: PackedStore(pack=tissue_pack_dir(t), tx_index=BACKENDS[args.emb_backend])
              for t in tissues}
    biotype = load_biotype()

    n_train = sum(len(ids) for _, ids in train_split)
    n_val = sum(len(ids) for _, ids in val_split)
    print(f"device={device}  holdout={holdout}  train_tissues={resolved_train}", file=sys.stderr)
    print(f"train_tx={n_train:,} (over {len(resolved_train)} tissues)  val_tx={n_val:,}  "
          f"test_tx={len(test_tx):,}  min_signal={args.min_train_signal}", file=sys.stderr)
    (out / "args.json").write_text(json.dumps(
        {**vars(args), "resolved_train_tissues": resolved_train,
         "n_train_tx": n_train, "n_val_tx": n_val, "n_test_tx": len(test_tx)}, indent=2))

    train_loader, train_samp = make_multi_loader(
        stores, train_split, args.budget, True, args.num_workers, seed=args.seed,
        input_mode=args.input_mode, use_orf=args.use_orf_track, cov_norm=args.cov_norm)
    val_loader, _ = make_multi_loader(
        stores, val_split, args.budget, False, args.num_workers,
        input_mode=args.input_mode, use_orf=args.use_orf_track, cov_norm=args.cov_norm)
    test_loader, _ = make_multi_loader(
        stores, [(holdout, test_tx)], args.budget, False, args.num_workers,
        input_mode=args.input_mode, use_orf=args.use_orf_track, cov_norm=args.cov_norm)

    d_emb = stores[holdout].load_emb(test_tx[0]).shape[1]
    extra_in = stores[holdout].orf_channels if args.use_orf_track else 0
    if args.use_orf_track and stores[holdout].orf_track is None:
        print("WARNING: --use_orf_track set but no ORF track loaded (RIBO_ORF_TRACK?)",
              file=sys.stderr)
    model = RiboSignalModel(d_emb=d_emb, channels=args.channels, n_blocks=args.n_blocks,
                            dropout=args.dropout, extra_in=extra_in,
                            n_attn_layers=args.n_attn_layers, n_heads=args.n_heads,
                            mixer=args.mixer, mamba_d_state=args.mamba_d_state,
                            mamba_d_conv=args.mamba_d_conv, mamba_expand=args.mamba_expand,
                            learn_start_context=args.learn_start_context,
                            fm_to_mixer=args.fm_to_mixer).to(device)
    n_par = sum(p.numel() for p in model.parameters())
    print(f"model params: {n_par:,}  d_emb={d_emb}  extra_in={extra_in}  "
          f"attn={args.n_attn_layers}  mixer={args.mixer}  backend={args.emb_backend}"
          f"  fm_to_mixer={args.fm_to_mixer}", file=sys.stderr)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    steps_per_epoch = max(1, len(train_samp))

    def lr_at(epoch_frac):
        if epoch_frac < args.warmup:
            return epoch_frac / max(1e-9, args.warmup)
        prog = (epoch_frac - args.warmup) / max(1e-9, args.epochs - args.warmup)
        return 0.5 * (1 + math.cos(math.pi * min(1.0, prog)))

    # CDS map for the ORF-aware frame-0 signal (tx_id -> (cds_start, cds_stop), tx coords).
    # Additive: val_pearson still drives best.pt + early stop; frame-0 -> best_frame0.pt only.
    cds_map = load_cds_map()
    print(f"frame-0 selection: {len(cds_map) if cds_map else 0} CDS tx loaded", file=sys.stderr)

    history = []
    best = -1e9
    best_f0 = -1e9
    best_state = None
    bad = 0
    start_epoch = 0
    ckpt_path = out / "last.pt"
    if args.resume and ckpt_path.exists():
        ck = torch.load(ckpt_path, map_location=device)
        model.load_state_dict(ck["model"])
        opt.load_state_dict(ck["opt"])
        start_epoch = ck["epoch"] + 1
        best = ck["best"]
        best_f0 = ck.get("best_f0", -1e9)
        bad = ck["bad"]
        if (out / "history.json").exists():
            history = json.loads((out / "history.json").read_text())
        if (out / "best.pt").exists():
            best_state = torch.load(out / "best.pt", map_location=device)
        print(f"RESUMED from epoch {start_epoch} (best val_pearson={best:.4f}, bad={bad})",
              file=sys.stderr)
    for epoch in range(start_epoch, args.epochs):
        model.train()
        train_samp.set_epoch(epoch)
        t0 = time.time()
        run_loss = run_pnll = run_cnt = run_pk = 0.0
        nb = 0
        for bi, batch in enumerate(train_loader):
            frac = epoch + bi / steps_per_epoch
            for g in opt.param_groups:
                g["lr"] = args.lr * lr_at(frac)
            feats = batch["feats"].to(device, non_blocking=True)
            mask = batch["mask"].to(device, non_blocking=True)
            counts = batch["counts"].to(device, non_blocking=True)
            logits, pred_lc = model(feats, mask)
            pnll = profile_multinomial_nll(logits, counts, mask)
            cnt = count_mse(pred_lc, counts, mask)
            loss = pnll + args.count_weight * cnt
            if args.peakiness_weight > 0:
                pk = profile_entropy_gap(logits, counts, mask)
                loss = loss + args.peakiness_weight * pk
                run_pk += pk.item()
            opt.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            run_loss += loss.item()
            run_pnll += pnll.item()
            run_cnt += cnt.item()
            nb += 1
        val = evaluate(model, val_loader, device, biotype, cap=args.eval_cap, cds=cds_map)
        vp = (val.get("all", {}).get("pearson_median", float("nan"))
              if val.get("n") else float("nan"))
        f0 = (val.get("all", {}).get("frame0_pred_median") if val.get("n") else None)
        dt = time.time() - t0
        rec = {"epoch": epoch, "lr": args.lr * lr_at(epoch + 1),
               "train_loss": run_loss / nb, "train_pnll": run_pnll / nb,
               "train_cnt": run_cnt / nb, "train_pk": run_pk / nb,
               "val": val, "sec": round(dt, 1)}
        history.append(rec)
        pk_str = f" pk={run_pk/nb:.4f}" if args.peakiness_weight > 0 else ""
        f0_str = f" f0={f0:.4f}" if f0 is not None else ""
        print(f"[e{epoch:02d}] loss={run_loss/nb:.4f} pnll={run_pnll/nb:.4f}{pk_str} "
              f"val_pearson={vp:.4f}{f0_str} ({dt:.0f}s)", file=sys.stderr)
        (out / "history.json").write_text(json.dumps(history, indent=2))
        if vp > best:
            best = vp
            bad = 0
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            torch.save(best_state, out / "best.pt")
        else:
            bad += 1
        # ORF-aware checkpoint: best by CDS-anchored frame-0 (additive -- does NOT drive early
        # stop, which stays on val_pearson; lets us compare the two selections post-hoc).
        if f0 is not None and f0 > best_f0:
            best_f0 = f0
            torch.save({k: v.detach().cpu().clone() for k, v in model.state_dict().items()},
                       out / "best_frame0.pt")
        # per-epoch full checkpoint so a mid-run kill (e.g. a filesystem-quota write failure,
        # which is how the first pilot died) resumes cleanly with --resume
        torch.save({"model": model.state_dict(), "opt": opt.state_dict(),
                    "epoch": epoch, "best": best, "best_f0": best_f0, "bad": bad}, ckpt_path)
        if bad >= args.patience:
            print(f"early stop at epoch {epoch} (best val_pearson={best:.4f})", file=sys.stderr)
            break

    if best_state is not None:
        model.load_state_dict(best_state)
    test = evaluate(model, test_loader, device, biotype, cap=None, cds=cds_map)
    (out / "test_metrics.json").write_text(json.dumps(test, indent=2))
    print(f"TEST ({holdout}):", json.dumps(test.get("all", {}), indent=2), file=sys.stderr)
    print(f"wrote {out}/test_metrics.json  best_val_pearson={best:.4f}", file=sys.stderr)


if __name__ == "__main__":
    main()
