#!/usr/bin/env python3
"""Train the dual-head dilated-CNN for per-nt Ribo-seq P-site profile prediction.

Splits are gene-disjoint by chromosome (fold 0 test, fold 1 val, folds 2-4 train).
Loss = multinomial profile NLL + count_weight * count MSE (methods.md section 5B).
Evaluation is scale-invariant per-transcript: Pearson/Spearman between the predicted
position distribution and the observed P-site profile, plus 3-nt periodicity recovery,
reported separately for protein_coding and lncRNA. Early stopping on val profile Pearson.

Runs inside the CUDA SIF (torch + numpy only). Example:
  singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 \
    -B /private/groups/carpenterlab/emalekos/RNAZoo_meta:/work $SIF \
    python3 .../train.py --epochs 40 --out results/fold0_rinalmo
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dataset import (  # noqa: E402
    BACKENDS,
    PACK,
    PackedStore,
    RiboDataset,
    TokenBudgetSampler,
    collate_pad,
    load_split,
)
from model import RiboSignalModel, count_mse, profile_multinomial_nll  # noqa: E402

MIN_SIGNAL = 50   # transcripts with fewer total P-sites are skipped in eval metrics


def load_biotype():
    d = {}
    with (PACK / "pack_meta.tsv").open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[c["tx_id"]]] = f[c["transcript_type"]]
    return d


def _rank(x):
    r = np.empty_like(x)
    r[np.argsort(x, kind="stable")] = np.arange(len(x))
    return r


def pearson(a, b):
    a = a - a.mean()
    b = b - b.mean()
    d = math.sqrt(float(a @ a) * float(b @ b))
    return float(a @ b) / d if d > 0 else float("nan")


def spearman(a, b):
    return pearson(_rank(a).astype(np.float64), _rank(b).astype(np.float64))


def periodicity(profile):
    """Frame (lag 3/6/9/12) minus off-frame (1/2/4/5/7/8) autocorrelation of a profile."""
    x = profile - profile.mean()
    den = float(x @ x)
    if den <= 0 or len(x) <= 13:
        return float("nan")
    ac = [float(x[:-lag] @ x[lag:]) / den for lag in range(1, 13)]
    frame = np.mean([ac[2], ac[5], ac[8], ac[11]])
    off = np.mean([ac[0], ac[1], ac[3], ac[4], ac[6], ac[7]])
    return float(frame - off)


def frame0_frac(profile, ts, te):
    """Fraction of a profile's mass in frame 0 within CDS [ts,te) (ts = CDS start -> frame 0):
    the ORF-aware periodicity that tracks ORF-call quality where val_pearson saturates."""
    seg = profile[ts:te]
    s = float(seg.sum())
    if seg.size < 3 or s <= 0:
        return float("nan")
    return float(seg[np.arange(seg.size) % 3 == 0].sum()) / s


@torch.no_grad()
def evaluate(model, loader, device, biotype, cap=None, cds=None):
    model.eval()
    rows = []
    seen = 0
    for batch in loader:
        feats = batch["feats"].to(device)
        mask = batch["mask"].to(device)
        logits, pred_lc = model(feats, mask)
        logp = torch.log_softmax(logits.masked_fill(~mask, -1e9), dim=1)
        p = torch.exp(logp).cpu().numpy()
        counts = batch["counts"].numpy()
        for i, tx in enumerate(batch["tx"]):
            L = batch["lengths"][i]
            c = counts[i, :L]
            N = c.sum()
            if N < MIN_SIGNAL:
                continue
            pi = p[i, :L]
            q = c / N
            f0 = float("nan")   # CDS-anchored predicted frame-0 fraction (ORF-aware signal)
            if cds is not None:
                cc = cds.get(tx)
                if cc is not None and cc[0] < min(cc[1], L):
                    f0 = frame0_frac(pi, cc[0], min(cc[1], L))
            rows.append((biotype.get(tx, "unknown"),
                         pearson(pi, q), spearman(pi, q),
                         periodicity(pi), periodicity(q), f0))
        seen += len(batch["tx"])
        if cap and seen >= cap:
            break
    if not rows:
        return {"n": 0}
    arr = rows
    def agg(sel):
        sub = [r for r in arr if sel(r)]
        if not sub:
            return None
        f0s = [r[5] for r in sub if not math.isnan(r[5])]
        return {
            "n": len(sub),
            "pearson_median": float(np.median([r[1] for r in sub])),
            "spearman_median": float(np.median([r[2] for r in sub])),
            "period_pred_median": float(np.nanmedian([r[3] for r in sub])),
            "period_obs_median": float(np.nanmedian([r[4] for r in sub])),
            "frame0_pred_median": float(np.median(f0s)) if f0s else None,
            "n_frame0": len(f0s),
        }
    return {
        "n": len(arr),
        "all": agg(lambda r: True),
        "protein_coding": agg(lambda r: r[0] == "protein_coding"),
        "lncRNA": agg(lambda r: r[0] == "lncRNA"),
    }


def make_loader(store, tx_ids, budget, shuffle, workers, seed=0, input_mode="both",
                use_orf=False, cov_norm="raw"):
    ds = RiboDataset(store, tx_ids, input_mode=input_mode, use_orf=use_orf, cov_norm=cov_norm)
    lens = [store.length_of(t) for t in tx_ids]
    samp = TokenBudgetSampler(lens, budget=budget, shuffle=shuffle, seed=seed)
    loader = DataLoader(ds, batch_sampler=samp, collate_fn=collate_pad,
                        num_workers=workers, pin_memory=torch.cuda.is_available(),
                        persistent_workers=workers > 0)
    return loader, samp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--weight_decay", type=float, default=1e-2)
    ap.add_argument("--channels", type=int, default=256)
    ap.add_argument("--n_blocks", type=int, default=10)
    ap.add_argument("--dropout", type=float, default=0.1)
    ap.add_argument("--budget", type=int, default=24000)
    ap.add_argument("--count_weight", type=float, default=0.1)
    ap.add_argument("--warmup", type=int, default=2)
    ap.add_argument("--patience", type=int, default=6)
    ap.add_argument("--eval_cap", type=int, default=2000)
    ap.add_argument("--num_workers", type=int, default=6)
    ap.add_argument("--max_train", type=int, default=0, help="cap train tx for a smoke")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--test_fold", type=int, default=0,
                    help="held-out test fold (val = next fold, train = the rest)")
    ap.add_argument("--input_mode", default="both", choices=["both", "emb", "cov"],
                    help="ablation: both inputs, embedding-only, or coverage-only")
    ap.add_argument("--cov_norm", default="global_mean", choices=["raw", "global_mean"],
                    help="RNA-seq coverage norm: raw = log1p(pooled depth) (legacy, carries "
                         "absolute depth); global_mean = divide by the dataset global mean depth "
                         "first -> depth-normalized, transferable across datasets")
    ap.add_argument("--emb_backend", default="rinalmo", choices=list(BACKENDS),
                    help="which FM per-token embedding feeds the input: rinalmo (1280-d), "
                         "orthrus 4-track (512-d), or concat (1792-d)")
    ap.add_argument("--use_orf_track", action="store_true",
                    help="append the sequence-derived ORF-candidate track (data/packed/"
                         "orf_track.npy) to the input; annotation-free, works for uORF/dORF/"
                         "lncRNA ORFs")
    ap.add_argument("--n_attn_layers", type=int, default=0,
                    help="self-attention layers after the conv body (0 = pure CNN) for "
                         "full-transcript context")
    ap.add_argument("--n_heads", type=int, default=8)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = args.device if (args.device != "cuda" or torch.cuda.is_available()) else "cpu"
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "args.json").write_text(json.dumps(vars(args), indent=2))

    store = PackedStore(tx_index=BACKENDS[args.emb_backend])
    biotype = load_biotype()
    train_ids, val_ids, test_ids = load_split(store, test_fold=args.test_fold)
    if args.max_train:
        train_ids = train_ids[:args.max_train]
        val_ids = val_ids[:max(200, args.eval_cap)]
    print(f"device={device}  train={len(train_ids):,} val={len(val_ids):,} "
          f"test={len(test_ids):,}", file=sys.stderr)

    train_loader, train_samp = make_loader(store, train_ids, args.budget, True,
                                           args.num_workers, seed=args.seed,
                                           input_mode=args.input_mode,
                                           use_orf=args.use_orf_track, cov_norm=args.cov_norm)
    val_loader, _ = make_loader(store, val_ids, args.budget, False, args.num_workers,
                                input_mode=args.input_mode, use_orf=args.use_orf_track,
                                cov_norm=args.cov_norm)
    test_loader, _ = make_loader(store, test_ids, args.budget, False, args.num_workers,
                                 input_mode=args.input_mode, use_orf=args.use_orf_track,
                                 cov_norm=args.cov_norm)

    d_emb = store.load_emb(train_ids[0]).shape[1]
    extra_in = store.orf_channels if args.use_orf_track else 0
    if args.use_orf_track and store.orf_track is None:
        print("WARNING: --use_orf_track set but data/packed/orf_track.npy absent; "
              "running without it", file=sys.stderr)
    model = RiboSignalModel(d_emb=d_emb, channels=args.channels, n_blocks=args.n_blocks,
                            dropout=args.dropout, extra_in=extra_in,
                            n_attn_layers=args.n_attn_layers, n_heads=args.n_heads).to(device)
    n_par = sum(p.numel() for p in model.parameters())
    print(f"model params: {n_par:,}  d_emb={d_emb}  extra_in={extra_in}  "
          f"attn={args.n_attn_layers}  backend={args.emb_backend}", file=sys.stderr)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr,
                            weight_decay=args.weight_decay)
    steps_per_epoch = max(1, len(train_samp))

    def lr_at(epoch_frac):
        if epoch_frac < args.warmup:
            return epoch_frac / max(1e-9, args.warmup)
        prog = (epoch_frac - args.warmup) / max(1e-9, args.epochs - args.warmup)
        return 0.5 * (1 + math.cos(math.pi * min(1.0, prog)))

    history = []
    best = -1e9
    best_state = None
    bad = 0
    for epoch in range(args.epochs):
        model.train()
        train_samp.set_epoch(epoch)
        t0 = time.time()
        run_loss = run_pnll = run_cnt = 0.0
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
            opt.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            run_loss += loss.item()
            run_pnll += pnll.item()
            run_cnt += cnt.item()
            nb += 1
        val = evaluate(model, val_loader, device, biotype, cap=args.eval_cap)
        vp = (val.get("all", {}).get("pearson_median", float("nan"))
              if val.get("n") else float("nan"))
        dt = time.time() - t0
        rec = {"epoch": epoch, "lr": args.lr * lr_at(epoch + 1),
               "train_loss": run_loss / nb, "train_pnll": run_pnll / nb,
               "train_cnt": run_cnt / nb, "val": val, "sec": round(dt, 1)}
        history.append(rec)
        print(f"[e{epoch:02d}] loss={run_loss/nb:.4f} pnll={run_pnll/nb:.4f} "
              f"val_pearson={vp:.4f} ({dt:.0f}s)", file=sys.stderr)
        (out / "history.json").write_text(json.dumps(history, indent=2))
        if vp > best:
            best = vp
            bad = 0
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            torch.save(best_state, out / "best.pt")
        else:
            bad += 1
            if bad >= args.patience:
                print(f"early stop at epoch {epoch} (best val_pearson={best:.4f})",
                      file=sys.stderr)
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    test = evaluate(model, test_loader, device, biotype, cap=None)
    (out / "test_metrics.json").write_text(json.dumps(test, indent=2))
    print("TEST:", json.dumps(test.get("all", {}), indent=2), file=sys.stderr)
    print(f"wrote {out}/test_metrics.json  best_val_pearson={best:.4f}", file=sys.stderr)


if __name__ == "__main__":
    main()
