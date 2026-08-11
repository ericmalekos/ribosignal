#!/usr/bin/env python3
"""Phase 3.1: sequence saliency aggregated over thousands of ORFs, by class and by position.

Figure C10 did this for THREE exemplar transcripts and found a suggestive asymmetry: the start codon
ranked 3rd in saliency for a uORF and 10th for a lncRNA ORF, but 292nd for a canonical CDS. With n=3
that is an anecdote. This runs the same gradient over thousands of ORFs and aggregates, so the
asymmetry either survives as a distribution or does not.

METHOD (identical gradient to C10, so the two are comparable)
  target = sum of predicted profile log-probability over the ORF's IN-FRAME (frame-0) positions
  saliency = |d target / d one-hot|, summed over the 4 nucleotide channels, per position
  The coverage channel's gradient is recorded separately: it separates "reads the sequence" from
  "reads the RNA-seq", which is the distinction the internal-ORF failure hinges on.

Per ORF we record:
  * start_rank        rank of the start codon's saliency among all positions (1 = highest)
  * in_orf_ratio      mean saliency inside the ORF / mean over the transcript
  * seq_vs_cov        total sequence saliency / (sequence + coverage) saliency
  * a window of saliency around the start codon, frame-resolved, for the positional aggregate

SUBSTRATE MATTERS HERE. Mouse is primary because its transcripts were never in training; human
Hepatocytes is retained as the deliberately CONTAMINATED arm, since 99.7% of its test transcripts
were seen with byte-identical sequence and gradient exposure is worth +0.068 profile correlation
(see docs/attribution_experiments_plan.md and the memorisation page). A feature present on BOTH
substrates is a learned rule; one present only on human is a memorisation candidate.
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys
from collections import defaultdict

import numpy as np

NEW = pathlib.Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model")
sys.path.insert(0, str(NEW / "scripts"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="model run dir containing best.pt + args.json")
    ap.add_argument("--pack", required=True)
    ap.add_argument("--orf-track", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--collapsed", required=True,
                    help="a *_collapsed.txt giving ORF class + coordinates")
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-class", type=int, default=300, help="ORFs sampled per class")
    ap.add_argument("--window", type=int, default=60, help="+/- nt around the start codon")
    ap.add_argument("--seed", type=int, default=0)
    # mamba_ssm dispatches to causal_conv1d_cuda and raises "Expected x.is_cuda() to be true" on CPU
    # tensors, even when the module imports fine inside the GPU image. attn has no such constraint.
    ap.add_argument("--device", default="cpu", choices=["cpu", "cuda"])
    a = ap.parse_args()

    os.environ["RIBO_ORF_TRACK"] = a.orf_track
    os.environ["RIBO_ONEHOT_FASTA"] = a.fasta
    import torch
    from dataset import PackedStore, RiboDataset, BACKENDS
    from model import RiboSignalModel

    run = pathlib.Path(a.run)
    cfg = json.load(open(run / "args.json")) if (run / "args.json").exists() else {}
    C = dict(channels=cfg.get("channels", 256), n_blocks=cfg.get("n_blocks", 10),
             n_attn_layers=cfg.get("n_attn_layers", 2), n_heads=cfg.get("n_heads", 8),
             input_mode=cfg.get("input_mode", "both"), cov_norm=cfg.get("cov_norm", "global_mean"),
             use_orf=cfg.get("use_orf_track", True), mixer=cfg.get("mixer", "attn"),
             learn_start_context=cfg.get("learn_start_context", False))

    # ---- ORFs to attribute, sampled per class from a collapsed call table ----
    rng = np.random.default_rng(a.seed)
    by_class = defaultdict(list)
    with open(a.collapsed) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(hdr)}
        need = ("ORF_type", "transcript_id", "ORF_tstart", "ORF_tstop")
        if not all(k in ci for k in need):
            print(f"  collapsed table lacks {need}; has {hdr[:12]}", file=sys.stderr)
            return 1
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                by_class[f[ci["ORF_type"]]].append(
                    (f[ci["transcript_id"]], int(f[ci["ORF_tstart"]]), int(f[ci["ORF_tstop"]])))
            except (ValueError, IndexError):
                continue
    picks = []
    for cls, items in by_class.items():
        idx = rng.permutation(len(items))[:a.per_class]
        picks += [(cls, *items[i]) for i in idx]
    print(f"  {a.label}: {len(picks)} ORFs over {len(by_class)} classes "
          f"({', '.join(f'{k}:{min(len(v), a.per_class)}' for k, v in sorted(by_class.items()))})")

    store = PackedStore(pack=pathlib.Path(a.pack), tx_index=BACKENDS["onehot"])
    keep = [p for p in picks if p[1] in store.row]
    print(f"  {len(keep)} of {len(picks)} ORFs on transcripts present in the pack")
    txs = sorted({p[1] for p in keep})
    ds = RiboDataset(store, txs, input_mode=C["input_mode"], use_orf=C["use_orf"],
                     cov_norm=C["cov_norm"])
    row = {t: i for i, t in enumerate(txs)}
    d_emb = store.load_emb(txs[0]).shape[1]
    model = RiboSignalModel(d_emb=d_emb, channels=C["channels"], n_blocks=C["n_blocks"],
                            dropout=0.0, extra_in=store.orf_channels if C["use_orf"] else 0,
                            n_attn_layers=C["n_attn_layers"], n_heads=C["n_heads"],
                            mixer=C["mixer"], learn_start_context=C["learn_start_context"])
    model.load_state_dict(torch.load(run / "best.pt", map_location="cpu"))
    model.eval()
    dev = torch.device(a.device)
    model.to(dev)

    W = a.window
    per_orf, prof = [], defaultdict(lambda: np.zeros(2 * W + 1))
    prof_n = defaultdict(int)
    cache = {}
    for n, (cls, tx, ts, te) in enumerate(keep):
        if te <= ts or te - ts < 30:
            continue
        if tx not in cache:
            cache.clear()                       # one transcript resident at a time; keeps RAM flat
            cache[tx] = ds[row[tx]]
        item = cache[tx]
        feats = item["feats"].unsqueeze(0).clone().to(dev).requires_grad_(True)
        L = feats.shape[1]
        if te > L:
            continue
        mask = torch.ones(1, L, dtype=torch.bool, device=dev)
        logits, _ = model(feats, mask)
        logp = torch.log_softmax(logits[0], dim=0)
        fr0 = torch.arange(ts, te, 3, device=dev)
        target = logp[fr0].sum()
        model.zero_grad(set_to_none=True)
        if feats.grad is not None:
            feats.grad = None
        target.backward()
        g = feats.grad[0].detach().cpu().numpy()
        sal = np.abs(g[:, :4]).sum(axis=1)                # sequence saliency per nt
        cov = float(np.abs(g[:, -1]).sum())               # coverage-channel saliency, whole tx
        s_tot = float(sal.sum())
        per_orf.append(dict(
            label=a.label, orf_class=cls, tx=tx, L=int(L), orf_len=int(te - ts),
            start_rank=int((sal >= sal[ts]).sum()),
            in_orf_ratio=float(sal[ts:te].mean() / max(sal.mean(), 1e-12)),
            seq_frac=float(s_tot / max(s_tot + cov, 1e-12))))
        lo, hi = ts - W, ts + W + 1
        if lo >= 0 and hi <= L:
            m = sal.mean()
            if m > 0:
                prof[cls] += sal[lo:hi] / m
                prof_n[cls] += 1
        if (n + 1) % 200 == 0:
            print(f"    {n+1}/{len(keep)} ORFs", flush=True)

    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("label\torf_class\ttx\tL\torf_len\tstart_rank\tin_orf_ratio\tseq_frac\n")
        for r in per_orf:
            fh.write(f"{r['label']}\t{r['orf_class']}\t{r['tx']}\t{r['L']}\t{r['orf_len']}\t"
                     f"{r['start_rank']}\t{r['in_orf_ratio']:.6f}\t{r['seq_frac']:.6f}\n")
    pf = out.with_name(out.stem + "_startprofile.tsv")
    with open(pf, "w") as fh:
        fh.write("label\torf_class\toffset\tmean_rel_saliency\tn\n")
        for cls, v in prof.items():
            if prof_n[cls]:
                for i, x in enumerate(v / prof_n[cls]):
                    fh.write(f"{a.label}\t{cls}\t{i - W}\t{x:.6f}\t{prof_n[cls]}\n")
    print(f"  wrote {out.name} ({len(per_orf)} ORFs) + {pf.name}")
    for cls in sorted({r["orf_class"] for r in per_orf}):
        v = [r for r in per_orf if r["orf_class"] == cls]
        print(f"    {cls:<14} n={len(v):<5} median start_rank={np.median([x['start_rank'] for x in v]):>7.0f}"
              f"  in_orf_ratio={np.median([x['in_orf_ratio'] for x in v]):.2f}"
              f"  seq_frac={np.median([x['seq_frac'] for x in v]):.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
