#!/usr/bin/env python3
"""Qualitative + summary figures for a trained fold-0 Ribo-seq signal model.

Loads a run's best checkpoint, predicts per-nt P-site profiles on the held-out test
chromosomes, and makes a figure with (a) observed vs predicted profile for a few example
transcripts with a periodicity zoom, and (b) a summary of per-transcript Pearson and
observed-vs-predicted periodicity, split by biotype. Runs on CPU (cas12a) or in the SIF.

Usage:
  plot_predictions.py --run results/fold0_rinalmo [--n_examples 3] [--eval_cap 3000]
Writes:
  <run>/figures/fold0_predictions.png
  <run>/figures/FIGURE_DATA_INPUTS.md
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import torch  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dataset import PackedStore, RiboDataset, load_split  # noqa: E402
from model import RiboSignalModel  # noqa: E402
from train import load_biotype, pearson, periodicity  # noqa: E402

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")


@torch.no_grad()
def predict_one(model, ds, k, device):
    it = ds[k]
    feats = it["feats"].unsqueeze(0).to(device)
    mask = torch.ones(1, it["length"], dtype=torch.bool, device=device)
    logits, _ = model(feats, mask)
    p = torch.softmax(logits[0], dim=0).cpu().numpy()
    c = it["counts"].numpy()
    return c, p, it["tx"], it["length"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--n_examples", type=int, default=3)
    ap.add_argument("--eval_cap", type=int, default=3000)
    ap.add_argument("--device", default="cpu")
    args = ap.parse_args()
    run = Path(args.run)
    if not run.is_absolute():
        run = NEW / run
    figdir = run / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    cfg = json.loads((run / "args.json").read_text())
    store = PackedStore()
    biotype = load_biotype()
    _, _, test_ids = load_split(store)
    d_emb = np.load(store.emb_path[test_ids[0]]).shape[1]
    # Read the full architecture out of args.json, as every other loader does. This previously passed
    # only channels/n_blocks, so load_state_dict raised on ANY orf_v2 checkpoint (it missed the
    # ORF-track input width, the mixer stack, and the learned-Kozak gate). Pre-existing bug, fixed here.
    model = RiboSignalModel(d_emb=d_emb, channels=cfg["channels"],
                            n_blocks=cfg["n_blocks"], dropout=0.0,
                            extra_in=(store.orf_channels
                                      if getattr(store, "orf_track", None) is not None else 0),
                            n_attn_layers=cfg.get("n_attn_layers", 0),
                            n_heads=cfg.get("n_heads", 8),
                            mixer=cfg.get("mixer", "transformer"),
                            mamba_d_state=cfg.get("mamba_d_state", 16),
                            mamba_d_conv=cfg.get("mamba_d_conv", 4),
                            mamba_expand=cfg.get("mamba_expand", 2),
                            learn_start_context=cfg.get("learn_start_context", False),
                            fm_to_mixer=cfg.get("fm_to_mixer", False))
    model.load_state_dict(torch.load(run / "best.pt", map_location=args.device))
    model.eval().to(args.device)

    ds = RiboDataset(store, test_ids)
    # scan a capped subset: collect per-tx metrics + remember well-predicted pc examples
    stats = []   # (tx, biotype, pearson, per_pred, per_obs, total)
    examples = []
    n = min(args.eval_cap, len(test_ids))
    for k in range(n):
        c, p, tx, L = predict_one(model, ds, k, args.device)
        N = c.sum()
        if N < 50:
            continue
        q = c / N
        r = pearson(p, q)
        bt = biotype.get(tx, "unknown")
        stats.append((tx, bt, r, periodicity(p), periodicity(q), float(N)))
        if bt == "protein_coding" and 400 < L < 4000:
            examples.append((r, k, tx, c, p, L))
    examples.sort(reverse=True)  # best-predicted pc first
    examples = examples[:args.n_examples]

    bt = np.array([s[1] for s in stats])
    rr = np.array([s[2] for s in stats])
    pp = np.array([s[3] for s in stats])
    oo = np.array([s[4] for s in stats])
    is_pc = bt == "protein_coding"
    is_lnc = bt == "lncRNA"

    nrows = len(examples) + 1
    fig, axes = plt.subplots(nrows, 2, figsize=(13, 2.6 * nrows))
    if nrows == 1:
        axes = axes[None, :]

    for i, (r, _k, tx, c, p, L) in enumerate(examples):
        N = c.sum()
        q = c / N
        ax = axes[i, 0]
        ax.plot(q, lw=0.5, color="#1f77b4", label="observed")
        ax.plot(p, lw=0.5, color="#d62728", alpha=0.8, label="predicted")
        ax.set_title(f"{tx}  L={L}  Pearson={r:.2f}", fontsize=9)
        ax.set_ylabel("profile")
        ax.legend(fontsize=7, loc="upper right")
        # periodicity zoom: a 90-nt window at the peak observed density
        w = 90
        center = int(np.argmax(np.convolve(q, np.ones(30) / 30, mode="same")))
        a = max(0, center - w // 2)
        b = min(L, a + w)
        az = axes[i, 1]
        az.bar(np.arange(a, b), q[a:b], width=0.9, color="#1f77b4", label="observed")
        az.plot(np.arange(a, b), p[a:b], color="#d62728", lw=1.0, label="predicted")
        az.set_title(f"periodicity zoom [{a}:{b}]", fontsize=9)
        az.legend(fontsize=7, loc="upper right")

    # summary row
    axs = axes[-1, 0]
    axs.hist(rr[is_pc], bins=40, range=(-0.2, 1), alpha=0.6, label=f"pc (n={is_pc.sum()})",
             color="#2ca02c", density=True)
    if is_lnc.sum():
        axs.hist(rr[is_lnc], bins=40, range=(-0.2, 1), alpha=0.6,
                 label=f"lncRNA (n={is_lnc.sum()})", color="#9467bd", density=True)
    axs.axvline(np.median(rr[is_pc]), color="#2ca02c", ls="--", lw=1)
    axs.set_title("per-transcript profile Pearson (test)", fontsize=9)
    axs.set_xlabel("Pearson r")
    axs.legend(fontsize=7)

    axp = axes[-1, 1]
    axp.scatter(oo[is_pc], pp[is_pc], s=4, alpha=0.3, color="#2ca02c", label="pc")
    if is_lnc.sum():
        axp.scatter(oo[is_lnc], pp[is_lnc], s=4, alpha=0.3, color="#9467bd", label="lncRNA")
    lim = [min(oo.min(), pp.min()), max(oo.max(), pp.max())]
    axp.plot(lim, lim, color="k", lw=0.6, ls=":")
    axp.set_title("periodicity: observed vs predicted", fontsize=9)
    axp.set_xlabel("observed")
    axp.set_ylabel("predicted")
    axp.legend(fontsize=7)

    fig.tight_layout()
    out_png = figdir / "fold0_predictions.png"
    fig.savefig(out_png, dpi=140)
    print(f"wrote {out_png}  (test tx scored: {len(stats)})", file=sys.stderr)

    med_pc = float(np.median(rr[is_pc])) if is_pc.sum() else float("nan")
    med_lnc = float(np.median(rr[is_lnc])) if is_lnc.sum() else float("nan")
    (figdir / "FIGURE_DATA_INPUTS.md").write_text(
        "# fold0_predictions.png data inputs\n\n"
        f"Generated by `scripts/plot_predictions.py --run {args.run}`.\n\n"
        "Every panel derives from the trained model applied to the held-out test\n"
        "chromosomes (fold 0). Inputs:\n\n"
        f"- Checkpoint: `{run}/best.pt` (config `{run}/args.json`).\n"
        "- Model input features per transcript: RiNALMo per-token embeddings\n"
        "  (`data/rinalmo_token_emb/`, via `data/rinalmo_token_emb/tx_index.tsv`) plus\n"
        "  `log1p` of the pooled RNAseq coverage from `data/packed/coverage.npy`.\n"
        "- Observed P-site profile (the label): `data/packed/target_counts.npy`\n"
        "  (packed from `data/target/Fibroblast_psites_pooled.hd5`).\n"
        "- Test transcript set: fold 0 of `data/splits/fibroblast_chrom_kfold.json`\n"
        "  (chromosomes chr1/7/8/22, gene-disjoint from train/val).\n"
        "- Biotype labels: `data/packed/pack_meta.tsv` (transcript_type).\n\n"
        "Example rows: the best-predicted protein_coding test transcripts with\n"
        "400 < L < 4000. Summary row: all scored test transcripts with >= 50 total\n"
        "P-sites.\n\n"
        f"Headline: median per-transcript profile Pearson pc={med_pc:.3f}, "
        f"lncRNA={med_lnc:.3f}.\n"
    )
    print(f"wrote {figdir}/FIGURE_DATA_INPUTS.md", file=sys.stderr)


if __name__ == "__main__":
    main()
