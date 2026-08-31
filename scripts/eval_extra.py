#!/usr/bin/env python3
"""Extra evaluation for a trained run: magnitude (count-head) accuracy + a posture-B proxy.

The main test metric (profile Pearson) scores shape only. This adds:

  1. Magnitude metric -- per-transcript correlation between the count head's predicted
     log-count and the observed log1p(total P-sites), across test transcripts. This is
     where the RNAseq coverage input should earn its keep (coverage feeds the count head),
     so it completes the ablation picture that profile Pearson misses.

  2. Posture-B proxy -- profile Pearson restricted to test transcripts whose gene has a
     single annotated isoform. For those genes there is no within-gene isoform multimapping,
     so posture A (count every isoform) and posture B (primary only) are identical. If the
     single-isoform Pearson matches the overall Pearson, the ~18x multimap inflation is not
     distorting the shape result -- a cheap check that needs no target rebuild.

Loads the run's args.json and the checkpoint named by --ckpt (default best.pt), evaluates on
that run's held-out test fold, writes
<run>/extra_metrics.json. Reuses dataset.py + model.py + train.py helpers.
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
    load_split,
    loto_split,
    tissue_pack_dir,
)
from model import RiboSignalModel  # noqa: E402
from train import load_biotype, pearson, spearman  # noqa: E402

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
TX2B = NEW / "data" / "tx2biotype.tsv"
TX2CDS = NEW / "data" / "tx2cds.tsv"
MIN_SIGNAL = 50
# uORF-aware eval: uORFs (5'UTR) are the most common alternative ORF, but the whole-transcript
# profile Pearson is CDS-dominated. Score the profile Pearson restricted to the 5'UTR (uORF) and
# 3'UTR (dORF) windows, among pc transcripts whose window carries real ribosome signal (a
# translated uORF/dORF). Windows come from the annotation ONLY at eval time; training is
# annotation-free (the model never sees the CDS boundary).
MIN_UTR_LEN = 30       # need enough positions to correlate within the window
MIN_UTR_SIGNAL = 20    # P-sites in the window to call it a translated uORF/dORF and score it


def load_tx2cds():
    """tx_id -> (utr5_len, cds_len, utr3_len, has_start_codon); coding transcripts only."""
    out = {}
    if not TX2CDS.exists():
        return out
    with TX2CDS.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {k: hdr.index(k) for k in ("tx_id", "utr5_len", "cds_len", "utr3_len",
                                        "has_start_codon")}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[ix["tx_id"]]] = (int(f[ix["utr5_len"]]), int(f[ix["cds_len"]]),
                                   int(f[ix["utr3_len"]]), int(f[ix["has_start_codon"]]))
    return out


def gene_isoform_counts():
    """tx_id -> gene_id, and gene_id -> number of annotated transcripts (full annotation)."""
    tx2gene = {}
    gene_n = {}
    with TX2B.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ti, gi = hdr.index("tx_id"), hdr.index("gene_id")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            tx2gene[f[ti]] = f[gi]
            gene_n[f[gi]] = gene_n.get(f[gi], 0) + 1
    return tx2gene, gene_n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--eval_cap", type=int, default=0, help="0 = all test tx")
    ap.add_argument("--budget", type=int, default=24000)
    ap.add_argument("--device", default="cuda")
    # Checkpoint selection is the point of the audit comparison: train_loto.py saves best.pt on
    # val pearson_median and best_frame0.pt on the CDS-anchored frame-0 fraction, and the two are
    # different epochs in all 30 runs that have both. Being able to evaluate either on the SAME
    # held-out fold is what makes the criteria comparable.
    ap.add_argument("--ckpt", default="best.pt",
                    help="checkpoint filename inside --run (best.pt | best_frame0.pt)")
    ap.add_argument("--out-suffix", default="",
                    help="appended to the output filenames so two checkpoints do not overwrite")
    args = ap.parse_args()
    run = Path(args.run)
    if not run.is_absolute():
        run = NEW / run
    cfg = json.loads((run / "args.json").read_text())
    device = args.device if (args.device != "cuda" or torch.cuda.is_available()) else "cpu"

    backend = cfg.get("emb_backend", "rinalmo")
    biotype = load_biotype()
    tx2gene, gene_n = gene_isoform_counts()
    tx2cds = load_tx2cds()
    if cfg.get("holdout"):
        # LOTO run: test on the held-out tissue's own pack (its target + coverage + depth mean)
        holdout = cfg["holdout"]
        store = PackedStore(pack=tissue_pack_dir(holdout), tx_index=BACKENDS[backend])
        _, _, _, test_ids = loto_split(
            holdout, train_tissues=cfg.get("resolved_train_tissues"),
            min_signal=cfg.get("min_train_signal", 50), val_fold=cfg.get("val_fold", 0))
    else:
        store = PackedStore(tx_index=BACKENDS[backend])
        _, _, test_ids = load_split(store, test_fold=cfg.get("test_fold", 0))

    use_orf = cfg.get("use_orf_track", False)
    ds = RiboDataset(store, test_ids, input_mode=cfg.get("input_mode", "both"),
                     use_orf=use_orf, cov_norm=cfg.get("cov_norm", "raw"))
    lens = [store.length_of(t) for t in test_ids]
    # Cap the eval batch budget at the run's own training budget: attention runs train at a
    # lower budget (O(L^2) memory), and evaluating them at the 24000 default OOMs the encoder.
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
    model.load_state_dict(torch.load(run / args.ckpt, map_location=device))
    model.eval()

    rows = []   # (biotype, single_iso, profile_pearson, pred_logcount, obs_log1p_N)
    seen = 0
    with torch.no_grad():
        for batch in torch.utils.data.DataLoader(ds, batch_sampler=samp,
                                                 collate_fn=collate_pad, num_workers=4):
            feats = batch["feats"].to(device)
            mask = batch["mask"].to(device)
            logits, pred_lc = model(feats, mask)
            logp = torch.log_softmax(logits.masked_fill(~mask, -1e9), dim=1)
            p = torch.exp(logp).cpu().numpy()
            plc = pred_lc.cpu().numpy()
            counts = batch["counts"].numpy()
            for i, tx in enumerate(batch["tx"]):
                L = batch["lengths"][i]
                c = counts[i, :L]
                N = c.sum()
                if N < MIN_SIGNAL:
                    continue
                q = c / N
                g = tx2gene.get(tx)
                single = 1 if (g is not None and gene_n.get(g, 99) == 1) else 0
                # uORF (5'UTR) / dORF (3'UTR) profile Pearson: pc transcripts with a 5'-complete
                # annotated CDS, scored only where the window carries a translated ORF's signal.
                u5p = u3p = float("nan")
                u5n = u3n = 0
                cds = tx2cds.get(tx)
                if cds is not None:
                    u5len, _cdslen, u3len, has_sc = cds
                    if has_sc and MIN_UTR_LEN <= u5len <= L:
                        cu = c[:u5len]
                        u5n = int(cu.sum())
                        if u5n >= MIN_UTR_SIGNAL:
                            u5p = pearson(p[i, :u5len], cu)
                    if MIN_UTR_LEN <= u3len <= L:
                        cd = c[L - u3len:]
                        u3n = int(cd.sum())
                        if u3n >= MIN_UTR_SIGNAL:
                            u3p = pearson(p[i, L - u3len:L], cd)
                # tx..int(N) appended after r[0..4] so block()'s indices are unchanged; the
                # trailing fields feed the pertx.tsv dump + the uORF/dORF blocks.
                rows.append((biotype.get(tx, "unknown"), single,
                             pearson(p[i, :L], q), float(plc[i]), float(np.log1p(N)),
                             tx, int(N), u5p, u5n, u3p, u3n))
            seen += len(batch["tx"])
            if args.eval_cap and seen >= args.eval_cap:
                break

    def block(sel, label):
        sub = [r for r in rows if sel(r)]
        if len(sub) < 8:
            return {"n": len(sub)}
        prof = np.array([r[2] for r in sub])
        plc = np.array([r[3] for r in sub])
        obs = np.array([r[4] for r in sub])
        return {
            "label": label, "n": len(sub),
            "profile_pearson_median": float(np.median(prof)),
            "count_pearson": pearson(plc, obs),      # magnitude accuracy across tx
            "count_spearman": spearman(plc, obs),
        }

    pc = lambda r: r[0] == "protein_coding"  # noqa: E731

    def utr_block(idx, label):
        """Median profile Pearson within a UTR window over pc tx that had a scorable window."""
        vals = [r[idx] for r in rows if r[0] == "protein_coding" and r[idx] == r[idx]]
        if len(vals) < 8:
            return {"label": label, "n": len(vals)}
        return {"label": label, "n": len(vals),
                "profile_pearson_median": float(np.median(vals))}

    out = {
        "run": str(run.name),
        "emb_backend": backend,
        "input_mode": cfg.get("input_mode", "both"),
        "test_fold": cfg.get("test_fold", 0),
        "n_scored": len(rows),
        "protein_coding": block(pc, "protein_coding"),
        "lncRNA": block(lambda r: r[0] == "lncRNA", "lncRNA"),
        "pc_single_isoform_gene": block(lambda r: pc(r) and r[1] == 1,
                                        "pc single-isoform gene"),
        "pc_multi_isoform_gene": block(lambda r: pc(r) and r[1] == 0,
                                       "pc multi-isoform gene"),
        "uorf_5utr": utr_block(7, "pc 5'UTR (uORF)"),
        "dorf_3utr": utr_block(9, "pc 3'UTR (dORF)"),
    }
    out["checkpoint"] = args.ckpt
    (run / f"extra_metrics{args.out_suffix}.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2), file=sys.stderr)
    print(f"wrote {run / ('extra_metrics' + args.out_suffix + '.json')}", file=sys.stderr)

    # Per-transcript dump for the multifold pooled aggregation: pooling every fold's test tx
    # gives a lncRNA median over ~1,239 transcripts rather than ~264 per fold. Same MIN_SIGNAL
    # filter as train.py's evaluate(), so the pooled set matches the reported per-fold medians.
    def fx(x):
        return f"{x:.6f}" if isinstance(x, float) and x == x else "nan"

    with (run / f"pertx{args.out_suffix}.tsv").open("w") as fh:
        fh.write("tx_id\tbiotype\tn_psites\tprofile_pearson\t"
                 "utr5_pearson\tutr5_psites\tutr3_pearson\tutr3_psites\n")
        for r in rows:
            fh.write(f"{r[5]}\t{r[0]}\t{r[6]}\t{r[2]:.6f}\t"
                     f"{fx(r[7])}\t{r[8]}\t{fx(r[9])}\t{r[10]}\n")
    print(f"wrote {run / ('pertx' + args.out_suffix + '.tsv')} ({len(rows)} tx)", file=sys.stderr)


if __name__ == "__main__":
    main()
