#!/usr/bin/env python3
"""Assemble the cross-study / cross-species held-out comparison table across backends:
profile Pearson (shape) + localization non-canonical AUROC + RiboCode drop-in F1.

Profile Pearson is computed here directly from dump_pred_profiles.py's pred_profiles.npz
(per-tx predicted softmax profile vs observed pooled P-site counts; Pearson is scale-invariant,
so pred-vs-counts == pred-vs-(counts/sum)), split protein_coding vs lncRNA. Localization and
drop-in are read from the eval JSONs. Mirrors the LOTO comparison table for the held-out sets.

Usage: assemble_heldout.py --dataset <name> [--backends onehot orthrus rinalmo]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")


def load_biotype(dataset):
    """tx_id -> 'protein_coding' | 'lncRNA' | other. Human from tx2biotype; mouse from the
    universe tsv."""
    out = {}
    if dataset.startswith("mouse"):
        p = NEW / "data" / "heldout_refs" / "mouse_wang_universe.tsv"
        with p.open() as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ti, bi = hdr.index("tx_id"), hdr.index("biotype")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                out[f[ti]] = f[bi]
    else:
        p = NEW / "data" / "tx2biotype.tsv"
        with p.open() as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ti, bi = hdr.index("tx_id"), hdr.index("transcript_type")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                out[f[ti]] = f[bi]
    return out


def pearson(a, b):
    a = a.astype(np.float64)
    b = b.astype(np.float64)
    a -= a.mean()
    b -= b.mean()
    na, nb = np.sqrt((a * a).sum()), np.sqrt((b * b).sum())
    if na == 0 or nb == 0:
        return np.nan
    return float((a * b).sum() / (na * nb))


def profile_pearson(npz_path, biotype, min_signal=50):
    d = np.load(npz_path, allow_pickle=True)
    tx = d["tx_ids"]
    lens = d["lengths"]
    pred = d["pred_flat"]
    obs = d["obs_flat"]
    off = np.zeros(len(lens) + 1, dtype=np.int64)
    off[1:] = np.cumsum(lens)
    pc, lnc = [], []
    for k, t in enumerate(tx):
        s, e = off[k], off[k + 1]
        o = obs[s:e]
        if o.sum() < min_signal:
            continue
        r = pearson(pred[s:e], o)
        if r != r:
            continue
        bt = biotype.get(str(t), "")
        if bt == "protein_coding":
            pc.append(r)
        elif bt == "lncRNA":
            lnc.append(r)
    med = lambda v: float(np.median(v)) if v else float("nan")  # noqa: E731
    return med(pc), med(lnc), len(pc), len(lnc)


def loc_metrics(run_out, dataset):
    p = run_out / f"localization_metrics_{dataset}.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text())["discrimination"]
    nc = d["noncanonical"]
    return (nc["auroc_pred_frame0_lengthctrl"],
            nc["auroc_obs_frame0_ceiling_lengthctrl"],
            d["novel"]["auroc_pred_frame0_lengthctrl"])


def dropin_metrics(run_out):
    p = run_out / "dropin" / "dropin_metrics.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    sa = d.get("pred_obsdepth_vs_real", {})
    nv = sa.get("recall_per_type", {}).get("novel", {})
    return (sa.get("precision"), sa.get("recall"), sa.get("f1"),
            nv.get("recall"), nv.get("recovered"), nv.get("n_real"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--backends", nargs="+", default=["onehot", "orthrus", "rinalmo"])
    args = ap.parse_args()
    biotype = load_biotype(args.dataset)
    root = NEW / "results" / "heldout" / args.dataset

    print(f"\n=== held-out {args.dataset}: backend comparison (orf_v2_attn) ===\n")
    print(f"{'backend':8} {'profile pc':>11} {'profile lnc':>12} {'n_pc':>7} {'n_lnc':>6} | "
          f"{'noncanon AUROC':>15} {'% ceil':>7} {'novel AUROC':>12} | "
          f"{'drop-in F1':>10} {'prec':>6} {'recall':>7} {'novel recall':>13}")
    for bk in args.backends:
        run_out = root / bk
        npz = run_out / "dropin" / "pred_profiles.npz"
        cells = [f"{bk:8}"]
        if npz.exists():
            pc, lnc, npc, nlnc = profile_pearson(npz, biotype)
            cells.append(f"{pc:>11.4f} {lnc:>12.4f} {npc:>7} {nlnc:>6}")
        else:
            cells.append(f"{'(no npz)':>11} {'':>12} {'':>7} {'':>6}")
        lm = loc_metrics(run_out, args.dataset)
        if lm:
            cells.append(f"| {lm[0]:>15.3f} {100 * lm[0] / lm[1]:>6.1f}% {lm[2]:>12.3f}")
        else:
            cells.append(f"| {'(pending)':>15} {'':>7} {'':>12}")
        dm = dropin_metrics(run_out)
        if dm and dm[2] is not None:
            nr = f"{dm[3]:.3f} ({dm[4]}/{dm[5]})" if dm[3] is not None else "NA"
            cells.append(f"| {dm[2]:>10.3f} {dm[0]:>6.3f} {dm[1]:>7.3f} {nr:>13}")
        else:
            cells.append(f"| {'(pending)':>10} {'':>6} {'':>7} {'':>13}")
        print(" ".join(cells))
    print()


if __name__ == "__main__":
    main()
