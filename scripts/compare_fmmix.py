#!/usr/bin/env python3
"""Compare the FM-embedding sweep against its direct-to-mixer (fm_to_mixer) twin.

Question: the FM sweep (Task 41) found NO lift from RiNALMo / Orthrus / HydraRNA embeddings over a 4-d
one-hot input. But by construction the embedding reaches the global mixer only after being squeezed to
`channels` by a SHARED 1x1 in_proj -- jointly with the ORF track and coverage -- and then rewritten by
10 dilated conv blocks. For RiNALMo that is 1280 -> 256 in a single 1x1 conv before any spatial
processing, a 5x bottleneck the 4-d one-hot backend never pays. So the sweep confounded "does the FM
carry useful signal" with "does the FM survive this fusion".

The `--fm_to_mixer` arm adds a dedicated projection of the raw embedding straight onto the mixer input,
gated by a scalar initialised to ZERO. Training therefore starts numerically identical to the baseline
and must actively open the path, which makes the learned |gate| a direct readout of how much the model
wants the un-convolved embedding. A gate that stays near zero is an interpretable NEGATIVE result: the
model was handed the embedding and declined it.

cas12a env for the metrics; torch is only needed for --gate (reads best.pt).
"""
import argparse
import json
import os

NEW = "/private/groups/carpenterlab/emalekos/RNAZoo_meta/RNAZoo/experiments/riboseq_signal_model"
DIMS = {"onehot": 4, "rinalmo": 1280, "orthrus": 512, "hydrarna": 1024}


def run_dir(backend, fmmix):
    return (f"{NEW}/results/loto/orf_v2_attn_{backend}_noBrain_nokozak_mm1_holdout_Hepatocytes"
            + ("_fmmix" if fmmix else ""))


def metrics(d):
    f = f"{d}/test_metrics.json"
    if not os.path.exists(f):
        return None
    m = json.load(open(f))
    return {"pc": m["protein_coding"]["pearson_median"],
            "per": m["protein_coding"]["period_pred_median"],
            "lnc": m["lncRNA"]["pearson_median"]}


def gate(d):
    """Learned scalar gate on the direct FM->mixer path (None for the baseline arm)."""
    f = f"{d}/best.pt"
    if not os.path.exists(f):
        return None
    try:
        import torch
        sd = torch.load(f, map_location="cpu")
    except Exception:
        return None
    for k in ("fm_gate", "module.fm_gate"):
        if k in sd:
            return float(sd[k].reshape(-1)[0])
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gate", action="store_true", help="also read the learned gate from best.pt (needs torch)")
    a = ap.parse_args()

    base_oh = metrics(run_dir("onehot", False))
    ref = base_oh["pc"] if base_oh else None
    print(f"Hepatocytes hold-out, n=33,918 (32,869 pc + 1,049 lncRNA). "
          f"one-hot reference pc Pearson = {ref:.4f}\n" if ref else "")
    hdr = f"  {'backend':<10}{'d_emb':>6}{'arm':>10}{'pc Pearson':>12}{'vs 1hot':>9}{'vs base':>9}{'period':>9}{'lncRNA':>9}"
    if a.gate:
        hdr += f"{'gate':>9}"
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    for b in ("onehot", "rinalmo", "orthrus", "hydrarna"):
        base = metrics(run_dir(b, False))
        for fmmix in (False, True):
            if b == "onehot" and fmmix:
                continue
            d = run_dir(b, fmmix)
            m = metrics(d)
            arm = "fm->mixer" if fmmix else "baseline"
            if m is None:
                print(f"  {b:<10}{DIMS[b]:>6}{arm:>10}{'(pending)':>12}")
                continue
            v1 = f"{m['pc'] - ref:+.4f}" if ref is not None and b != "onehot" else "--"
            v2 = f"{m['pc'] - base['pc']:+.4f}" if (fmmix and base) else "--"
            line = (f"  {b:<10}{DIMS[b]:>6}{arm:>10}{m['pc']:>12.4f}{v1:>9}{v2:>9}"
                    f"{m['per']:>9.4f}{m['lnc']:>9.4f}")
            if a.gate:
                g = gate(d)
                line += f"{(f'{abs(g):.4f}' if g is not None else '--'):>9}"
            print(line)
    print("\n  'vs base' is the controlled A/B: SAME backend, SAME data, SAME hyperparameters, "
          "only the\n  fusion point differs. 'gate' near 0 => the model declined the direct embedding path.")


if __name__ == "__main__":
    main()
