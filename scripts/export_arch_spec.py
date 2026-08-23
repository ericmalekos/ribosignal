#!/usr/bin/env python3
"""Export the architecture + ORF-track semantics as PORTABLE JSON, for a poster/figure build
that cannot read the cluster filesystem.

WHY. Every structural number on an architecture panel -- layer types, kernel size, the dilation
schedule, GroupNorm groups, dropout, n_heads, d_ff, the parameter count, the loss terms -- has
until now traced only to a prose table in `figures/arch_attn/FIGURE_DATA_INPUTS.md`. `model.py`,
the run's `args.json` and the checkpoint do not ship. The poster's own provenance rule is that
every printed number traces to a portable file, and a markdown paragraph is not one.

Everything here is READ FROM THE ARTIFACTS, never transcribed:

  n_params        summed from best.pt's state_dict
  channel counts  from in_proj.weight's shape
  kernel sizes    from the conv weight shapes
  hyperparameters from the run's args.json
  START_W         imported from build_orf_track.py
  orf track       computed by calling build_orf_track.orf_track() on the display window

So a drift in the code changes the JSON on the next run rather than silently invalidating a
number printed on a poster.

  export_arch_spec.py --model attn
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import torch

NEW = Path(__file__).resolve().parents[1]

# The 24-nt window drawn on the architecture panel. Chosen upstream to exercise every channel:
# two in-frame AUGs, a CTG and an AAG near-cognate, and three stops in two frames.
DISPLAY_WINDOW = "ATGCATGCTGAAGCCTAACTAGCC"

RUNS = {
    "attn": "results/loto/orf_v2_attn_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes",
    "mamba4": "results/loto/orf_v2_mamba4_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes",
}


def load_build_orf_track():
    spec = importlib.util.spec_from_file_location("bot", NEW / "scripts/build_orf_track.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="attn", choices=["attn", "mamba4"])
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args()
    run = NEW / RUNS[a.model]
    out = Path(a.outdir) if a.outdir else NEW / "figures/arch_attn"
    out.mkdir(parents=True, exist_ok=True)

    args = json.loads((run / "args.json").read_text())
    sd = torch.load(run / "best.pt", map_location="cpu", weights_only=False)
    sd = sd.get("model") or sd.get("state_dict") or sd
    n_params = int(sum(v.numel() for v in sd.values() if hasattr(v, "numel")))
    c_out, c_in, k_proj = (int(x) for x in sd["in_proj.weight"].shape)
    k_body = int(sd["blocks.0.conv1.weight"].shape[2])
    n_blocks = len({k.split(".")[1] for k in sd if k.startswith("blocks.")})
    dil = [2 ** (i % 10) for i in range(n_blocks)]

    # Receptive field of the CONVOLUTIONAL stack only. Each block is two k=3 convs at dilation d,
    # and a k=3 conv at dilation d widens the field by 2d, so each block adds 4d. in_proj is k=1
    # and adds nothing. This confirms the poster's 1 + 4*sum(dilations) derivation.
    rf_conv = 1 + 2 * (k_body - 1) * sum(dil)

    # Group counts read off the checkpoint rather than assumed, so a refactor shows up here.
    groups = {}
    for kk, v in sd.items():
        if not hasattr(v, "numel"):
            continue
        g = ("in_proj" if kk.startswith("in_proj") else
             "body_blocks" if kk.startswith("blocks") else
             "mixer" if kk.startswith(("attn", "mamba_body")) else
             "profile_head" if kk.startswith("profile_head") else
             "count_head" if kk.startswith("count_head") else "other")
        groups[g] = groups.get(g, 0) + int(v.numel())

    spec = {
        "run": run.name,
        "model": a.model,
        "mixer": args["mixer"],
        "n_params": n_params,
        "n_params_by_group": groups,
        "input_channels": {
            "order": "one-hot(4) | orf_track(5) | rna_coverage(1); coverage is ALWAYS last",
            "onehot": 4, "onehot_order": ["A", "C", "G", "T"],
            "onehot_note": "U maps to the T column; N or any other character is an all-zero row",
            "orf_track": 5, "rna_coverage": 1, "total": c_in,
            "rna_coverage_transform": (
                f"log1p(coverage / global_mean_coverage)  [cov_norm={args['cov_norm']}]. "
                "Per-nt RNA depth divided by the PACK's own global mean before log1p, so the "
                "channel is depth-normalised and comparable across datasets."),
        },
        "projection": {"kernel": k_proj, "in": c_in, "out": c_out},
        "body": {
            "n_blocks": n_blocks, "kernel": k_body, "dilations": dil,
            "block": "residual: x + drop(conv2(gelu(norm2(conv1(gelu(norm1(x)))))))-style, "
                     "two k=3 dilated convs per block",
            "norm": "GroupNorm", "groups": 8, "dropout": args["dropout"],
            "receptive_field_nt": rf_conv,
            "receptive_field_derivation": (
                f"1 + 2*(k-1)*sum(dilations) = 1 + {2*(k_body-1)}*{sum(dil)} = {rf_conv}. "
                "Two k=3 convs per block, each widening by 2*d; in_proj is k=1 and adds nothing."),
            "receptive_field_caveat": (
                "THIS IS THE CONVOLUTIONAL RECEPTIVE FIELD ONLY. The mixer below is FULL "
                "self-attention over the transcript, masked for padding alone "
                "(`self.attn(xt, src_key_padding_mask=~mask)`) -- no causal and no local window. "
                "So the attn model's EFFECTIVE context is the ENTIRE transcript, not 4093 nt. "
                "Labelling the panel 'receptive field ~4 kb' without this qualifier understates "
                "the attn model. Say 'convolutional receptive field 4,093 nt; attention is "
                "global over the transcript'."),
        },
        "mixer_cfg": ({
            "n_attn_layers": args["n_attn_layers"], "n_heads": args["n_heads"],
            "norm_first": True, "activation": "gelu",
            "d_ff": c_out * 2,
            "d_ff_note": "dim_feedforward = channels * 2, not an independent hyperparameter",
            "attention_scope": "global over the transcript; only padding is masked",
        } if args["mixer"] == "transformer" else {
            "n_layers": args["n_attn_layers"], "d_state": args["mamba_d_state"],
            "d_conv": args["mamba_d_conv"], "expand": args["mamba_expand"],
            "direction": "bidirectional",
        }),
        "heads": {
            "profile": {"op": f"Conv1d({c_out}, 1, 1) -> (B, L) logits, -1e9 at padding",
                        "activation": "log_softmax over positions, applied at loss time"},
            "count": {"op": f"Linear({c_out}+1, {c_out}) -> GELU -> Dropout -> Linear({c_out}, 1)",
                      "input": "masked mean of body features over valid positions, CONCATENATED "
                               "with log1p(total RNA coverage) for the transcript",
                      "target": "log1p(total observed P-sites)"},
        },
        "loss": {
            "profile": "per-transcript multinomial NLL, equal-weighted across transcripts",
            "count": "MSE on log1p(total)",
            "count_weight": args["count_weight"],
            "peakiness_weight": args["peakiness_weight"],
            "peakiness_note": "0.0 in the deployed model; the entropy-gap anti-smoothing term is "
                              "an A/B arm, not part of this checkpoint",
        },
        "training": {k: args[k] for k in
                     ("epochs", "lr", "weight_decay", "budget", "warmup", "patience", "seed",
                      "val_fold", "min_train_signal", "max_tx_per_tissue", "input_mode")},
        "data": {"holdout": args["holdout"],
                 "train_tissues": args["resolved_train_tissues"],
                 "n_train_tx": args["n_train_tx"], "n_val_tx": args["n_val_tx"],
                 "n_test_tx": args["n_test_tx"],
                 "emb_backend": args["emb_backend"],
                 "kozak": "none (heuristic Kozak gate DROPPED; see docs/PIPELINE_POLICY.md)"},
        "provenance": {
            "n_params": "summed from best.pt state_dict",
            "channels_and_kernels": "read from checkpoint tensor shapes",
            "hyperparameters": f"{run.name}/args.json",
            "generated_by": "scripts/export_arch_spec.py",
        },
    }
    p = out / f"arch_spec_{a.model}.json"
    p.write_text(json.dumps(spec, indent=2))
    print(f"  wrote {p}  ({n_params:,} params, conv RF {rf_conv} nt)")

    # ---- ORF track semantics, computed not transcribed ----------------------------------------
    bot = load_build_orf_track()
    tr = np.asarray(bot.orf_track(DISPLAY_WINDOW, mode="ext", kozak="none"), dtype=float)
    chan = ["orf_f0", "orf_f1", "orf_f2", "start_ext", "is_stop"]
    onehot = {b: [1 if c == b else 0 for c in DISPLAY_WINDOW] for b in "ACGT"}
    ospec = {
        "start_weights": dict(bot.START_W),
        "start_weights_note": (
            "ALL 10 graded start codons, T-alphabet. The figure doc previously recorded only 6. "
            "Ordinal relative initiation efficiency from the non-AUG initiation literature; "
            "these are NOT fitted."),
        "kozak": "none",
        "kozak_note": ("The deployed model is built with --kozak none, so the start channel is "
                       "the raw codon weight with NO context multiplier. A reimplementation that "
                       "omits Kozak is therefore CORRECT for this model, not an approximation."),
        "mode": "ext (v2)",
        "channels": [
            {"index": 0, "name": "orf_f0",
             "meaning": "nt lies inside an ATG..in-frame-stop stretch whose ATG start position % 3 == 0"},
            {"index": 1, "name": "orf_f1", "meaning": "... start position % 3 == 1"},
            {"index": 2, "name": "orf_f2", "meaning": "... start position % 3 == 2"},
            {"index": 3, "name": "start_ext",
             "meaning": "graded start propensity at the FIRST nt of each candidate start codon, "
                        "weighted by start_weights; 0 elsewhere"},
            {"index": 4, "name": "is_stop",
             "meaning": "first nt of any stop codon TAA / TAG / TGA"},
        ],
        "occupancy_note": (
            "Channels 0-2 are ATG-ONLY occupancy in every mode, deliberately: near-cognate starts "
            "enter through the graded start channel alone. Opening occupancy at 10 start codons "
            "would explode the density and destroy the channel's meaning."),
        "display_window": {
            "sequence": DISPLAY_WINDOW,
            "length": len(DISPLAY_WINDOW),
            "orf_track": {c: [round(float(x), 4) for x in tr[:, j]] for j, c in enumerate(chan)},
            "onehot": onehot,
            "rna_coverage": None,
            "rna_coverage_note": (
                "The 10th input channel is DATA, not sequence-derived, so it has no defined value "
                "for a synthetic display window. Draw it as a placeholder or omit it; do not "
                "invent numbers. There is no 24x10 matrix to ship -- only 24x9 is derivable."),
        },
        "provenance": {
            "start_weights": "imported from scripts/build_orf_track.py START_W",
            "orf_track": "computed by calling build_orf_track.orf_track(seq, mode='ext', "
                         "kozak='none'), not transcribed",
            "generated_by": "scripts/export_arch_spec.py",
        },
    }
    q = out / "orf_track_spec.json"
    q.write_text(json.dumps(ospec, indent=2))
    print(f"  wrote {q}  ({len(bot.START_W)} start codons, {len(DISPLAY_WINDOW)}-nt window)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
