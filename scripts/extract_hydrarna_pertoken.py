#!/usr/bin/env python3
"""HydraRNA per-token embeddings for the riboseq Fibroblast universe -- VARIABLE length: one
<tx_id>_tokens.npy of shape (L, 1024) float16 per transcript (matches the rinalmo/orthrus stores the
Ribo-seq model reads via a tx_index). Universe is pre-capped at <=9996 nt (within HydraRNA's ~10K
limit), so no truncation. Resume-safe (skips a tx whose .npy exists). Processes a [--start,--end)
transcript range for array chunking. GPU-only (Mamba selective-scan). Runs inside the HydraRNA SIF.

Adapted from multirm_fm_probe/extract_pertoken_hydrarna.py (fixed memmap -> per-tx .npy).
"""
import argparse
import os
import sys

import numpy as np
import torch
from fairseq import checkpoint_utils, data, options, tasks

WEIGHTS = "/opt/hydrarna_weights/HydraRNA_model.pt"
DICT = "/opt/hydrarna_weights/dict"
EMB_DIM = 1024


def parse_fasta(path):
    h, p = None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if h is not None:
                yield h, "".join(p)
            h = line[1:].split()[0]
            p = []
        elif line:
            p.append(line.strip())
    if h is not None:
        yield h, "".join(p)


def load_model():
    parser = options.get_generation_parser(default_task="masked_lm_span")
    a = options.parse_args_and_arch(parser, [DICT])
    task = tasks.setup_task(a)
    print(f"| loading HydraRNA from {WEIGHTS}", file=sys.stderr)
    models, _ = checkpoint_utils.load_model_ensemble([WEIGHTS], task=task)
    return models[0], task


def embed_pertoken(model, task, seq, device):
    tokenized = "<s> " + " ".join(list(seq))
    tokens = task.source_dictionary.encode_line(tokenized, add_if_not_exist=False)
    batch = data.monolingual_dataset.collate(
        samples=[{"id": -1, "source": tokens, "target": tokens}],
        pad_idx=task.source_dictionary.pad(),
        eos_idx=task.source_dictionary.eos(),
    )
    src = batch["net_input"]["src_tokens"].to(device)
    with torch.no_grad():
        features = model.encoder.extract_features(src_tokens=src)
    return features[0][0, 1:-1, :].float().cpu().numpy()  # (L, 1024), BOS+EOS stripped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="universe FASTA (variable length)")
    ap.add_argument("-o", "--output", required=True, help="output dir for <tx_id>_tokens.npy")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=-1)
    ap.add_argument("--log-every", type=int, default=100)
    a = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        sys.exit("ERROR: HydraRNA requires CUDA (Mamba selective-scan)")
    os.makedirs(a.output, exist_ok=True)
    model, task = load_model()
    model = model.to(device).half().eval()

    seqs = list(parse_fasta(a.input))
    end = len(seqs) if a.end < 0 else min(a.end, len(seqs))
    print(f"range [{a.start},{end}) of {len(seqs)} tx -> {a.output}", file=sys.stderr)
    done = skip = bad = 0
    for i in range(a.start, end):
        txid, seq = seqs[i]
        out = os.path.join(a.output, f"{txid}_tokens.npy")
        if os.path.exists(out):
            skip += 1
            continue
        seq = seq.upper().replace("T", "U")
        emb = embed_pertoken(model, task, seq, device)          # (L, 1024)
        if emb.shape[0] != len(seq) or emb.shape[1] != EMB_DIM:
            print(f"  SHAPE-MISMATCH {txid}: emb {emb.shape} vs seq {len(seq)}", file=sys.stderr)
            bad += 1
            continue
        tmp = out + ".tmp"
        # file handle so np.save keeps the exact name (it re-appends .npy to a str path)
        with open(tmp, "wb") as fh:
            np.save(fh, emb.astype(np.float16))
        os.replace(tmp, out)
        done += 1
        if done % a.log_every == 0:
            print(f"  {i + 1}/{end} written={done} skip={skip}", file=sys.stderr)
            sys.stderr.flush()
    print(f"chunk [{a.start},{end}): written={done} skip={skip} bad={bad}", file=sys.stderr)


if __name__ == "__main__":
    main()
