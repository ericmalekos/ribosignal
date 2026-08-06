#!/usr/bin/env python3
"""Orthrus 4-track PER-TOKEN embedding extraction for the Fibroblast universe.

The averaged extractor (Embeddings/orthrus6t/scripts/extract_orthrus4t.py) calls
model.representation(...) which mean-pools over the length axis and returns (1, 512).
The per-nt Ribo-seq signal model needs the PRE-pool per-token tensor instead:
model.forward(x, channel_last=True) returns (1, L, 512) with one row per nucleotide and
NO CLS/EOS special tokens (unlike RiNALMo, which prepends CLS + appends EOS and must be
sliced 1:L+1). Orthrus one-hot input therefore aligns 1:1 to the transcript nt axis, the
P-site target, the RNAseq coverage, and the RiNALMo per-token embeddings.

Mirrors extract_lncrna_rinalmo.py's operational contract:
  - writes one <tx_id>_tokens.npy of shape (L, 512) per transcript
  - skip-if-exists so a resumed job recovers from a slurm timeout without recompute
  - length-sorted ascending iteration so any OOM only aborts the long tail
  - per-record OOM guard -> skipped.tsv, continue
  - --dtype float16 storage (compute is float32); model output is float32

Orthrus 4-track is sequence-only (A/C/G/T one-hot, U->T), so it needs no annotation and
is the fair FM-vs-FM counterpart to RiNALMo. The 6-track variant additionally requires
per-nt CDS + 5' splice indicator channels built from a GTF, which for Ribo-seq would leak
the very CDS location the model should learn, so it is intentionally not used here.

Requires CUDA (mamba-ssm selective-scan kernel is CUDA-only).
"""
import argparse
import os
import sys
import time
from collections.abc import Iterator

import numpy as np
import torch


def iter_fasta(path: str) -> Iterator[tuple[str, str]]:
    header = None
    chunks: list = []
    with open(path) as fh:
        for line in fh:
            if not line:
                continue
            if line[0] == ">":
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)


def load_sequences(fasta: str, max_length: int) -> list:
    items = []
    n_over = 0
    for header, seq in iter_fasta(fasta):
        tx_id = header.split()[0]
        seq = seq.upper()
        if not seq:
            continue
        if len(seq) > max_length:
            n_over += 1
            continue
        items.append((tx_id, seq))
    print(f"  loaded {len(items):,} sequences (dropped {n_over:,} over "
          f"max_length={max_length})", file=sys.stderr)
    return items


def seq_to_oh(seq: str) -> np.ndarray:
    """(L, 4) float32 one-hot, order [A, C, G, T]; U->T; unknown/N -> all-zero row."""
    L = len(seq)
    out = np.zeros((L, 4), dtype=np.float32)
    idx = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3}
    for i, b in enumerate(seq):
        j = idx.get(b)
        if j is not None:
            out[i, j] = 1
    return out


def safe_path(out_dir: str, tx_id: str) -> str:
    safe = tx_id.replace("/", "_").replace(" ", "_")
    return os.path.join(out_dir, f"{safe}_tokens.npy")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output", required=True, help="output directory")
    ap.add_argument("--weights", default="/weights/orthrus-large-4-track")
    ap.add_argument("--max-length", type=int, required=True,
                    help="hard upper-length cap; longer records are skipped")
    ap.add_argument("--checkpoint-every", type=int, default=500,
                    help="progress-print cadence (per-token files are written as they go)")
    ap.add_argument("--dtype", choices=["float32", "float16"], default="float16",
                    help="per-token output dtype")
    ap.add_argument("--limit", type=int, default=None,
                    help="optional: process only the first N records after sorting (smoke)")
    args = ap.parse_args()

    os.makedirs(args.output, exist_ok=True)
    skipped_path = os.path.join(args.output, "skipped.tsv")
    new_skipped = not os.path.exists(skipped_path)
    skipped_fh = open(skipped_path, "a")
    if new_skipped:
        skipped_fh.write("tx_id\tlength\treason\n")

    if not torch.cuda.is_available():
        print("ERROR: CUDA not available (mamba-ssm requires CUDA)", file=sys.stderr)
        return 1
    dev = torch.device("cuda")
    name = torch.cuda.get_device_name(dev)
    total_gb = torch.cuda.get_device_properties(dev).total_memory / 1e9
    print(f"GPU: {name} ({total_gb:.1f} GB total)", file=sys.stderr)

    print(f"Loading Orthrus 4-track from {args.weights}...", file=sys.stderr)
    t0 = time.time()
    from transformers import AutoModel
    model = AutoModel.from_pretrained(args.weights, trust_remote_code=True).to(dev).eval()
    print(f"  model load took {time.time()-t0:.1f}s", file=sys.stderr)

    items = load_sequences(args.input, args.max_length)
    items.sort(key=lambda kv: len(kv[1]))
    if args.limit is not None:
        items = items[:args.limit]
        print(f"  --limit applied; processing {len(items):,} records", file=sys.stderr)

    np_dtype = np.float16 if args.dtype == "float16" else np.float32
    n_done = 0
    n_skip_existing = 0
    n_oom = 0
    t_start = time.time()

    for tx_id, seq in items:
        per_token_path = safe_path(args.output, tx_id)
        if os.path.exists(per_token_path):
            n_skip_existing += 1
            continue

        L = len(seq)
        torch.cuda.empty_cache()
        try:
            oh = seq_to_oh(seq)
            x = torch.from_numpy(oh).to(dev).unsqueeze(0)         # (1, L, 4)
            with torch.no_grad():
                h = model.forward(x, channel_last=True)           # (1, L, 512) per-token
            emb = h.squeeze(0).float().cpu().numpy()              # (L, 512), no CLS/EOS
            del x, h
        except torch.cuda.OutOfMemoryError:
            skipped_fh.write(f"{tx_id}\t{L}\tOOM\n")
            skipped_fh.flush()
            torch.cuda.empty_cache()
            n_oom += 1
            print(f"  OOM @ tx={tx_id} L={L}; skipping", file=sys.stderr)
            continue
        except Exception as e:
            skipped_fh.write(f"{tx_id}\t{L}\tERR:{type(e).__name__}\n")
            skipped_fh.flush()
            torch.cuda.empty_cache()
            print(f"  ERROR {type(e).__name__} @ tx={tx_id} L={L}: {e}", file=sys.stderr)
            continue

        if emb.shape[0] != L:
            skipped_fh.write(f"{tx_id}\t{L}\tLEN_MISMATCH:{emb.shape[0]}\n")
            skipped_fh.flush()
            print(f"  LEN MISMATCH @ tx={tx_id}: emb L={emb.shape[0]} != {L}",
                  file=sys.stderr)
            continue

        tmp_path = per_token_path[:-4] + ".partial.npy"
        np.save(tmp_path, emb.astype(np_dtype, copy=False))
        os.replace(tmp_path, per_token_path)
        n_done += 1

        if n_done % args.checkpoint_every == 0:
            elapsed = time.time() - t_start
            rate = n_done / elapsed if elapsed > 0 else 0
            print(f"  [{n_done:,} new / {len(items):,}] last_L={L} "
                  f"rate={rate:.2f} seq/s elapsed={elapsed/60:.1f}m", file=sys.stderr)

    skipped_fh.close()
    elapsed = time.time() - t_start
    print(f"\nDONE: wrote {n_done:,} new, skipped {n_skip_existing:,} pre-existing, "
          f"{n_oom} OOM in {elapsed/60:.1f}m", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
