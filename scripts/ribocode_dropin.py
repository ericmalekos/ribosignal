#!/usr/bin/env python3
"""RiboCode drop-in test (Task 14): run RiboCode's ORF caller (detectORF.main) on a SUBSTITUTED
per-nt P-site density -- the model's PREDICTED profile instead of experimental Ribo-seq -- and see
how well the resulting ORF calls track the calls RiboCode makes on the real data.

Runs in the `ribocode` conda env (RiboCode v1.2.15). Bypasses BAM reading entirely: RiboCode's
density object `tpsites_sum` is just a {transcript_id: np.ndarray} dict (0-based, length == tx
length), so we build it from dump_pred_profiles.py's npz and pass it straight to detectORF.main with
the SAME prebuilt annotation DB and default parameters the official per-tissue run used (start ATG,
min_AA 5, pval 0.05, longest_orf, stouffer none, fdr_bh).

Three density variants (identical caller + params; ONLY the density array differs):
  real            real pooled P-site counts (reference; reproduces official <Tissue>_collapsed.txt)
  pred_obsdepth   predicted SHAPE * real per-tx total, rounded to int  (isolates profile shape)
  pred_preddepth  predicted SHAPE * count-head predicted total, rounded  (standalone, no real data)

Rounding to int is deliberate: a smooth density has almost no exact zeros, so it would pass
RiboCode's per-codon coverage/tie logic more easily than real integer counts and look artificially
more periodic. Scaling to realistic depth is also required -- RiboCode hard-gates per-tx sum>=5,
per-ORF frame0 sum>=5, and >=5 nonzero frame0 codons (detectORF.py:306).

Usage (one variant per invocation, so they parallelize as an sbatch array):
  ribocode_dropin.py --profiles <pred_profiles.hd5> --annot <ribocode_annot dir> \
      --variant {real,pred_obsdepth,pred_preddepth} --out <out_dir> [--min_aa 5] [--pval 0.05]
Writes <out_dir>/<variant>.txt and <out_dir>/<variant>_collapsed.txt (RiboCode's native output).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from RiboCode import detectORF
from RiboCode.prepare_transcripts import load_transcripts_pickle


def build_density(variant, prof, obs, ptot, floor_mult=0.0, pred_scale=1.0, rng=None, poisson=False):
    """Per-nt density array for one transcript under the chosen variant, as float32 (rounded to int
    values for the pred variants so zeros/ties mirror real integer counts).

    floor_mult > 0 applies a PER-POSITION density floor to the predicted profile before scaling: a
    position whose predicted probability is below floor_mult / L (floor_mult x the uniform
    background 1/L) is zeroed. This removes the softmax's diffuse low-probability leak into the
    3'UTR (the source of the spurious dORF over-calls) while leaving the concentrated CDS/uORF peaks
    intact. The profile is NOT renormalized after flooring, so surviving (CDS) counts are preserved
    and the zeroed positions drop out. Pred variants only; real counts are untouched."""
    if variant == "real":
        return obs.astype(np.float32)
    p = prof.astype(np.float64)
    if floor_mult > 0 and p.size:
        p = np.where(p >= floor_mult / p.size, p, 0.0)
    if variant == "pred_obsdepth":
        depth = float(obs.sum())
    elif variant == "pred_preddepth":
        depth = float(ptot) * pred_scale   # effective-depth knob for the standalone variant (Check 2)
    else:
        raise ValueError(variant)
    lam = np.maximum(p * depth, 0.0)
    if poisson and rng is not None:
        return rng.poisson(lam).astype(np.float32)   # noise injection (Check 4b)
    return np.rint(lam).astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profiles", required=True, help="pred_profiles.npz from dump_pred_profiles")
    ap.add_argument("--annot", required=True, help="RiboCode prepare_transcripts annotation dir")
    ap.add_argument("--variant", required=True,
                    choices=["real", "pred_obsdepth", "pred_preddepth"])
    ap.add_argument("--out", required=True, help="output dir; writes <variant>[_collapsed].txt")
    ap.add_argument("--min_aa", type=int, default=5)
    ap.add_argument("--pval", type=float, default=0.05)
    ap.add_argument("--floor_mult", type=float, default=0.0,
                    help="per-position density floor: zero predicted positions below floor_mult/L "
                         "(floor_mult x uniform). Removes the 3'UTR softmax leak (spurious dORFs). "
                         "0 = off. Output tagged <variant>_floor<mult> when > 0.")
    ap.add_argument("--subsample", type=float, default=1.0,
                    help="binomially thin the REAL per-nt P-site counts to this fraction of depth "
                         "(models a shallower Ribo-seq experiment; each footprint kept i.i.d. w.p. f). "
                         "real variant only; 1.0 = off. Output tagged <variant>_depth<f>_seed<s>.")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed for --subsample / --pred_poisson.")
    ap.add_argument("--pred_scale", type=float, default=1.0,
                    help="scale the PREDICTED effective depth for pred_preddepth (the standalone "
                         "count-calibration / CDS-anchored depth knob, Check 2). 1.0 = model's own "
                         "predicted total; <1 stricter (fewer calls), >1 more permissive.")
    ap.add_argument("--pred_poisson", action="store_true",
                    help="sample the predicted density as Poisson(rate) instead of rounding "
                         "(injects realistic detection noise into the smooth profile; Check 4b).")
    ap.add_argument("--alt_start_codons", default="",
                    help="comma-separated alternative start codons (e.g. CTG) enabling non-AUG ORF "
                         "calling via RiboCode's ALTERNATIVE_START_CODON_LIST (== the RiboCode -A "
                         "flag). RiboCode uses these as a FALLBACK per common-stop: an alt-start ORF "
                         "is called only where that in-frame stop has no ATG, so alt-start calls are "
                         "DISJOINT from the ATG calls and do NOT reclassify canonical ORFs. Default "
                         "empty = ATG only, matching the official per-tissue runs. Filenames are NOT "
                         "tagged, so write alt-start runs to a SEPARATE --out dir (e.g. dropin_ctg).")
    args = ap.parse_args()
    alt_list = [c.strip().upper() for c in args.alt_start_codons.split(",") if c.strip()] or None

    annot = Path(args.annot)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print(f"loading annotation pickle {annot/'transcripts.pickle'} ...", file=sys.stderr)
    gene_dict, transcript_dict = load_transcripts_pickle(str(annot / "transcripts.pickle"))
    print(f"  {len(transcript_dict):,} transcripts, {len(gene_dict):,} genes", file=sys.stderr)

    # RiboCode iterates EVERY transcript_dict key and does tpsites_sum[tid].sum(); init all to zeros
    # (length-matched) so unpredicted transcripts are harmlessly skipped by the sum>=5 gate.
    tpsites_sum = {tid: np.zeros(t.length, dtype=np.float32)
                   for tid, t in transcript_dict.items()}

    npz = np.load(args.profiles, allow_pickle=False)
    tx_ids = npz["tx_ids"]
    lengths = npz["lengths"]
    pred_flat = npz["pred_flat"]
    obs_flat = npz["obs_flat"]
    pred_total = npz["pred_total"]
    off = np.concatenate([[0], np.cumsum(lengths)])
    n_ok = n_missing = n_lenmm = 0
    total_depth = 0.0
    rng = np.random.default_rng(args.seed)
    subsample = args.subsample if args.variant == "real" else 1.0  # thinning is real-variant only
    for j, tx in enumerate(tx_ids):
        t = transcript_dict.get(str(tx))
        if t is None:
            n_missing += 1
            continue
        prof = pred_flat[off[j]:off[j + 1]]
        if len(prof) != t.length:
            n_lenmm += 1
            continue
        obs = obs_flat[off[j]:off[j + 1]]
        if subsample < 1.0:
            obs = rng.binomial(obs.astype(np.int64), subsample)
        dens = build_density(args.variant, prof, obs, pred_total[j], args.floor_mult,
                             pred_scale=args.pred_scale, rng=rng, poisson=args.pred_poisson)
        tpsites_sum[str(tx)] = dens
        total_depth += float(dens.sum())
        n_ok += 1
    print(f"variant={args.variant} floor_mult={args.floor_mult} "
          f"alt_start_codons={alt_list}: {n_ok:,} tx injected, "
          f"{n_missing:,} not in annot DB, {n_lenmm:,} length-mismatch, "
          f"total P-sites={total_depth:,.0f}", file=sys.stderr)

    tag = args.variant if args.floor_mult <= 0 else f"{args.variant}_floor{args.floor_mult:g}"
    if subsample < 1.0:
        tag = f"{tag}_depth{subsample:g}_seed{args.seed}"
    outname = str(out / tag)
    detectORF.main(
        gene_dict=gene_dict, transcript_dict=transcript_dict, annot_dir=str(annot),
        tpsites_sum=tpsites_sum, total_psites_number=int(total_depth),
        pval_cutoff=args.pval, only_longest_orf=True,
        START_CODON=["ATG"], ALTERNATIVE_START_CODON_LIST=alt_list,
        STOP_CODON_LIST=["TAA", "TAG", "TGA"], MIN_AA_LENGTH=args.min_aa,
        outname=outname, output_gtf=False, output_bed=False,
        dependence_test=False, stouffer_adj="none", pval_adj="fdr_bh")
    print(f"wrote {outname}_collapsed.txt", file=sys.stderr)


if __name__ == "__main__":
    main()
