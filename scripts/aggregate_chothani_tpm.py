#!/usr/bin/env python3
"""Pool the per-run alignment-mode salmon quants into per-tissue TPM, and measure what the missing
decoys cost by comparing Fibroblast against the surviving decoy-aware table.

WHY FIBROBLAST IS THE CONTROL. Its decoy-aware mean TPM (`data/tpm/fibroblast_salmon_mean_tpm.tsv`)
is the ONLY Chothani TPM that survived the 2026-08-15 loss. Quantifying the same 32 libraries in
decoy-free alignment mode and comparing gives a direct read on whether the seven RECOVERED tissues
can be trusted for isoform ranking, or whether a 30-50 hour refetch is actually required.

  aggregate_chothani_tpm.py [--outdir ...]
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

NEW = Path(__file__).resolve().parents[1]
ALN = NEW / "data/tpm/chothani_alnmode"
DECOY_FIB = NEW / "data/tpm/fibroblast_salmon_mean_tpm.tsv"


def read_quant(p):
    """tx_id -> TPM, normalising the two DIFFERENT id conventions in play.

    The decoy-aware GENCODE index names targets with the FULL pipe-delimited FASTA header
    (`ENST...|ENSG...|OTTHUMG...|NAME-201|NAME|len|UTR5:..|CDS:..|`), while the alignment-mode
    run used RiboCode's `transcripts_sequence.fa`, whose names are bare versioned ENST. Keying on
    the bare id against a decoy quant silently yields ZERO for every transcript -- it produced a
    GTF2I total of 0.0 TPM before this was caught. Split on '|' and take field 0 either way."""
    with open(p) as fh:
        return {r["Name"].split("|")[0]: float(r["TPM"])
                for r in csv.DictReader(fh, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=ALN)
    a = ap.parse_args()

    per_tissue = defaultdict(list)
    for q in sorted(ALN.rglob("quant.sf")):
        per_tissue[q.parent.parent.name].append(q)
    if not per_tissue:
        raise SystemExit("no quant.sf under " + str(ALN))

    means, summary = {}, {}
    for tis, qs in sorted(per_tissue.items()):
        acc = defaultdict(float)
        for q in qs:
            for k, v in read_quant(q).items():
                acc[k] += v
        m = {k: v / len(qs) for k, v in acc.items()}
        means[tis] = m
        summary[tis] = {"n_runs": len(qs), "n_tx": len(m),
                        "n_tpm_ge1": sum(1 for v in m.values() if v >= 1)}
        out = a.outdir / f"{tis}_mean_tpm.tsv"
        with open(out, "w") as fh:
            fh.write("tx_id\tmean_tpm\tn_runs\n")
            for k in sorted(m):
                fh.write(f"{k}\t{m[k]:.6f}\t{len(qs)}\n")
        print(f"  {tis:<12} {len(qs):>3} runs -> {out.name}")

    # ---- the control: decoy-aware vs decoy-free on the SAME 32 Fibroblast libraries -------------
    cmp_out = None
    if "Fibroblast" in means and DECOY_FIB.exists():
        with open(DECOY_FIB) as fh:
            decoy = {r["tx_id"]: float(r["mean_tpm"]) for r in csv.DictReader(fh, delimiter="\t")}
        aln = means["Fibroblast"]
        shared = sorted(set(decoy) & set(aln))
        d = np.array([decoy[t] for t in shared])
        v = np.array([aln[t] for t in shared])
        both = (d > 0) | (v > 0)
        rho = spearmanr(d[both], v[both]).statistic
        ld, lv = np.log10(d[both] + 0.01), np.log10(v[both] + 0.01)
        cmp_out = {
            "n_shared_tx": len(shared),
            "n_either_expressed": int(both.sum()),
            "spearman_all": round(float(rho), 4),
            "spearman_tpm_ge1_in_decoy": round(float(
                spearmanr(d[d >= 1], v[d >= 1]).statistic), 4) if (d >= 1).sum() > 10 else None,
            "pearson_log10": round(float(np.corrcoef(ld, lv)[0, 1]), 4),
            "median_log2_ratio": round(float(np.median(
                np.log2(v[(d >= 1) & (v > 0)] / d[(d >= 1) & (v > 0)]))), 4),
            "n_decoy_ge1": int((d >= 1).sum()), "n_aln_ge1": int((v >= 1).sum()),
            "note": "decoy-aware (surviving) vs decoy-free alignment mode, SAME 32 libraries",
        }
        print("\n  Fibroblast decoy-aware vs decoy-free:")
        for k, val in cmp_out.items():
            print(f"    {k:<28} {val}")

    (a.outdir / "aggregate_summary.json").write_text(json.dumps(
        {"per_tissue": summary, "fibroblast_decoy_vs_alnmode": cmp_out,
         "method_caveat": "alignment-based salmon against transcriptome only; the ORIGINAL "
                          "chothani quant was decoy-aware selective alignment from FASTQ. Not a "
                          "drop-in replacement -- see the Fibroblast comparison for the cost."},
        indent=2))
    print(f"\n  wrote {a.outdir}/aggregate_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
