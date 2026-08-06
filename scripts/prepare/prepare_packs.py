#!/usr/bin/env python3
"""Batch pack builder: one pack per group from a samplesheet (nf-core style).

Samplesheet is a CSV/TSV with a header and columns:  group,assay,path
  group : pack/group name (rows with the same group are pooled into one pack)
  assay : "ribo" or "rna"
  path  : a STAR toTranscriptome BAM, OR a pre-made *_psites.hd5 / *_coverage.hd5

For each group, this routes to prepare_from_bams.py (if any path is a .bam) or prepare_pack.py
(if all inputs are hd5), forwarding the shared universe/reference/tool flags, and writes the
pack to <outdir>/<group>. All paths are CLI-driven; nothing is hardcoded.

Example
  prepare_packs.py --samplesheet samples.csv \
    --ref-pack data/packed --gtf gencode.v49.annotation.gtf --species human \
    --outdir packs/ --annot /path/ribocode_annot --ncrna-tx nc.txt --tx-to-gene t2g.tsv \
    --ribocode-bin $RIBOCODE_BIN --pysam-python /path/riboseq/bin/python
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_samplesheet(path):
    """-> {group: {"ribo": [paths], "rna": [paths]}}. Auto-detects CSV vs TSV by the header line."""
    text = Path(path).read_text().splitlines()
    if not text:
        raise SystemExit(f"empty samplesheet: {path}")
    delim = "\t" if "\t" in text[0] else ","
    groups = defaultdict(lambda: {"ribo": [], "rna": []})
    for row in csv.DictReader(text, delimiter=delim):
        g = row["group"].strip()
        assay = row["assay"].strip().lower()
        p = row["path"].strip()
        if assay not in ("ribo", "rna"):
            raise SystemExit(f"assay must be ribo|rna, got '{assay}' for group {g}")
        groups[g][assay].append(p)
    return groups


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--outdir", required=True, help="one pack written per group under here")
    uni = ap.add_mutually_exclusive_group(required=True)
    uni.add_argument("--universe-tx")
    uni.add_argument("--ref-pack")
    ap.add_argument("--species", default="human")
    ap.add_argument("--fasta")
    ap.add_argument("--build-orf-track", action="store_true")
    ap.add_argument("--kozak", default="none", choices=["none", "heuristic", "pwm"])
    ap.add_argument("--tx2biotype")
    ap.add_argument("--gtf")
    ap.add_argument("--coverage-only", action="store_true")
    # ribo-processing refs/tools (forwarded to prepare_from_bams when inputs are BAMs)
    ap.add_argument("--annot")
    ap.add_argument("--ncrna-tx")
    ap.add_argument("--tx-to-gene")
    ap.add_argument("--ribocode-bin")
    ap.add_argument("--pysam-python")
    ap.add_argument("--samtools")
    ap.add_argument("--dry-run", action="store_true")
    return ap.parse_args()


def shared_pack_flags(args, out):
    f = ["--group", None, "--species", args.species, "--out", out]  # group filled per call
    f += (["--ref-pack", args.ref_pack] if args.ref_pack else ["--universe-tx", args.universe_tx])
    if args.build_orf_track:
        f += ["--build-orf-track", "--fasta", args.fasta, "--kozak", args.kozak]
    if args.tx2biotype:
        f += ["--tx2biotype", args.tx2biotype]
    elif args.gtf:
        f += ["--gtf", args.gtf]
    return f


def main():
    args = parse_args()
    groups = read_samplesheet(args.samplesheet)
    outdir = Path(args.outdir)
    print(f"{len(groups)} group(s): {', '.join(groups)}", file=sys.stderr)

    for g, streams in groups.items():
        out = outdir / g
        ribo, rna = streams["ribo"], streams["rna"]
        has_bam = any(p.endswith(".bam") for p in ribo + rna)
        flags = shared_pack_flags(args, str(out))
        flags[1] = g  # fill group

        if has_bam:
            cmd = [sys.executable, str(HERE / "prepare_from_bams.py")]
            ribo_bam = [p for p in ribo if p.endswith(".bam")]
            ribo_hd5 = [p for p in ribo if not p.endswith(".bam")]
            rna_bam = [p for p in rna if p.endswith(".bam")]
            rna_hd5 = [p for p in rna if not p.endswith(".bam")]
            if ribo_bam:
                cmd += ["--ribo-bam", *ribo_bam]
            if ribo_hd5:
                cmd += ["--ribo-psites", *ribo_hd5]
            if rna_bam:
                cmd += ["--rna-bam", *rna_bam]
            if rna_hd5:
                cmd += ["--rna-coverage", *rna_hd5]
            fwd = (("--annot", args.annot), ("--ncrna-tx", args.ncrna_tx),
                   ("--tx-to-gene", args.tx_to_gene), ("--ribocode-bin", args.ribocode_bin),
                   ("--pysam-python", args.pysam_python), ("--samtools", args.samtools))
            for opt, val in fwd:
                if val:
                    cmd += [opt, val]
            if args.coverage_only:
                cmd.append("--coverage-only")
            if args.dry_run:
                cmd.append("--dry-run")
        else:
            cmd = [sys.executable, str(HERE / "prepare_pack.py"), "--rna-coverage", *rna]
            if args.coverage_only:
                cmd.append("--coverage-only")
            elif ribo:
                cmd += ["--ribo-psites", *ribo]
        cmd += flags

        print(f"\n=== group {g}: ribo={len(ribo)} rna={len(rna)} "
              f"({'BAM' if has_bam else 'hd5'}) -> {out} ===", file=sys.stderr)
        print("+ " + " ".join(cmd), file=sys.stderr)
        # hd5-only path (prepare_pack) has no --dry-run flag, so skip it entirely on dry-run;
        # the BAM path (prepare_from_bams) already got --dry-run appended and prints its own plan.
        if not (args.dry_run and not has_bam):
            subprocess.run(cmd, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
