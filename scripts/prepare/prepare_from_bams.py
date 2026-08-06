#!/usr/bin/env python3
"""Build a model-input pack straight from STAR **transcriptome** BAM paths (the BYO-BAM entry).

This is the "a future user passes paths to their STAR alignments on the CLI" front-end. It
orchestrates the validated per-BAM steps, then calls prepare_pack.py:

  Ribo-seq BAMs -> [filter_tx_heldout.py: drop ncRNA + cross-gene] -> samtools index
                -> RiboCode metaplots (auto per-read-length P-site offsets) -> RiboCode -g
                -> per-nt *_psites.hd5
  RNA-seq  BAMs -> rnaseq_coverage.py (auto strand/libtype, PE/SE) -> *_coverage.hd5
  -> prepare_pack.py (pool onto --universe-tx | --ref-pack, write the 7-file pack + optional
     no-Kozak ORF track)

No hardcoded project paths, tissue names, or sample counts. External tools are resolved from
CLI flags or environment ($RIBOCODE_BIN, $RIBO_PYSAM_PYTHON, $RIBO_SAMTOOLS), never baked in.
The user's input BAMs are copied into --workdir before the in-place ncRNA filter, so inputs are
never mutated. Requirement: the BAM @SQ reference names must be versioned transcript ids from
the SAME GENCODE annotation as the universe (prepare_pack fails loudly otherwise).

It does NOT align (bring transcriptome BAMs). Pass --ribo-psites / --rna-coverage instead of
--ribo-bam / --rna-bam to skip a stage whose hd5 you already have. --coverage-only (no Ribo)
gives a predict-only pack.

Example
  prepare_from_bams.py \
    --ribo-bam sampleA.toTranscriptome.bam sampleB.toTranscriptome.bam \
    --rna-bam  rnaA.toTranscriptome.bam \
    --annot /path/ribocode_annot --ncrna-tx ncrna_tx.txt --tx-to-gene tx2gene.tsv \
    --universe-tx universe_tx.txt --fasta universe.fa --build-orf-track \
    --group mycellline --species human --out packs/mycellline \
    --ribocode-bin $RIBOCODE_BIN --pysam-python /path/riboseq/bin/python
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import project_root  # noqa: E402


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    # inputs (BAM to run the stage, OR hd5 to skip it)
    ap.add_argument("--ribo-bam", nargs="*", default=[], help="STAR toTranscriptome Ribo BAMs")
    ap.add_argument("--rna-bam", nargs="*", default=[], help="STAR toTranscriptome RNA BAMs")
    ap.add_argument("--ribo-psites", nargs="*", default=[], help="pre-made *_psites.hd5 (skip)")
    ap.add_argument("--rna-coverage", nargs="*", default=[], help="pre-made *_coverage.hd5 (skip)")
    ap.add_argument("--coverage-only", action="store_true", help="predict-only (no Ribo target)")
    # universe / output (forwarded to prepare_pack)
    uni = ap.add_mutually_exclusive_group(required=True)
    uni.add_argument("--universe-tx")
    uni.add_argument("--ref-pack")
    ap.add_argument("--group", required=True)
    ap.add_argument("--species", default="human")
    ap.add_argument("--out", required=True)
    ap.add_argument("--fasta", help="universe FASTA (for --build-orf-track)")
    ap.add_argument("--build-orf-track", action="store_true")
    ap.add_argument("--kozak", default="none", choices=["none", "heuristic", "pwm"])
    ap.add_argument("--tx2biotype")
    ap.add_argument("--gtf")
    # ribo-processing refs
    ap.add_argument("--annot", help="RiboCode prepared annotation dir (required with --ribo-bam)")
    ap.add_argument("--ncrna-tx", help="ncRNA tx blocklist (filter_tx_heldout.py); optional")
    ap.add_argument("--tx-to-gene", help="tx->gene TSV for cross-gene drop; optional")
    # external tools (env-defaulted, never hardcoded)
    ap.add_argument("--ribocode-bin", default=os.environ.get("RIBOCODE_BIN"),
                    help="dir with metaplots + RiboCode (or $RIBOCODE_BIN)")
    ap.add_argument("--pysam-python", default=os.environ.get("RIBO_PYSAM_PYTHON", sys.executable),
                    help="python with pysam for filter_tx + coverage (or $RIBO_PYSAM_PYTHON)")
    ap.add_argument("--samtools", default=os.environ.get("RIBO_SAMTOOLS", "samtools"))
    ap.add_argument("--workdir", help="scratch for intermediates (default <out>/_work)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--dry-run", action="store_true", help="print the plan; run nothing")
    return ap.parse_args()


def main():
    args = parse_args()
    root = project_root()
    scripts = root / "scripts"
    out_dir = Path(args.out)
    work = Path(args.workdir) if args.workdir else out_dir / "_work"
    psites_dir = work / "psites"
    cov_dir = work / "coverage"
    dry = args.dry_run
    cmds = []

    def run(cmd, cwd=None):
        cmd = [str(c) for c in cmd]
        print("+ " + " ".join(cmd) + (f"   (cwd={cwd})" if cwd else ""), file=sys.stderr)
        cmds.append(cmd)
        if not dry:
            subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)

    # ---- validation ----------------------------------------------------------
    if args.ribo_bam and not args.annot:
        print("ERROR: --annot (RiboCode annotation) required with --ribo-bam", file=sys.stderr)
        return 2
    if args.ribo_bam and not args.ribocode_bin:
        print("ERROR: --ribocode-bin or $RIBOCODE_BIN required with --ribo-bam", file=sys.stderr)
        return 2
    if not (args.rna_bam or args.rna_coverage):
        print("ERROR: provide --rna-bam or --rna-coverage", file=sys.stderr)
        return 2
    if not dry:
        for d in (psites_dir, cov_dir):
            d.mkdir(parents=True, exist_ok=True)

    # ---- Ribo side: BAMs -> psites hd5 --------------------------------------
    ribo_psites = list(args.ribo_psites)
    if args.ribo_bam:
        bam_list = psites_dir / f"{args.group}_bam_list.txt"
        wbams = []
        for bam in args.ribo_bam:
            wbam = work / "ribo_bam" / Path(bam).name
            if not dry:
                wbam.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(bam, wbam)
            else:
                print(f"+ cp {bam} {wbam}", file=sys.stderr)
            if args.ncrna_tx and args.tx_to_gene:
                run([args.pysam_python, scripts / "heldout" / "filter_tx_heldout.py",
                     wbam, args.ncrna_tx, args.tx_to_gene])
            run([args.samtools, "index", wbam])
            wbams.append(wbam)
        if not dry:
            bam_list.write_text("\n".join(str(b) for b in wbams) + "\n")
        else:
            print(f"+ write {bam_list} ({len(wbams)} BAMs)", file=sys.stderr)
        run([Path(args.ribocode_bin) / "metaplots", "-a", args.annot,
             "-i", bam_list, "-o", args.group], cwd=psites_dir)
        config = psites_dir / f"{args.group}_pre_config.txt"
        run([Path(args.ribocode_bin) / "RiboCode", "-a", args.annot,
             "-c", config, "-l", "no", "-g", "-o", args.group], cwd=psites_dir)
        ribo_psites = [str(psites_dir)]

    # ---- RNA side: BAMs -> coverage hd5 -------------------------------------
    rna_coverage = list(args.rna_coverage)
    if args.rna_bam:
        for bam in args.rna_bam:
            name = Path(bam).name.split(".")[0]
            out_cov = cov_dir / f"{name}_coverage.hd5"
            run([args.pysam_python, scripts / "rnaseq_coverage.py", bam, out_cov, name])
        rna_coverage = [str(cov_dir)]

    # ---- pack ---------------------------------------------------------------
    cmd = [sys.executable, scripts / "prepare" / "prepare_pack.py",
           "--rna-coverage", *rna_coverage,
           "--group", args.group, "--species", args.species, "--out", out_dir]
    if args.coverage_only:
        cmd.append("--coverage-only")
    else:
        cmd += ["--ribo-psites", *ribo_psites]
    cmd += (["--ref-pack", args.ref_pack] if args.ref_pack else ["--universe-tx", args.universe_tx])
    if args.build_orf_track:
        cmd += ["--build-orf-track", "--fasta", args.fasta, "--kozak", args.kozak]
    if args.tx2biotype:
        cmd += ["--tx2biotype", args.tx2biotype]
    elif args.gtf:
        cmd += ["--gtf", args.gtf]
    run(cmd)

    if dry:
        print(f"\n# dry-run: {len(cmds)} commands planned, nothing executed.", file=sys.stderr)
    else:
        print(f"\n# pack built at {out_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
