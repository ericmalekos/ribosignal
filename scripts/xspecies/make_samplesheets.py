#!/usr/bin/env python3
"""Build alignment samplesheets by joining the download manifest to the MEASURED
adapter verdicts.

The adapter column is never guessed and never set per dataset: it comes from
`adapter_probe.tsv`, which measured each run's complete FASTQ. That matters because
`primate_rm_ribo` genuinely contains two regimes (4 runs TruSeq at 51 nt, 3 runs
N-padded at 35 nt), and because the worm Ribo arm uses the Illumina small-RNA adapter
`TGGAATTCTCGGGTGCCAAGG` where the standing cluster rule would have said TruSeq and
trimmed nothing.

A run with no probe row is REFUSED rather than defaulted. Defaulting is exactly how a
library gets aligned with the wrong adapter and produces a near-empty BAM with no error.

Outputs, under data/external/xspecies/samplesheets/:
  ribo_<species>.tsv   dataset run fastq adapter species umi cutadapt_extra
                       -> consumed by scripts/riboseq_align.sbatch, unchanged
  rna_<species>.tsv    dataset run r1 r2 adapter species layout
                       -> consumed by scripts/xspecies/align_rna_xspecies.sbatch

Species naming: the cross-species arm uses `human_refseq` for the human runs, not the
GENCODE `human`, so all nine species are annotated by one source. The existing
GENCODE-based `human_ruizorera` held-out arm is untouched.

Usage: python3 scripts/xspecies/make_samplesheets.py [--species X ...] [--check]
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
XDIR = NEW / "data" / "external" / "xspecies"
FQROOT = XDIR / "fastq"
OUTDIR = XDIR / "samplesheets"

# manifest species label -> the species key used by species.tsv / riboseq_align.sbatch
SPECIES_KEY = {
    "human": "human_refseq",
    "gorilla": "gorilla", "chimp": "chimp", "macaque": "macaque",
    "zebrafish": "zebrafish", "celegans": "celegans",
    "fly": "fly", "yeast": "yeast",
}


def load_probe() -> dict[str, dict]:
    p = XDIR / "adapter_probe.tsv"
    if not p.exists():
        sys.exit(f"missing {p}; run scripts/xspecies/probe_adapters.sh first")
    return {r["run"]: r for r in csv.DictReader(p.open(), delimiter="\t")}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", nargs="*", default=None)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    probe = load_probe()
    ov_path = XDIR / "adapter_overrides.tsv"
    overrides: dict[str, str] = {}
    if ov_path.exists():
        for line in ov_path.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) >= 2 and f[0] != "label":
                overrides[f[0]] = f[1]
    n_overridden: set[str] = set()
    needs_adapter: list[str] = []
    rows = list(csv.DictReader((XDIR / "download_manifest.tsv").open(), delimiter="\t"))
    OUTDIR.mkdir(parents=True, exist_ok=True)

    ribo: dict[str, list] = defaultdict(list)
    rna: dict[str, list] = defaultdict(list)
    missing_probe: list[str] = []
    missing_fastq: list[str] = []

    for r in rows:
        sp = SPECIES_KEY[r["species"]]
        if args.species and sp not in args.species:
            continue
        run, label = r["run_accession"], r["label"]
        files = [f for f in r["fastq_ftp"].split(";") if f]
        if len(files) > 1:
            files = [f for f in files if f.endswith(("_1.fastq.gz", "_2.fastq.gz"))]
        paths = [FQROOT / label / Path(f).name for f in files]
        absent = [str(p) for p in paths if not p.exists()]
        if absent:
            missing_fastq += absent
            continue

        pr = probe.get(run)
        if pr is None:
            missing_probe.append(f"{run} ({label})")
            continue
        adapter = pr["cutadapt_arg"]
        if label in overrides:
            adapter = overrides[label]
            n_overridden.add(label)
        elif adapter == "UNRESOLVED" or pr.get("verdict") == "NEEDS_ADAPTER":
            # The probe flagged that cutadapt would discard essentially every read here.
            # Refuse rather than emit a row that produces an empty BAM.
            needs_adapter.append(f"{run} ({label})")
            continue

        if r["assay"] == "ribo":
            # riboseq_align.sbatch columns: dataset run fastq adapter species umi cutadapt_extra
            ribo[sp].append([label, run, str(paths[0]), adapter, sp, "none", ""])
        else:
            r1 = str(paths[0])
            # "-" not "": an empty field between two tabs collapses under bash's
            # IFS-whitespace rule and shifts every subsequent column left.
            r2 = str(paths[1]) if len(paths) > 1 else "-"
            layout = "PAIRED" if r2 != "-" else "SINGLE"
            rna[sp].append([label, run, r1, r2, adapter, sp, layout])

    if needs_adapter:
        print(f"REFUSING: {len(needs_adapter)} run(s) are flagged NEEDS_ADAPTER by the probe -- "
              f"cutadapt would discard essentially every read.", file=sys.stderr)
        for m in needs_adapter[:8]:
            print(f"  {m}", file=sys.stderr)
        print("Determine the real adapter from the 3' end of the reads and add a row to\n"
              "data/external/xspecies/adapter_overrides.tsv with the measurement.", file=sys.stderr)
        return 1
    if missing_probe:
        print(f"REFUSING: {len(missing_probe)} run(s) have no adapter probe row.", file=sys.stderr)
        for m in missing_probe[:10]:
            print(f"  {m}", file=sys.stderr)
        print("Run scripts/xspecies/probe_adapters.sh for those labels. Not defaulting: a wrong\n"
              "adapter yields a near-empty BAM rather than an error.", file=sys.stderr)
        return 1
    if missing_fastq:
        print(f"NOTE: {len(missing_fastq)} FASTQ still absent or unverified, those runs skipped:",
              file=sys.stderr)
        for m in missing_fastq[:6]:
            print(f"  {Path(m).name}", file=sys.stderr)

    total = 0
    for sp in sorted(set(list(ribo) + list(rna))):
        if ribo[sp]:
            f = OUTDIR / f"ribo_{sp}.tsv"
            if not args.check:
                with f.open("w", newline="") as fh:
                    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                    w.writerow(["dataset", "run", "fastq", "adapter", "species", "umi", "cutadapt_extra"])
                    w.writerows(sorted(ribo[sp]))
            adapters = sorted({x[3] for x in ribo[sp]})
            print(f"  ribo_{sp:<14} {len(ribo[sp]):3d} runs   adapters: {', '.join(adapters)}")
            total += len(ribo[sp])
        if rna[sp]:
            f = OUTDIR / f"rna_{sp}.tsv"
            if not args.check:
                with f.open("w", newline="") as fh:
                    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                    w.writerow(["dataset", "run", "r1", "r2", "adapter", "species", "layout"])
                    w.writerows(sorted(rna[sp]))
            lay = sorted({x[6] for x in rna[sp]})
            print(f"  rna_{sp:<15} {len(rna[sp]):3d} runs   layout: {', '.join(lay)}")
            total += len(rna[sp])
    for l in sorted(n_overridden):
        print(f"  NOTE: {l} uses an adapter OVERRIDE from adapter_overrides.tsv")
    print(f"\n  {total} runs across {len(set(list(ribo)+list(rna)))} species -> {OUTDIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
