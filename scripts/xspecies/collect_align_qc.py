#!/usr/bin/env python3
"""Collect per-run alignment QC for the cross-species arm into one table.

Every number here comes from an artifact the aligner already wrote (cutadapt log,
STAR Log.final.out, the filter's own accounting), so this is a reader, not a
re-computation.

WHAT TO LOOK AT, and why each column earns its place:

  pct_trimmed_kept  cutadapt survival. A collapse here means the adapter was wrong.
                    An untrimmed Janich library once mapped 0.03% and produced a ~5 KB
                    BAM stub rather than an error.
  avg_len           STAR's average input read length. For a Ribo arm this should be
                    footprint-sized (~26-32 nt). If it is 50+, trimming did not happen.
  pct_unique        uniquely mapped. LOW IS NOT AUTOMATICALLY BAD for Ribo: with
                    --outFilterMultimapNmax 1 an rRNA-heavy library sends nearly
                    everything to "too many loci". Read it together with pct_multi.
  pct_multi         "mapped to too many loci". High + low unique = undepleted rRNA,
                    which is a property of the library, not a pipeline fault.
  pct_short         "unmapped: too short". THIS is the column that indicates a real
                    alignment problem: untrimmed adapter shows up here, not in pct_multi.

Usage: python3 scripts/xspecies/collect_align_qc.py [--tsv <out.tsv>]
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
ARMS = [
    ("ribo", "mm1", NEW / "data" / "xspecies_ribo_bam_mm1"),
    ("ribo", "mm25", NEW / "data" / "xspecies_ribo_bam_mm25"),
    ("rna", "mm10", NEW / "data" / "xspecies_rna_bam_mm10"),
    ("rna", "mm1", NEW / "data" / "xspecies_rna_bam_mm1"),
]


def star_metrics(p: Path) -> dict:
    out = {}
    if not p.exists():
        return out
    for line in p.read_text(errors="replace").splitlines():
        if "|" not in line:
            continue
        k, v = [x.strip() for x in line.split("|", 1)]
        out[k] = v
    def num(k, cast=float):
        v = out.get(k, "").rstrip("%")
        try:
            return cast(v)
        except Exception:
            return None
    return {
        "n_input": num("Number of input reads", int),
        "avg_len": num("Average input read length", int),
        "pct_unique": num("Uniquely mapped reads %"),
        "pct_multi_many": num("% of reads mapped to too many loci"),
        "pct_multi": num("% of reads mapped to multiple loci"),
        "pct_short": num("% of reads unmapped: too short"),
    }


def cutadapt_metrics(p: Path) -> dict:
    if not p.exists():
        return {}
    t = p.read_text(errors="replace")
    def grab(pat):
        m = re.search(pat, t)
        return int(m.group(1).replace(",", "")) if m else None
    total = grab(r"Total reads processed:\s+([\d,]+)") or \
            grab(r"Total read pairs processed:\s+([\d,]+)")
    kept = grab(r"Reads written \(passing filters\):\s+([\d,]+)") or \
           grab(r"Pairs written \(passing filters\):\s+([\d,]+)")
    with_ad = grab(r"Reads with adapters:\s+([\d,]+)")
    return {"n_raw": total, "n_kept": kept,
            "pct_trimmed_kept": round(100 * kept / total, 1) if total and kept else None,
            "pct_with_adapter": round(100 * with_ad / total, 1) if total and with_ad else None}


def filter_metrics(p: Path) -> dict:
    """filter_tx_heldout accounting, printed into the aligner's stdout log."""
    if not p.exists():
        return {}
    m = re.search(r"records total=([\d,]+) drop_ncRNA=([\d,]+) drop_xgene=([\d,]+)", p.read_text(errors="replace"))
    if not m:
        return {}
    tot, nc, xg = (int(x.replace(",", "")) for x in m.groups())
    return {"filt_total": tot, "filt_drop_ncrna": nc, "filt_drop_xgene": xg,
            "filt_pct_dropped": round(100 * (nc + xg) / tot, 2) if tot else None}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", type=Path, default=NEW / "results" / "xspecies_align_qc.tsv")
    args = ap.parse_args()

    rows = []
    for assay, posture, root in ARMS:
        if not root.exists():
            continue
        for logf in sorted(root.rglob("*.STAR.Log.final.out")):
            run = logf.name.replace(".STAR.Log.final.out", "")
            ds = logf.parent.name
            rec = {"dataset": ds, "run": run, "assay": assay, "posture": posture}
            rec.update(cutadapt_metrics(logf.parent / f"{run}.cutadapt.log"))
            rec.update(star_metrics(logf))
            rows.append(rec)

    if not rows:
        print("no alignment logs found yet", file=sys.stderr)
        return 1

    rows.sort(key=lambda r: (r["assay"], r["dataset"], r["run"]))
    hdr = f"  {'dataset':<24}{'run':<13}{'post':<6}{'raw':>11}{'kept%':>7}{'len':>5}{'uniq%':>7}{'many%':>7}{'short%':>7}"
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    for r in rows:
        print(f"  {r['dataset']:<24}{r['run']:<13}{r['posture']:<6}"
              f"{(r.get('n_raw') or 0):>11,}{(r.get('pct_trimmed_kept') or 0):>7.1f}"
              f"{(r.get('avg_len') or 0):>5}{(r.get('pct_unique') or 0):>7.2f}"
              f"{(r.get('pct_multi_many') or 0):>7.2f}{(r.get('pct_short') or 0):>7.2f}")

    print("\n  flags:")
    flagged = 0
    for r in rows:
        why = []
        if (r.get("pct_short") or 0) > 25:
            why.append(f"{r['pct_short']:.0f}% unmapped-too-short (adapter or length window?)")
        if (r.get("pct_trimmed_kept") or 100) < 40:
            why.append(f"only {r['pct_trimmed_kept']:.0f}% survived trimming")
        if r["assay"] == "ribo" and (r.get("avg_len") or 0) > 45:
            why.append(f"avg read {r['avg_len']} nt is not footprint-sized")
        if (r.get("pct_unique") or 0) < 5 and (r.get("pct_multi_many") or 0) < 50:
            why.append(f"only {r['pct_unique']:.1f}% unique and NOT explained by multimapping")
        if why:
            flagged += 1
            print(f"    {r['dataset']}/{r['run']}: " + "; ".join(why))
    if not flagged:
        print("    none")

    args.tsv.parent.mkdir(parents=True, exist_ok=True)
    cols = ["dataset", "run", "assay", "posture", "n_raw", "n_kept", "pct_trimmed_kept",
            "pct_with_adapter", "n_input", "avg_len", "pct_unique", "pct_multi",
            "pct_multi_many", "pct_short"]
    with args.tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"\n  {len(rows)} runs -> {args.tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
