#!/usr/bin/env python3
"""Write the cross-species dataset registry: one row per sequencing run, joining every
piece of provenance and QC produced along the way.

Sources, all of them artifacts the pipeline already wrote:
  download_manifest.tsv   what was selected and why (species, assay, tissue, study, layout)
  adapter_probe.tsv       the MEASURED adapter verdict per run
  adapter_overrides.tsv   where the probe was overridden, and on what evidence
  sra_provenance.tsv      which files came from NCBI SRA rather than ENA, and how they
                          were validated (an SRA file cannot match an ENA md5)
  xspecies_align_qc.tsv   per-run STAR and cutadapt metrics
  on disk                 whether the BAM, coverage hd5 and P-sites actually exist

The point of the registry is that every run can be traced from accession to artifact
without reading a log, and that the runs whose provenance DIFFERS from the norm (SRA
sourced, adapter overridden, dropped species) are visible rather than buried.

Usage: python3 scripts/xspecies/build_dataset_registry.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
X = NEW / "data" / "external" / "xspecies"
MD = NEW / "docs" / "XSPECIES_DATASET_REGISTRY.md"
TSV = NEW / "data" / "xspecies_dataset_registry.tsv"

DROPPED = {"fly_gse99920_ribo", "fly_gse99920_rna"}
DROP_REASON = ("RNase T1 footprinting: no 3-nt periodicity (15.2/56.2/28.7 at the dominant "
               "read length vs 93.9/2.6/3.6 for yeast). Alignments retained, modelling arm dropped.")


def load(path: Path, key: str) -> dict[str, dict]:
    if not path.exists():
        return {}
    with path.open() as fh:
        rows = [r for r in csv.DictReader(
            (l for l in fh if not l.startswith("#")), delimiter="\t")]
    return {r[key]: r for r in rows if r.get(key)}


def main() -> int:
    man = load(X / "download_manifest.tsv", "run_accession")
    probe = load(X / "adapter_probe.tsv", "run")
    qc = {}
    qp = NEW / "results" / "xspecies_align_qc.tsv"
    if qp.exists():
        for r in csv.DictReader(qp.open(), delimiter="\t"):
            qc.setdefault(r["run"], {})[r["posture"]] = r
    sra = set()
    sp = X / "sra_provenance.tsv"
    if sp.exists():
        for r in csv.DictReader(sp.open(), delimiter="\t"):
            sra.add(r["run"])
    overrides = {}
    op = X / "adapter_overrides.tsv"
    if op.exists():
        for line in op.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if len(f) >= 3 and f[0] != "label":
                overrides[f[0]] = f[2]

    rows = []
    for run, m in sorted(man.items(), key=lambda kv: (kv[1]["species"], kv[1]["assay"], kv[0])):
        label, assay = m["label"], m["assay"]
        pr = probe.get(run, {})
        q1 = qc.get(run, {}).get("mm1") or qc.get(run, {}).get("mm10") or {}
        bamdir = (NEW / "data" / ("xspecies_ribo_bam_mm1" if assay == "ribo"
                                  else "xspecies_rna_bam_mm10") / label)
        bam = bamdir / f"{run}.Aligned.toTranscriptome.out.bam"
        # Every posture this run exists in, and whether each is local or archived. A path that
        # is neither is genuinely absent, which is different from archived and must not read
        # the same: `-` in the table means the posture was never built for this arm.
        extra = {}
        for col, root, fname in (
            ("bam_mm25", "xspecies_ribo_bam_mm25", f"{run}.Aligned.toTranscriptome.out.bam"),
            ("bam_genome", "xspecies_ribo_bam_mm1_genome", f"{run}.genome.bam"),
            ("bigwig", "xspecies_ribo_bigwig_genome", f"{run}.bw"),
        ):
            if assay != "ribo":
                extra[col] = "-"
                continue
            d = NEW / "data" / root / label
            extra[col] = ("yes" if (d / fname).exists()
                          else "ARCHIVED" if (d.parent / f"{label}.ARCHIVED.md").exists()
                          or (d.parent.parent / f"{root}.ARCHIVED.md").exists()
                          else "-")
        psi25 = (NEW / "data" / "xspecies_psites_mm25" / label /
                 f"{run}.Aligned.toTranscriptome.out_psites.hd5")
        extra["psites_mm25"] = "yes" if psi25.exists() else "-"
        cov = NEW / "data" / "xspecies_rna_coverage_mm10" / label / f"{run}_coverage.hd5"
        psi = NEW / "data" / "xspecies_psites" / label / f"{run}.Aligned.toTranscriptome.out_psites.hd5"
        rows.append({
            "species": m["species"], "assay": assay, "tissue": m["tissue"],
            "label": label, "run": run, "project": m["project"],
            "sample_alias": m["sample_alias"], "layout": m["layout"],
            "gib": m["gib"], "read_len_est": m["read_len_est"],
            "source": "NCBI_SRA" if run in sra else "ENA",
            "adapter_measured": pr.get("best_adapter", ""),
            "adapter_pct": pr.get("best_pct", ""),
            "adapter_used": overrides.get(label) and "OVERRIDE" or pr.get("cutadapt_arg", ""),
            "adapter_override_reason": overrides.get(label, ""),
            "pct_trimmed_kept": q1.get("pct_trimmed_kept", ""),
            "avg_len": q1.get("avg_len", ""),
            "pct_unique": q1.get("pct_unique", ""),
            "pct_multi_many": q1.get("pct_multi_many", ""),
            "bam": "yes" if bam.exists() else ("ARCHIVED" if bamdir.with_suffix(".ARCHIVED.md").exists() or (bamdir.parent / f"{label}.ARCHIVED.md").exists() else "no"),
            "coverage_hd5": "yes" if cov.exists() else "-",
            "psites_hd5": "yes" if psi.exists() else "-",
            **extra,
            "status": "DROPPED" if label in DROPPED else "active",
            "note": DROP_REASON if label in DROPPED else "",
        })

    cols = list(rows[0].keys())
    with TSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    # Markdown summary, per arm rather than per run: 144 rows is a table nobody reads.
    from collections import defaultdict
    arms = defaultdict(list)
    for r in rows:
        arms[(r["species"], r["assay"], r["label"])].append(r)

    L = ["# Cross-species dataset registry", "",
         "Generated by `scripts/xspecies/build_dataset_registry.py`. Per-run detail, including the",
         "measured adapter and per-run STAR metrics, is in `data/xspecies_dataset_registry.tsv`.",
         "", "## Arms", "",
         "| species | assay | label | runs | GiB | source | adapter | kept% | uniq% | BAM | cov | psites | mm25 | genome bw | status |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for (spc, assay, label), rs in sorted(arms.items()):
        gib = sum(float(r["gib"]) for r in rs)
        srcs = sorted({r["source"] for r in rs})
        ad = sorted({r["adapter_used"] for r in rs if r["adapter_used"]})
        def mean(k):
            v = [float(r[k]) for r in rs if r.get(k)]
            return f"{sum(v)/len(v):.1f}" if v else "-"
        L.append(f"| {spc} | {assay} | `{label}` | {len(rs)} | {gib:.1f} | {'+'.join(srcs)} | "
                 f"{', '.join(ad) or '-'} | {mean('pct_trimmed_kept')} | {mean('pct_unique')} | "
                 f"{sum(r['bam']!='no' for r in rs)}/{len(rs)} | "
                 f"{sum(r['coverage_hd5']=='yes' for r in rs)}/{len(rs)} | "
                 f"{sum(r['psites_hd5']=='yes' for r in rs)}/{len(rs)} | "
                 f"{sum(r['bam_mm25']!='-' for r in rs)}/{len(rs)} | "
                 f"{sum(r['bigwig']!='-' for r in rs)}/{len(rs)} | {rs[0]['status']} |")

    nsra = sum(1 for r in rows if r["source"] == "NCBI_SRA")
    novr = sorted({r["label"] for r in rows if r["adapter_override_reason"]})
    L += ["", "## Runs whose provenance differs from the norm", "",
          f"**{nsra} run(s) came from NCBI SRA, not ENA.** An SRA-derived FASTQ cannot reproduce",
          "ENA's published md5 (fasterq-dump regenerates the file, so gzip framing and read order",
          "differ), so these were validated on read count against ENA's `read_count` and, where a",
          "mate was already ENA-verified, on read-ID identity. See `sra_provenance.tsv`.",
          "`verify_and_repair.sh` skips them rather than deleting them as corrupt.", ""]
    if novr:
        L += ["**Adapter overrides:**", ""]
        for l in novr:
            r = next(x for x in rows if x["label"] == l)
            L += [f"- `{l}`: {r['adapter_override_reason']}", ""]
    L += ["## Postures", "",
          "Each Ribo run exists in up to four forms. `yes` = on disk, `ARCHIVED` = in the warm",
          "archive with a stub in place (see `data/ARCHIVE_INDEX.tsv`), `-` = never built.", "",
          "- **mm1 transcriptome** (`data/xspecies_ribo_bam_mm1/`): the ONLY posture used as a",
          "  model target, per the project-wide Ribo-seq rule.",
          "- **mm25 transcriptome** (`data/xspecies_ribo_bam_mm25/`): the multimap DIAGNOSTIC.",
          "  Call sets in `data/xspecies_psites_mm25/`, comparison in",
          "  `results/xspecies_mm1_vs_mm25_calls.tsv`. Six of seven comparable arms sit at",
          "  1.000-1.013x the mm1 call count. NOT a model target; must not feed a pack.",
          "- **mm1 genome** (`data/xspecies_ribo_bam_mm1_genome/`) and its bigwigs",
          "  (`data/xspecies_ribo_bigwig_genome/`, `--binSize 1`, unnormalized): browser tracks.",
          "  Same posture as mm1 transcriptome; the ncRNA drop is by genomic interval rather than",
          "  by transcript id. Contigs are RefSeq-style, so load",
          "  `genomes/xspecies_refs/<species>/genome.fna` as the IGV reference.", "",
          "## Dropped", "",
          f"- `fly_gse99920_*`: {DROP_REASON} See methods.md section G.", ""]
    MD.write_text("\n".join(L) + "\n")
    print(f"  {len(rows)} runs -> {TSV}")
    print(f"  summary -> {MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
