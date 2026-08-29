#!/usr/bin/env python3
"""Apply the cross-species run-selection rules and emit one download manifest.

Reads the full ENA run tables written by build_run_tables.sh and writes
data/external/xspecies/download_manifest.tsv, which is the single authoritative
list of what gets fetched. Every exclusion is stated as a rule here rather than
as a hand-typed accession list, so the reasoning survives and re-running after an
ENA metadata correction picks the change up.

Two things this deliberately does NOT do:

  * It does not trust `library_strategy` to tell Ribo from RNA. That vocabulary
    term postdates most of these submissions: zero worm/fly/yeast/zebrafish runs
    are labelled Ribo-Seq, and they hide under RNA-Seq or OTHER. Assay is decided
    per project by the rule that is actually reliable for that project
    (scientific_name + library_strategy for PRJEB65856, sample_title everywhere
    else), and the chosen basis is recorded in the `assay_basis` column.

  * It does not treat the read length it computes as ground truth. ENA is
    inconsistent about whether read_count counts spots or mates, so
    base_count/read_count is 2x the read length on some projects and 1x on
    others. `read_len_est` applies a heuristic and `len_basis` says which branch
    fired. The real length and the real adapter status are measured off the
    FASTQ in Phase 4; this column only orders the work and sizes the indexes.

Usage: python3 scripts/xspecies/select_runs.py [--check]
  --check  validate against the expected per-species counts and exit nonzero on
           any mismatch, without rewriting the manifest.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
XDIR = NEW / "data" / "external" / "xspecies"
MANIFEST = XDIR / "download_manifest.tsv"

# The 8 shared timepoints between the worm Ribo (contDevB) and RNA (contDevA)
# timecourses. The arms are replicate timecourses, matched by stage rather than
# by lysate -- recorded in the registry so no later analysis assumes paired samples.
WORM_HOURS = (22, 24, 26, 28, 30, 32, 34, 36)
WORM_RIBO_RE = re.compile(r"^Footprint_contDevB_N2_(\d+)hr$")
WORM_RNA_RE = re.compile(r"^RNASeq_polyA_contDevA_N2_(\d+)hr$")

# sample_alias on PRJEB65856 is <species>_<tissue>_<individual>_<mR|Ri>.
PRIMATE_ALIAS_RE = re.compile(r"^(hs|pt|gg|rm)_(CM|lv|iPS)_(.+)_(mR|Ri)$")
SPECIES_FULL = {"hs": "human", "pt": "chimp", "gg": "gorilla", "rm": "macaque"}
TISSUE_FULL = {"CM": "ipsc_cm", "lv": "left_ventricle", "iPS": "ipsc"}

# Download order: smallest first, so a problem in the fetcher surfaces cheaply.
ORDER = [
    "ruizorera_hsCM_ribo",
    "fly_gse99920_rna", "fly_gse99920_ribo",
    "yeast_gse173654_rna", "yeast_gse173654_ribo",
    "worm_gse52861_rna", "worm_gse52905_ribo",
    "primate_gg_rna", "primate_gg_ribo",
    "primate_rm_rna", "primate_rm_ribo",
    "primate_pt_rna", "primate_pt_ribo",
    "zf_gse70549_rna", "zf_gse46512_ribo",
]

# Expected counts, asserted by --check. A silent drift here is exactly the kind
# of thing that produces a species with a missing replicate nobody notices.
EXPECT = {
    ("human", "ribo"): 5,
    ("gorilla", "rna"): 3, ("gorilla", "ribo"): 3,
    ("chimp", "rna"): 8, ("chimp", "ribo"): 8,
    ("macaque", "rna"): 7, ("macaque", "ribo"): 7,
    ("zebrafish", "rna"): 24, ("zebrafish", "ribo"): 8,
    ("celegans", "rna"): 10, ("celegans", "ribo"): 11,
    ("yeast", "rna"): 16, ("yeast", "ribo"): 16,
    ("fly", "rna"): 9, ("fly", "ribo"): 9,
}


def read_table(acc: str) -> list[dict]:
    path = XDIR / f"{acc}_runs.tsv"
    if not path.exists():
        sys.exit(f"missing run table {path}; run scripts/xspecies/build_run_tables.sh first")
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def gib(row: dict) -> float:
    return sum(int(b) for b in row["fastq_bytes"].split(";") if b) / 1073741824


def read_len(row: dict) -> tuple[int, str]:
    """Estimated per-mate read length, and which branch produced it.

    ENA is not consistent about whether read_count counts spots or mates, so
    base_count/read_count is 2x the read length on some paired projects and 1x
    on others. Anything over 120 on a paired run is taken as a spot total.
    """
    rc, bc = int(row["read_count"] or 0), int(row["base_count"] or 0)
    if rc <= 0:
        return 0, "unknown"
    raw = bc / rc
    if row["library_layout"] == "PAIRED" and raw > 120:
        return round(raw / 2), "paired_spot_halved"
    return round(raw), "direct"


def n_files(row: dict) -> int:
    return len([f for f in row["fastq_ftp"].split(";") if f])


def emit(rows: list[dict], **kw) -> list[dict]:
    out = []
    for r in rows:
        rl, basis = read_len(r)
        rec = dict(kw)
        rec.update({
            "run_accession": r["run_accession"],
            "study_accession": r["study_accession"],
            "sample_alias": r["sample_alias"] or r["sample_title"],
            "layout": r["library_layout"],
            "n_fastq_files": n_files(r),
            "read_len_est": rl,
            "len_basis": basis,
            "gib": f"{gib(r):.2f}",
            "library_strategy": r["library_strategy"],
            "library_selection": r["library_selection"],
            "fastq_ftp": r["fastq_ftp"],
            "fastq_md5": r["fastq_md5"],
        })
        out.append(rec)
    return out


def select() -> list[dict]:
    picked: list[dict] = []

    # ---- PRJEB65856 Ruiz-Orera -------------------------------------------
    # Keep every run whose tissue is iPSC-CM or left ventricle; drop the 12
    # undifferentiated-iPSC runs, which are RNA-only with no Ribo partner.
    # Human RNA is already on disk byte-verified in
    # data/external/heldout_rna_mm10/ruizorera_hsCM_rna/, so only the human Ribo
    # arm is fetched: it exists locally as BAMs and P-sites but has no FASTQ
    # anywhere, which would leave human as the one iPSC-CM species with no mm25
    # diagnostic and no option to re-trim.
    for r in read_table("PRJEB65856"):
        m = PRIMATE_ALIAS_RE.match(r["sample_alias"] or "")
        if not m:
            sys.exit(f"PRJEB65856 alias did not parse: {r['sample_alias']!r} ({r['run_accession']})")
        sp, tissue, _indiv, assay_sfx = m.groups()
        if tissue == "iPS":
            continue
        assay = "ribo" if assay_sfx == "Ri" else "rna"
        species = SPECIES_FULL[sp]
        if species == "human" and assay == "rna":
            continue
        label = ("ruizorera_hsCM_ribo" if species == "human"
                 else f"primate_{sp}_{assay}")
        picked += emit([r], label=label, species=species, assay=assay,
                       tissue=TISSUE_FULL[tissue], project="PRJEB65856",
                       assay_basis="sample_alias suffix mR|Ri + library_strategy agree")

    # ---- zebrafish: GSE46512 Ribo + GSE70549 RNA, matched stage-for-stage --
    # GSE70549's own summary says it is "sample matched to a previous ribosome
    # profiling dataset GSE46512". Both arms are labelled RNA-Seq; the Ribo arm
    # is identified by library_selection="size fractionation".
    picked += emit(read_table("PRJNA200706"), label="zf_gse46512_ribo",
                   species="zebrafish", assay="ribo", tissue="embryo",
                   project="PRJNA200706",
                   assay_basis='library_selection="size fractionation" + RPF-Seq sample_title')
    picked += emit(read_table("PRJNA288987"), label="zf_gse70549_rna",
                   species="zebrafish", assay="rna", tissue="embryo",
                   project="PRJNA288987", assay_basis="whole project is the RNA arm")

    # ---- C. elegans: the 8 shared timepoints of the two timecourses ---------
    for tbl, rx, label, assay, proj in (
        ("PRJNA230441", WORM_RIBO_RE, "worm_gse52905_ribo", "ribo", "PRJNA230441"),
        ("PRJNA230374", WORM_RNA_RE, "worm_gse52861_rna", "rna", "PRJNA230374"),
    ):
        keep = [r for r in read_table(tbl)
                if (m := rx.match(r["sample_title"] or "")) and int(m.group(1)) in WORM_HOURS]
        picked += emit(keep, label=label, species="celegans", assay=assay,
                       tissue="larval_timecourse", project=proj,
                       assay_basis="sample_title timecourse+assay prefix")

    # ---- yeast: whole project, 8 strains x 2 reps, both arms ---------------
    for r in read_table("PRJNA726548"):
        t = r["sample_title"] or ""
        if t.startswith("Ribo_"):
            assay, label = "ribo", "yeast_gse173654_ribo"
        elif t.startswith("RNA_"):
            assay, label = "rna", "yeast_gse173654_rna"
        else:
            sys.exit(f"PRJNA726548 unexpected sample_title {t!r} ({r['run_accession']})")
        picked += emit([r], label=label, species="yeast", assay=assay,
                       tissue="whole_cell", project="PRJNA726548",
                       assay_basis="sample_title RNA_|Ribo_ prefix")

    # ---- fly: 3 genotypes x 3 reps, both arms; drop the 3 TRAP runs --------
    # TRAP is a different assay (affinity-tagged ribosome pulldown), not
    # footprinting, so it cannot serve as a per-nt Ribo target.
    for r in read_table("PRJNA390134"):
        t = r["sample_title"] or ""
        if "translational profiling (TRAP)" in t:
            continue
        if "ribosome profiling" in t:
            assay, label = "ribo", "fly_gse99920_ribo"
        elif "transcriptional profiling" in t:
            assay, label = "rna", "fly_gse99920_rna"
        else:
            sys.exit(f"PRJNA390134 unexpected sample_title {t!r} ({r['run_accession']})")
        picked += emit([r], label=label, species="fly", assay=assay,
                       tissue="larval_muscle", project="PRJNA390134",
                       assay_basis="sample_title assay phrase")

    picked.sort(key=lambda r: (ORDER.index(r["label"]) if r["label"] in ORDER else 99,
                               r["run_accession"]))
    return picked


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    rows = select()

    counts: dict[tuple[str, str], int] = {}
    for r in rows:
        counts[(r["species"], r["assay"])] = counts.get((r["species"], r["assay"]), 0) + 1

    bad = []
    for k, want in EXPECT.items():
        got = counts.get(k, 0)
        if got != want:
            bad.append(f"  {k[0]:10s} {k[1]:4s}  expected {want}, got {got}")
    for k in counts:
        if k not in EXPECT:
            bad.append(f"  {k[0]:10s} {k[1]:4s}  unexpected group, got {counts[k]}")

    total_gib = sum(float(r["gib"]) for r in rows)
    print(f"{len(rows)} runs selected, {total_gib:.1f} GiB")
    for (sp, assay), n in sorted(counts.items()):
        g = sum(float(r["gib"]) for r in rows if r["species"] == sp and r["assay"] == assay)
        print(f"  {sp:10s} {assay:4s} {n:3d} runs  {g:7.1f} GiB")

    # ENA sometimes lists an extra merged SRR.fastq.gz alongside _1/_2 on paired
    # runs. fetch_heldout_rna.sh would download all three, wasting ~50% of the
    # bytes, so flag it here rather than discovering it as an inflated transfer.
    tri = [r for r in rows if r["layout"] == "PAIRED" and r["n_fastq_files"] != 2]
    single = [r for r in rows if r["layout"] == "SINGLE" and r["n_fastq_files"] != 1]
    if tri:
        print(f"\nWARNING: {len(tri)} PAIRED runs do not have exactly 2 fastq files "
              f"(redundant merged file): {[r['run_accession'] for r in tri][:8]}")
    if single:
        print(f"WARNING: {len(single)} SINGLE runs do not have exactly 1 fastq file: "
              f"{[r['run_accession'] for r in single][:8]}")
    nomd5 = [r for r in rows
             if len(r["fastq_md5"].split(";")) != r["n_fastq_files"] or not r["fastq_md5"]]
    if nomd5:
        print(f"WARNING: {len(nomd5)} runs have an md5 count that does not match the file count")

    if bad:
        print("\nCOUNT MISMATCH:")
        print("\n".join(bad))
        return 1

    if args.check:
        print("\ncheck OK")
        return 0

    cols = ["label", "species", "assay", "tissue", "project", "run_accession",
            "study_accession", "sample_alias", "layout", "n_fastq_files",
            "read_len_est", "len_basis", "gib", "library_strategy",
            "library_selection", "assay_basis", "fastq_ftp", "fastq_md5"]
    # lineterminator="\n" is not cosmetic: csv defaults to "\r\n", which makes the
    # LAST column of every row end in a carriage return. The downstream shell
    # fetcher compares fastq_md5 by string equality, so a stray \r would fail
    # every checksum on every file. Caught by a dry run; do not remove.
    with MANIFEST.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {MANIFEST}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
