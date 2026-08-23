#!/usr/bin/env python3
"""One row per Ribo-seq dataset the project touches: depth, frame-0, role, recipe status.

WHY. The poster's dataset-census panel (2a) needs depth and frame-0 for every dataset, and those
numbers currently live in five places -- `P3_values.json`, `E3_values.json`, `riboseq_qc_values.json`,
per-pack `target_counts.npy`, and prose in results.md. Four rows had no depth at all and two had no
frame-0. Assembling that by hand across five sources is how a poster ends up mixing quantities.

RNA COLUMNS (added 2026-08-16). The model is fed sequence + RNA-seq, so a census showing only the
target modality does not show what the model consumes.

  n_rna_libraries      count of RNA coverage files pooled into the pack
  rna_coverage_total   summed RNA coverage over the pack's own universe -- a RAW TOTAL, verified
                       equal to coverage.npy.sum() on int32 per-nt counts. NOT normalised, NOT
                       depth-corrected. The pack's `global_mean_coverage` is exactly
                       rna_coverage_total / sum_L, i.e. the per-nt mean of this same quantity.
  rna_source           where the count came from, because it is NOT uniformly reliable (below)

WHY n_rna_libraries IS NOT SIMPLY coverage_norm.json's n_rna_samples. Every `_canon` pack reports
n_rna_samples = 0 while holding real coverage -- they were built by REUSING their union twin's
coverage verbatim (correct: the final recipe leaves RNA untouched) and the rebuild never carried the
count forward. packed_canon reports 0 against 1.2 TRILLION coverage. Shipping that would read as
"no RNA input", which is false and is exactly the psites_at_cds class of error.

Resolution order, most authoritative first:
  1. provenance.json `rna_coverage_inputs`  -- the actual list of files pooled
  2. coverage_norm.json `n_rna_samples`     -- when > 0
  3. the union twin's count, ACCEPTED ONLY IF coverage_total matches byte-for-byte, which proves
     the coverage really is the same pooled data
Anything unresolved is left EMPTY, never 0.

TWO QUANTITIES THAT LOOK ALIKE AND ARE NOT:
  * `psites` -- total per-nt P-sites over the pack's own universe. THE LIBRARY DEPTH. Read from
    target_counts.npy, the same quantity P3 reports.
  * `f0_frac` -- pooled frame-0 share at ANNOTATED START CODONS, from RiboCode metaplots. A quality
    ratio, not a depth. The QC panel's `psites_at_start_window` is ~0.5% of `psites` and must never
    be substituted for it (its ratio to real depth swings 0.33-0.65% across the Chothani tissues,
    enough to invert the ES/Fat ordering).

RECIPE STATUS IS CARRIED PER ROW, because the project mixes substrates and a census that hides that
invites cross-substrate comparison. `on_recipe=False` means the pack predates the 2026-08-12 final
recipe (no ncRNA + cross-gene filter), so its depth is inflated relative to on-recipe rows.

  python3 scripts/build_dataset_census.py
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys

import numpy as np

NEW = pathlib.Path(__file__).resolve().parents[1]

# (census label, QC label or None, pack dir or None, role, on_recipe, note)
# QC label is how the dataset appears in riboseq_qc_values.json; None means no metaplots table.
ROWS = [
    ("ES",                  "ES",                       "data/packed_canon_ES",          "train",   True,  ""),
    ("Fat",                 "Fat",                      "data/packed_canon_Fat",         "train",   True,  ""),
    ("Fibroblast",          "Fibroblast",               "data/packed_canon",             "train",   True,  ""),
    ("HA_EC",               "HA_EC",                    "data/packed_canon_HA_EC",       "train",   True,  ""),
    ("HCAEC",               "HCAEC",                    "data/packed_canon_HCAEC",       "train",   True,  ""),
    ("Hepatocytes",         "Hepatocytes",              "data/packed_canon_Hepatocytes", "holdout", True,  "LOTO holdout, not trained on"),
    ("HUVEC",               "HUVEC",                    "data/packed_canon_HUVEC",       "train",   True,  ""),
    ("VSMC",                "VSMC",                     "data/packed_canon_VSMC",        "train",   True,  ""),
    ("Brain",               "Brain (dropped)",          "data/packed_Brain",             "dropped", False,
     "EXCLUDED on period_obs 0.044. Depth is OFF-RECIPE (pack 2026-07-11); frame-0 is on-recipe."),
    ("iPSC-CM (Ruiz-Orera)", "Ruiz-Orera (human)",      "data/packed_heldout_human_ruizorera", "heldout", False,
     "RECIPE STATUS UNKNOWN, not confirmed off-recipe: BAMs + source hd5 deleted and no align script "
     "survives, so nothing is left to check. Pack predates the 2026-08-12 adoption. FASTQs "
     "re-downloaded from ENA 2026-08-15 (md5-verified) to settle it by rebuilding. Task #92."),
    ("THP-1 GSE208041",     "THP-1 GSE208041",          "data/packed_heldout_human_gse208041_canon", "heldout", True, ""),
    ("CAR-T GSE304796",     "CAR-T GSE304796",          "data/packed_heldout_human_cart_canon", "heldout", True, ""),
    ("THP-1 GSE39561",      None,                       "data/packed_heldout_human_gse39561_canon", "dropped", True,
     "EXCLUDED: zero periodic read lengths. 0 P-sites IS the measurement. metaplots wrote no table, "
     "so frame-0 comes from prose (40.7%, results.md:2813) -- same quantity as f0_frac."),
    ("Wang liver (mouse)",  "Wang liver (mouse)",       "data/liver3x3/pool_wang_canon", "heldout", True,  ""),
    ("GSE120762 NT (mouse)", "GSE120762 NT (mouse)",    "data/packed_heldout_mouse_gse120762_nt_canon", "heldout", True, ""),
    ("GSE120762 LPS (mouse)", "GSE120762 LPS (mouse)",  "data/packed_heldout_mouse_gse120762_lps_canon", "heldout", True, ""),
    ("GSE155087 T-cell (mouse)", "GSE155087 T-cell (mouse)", "data/packed_heldout_mouse_gse155087_tcell_canon", "heldout", True, ""),
    ("B721.221",            "B721.221",                 None,                            "pgx",     True,
     "ON-RECIPE (align_b721_ribo.sbatch: EndToEnd + mm1 + filter_tx_heldout; 7/7 BAMs carry the "
     "filter's @PG). A 2026-08-15 audit wrongly called this off-recipe by scraping logs/. "
     "Depth is the MEASURED reference (Ouspenskaia), not a pack: "
     "packed_heldout_human_b721 is COVERAGE-ONLY BY DESIGN -- Ribo-free input for standalone "
     "prediction, target_counts deliberately 0. Do not read that 0 as a measurement."),
    ("HBL-1",               "HBL-1",                    None,                            "pgx",     True,  "no pack in this project"),
]
# Frame-0 values with no metaplots table, taken from prose. Source recorded so it is auditable.
PROSE_F0 = {"THP-1 GSE39561": (0.407, "results.md:2813 (pooled frame concentration at annotated CDS starts)")}

# Depth that does not come from a pack's target_counts. Both entries exist because a bare 0 would be
# read as "no signal", and only ONE of the project's two zeros means that:
#   GSE39561  psites = 0  -> A MEASUREMENT. Zero periodic read lengths; the library has no usable
#                            Ribo-seq signal. Recorded as 0 above, correctly.
#   B721.221  pack has 0  -> NOT a measurement. Its pack is coverage-only by design (the Ribo-free
#                            input for standalone prediction). Its real depth is the measured
#                            reference it is SCORED AGAINST, which lives outside the pack.
EXTERNAL_DEPTH = {
    "B721.221": (327_000_000, "Ouspenskaia measured Ribo-seq footprints, the A4 reference "
                              "(figures/A4_b721_ground_truth/FIGURE_DATA_INPUTS.md); the pack itself "
                              "is coverage-only")}


# Packs whose RNA facts must be read from a pack OTHER than the one in the ROWS table. B721's ROWS
# entry has pack=None because its DEPTH is the external measured reference, but its pack exists and
# is coverage-only -- RNA is precisely what it holds, so the columns would otherwise be blank for the
# one dataset where RNA is the entire input.
RNA_PACK_OVERRIDE = {"B721.221": "data/packed_heldout_human_b721"}
# Counts recoverable from a run->tissue map when no pack records them. Brain's pack predates the
# provenance convention and has no union twin.
RNA_FROM_SRR_MAP = {"Brain": ("data/loto_rnaseq_srr_tissue.tsv", "Brain")}


def rna_facts(pack_rel):
    """-> (n_rna_libraries, rna_coverage_total, source). Empty count rather than a misleading 0."""
    if not pack_rel:
        return None, None, ""
    d = NEW / pack_rel
    total = None
    cn = d / "coverage_norm.json"
    cj = json.loads(cn.read_text()) if cn.exists() else {}
    total = cj.get("coverage_total")
    # 1. provenance: the list of files actually pooled
    pv = d / "provenance.json"
    if pv.exists():
        try:
            ri = json.loads(pv.read_text()).get("rna_coverage_inputs")
            if isinstance(ri, list) and ri:
                return len(ri), total, "provenance.json rna_coverage_inputs"
        except Exception:
            pass
    # 2. coverage_norm, when it is not the misleading 0
    n = cj.get("n_rna_samples")
    if n:
        return n, total, "coverage_norm.json n_rna_samples"
    # 3. the union twin, ONLY if coverage_total matches exactly
    twin = None
    if d.name.startswith("packed_canon"):
        twin = NEW / "data" / d.name.replace("packed_canon", "packed_union", 1)
    elif d.name.endswith("_canon"):
        twin = NEW / "data" / d.name[:-len("_canon")]
    if twin and (twin / "coverage_norm.json").exists():
        tj = json.loads((twin / "coverage_norm.json").read_text())
        if tj.get("n_rna_samples") and total is not None and tj.get("coverage_total") == total:
            return tj["n_rna_samples"], total, (
                f"{twin.name} (coverage_total identical, so the same pooled RNA)")
    return None, total, "UNRESOLVED -- left empty, not 0"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(NEW / "results/dataset_census"))
    a = ap.parse_args()
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)

    qc = {d["label"]: d for d in json.loads(
        (NEW / "figures/S_riboseq_qc/riboseq_qc_values.json").read_text())["datasets"]}

    rows = []
    for label, qlab, pack, role, on_recipe, note in ROWS:
        r = {"dataset": label, "role": role, "on_recipe": on_recipe, "note": note,
             "psites": None, "n_tx": None, "f0_frac": None, "f0_source": None,
             "n_libraries": None, "mode_read_length": None, "pack": pack or "",
             "psites_source": None,
             "n_rna_libraries": None, "rna_coverage_total": None, "rna_source": None,
             "accession_rna": None}
        if pack:
            f = NEW / pack / "target_counts.npy"
            if f.exists():
                r["psites"] = int(np.asarray(np.load(f, mmap_mode="r")).sum())
                r["n_tx"] = sum(1 for _ in open(NEW / pack / "tx_order.txt"))
        if label in EXTERNAL_DEPTH:
            r["psites"], r["psites_source"] = EXTERNAL_DEPTH[label]
        elif r["psites"] is not None:
            r["psites_source"] = f"{pack}/target_counts.npy"
        nr, rtot, rsrc = rna_facts(RNA_PACK_OVERRIDE.get(label, pack))
        if nr is None and label in RNA_FROM_SRR_MAP:
            rel, tis = RNA_FROM_SRR_MAP[label]
            f = NEW / rel
            if f.exists():
                n = sum(1 for ln in open(f) if ln.rstrip("\n").endswith("\t" + tis))
                if n:
                    nr, rsrc = n, f"{pathlib.Path(rel).name} (run->tissue map)"
        r["n_rna_libraries"], r["rna_coverage_total"], r["rna_source"] = nr, rtot, rsrc
        if qlab and qlab in qc:
            d = qc[qlab]
            r.update(f0_frac=d["f0_frac"], f0_source="metaplots (riboseq_qc_values.json)",
                     n_libraries=d["n_libraries"], mode_read_length=d["mode_read_length"])
        elif label in PROSE_F0:
            r["f0_frac"], r["f0_source"] = PROSE_F0[label]
        rows.append(r)

    # accession_rna from the registry: Chothani's RNA is GSE182372, a DIFFERENT series from the
    # Ribo-seq GSE182371, and citing the wrong one has already happened in this project.
    try:
        sys.path.insert(0, str(NEW / "scripts"))
        from dataset_labels import _rows as _reg_rows
        reg = {r["display"]: r for r in _reg_rows()}
        by_key = {r["key"]: r for r in _reg_rows()}
        ALIAS = {"iPSC-CM (Ruiz-Orera)": "human_ruizorera", "THP-1 GSE208041": "gse208041",
                 "CAR-T GSE304796": "cart", "THP-1 GSE39561": "gse39561", "B721.221": "B721",
                 "Wang liver (mouse)": "wang", "GSE120762 NT (mouse)": "gse120762_nt",
                 "GSE120762 LPS (mouse)": "gse120762_lps",
                 "GSE155087 T-cell (mouse)": "gse155087"}
        for r in rows:
            k = ALIAS.get(r["dataset"])
            if k is None and r["dataset"] in ("ES", "Fat", "Fibroblast", "HA_EC", "HCAEC",
                                              "Hepatocytes", "HUVEC", "VSMC", "Brain"):
                k = "chothani"
            src = by_key.get(k) if k else None
            if src:
                a = src.get("accession_rna", "").strip()
                r["accession_rna"] = a if a and a != "-" else ""
    except Exception as e:                                    # noqa: BLE001
        print(f"  WARNING: accession_rna not resolved: {e}")

    cols = ["dataset", "role", "on_recipe", "psites", "psites_source", "n_tx", "f0_frac",
            "f0_source", "n_libraries", "mode_read_length",
            "n_rna_libraries", "rna_coverage_total", "rna_source", "accession_rna",
            "pack", "note"]
    with open(out / "dataset_census.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r[k] is None else r[k]) for k in cols})
    (out / "dataset_census.json").write_text(json.dumps(
        {"quantities": {
            "psites": "total per-nt P-sites over the pack's own universe = LIBRARY DEPTH "
                      "(same quantity as P3_values.json psites)",
            "f0_frac": "pooled frame-0 share at ANNOTATED START CODONS from RiboCode metaplots "
                       "-- a quality ratio, NOT a depth",
            "n_rna_libraries": "count of RNA coverage files pooled into the pack. See rna_source "
                               "per row -- it is NOT uniformly from one place, and every _canon "
                               "pack reports n_rna_samples=0 in its own coverage_norm.json despite "
                               "holding real coverage (they reused their union twin's coverage "
                               "verbatim). Unresolved rows are EMPTY, never 0.",
            "rna_coverage_total": "RAW summed RNA coverage over the pack's own universe. VERIFIED "
                                  "equal to coverage.npy.sum() on int32 per-nt counts -- it is a "
                                  "TOTAL, not normalised and not depth-corrected, so it may be "
                                  "labelled as depth. The pack's global_mean_coverage is exactly "
                                  "this value / sum_L, i.e. the per-nt mean of the same quantity. "
                                  "Units are per-nt read depth summed over all transcript "
                                  "positions, so it is NOT comparable to psites (a count of "
                                  "footprints); plot them on separate axes.",
            "accession_rna": "the RNA-seq series, which for Chothani is GSE182372 -- a DIFFERENT "
                             "series from the Ribo-seq GSE182371. Citing the wrong one has already "
                             "happened in this project (see DATASET_LABELS.md).",
            "on_recipe": "False = pack predates the 2026-08-12 final recipe (no ncRNA + cross-gene "
                         "filter), so its depth is inflated relative to on-recipe rows"},
         "do_not_substitute": "riboseq_qc's psites_at_start_window is ~0.5% of psites and is NOT a "
                              "depth proxy; its ratio to real depth swings 0.33-0.65% across the "
                              "Chothani tissues, enough to invert the ES/Fat ordering",
         "rows": rows}, indent=2))

    print(f"  {'dataset':<26}{'role':<9}{'rec':>5}{'P-sites':>15}{'frame-0':>9}  note")
    for r in rows:
        ps = f"{r['psites']:,}" if r["psites"] is not None else "--"
        f0 = f"{r['f0_frac']*100:.1f}%" if r["f0_frac"] is not None else "--"
        print(f"  {r['dataset']:<26}{r['role']:<9}{'on' if r['on_recipe'] else 'OFF':>5}"
              f"{ps:>15}{f0:>9}  {r['note'][:44]}")
    miss = [r["dataset"] for r in rows if r["psites"] is None and r["role"] != "pgx"]
    if miss:
        print(f"\n  no depth available: {miss}")
    print(f"\n  wrote {out}/dataset_census.{{tsv,json}}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
