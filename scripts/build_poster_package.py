#!/usr/bin/env python3
"""Assemble a self-contained, downloadable poster package: figures + the data to remake them.

The point is that someone can download ONE directory, open the figures, edit a generator, re-run it,
and get the same figure back -- without the 20 TB of BAMs and packs behind the project. So each panel
ships its PDF/PNG, its `*_values.json` (the exact numbers plotted, written by the generator AFTER
drawing so they cannot drift), its TSV export, its `FIGURE_DATA_INPUTS.md` provenance note, and its
generator script.

WHAT IS NOT SELF-CONTAINED, AND WHY THAT IS STATED RATHER THAN HIDDEN. The generators read absolute
paths into the project (packs, search results, dropin dumps). Re-running one on a different machine
needs the project tree. What IS portable is the values JSON/TSV: every plotted number is there, so a
figure can be restyled, recombined or replotted from the package alone. The README says so per panel.

THE README IS NOT WRITTEN HERE, AND MUST NOT BE. This script starts with `shutil.rmtree(out)`, so
anything hand-written inside `poster_package/` is destroyed on every refresh -- which is what nearly
happened to the 18 KB assembly README. It now lives in the repo at `docs/POSTER_PACKAGE_README.md`
and is COPIED in, so a rebuild refreshes it instead of deleting it. Edit it there.

Run: python3 scripts/build_poster_package.py [--out <dir>] [--model attn]
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

NEW = Path(__file__).resolve().parents[1]
FIG = NEW / "figures"


def dataset_labels_md():
    """-> the old-label -> new-label cross-reference, rendered from the registry.

    Generated rather than transcribed: `data/dataset_registry.tsv` is the single source of truth, and
    a hand-copied table in the package would be one more thing to drift. Internal keys are what
    appear in directory names, pack paths and every values JSON; display labels are what a figure,
    table or poster should print.
    """
    sys.path.insert(0, str(NEW / "scripts"))
    from dataset_labels import CHOTHANI_TISSUES, _rows, display  # noqa: E402

    L = ["# Dataset identifiers: internal key -> display label",
         "",
         "The labels in this project grew as a mix of surnames (`janich`, `wang`), GEO accessions",
         "(`gse243134`), cell types (`Hepatocytes`, `HUVEC`) and product names (`cart`). That is not a",
         "style problem. The surname-labelled datasets were surname-labelled because **no accession was",
         "recorded anywhere**, and a surname breaks exactly where GSE243134 sits: a Genome Biology paper",
         "with multiple first authors. Two datasets are also **both THP-1** (GSE208041, GSE39561), so a",
         "cell-type label alone is ambiguous in this project's own data.",
         "",
         "**Display format:** `<species>_<tissue/cell type>_<primary accession>`. The accession is the",
         "unambiguous key, the species+tissue prefix is what a human scans for, and the pair resolves the",
         "two-THP-1 collision that neither component resolves alone.",
         "",
         "**Internal keys are deliberately NOT renamed.** They are load-bearing in directory names, pack",
         "paths, script arguments and every `*_values.json`. Mapping happens at render time",
         "(`scripts/dataset_labels.py`), which buys consistency without a rename whose blast radius is",
         "the whole project. **Use the display label on the poster; expect the internal key in filenames.**",
         "",
         "| internal key (in paths/JSONs) | display label (for the poster) | accession | species | tissue / cell type |",
         "|---|---|---|---|---|"]
    for r in sorted(_rows(), key=lambda r: (r["species"], r["display"])):
        acc = r["accession_ribo"] if r["accession_ribo"] not in ("", "-") else (r["sra"] or "-")
        L.append(f"| `{r['key']}` | **{r['display']}** | {acc} | {r['species']} | {r['tissue']} |")
    L += ["",
          "## Chothani (GSE182371) expands per tissue",
          "",
          "One study supplied 9 cell types, so its row carries a `<tissue>` placeholder and the tissue",
          "name is substituted. A bare Chothani tissue name resolves directly.",
          "",
          "| internal key | display label | role |",
          "|---|---|---|"]
    for t in CHOTHANI_TISSUES:
        role = ("**dropped from training** (low periodicity, period_obs 0.044)" if t == "Brain"
                else "**held out** (LOTO test)" if t == "Hepatocytes" else "training")
        L.append(f"| `{t}` | **{display(t)}** | {role} |")
    L += ["",
          "## Accession pitfalls that have already caused errors here",
          "",
          "- **Chothani is GSE182371** (Ribo-seq). `GSE182372` is its RNA-seq and `GSE182377` the",
          "  SuperSeries. A sibling project cites **GSE183670**, which is an unrelated lung",
          "  adenocarcinoma study -- do not copy that citation.",
          "- **GSE304796 (CAR-T)** is labelled `library_strategy: RNA-Seq` in GEO. It is Ribo-seq.",
          "- **PRJEB65856 (iPSC-CM)** is ENA, not GEO: a GSE-only naming scheme cannot express it.",
          "- **PXD008733 (B721.221)** combines TWO studies -- Ribo-seq from Ouspenskaia, RNA-seq from",
          "  Sarkizova -- which no single accession expresses.",
          "- **GSE120762** covers both the untreated and LPS macrophage arms; the treatment, not the",
          "  accession, distinguishes them.",
          "",
          "Generated by `scripts/build_poster_package.py` from `data/dataset_registry.tsv`.",
          "Re-run `python3 scripts/dataset_labels.py` to print the same table at the terminal.",
          ""]
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(NEW / "poster_package"))
    ap.add_argument("--model", default="attn", choices=["attn", "mamba4"],
                    help="which model the poster emphasises; recorded in the manifest")
    a = ap.parse_args()
    out = Path(a.out)
    if out.exists():
        shutil.rmtree(out)
    (out / "figures").mkdir(parents=True)

    manifest = {"emphasised_model": a.model, "figures": {}}
    # `_`-prefixed directories are RETIRED figures, kept in the repo for the record but not shipped.
    # S_rescore_substrate was the first: its claim was a tryptic-vs-HLA contrast and the tryptic
    # dataset (A549) was removed, so the panel has no second half to compare against.
    for d in sorted(p for p in FIG.iterdir() if p.is_dir() and not p.name.startswith("_")):
        gen = sorted(d.glob("make_*.py")) or sorted(d.glob("*.py"))
        arts = sorted(d.glob("*.pdf")) + sorted(d.glob("*.png"))
        # ALL json, not just *values*.json. The narrower glob silently dropped
        # E3_all_rna_cells.json -- a portable data file that simply did not have "values" in its
        # name -- while shipping its .tsv sibling, so the package looked complete. Any JSON a
        # generator writes into a figure dir is portable data by construction.
        vals = sorted(set(d.glob("*.json")) | set(d.glob("*.tsv")))
        docs = sorted(d.glob("FIGURE_DATA_INPUTS.md"))
        if not arts and not vals:
            continue
        dst = out / "figures" / d.name
        dst.mkdir(parents=True, exist_ok=True)
        n = 0
        for f in arts + vals + docs + gen:
            if f.is_file():
                shutil.copy2(f, dst / f.name)
                n += 1
        manifest["figures"][d.name] = {
            "outputs": [f.name for f in arts],
            "values": [f.name for f in vals],
            "generator": [f.name for f in gen],
            "has_provenance_doc": bool(docs),
            "n_files": n,
        }
        print(f"  {d.name:<28} {n:>3} files ({len(arts)} figure, {len(vals)} data)")

    # Project-level documents a reader needs to interpret any of it.
    # docs/POSTER_MSG_*.md are the REPLIES to the poster session's requests. They shipped only by
    # being pasted into a prompt, which meant a rsync could land new files with no statement of
    # what changed or which caveats attach to them. Globbed, not listed, so a new reply ships
    # automatically instead of being remembered.
    msgs = sorted(str(q.relative_to(NEW)) for q in (NEW / "docs").glob("POSTER_MSG_*.md"))
    for rel in (*msgs,
                "docs/POSTER_SESSION_HANDOFF.md",
                "docs/POSTER_LAYOUT.md", "docs/STATUS_CURRENT_VS_ARCHIVED.md",
                "docs/PIPELINE_POLICY.md", "figures/README.md", "results.md", "methods.md"):
        src = NEW / rel
        if src.exists():
            shutil.copy2(src, out / ("PROJECT_" + Path(rel).name))
            manifest.setdefault("project_docs", []).append("PROJECT_" + Path(rel).name)

    # The assembly README, from the repo -- NOT written here, so a rebuild refreshes it rather than
    # deleting a hand-written file (see module docstring).
    src = NEW / "docs/POSTER_PACKAGE_README.md"
    if src.exists():
        shutil.copy2(src, out / "README.md")
        manifest["readme_source"] = "docs/POSTER_PACKAGE_README.md"
        print("  README.md            <- docs/POSTER_PACKAGE_README.md")
    else:
        print("  WARNING: docs/POSTER_PACKAGE_README.md missing -- package ships with NO README")

    # Portable DATA that is not a figure. The builder copies figures/ and a fixed list of docs, so a
    # results/ artifact silently does not ship -- which is what happened to the held-out P-site
    # totals the poster's dataset-census panel needed.
    for rel in ("results/dataset_census/dataset_census.tsv",
                "results/dataset_census/dataset_census.json",
                "results/heldout_psite_totals.json",
                # B9 plots a subset of these strata; the file carries all of them, including the
                # per-model EXCLUSIVE rows that no figure draws. The poster needs the attn numbers.
                "results/merged_liver_ribocode/overcall/overcall_meta.json"):
        src = NEW / rel
        if src.exists():
            shutil.copy2(src, out / Path(rel).name)
            manifest.setdefault("data_files", []).append(Path(rel).name)
            print(f"  {Path(rel).name:<28} <- {rel}")

    # Dataset identifier cross-reference + the registry it is generated from.
    try:
        (out / "DATASET_LABELS.md").write_text(dataset_labels_md())
        shutil.copy2(NEW / "data/dataset_registry.tsv", out / "dataset_registry.tsv")
        manifest["dataset_labels"] = ["DATASET_LABELS.md", "dataset_registry.tsv"]
        print("  DATASET_LABELS.md    <- data/dataset_registry.tsv (generated)")
    except Exception as e:                                    # noqa: BLE001
        print(f"  WARNING: dataset label table not built: {e}")

    (out / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))
    tot = sum(v["n_files"] for v in manifest["figures"].values())
    print(f"\n  {len(manifest['figures'])} figure folders, {tot} files -> {out}")
    print(f"  project docs: {len(manifest.get('project_docs', []))}")
    return manifest


if __name__ == "__main__":
    main()
