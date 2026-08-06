#!/usr/bin/env python3
"""CLASS-AWARE rescoring: the correct way to combine MS2Rescore with cryptic-peptide FDR. Global mokapot
(db_comparison_rescored.md) optimizes canonical-dominated target/decoy separation and interleaves novel
targets/decoys, collapsing the novel-class FDR. Fix: re-train mokapot on the NOVEL CLASS ONLY (novel_t vs
novel_d) using the same MS2Rescore features (79 `rescoring:` cols), so the model learns what separates real
novel peptides from novel decoys. Then novel peptides at 1% mokapot FDR = the class-aware count. Run in the
ms2rescore env (mokapot 0.10 + pandas). Compares model vs null; reference: raw hyperscore 10/13, global
rescored 1/6."""
from __future__ import annotations

import ast
import glob
import sys
from pathlib import Path

import mokapot
import numpy as np
import pandas as pd

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")


def parse_accs(s):
    try:
        return tuple(ast.literal_eval(s))
    except Exception:
        return tuple(a.strip(" '\"[]") for a in str(s).split(",") if a.strip(" '\"[]"))


def novel_class(df):
    """keep rank-1 PSMs that are novel_t (all nuORF| targets) or novel_d (all REV_nuORF| decoys)."""
    df = df[df["rank"].astype(str).isin(["1", "1.0"])].copy()
    accs = df["protein_list"].map(parse_accs)
    is_dec = accs.map(lambda a: all(x.startswith("REV_") for x in a) and len(a) > 0)
    novel_t = accs.map(lambda a: len(a) > 0 and all(x.startswith("nuORF|") for x in a))
    novel_d = accs.map(lambda a: len(a) > 0 and all(x.startswith("REV_nuORF|") for x in a))
    keep = novel_t | novel_d
    out = df[keep].copy()
    out["is_target"] = novel_t[keep].values          # True = novel_t, False = novel_d
    return out


def classaware_qcount(db, fdr=0.01):
    files = sorted(glob.glob(str(NEW / f"proteogenomics/data/A549_pilot/rescore/{db}/*.psms.tsv")))
    df = pd.concat([pd.read_csv(f, sep="\t", low_memory=False) for f in files], ignore_index=True)
    nov = novel_class(df)
    feat = [c for c in nov.columns if c.startswith("rescoring:")]
    # numeric, drop all-NaN/constant features, fill remaining NaN
    X = nov[feat].apply(pd.to_numeric, errors="coerce")
    X = X.loc[:, X.notna().any() & (X.nunique() > 1)].fillna(0.0)
    specid = (nov["run"].astype(str) + "::" + nov["spectrum_id"].astype(str)).reset_index(drop=True)
    nov = pd.concat([nov[["peptidoform", "is_target"]].reset_index(drop=True),
                     specid.rename("specid"), X.reset_index(drop=True)], axis=1)
    nt, nd = int(nov.is_target.sum()), int((~nov.is_target).sum())
    print(f"[{db}] novel-class PSMs: {nt:,} target / {nd:,} decoy; {X.shape[1]} features; "
          f"{nov.specid.nunique():,} unique spectra", file=sys.stderr)
    ds = mokapot.LinearPsmDataset(
        psms=nov, target_column="is_target",
        spectrum_columns=["specid"], peptide_column="peptidoform",
        feature_columns=list(X.columns))
    # relaxed train_fdr: the novel class is too sparse to bootstrap at the default 1%
    for tfdr in (0.05, 0.1, 0.25):
        try:
            results, _ = mokapot.brew([ds], model=mokapot.PercolatorModel(train_fdr=tfdr))
            peps = results.peptides
            qcol = [c for c in peps.columns if "q-value" in c.lower()][0]
            print(f"[{db}] class-aware mokapot trained at train_fdr={tfdr}", file=sys.stderr)
            return int((peps[qcol] <= fdr).sum()), nt, nd
        except Exception as e:
            print(f"[{db}] train_fdr={tfdr} failed: {type(e).__name__}: {str(e)[:90]}", file=sys.stderr)
    return None, nt, nd   # untrainable -- novel class too sparse for rescoring


def main():
    dbsize = {}
    for db in ("model", "null"):
        cm = NEW / f"proteogenomics/data/A549_pilot/db_fresh/{db}_class_map.tsv"
        dbsize[db] = (sum(1 for _ in open(cm)) - 1) if cm.exists() else 0
    res = {}
    for db in ("model", "null"):
        res[db] = classaware_qcount(db)
    L = ["# A549 (fresh universe): model vs null -- CLASS-AWARE MS2Rescore (novel-class mokapot)\n",
         "Novel peptides at 1% FDR from mokapot re-trained on the NOVEL class only (novel_t vs novel_d) using "
         "the MS2Rescore features. Contrast: raw hyperscore 10/13 (2.31x); global-rescored 1/6 (0.50x).\n",
         "| DB | novel ORFs | novel peptides @1% (class-aware) | rate (pep/1k ORFs) | novel-class PSMs (t/d) |",
         "|----|---:|---:|---:|---:|"]
    for db in ("model", "null"):
        n_pep, nt, nd = res[db]
        n = dbsize[db]
        cell = "**untrainable**" if n_pep is None else f"**{n_pep:,}**"
        rate = "--" if n_pep is None else f"{1000*n_pep/n if n else 0:.2f}"
        L.append(f"| {db} | {n:,} | {cell} | {rate} | {nt:,}/{nd:,} |")
    mp = res["model"][0]; npn = res["null"][0]
    L += ["\n## Verdict (class-aware rescoring)\n"]
    if mp is None or npn is None:
        L.append("- **Class-aware rescoring is UNTRAINABLE**: the novel class is too sparse (~10-25 real "
                 "novel peptides among ~35k noise + ~26k decoys) for mokapot's semi-supervised bootstrap, "
                 "even at a relaxed train_fdr. Rescoring needs a reasonably populated target class; a "
                 "tryptic whole-proteome yields far too few cryptic-ORF peptides.")
    else:
        rm = 1000*mp/dbsize["model"] if dbsize["model"] else 0
        rn = 1000*npn/dbsize["null"] if dbsize["null"] else 0
        L += [f"- Model {mp:,} novel peptides ({rm:.2f}/1k); null {npn:,} ({rn:.2f}/1k).",
              f"- **model/null rate ratio: {rm/rn:.2f}x** (raw hyperscore 2.31x; global-rescored 0.50x)."]
    L.append("\n## Overall MS2Rescore verdict on A549 (tryptic whole-proteome)\n")
    L.append("| approach | model | null | model/null | note |")
    L.append("|---|---:|---:|---:|---|")
    L.append("| raw hyperscore class-FDR | 10 | 13 | 2.31x | the correct method here |")
    L.append("| global MS2Rescore + class-FDR | 1 | 6 | 0.50x | mokapot interleaves novel t/d |")
    ca = "untrainable" if (mp is None or npn is None) else f"{mp}/{npn}"
    L.append(f"| class-aware MS2Rescore | {ca.split('/')[0] if '/' in ca else ca} | "
             f"{ca.split('/')[1] if '/' in ca else ''} | -- | novel class too sparse to train |")
    out = NEW / "proteogenomics/data/A549_pilot/db_comparison_rescored_classaware.md"
    out.write_text("\n".join(L) + "\n")
    print("\n".join(L)); print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
