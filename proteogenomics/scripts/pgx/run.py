#!/usr/bin/env python3
"""End-to-end `pgx` orchestrator: predicted signal -> calibrated ORF calls -> databases -> search
-> FDR-controlled table. Every input is a CLI argument; no paths are baked in.

    universe   salmon quant.sf -> expressed universe (TPM >= 1, pc + lncRNA, no chrM, length cap)
    calibrate  Poisson/theta sweep, CDS-anchored          -> theta.json
    call       two-arm RiboCode calling on the prediction -> calls_<arm>.tsv
    extend     N-terminal extension scan (own entry class) -> extensions.tsv
    dbs        db_gencode / db_model_<arm> / db_null_atg / db_null_nc
    search     MSFragger, tryptic or nonspecific, hash-cached
    report     global FDR canonical + class-specific FDR novel -> table.md

Layout under --out:  universe/ calib/ calls/ extensions/ db/ search/<arm>/ report/

Every step is skip-if-done, so a rerun after a failure resumes rather than repeats. MSFragger is
the only step that is submitted rather than run inline (hours, 160 GB); the rest are minutes.
Run the whole thing on a compute node -- the calibration sweep alone is several RiboCode passes.

  python -m pgx.run --species mouse --label BMDM --profiles p.npz --salmon q1.sf q2.sf \
      --mzml '<dir>/BMDM_Rep*_F*.mzML' --enzyme tryptic --out <dir> --submit-searches
  python -m pgx.run ... --report-only          # after the searches land
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

from .refs import require, tool_path

DB_ARMS = ("gencode", "null_atg", "null_nc")     # + model_<arm>, added at runtime


def sh(cmd, **kw):
    print("+ " + " ".join(str(c) for c in cmd), file=sys.stderr)
    subprocess.run([str(c) for c in cmd], check=True, **kw)


def done(path):
    p = Path(path)
    return p.exists() and p.stat().st_size > 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True, choices=["human", "mouse"])
    ap.add_argument("--label", required=True, help="population / dataset name (row label)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--profiles", required=True, help="pred_profiles.npz from dump_pred_profiles")
    ap.add_argument("--salmon", nargs="*", default=[], help="quant.sf files (mean-TPM merged)")
    ap.add_argument("--universe-fa", default=None, help="prebuilt universe FASTA (skips --salmon)")
    ap.add_argument("--mzml", default=None, help="mzML directory or glob")
    ap.add_argument("--enzyme", default="tryptic", choices=["tryptic", "nonspecific"])
    ap.add_argument("--template", default=None,
                    help="MSFragger params template: either a path, or a bare basename resolved "
                         "against proteogenomics/msfragger/. Defaults to the FROZEN template for "
                         "the enzyme (pgx.search.DEFAULT_TEMPLATE). Set it when the instrument "
                         "differs from the default -- e.g. A549 is TMT-labelled, so it needs "
                         "fragger_a549_tmt_frozen.params rather than the label-free tryptic "
                         "default. Pass an unfrozen template only to reproduce a pre-2026-08-04 run.")
    # universe
    ap.add_argument("--min-tpm", type=float, default=1.0)
    ap.add_argument("--max-len", type=int, default=10000)
    ap.add_argument("--biotype-field", default="transcript_type",
                    choices=["transcript_type", "gene_type"],
                    help="which biotype the {protein_coding, lncRNA} filter applies to. "
                         "transcript_type (default) matches every existing universe in this "
                         "project; gene_type also admits NMD / retained_intron / TEC isoforms of "
                         "coding genes, where non-canonical ORFs are enriched. See pgx.universe.")
    # calling / calibration
    ap.add_argument("--cds-recall", type=float, default=0.90)
    ap.add_argument("--recall-mode", default="relative", choices=["relative", "absolute"])
    ap.add_argument("--thetas", default=None)
    ap.add_argument("--pval", type=float, default=0.05,
                    help="RiboCode per-ORF significance cutoff. This is the REAL significance "
                         "lever: RiboCode filters on it before BH, so --orf-qvalue at 0.05 is a "
                         "no-op (measured: BH moves p by at most 4e-4 on the surviving set).")
    ap.add_argument("--orf-qvalue", type=float, default=0.05)
    ap.add_argument("--min-aa-call", type=int, default=5, help="RiboCode MIN_AA_LENGTH")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--db-arms", default="standard,poisson",
                    help="which calling arms get their own model database and search. Default is "
                         "BOTH, so every report shows the uncalibrated arm (theta = 1, "
                         "deterministic rounding) next to the Poisson-calibrated one and the "
                         "benefit of calibration is always visible rather than asserted. Pass a "
                         "single arm to halve the MSFragger cost.")
    # extensions
    ap.add_argument("--start-codons", default="near_cognate")
    ap.add_argument("--ext-qvalue", type=float, default=0.05)
    ap.add_argument("--min-ext-aa", type=int, default=5)
    ap.add_argument("--min-nonzero", type=int, default=5)
    # databases
    ap.add_argument("--min-aa", type=int, default=7, help="minimum protein length in the DBs")
    ap.add_argument("--dbs", default="gencode,model,null_atg,null_nc")
    ap.add_argument("--null-starts", default="ATG")
    ap.add_argument("--nc-starts", default="near_cognate")
    ap.add_argument("--shared-search-root", default=None,
                    help="where to put the arms that do NOT depend on the model (gencode, "
                         "null_atg, null_nc). Those depend only on (species, universe, mzML, "
                         "enzyme), so pointing several model runs at one shared root lets the "
                         "content-addressed cache skip them outright. Worth setting: a null_nc "
                         "search over ~2.5M targets takes hours and is byte-identical across "
                         "checkpoints. The shared arms are symlinked back under <out>/search so "
                         "reporting still sees one directory.")
    # control
    ap.add_argument("--submit-searches", action="store_true")
    ap.add_argument("--skip-searches", action="store_true")
    ap.add_argument("--report-only", action="store_true")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    require(a.species)          # fail fast on any missing reference, not mid-pipeline
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    pkg = Path(__file__).resolve().parent            # .../pgx
    scripts = pkg.parent                             # .../proteogenomics/scripts
    py = str(tool_path("python"))
    rpy = str(tool_path("ribocode_python"))
    env = {**os.environ, "PYTHONPATH": f"{scripts}{os.pathsep}{os.environ.get('PYTHONPATH', '')}"}
    M = [py, "-m"]
    want_dbs = {x.strip() for x in a.dbs.split(",") if x.strip()}
    # null_nc is NOT usable with nonspecific digestion. Nonspecific emits every 8-14mer substring,
    # so the ~2.5-3M near-cognate null generates ~2.9 BILLION peptides against MSFragger's ~2 billion
    # hard cap and the search dies with "Too many peptides were generated". Measured on all three
    # immunopeptidome lines (2,857,862,951 peptides, HBL-1). Tryptic runs are unaffected -- the same
    # database size is fine there -- so this is an enzyme limit, not a database-size limit.
    # DECISION (2026-08-04): immunopeptidome runs use null_atg as the null. Splitting the database
    # was considered and declined; null_atg is the operative specificity control, and the cost is
    # that ncStart (peptides an AUG-only null cannot reach) is not measurable under nonspecific.
    if a.enzyme == "nonspecific" and "null_nc" in want_dbs:
        want_dbs.discard("null_nc")
        print("[pgx] dropping null_nc: nonspecific digestion of a ~2.5-3M-sequence near-cognate "
              "null exceeds MSFragger's ~2e9 peptide cap. Using null_atg as the null.",
              file=sys.stderr)
    db_arms = [x.strip() for x in a.db_arms.split(",") if x.strip()]
    bad = [x for x in db_arms if x not in ("standard", "poisson")]
    if bad:
        raise SystemExit(f"unknown calling arm(s) {bad}; known: standard, poisson")

    uni_fa = Path(a.universe_fa) if a.universe_fa else out / "universe" / f"{a.label}_universe.fa"
    calib = out / "calib"
    calls = out / "calls"
    exts = out / "extensions"
    dbd = out / "db"
    searches = out / "search"
    rep = out / "report"

    if not a.report_only:
        # 1. universe
        if not a.universe_fa and (a.force or not done(uni_fa)):
            if not a.salmon:
                raise SystemExit("need --salmon (or --universe-fa) to build the universe")
            sh(M + ["pgx.universe", "--salmon", *a.salmon, "--species", a.species,
                    "--out-prefix", out / "universe" / a.label,
                    "--min-tpm", a.min_tpm, "--max-len", a.max_len,
                    "--biotype-field", a.biotype_field], env=env, cwd=scripts)
        else:
            print(f"[skip] universe -> {uni_fa}", file=sys.stderr)

        # 2. calibrate (each theta is itself skip-if-done inside pgx.calibrate)
        cmd = M + ["pgx.calibrate", "--profiles", a.profiles, "--species", a.species,
                   "--out", calib, "--cds-recall", a.cds_recall, "--recall-mode", a.recall_mode,
                   "--orf-qvalue", a.orf_qvalue, "--min-aa", a.min_aa_call, "--pval", a.pval,
                   "--seed", a.seed, "--universe-tx", str(uni_fa).replace("_universe.fa",
                                                                         "_universe_tx.txt")]
        if a.thetas:
            cmd += ["--thetas", a.thetas]
        sh(cmd, env=env, cwd=scripts)
        # Each arm's own theta comes from call_summary.json below (standard = 1.0, poisson = theta*).
        theta_json = calib / "theta.json"

        # 3. two-arm calling
        if a.force or not all(done(calls / f"calls_{x}.tsv") for x in db_arms):
            sh(M + ["pgx.call_orfs", "--profiles", a.profiles, "--species", a.species,
                    "--out", calls, "--theta-json", theta_json, "--calib-dir", calib,
                    "--orf-qvalue", a.orf_qvalue, "--min-aa", a.min_aa_call,
                    "--pval", a.pval, "--seed", a.seed], env=env, cwd=scripts)
        else:
            print(f"[skip] calls -> {calls}", file=sys.stderr)

        # 4. N-terminal extensions, PER MODEL ARM. Each arm gets its own scan because the test
        # runs on that arm's density: theta* with Poisson sampling, or theta = 1 deterministic.
        # Reusing one arm's extensions for the other would silently mix stringencies.
        summary = json.loads((calls / "call_summary.json").read_text())
        for arm in db_arms:
            d = exts / arm
            if not a.force and done(d / "extensions.tsv"):
                print(f"[skip] extensions/{arm} -> {d}", file=sys.stderr)
                continue
            sh([rpy, "-m", "pgx.extensions", "--profiles", a.profiles,
                "--calls-collapsed", summary[arm]["collapsed"], "--fasta", uni_fa,
                "--species", a.species, "--theta", summary[arm]["theta"],
                "--start-codons", a.start_codons, "--out", d,
                "--ext-qvalue", a.ext_qvalue, "--orf-qvalue", a.orf_qvalue,
                "--min-ext-aa", a.min_ext_aa, "--min-nonzero", a.min_nonzero],
               env=env, cwd=scripts)

        # 5. databases. The model-independent ones are built once; each calling arm then gets its
        # own db_model_<arm> so the report always shows calibrated vs uncalibrated side by side.
        shared_dbs = ",".join(sorted(want_dbs - {"model"}))
        if shared_dbs and (a.force or not done(dbd / "db_summary.json")):
            sh(M + ["pgx.build_dbs", "--species", a.species, "--out", dbd,
                    "--universe-fa", uni_fa, "--dbs", shared_dbs, "--min-aa", a.min_aa,
                    "--null-starts", a.null_starts, "--nc-starts", a.nc_starts],
               env=env, cwd=scripts)
        elif shared_dbs:
            print(f"[skip] shared databases -> {dbd}", file=sys.stderr)
        if "model" in want_dbs:
            for arm in db_arms:
                if not a.force and done(dbd / f"db_model_{arm}.fasta"):
                    print(f"[skip] db_model_{arm}", file=sys.stderr)
                    continue
                sh(M + ["pgx.build_dbs", "--species", a.species, "--out", dbd,
                        "--calls", calls / f"calls_{arm}.tsv",
                        "--extensions", exts / arm / "extensions.tsv", "--arm", arm,
                        "--dbs", "model", "--min-aa", a.min_aa], env=env, cwd=scripts)

    arms = [x for x in DB_ARMS if x in want_dbs]
    if "model" in want_dbs:
        arms += [f"model_{x}" for x in db_arms]

    # Arms whose search depends only on (species, universe, mzML, enzyme) -- not on the model --
    # may live in a shared root so repeated model runs reuse them via the content-addressed cache.
    shared = Path(a.shared_search_root) / a.label if a.shared_search_root else None

    def search_dir(arm):
        return (shared / arm) if (shared and not arm.startswith("model_")) else (searches / arm)

    # 6. searches
    if not a.skip_searches and not a.report_only:
        if not a.mzml:
            raise SystemExit("need --mzml to search (or pass --skip-searches)")
        sb = pkg / "search.sbatch"
        for arm in arms:
            db = dbd / f"db_{arm}.fasta"
            if not db.exists():
                print(f"[warn] no {db}, skipping search", file=sys.stderr)
                continue
            cmd = ["sbatch", "--parsable", "--job-name", f"pgxmsf_{a.label}_{arm}", str(sb)]
            e = {**env, "DB": str(db), "MZML": a.mzml, "OUT": str(search_dir(arm)),
                 "ENZYME": a.enzyme, "PGX_SCRIPTS": str(scripts)}
            if a.template:
                e["TEMPLATE"] = a.template
            if a.submit_searches:
                jid = subprocess.run(cmd, env=e, capture_output=True, text=True, check=True)
                print(f"submitted {arm}: job {jid.stdout.strip()}", file=sys.stderr)
            else:
                print(f"  DB={db} MZML='{a.mzml}' OUT={search_dir(arm)} ENZYME={a.enzyme} "
                      f"PGX_SCRIPTS={scripts} sbatch {sb}")
        if not a.submit_searches:
            print("\n(searches NOT submitted; pass --submit-searches, then rerun with "
                  "--report-only)", file=sys.stderr)
            return

    # 7. report. Symlink any shared arm back under <out>/search so reporting sees one directory.
    searches.mkdir(parents=True, exist_ok=True)
    for arm in arms:
        src, dst = search_dir(arm), searches / arm
        if src != dst and src.exists() and not dst.exists():
            dst.symlink_to(src)
    # An arm counts as searched if EITHER the compact digest (rank1.tsv.gz, what search.sbatch
    # writes now) or raw MSFragger *.tsv is present. Globbing only "*.tsv" made this print
    # "no search output yet" for completed searches, because "rank1.tsv.gz" does not match it.
    def searched(arm):
        d = searches / arm
        return d.exists() and (any(d.glob("rank1.tsv.gz")) or any(d.glob("*.tsv")))

    have = [arm for arm in arms if searched(arm)]
    if not have:
        print("no search output yet; rerun with --report-only once the searches finish",
              file=sys.stderr)
        return
    missing = [arm for arm in arms if arm not in have]
    if missing:
        print(f"[pgx] WARNING: reporting WITHOUT {missing} -- those arms have no search output. "
              "Every table below omits them; that is a gap, not a zero.", file=sys.stderr)
    rep.mkdir(parents=True, exist_ok=True)
    sh(M + ["pgx.report", "--search-root", searches, "--db-dir", dbd, "--arms", ",".join(have),
            "--baseline", "gencode", "--label", a.label,
            "--out", rep / "table.md", "--json-out", rep / "table.json"], env=env, cwd=scripts)


if __name__ == "__main__":
    main()
