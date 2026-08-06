#!/usr/bin/env python3
"""CPAT / CPC2 coding-potential database arms -- the comparator locked decision D3 asks for.

WHY THIS EXISTS. Every existing comparison is model-vs-naive-enumeration. "The model beats a
3-frame null" is a weak claim: the null does no selection at all. The question a reviewer asks is
whether the model beats an ESTABLISHED coding-potential selector. CPAT and CPC2 are the standard
ones, so this builds them as database arms.

CRITICAL DESIGN POINT: these arms draw from the SAME candidate pool as `null_atg` -- the identical
`candidate_orfs` enumeration over the identical expressed universe, with the identical class filter
and min-length. The ONLY thing that differs is the selector:

    null_atg   keep every enumerated ORF
    cpat       keep ORFs with CPAT  coding probability >= --cpat-cut
    cpc2       keep ORFs whose CPC2 label is `coding`
    model      keep ORFs RiboCode calls translated on the PREDICTED profile

So a difference in discovery density is attributable to the selection rule and nothing else. Scoring
the TRANSCRIPTS instead (the usual CPC2 use) would be the wrong granularity -- a uORF on a coding
mRNA sits on a transcript any tool calls coding, which would make the comparison meaningless.

    python -m pgx.coding_potential --universe-fa <fa> --species human --out <dir> \
        [--tools cpat,cpc2] [--cpat-cut 0.364] [--min-aa 7] [--starts ATG]

Writes: orfs.nt.fa (the shared candidate pool), cpat.tsv / cpc2.txt (raw scorer output),
db_cpat.fasta / db_cpc2.fasta (+ class maps), and summary.json.

CPAT's default human coding cutoff is 0.364 (Wang et al. 2013, the value at which the human
logit model's TPR/FPR cross). It is exposed so the arm can be re-thresholded without re-scoring.
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

from . import rc_io
from .build_dbs import NOVEL_CLASSES, _add, load_canonical, write_db
from .refs import species_refs
from .seqtools import candidate_orfs, cds_relationship, iter_fasta, orf_class, translate

CPAT_DEFAULT_CUT = 0.364          # Wang 2013 human cutoff


def enumerate_pool(uni_fa, cds_of, biotype_of, min_aa, starts, canon_seqs):
    """The candidate pool -- byte-for-byte the same enumeration null_atg uses.

    canon_seqs is REQUIRED: build_dbs._add drops any ORF whose protein is identical to a canonical
    GENCODE protein, so the null's candidate set excludes them. Omitting that filter here made the
    pool a superset by 19 proteins on a 1,000-transcript test -- small, but it would mean CPAT/CPC2
    were scoring a different candidate set than the null they are compared against.
    """
    min_nt = (min_aa + 1) * 3
    pool = []                      # (orf_id, tx, a, e, codon, klass, nt_seq, protein)
    for tx, seq in iter_fasta(uni_fa):
        cds, bt = cds_of.get(tx), biotype_of.get(tx)
        for a, e, codon in candidate_orfs(seq, min_nt, starts):
            if e > len(seq):
                continue
            klass = orf_class(cds_relationship(a, e, cds), bt)
            if klass not in NOVEL_CLASSES:
                continue
            nt = seq[a:e].upper().replace("U", "T")
            prot = translate(nt[:-3])
            if len(prot) < min_aa or "*" in prot or prot in canon_seqs:
                continue
            pool.append((f"{tx}_{a}_{e}", tx, a, e, codon, klass, nt, prot))
    return pool


def _run(cmd, what):
    """Run a scorer, and on failure surface ITS stderr -- capture_output hides the real message."""
    for f in (cmd[2], cmd[4]) if what == "cpat" else ():
        if not Path(f).exists():
            raise SystemExit(f"{what}: missing input {f}\n"
                             "(models must be on the GROUP filesystem, not /data/tmp -- that is "
                             "node-local and invisible to compute nodes)")
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"{what} failed (exit {r.returncode}):\n"
                         f"--- stderr ---\n{r.stderr[-2000:]}\n--- stdout ---\n{r.stdout[-800:]}")
    return r


def run_cpat(nt_fa, out_prefix, hexamer, logit, py):
    _run([py, "-x", str(hexamer), "-d", str(logit), "-g", str(nt_fa),
          "-o", str(out_prefix), "--top-orf=1", "--min-orf=25"], "cpat")
    best = Path(str(out_prefix) + ".ORF_prob.best.tsv")
    scores = {}
    with best.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            # seq_ID is the ORF id we wrote; CPAT rescans it, so keep its coding probability
            scores[r["seq_ID"]] = float(r["Coding_prob"])
    return scores


def run_cpc2(nt_fa, out_prefix, py, script):
    _run([py, str(script), "-i", str(nt_fa), "-o", str(out_prefix)], "cpc2")
    res = {}
    with open(str(out_prefix) + ".txt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = r.get("#ID") or r.get("ID")
            res[key] = (float(r["coding_probability"]), r["label"])
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--universe-fa", required=True)
    ap.add_argument("--species", required=True, choices=["human", "mouse"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--tools", default="cpat,cpc2")
    ap.add_argument("--cpat-cut", type=float, default=CPAT_DEFAULT_CUT)
    ap.add_argument("--min-aa", type=int, default=7)
    ap.add_argument("--starts", default="ATG",
                    help="comma-separated start codons; must MATCH the null arm being compared to")
    E = "/private/groups/carpenterlab/emalekos/conda_envs"
    M = str(Path(__file__).resolve().parents[2] / "refs" / "cpat_models")  # GROUP fs;
    # NOT /data/tmp -- that is node-local, so a SLURM job on another node sees no models.
    ap.add_argument("--cpat-bin", default=f"{E}/quant/bin/cpat")
    ap.add_argument("--cpat-hexamer", default=f"{M}/Human_Hexamer.tsv")
    ap.add_argument("--cpat-logit", default=f"{M}/Human_logitModel.RData")
    ap.add_argument("--cpc2-python", default=f"{E}/cpc2/bin/python")
    ap.add_argument("--cpc2-script", default=f"{E}/cpc2/bin/CPC2.py")
    a = ap.parse_args()

    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    tools = {t.strip() for t in a.tools.split(",") if t.strip()}
    starts = tuple(x.strip().upper() for x in a.starts.split(",") if x.strip())

    # Same reference resolution build_dbs uses, so canonical base and class assignment match exactly.
    refs = species_refs(a.species)
    canon = load_canonical(refs["gencode_proteome"], a.min_aa)
    canon_seqs = set(canon)
    cds_of = rc_io.load_tx2cds(refs["tx2cds"])
    biotype_of = rc_io.load_biotype(refs["tx2biotype"])

    pool = enumerate_pool(a.universe_fa, cds_of, biotype_of, a.min_aa, starts, canon_seqs)
    print(f"candidate pool (same enumeration as null_atg): {len(pool):,} ORFs", file=sys.stderr)

    nt_fa = out / "orfs.nt.fa"
    with nt_fa.open("w") as fh:
        for oid, _tx, _a, _e, _c, _k, nt, _p in pool:
            fh.write(f">{oid}\n{nt}\n")

    by_id = {p[0]: p for p in pool}
    summary = {"n_pool": len(pool), "starts": list(starts), "min_aa": a.min_aa,
               "cpat_cut": a.cpat_cut, "arms": {}}

    if "cpat" in tools:
        scores = run_cpat(nt_fa, out / "cpat", a.cpat_hexamer, a.cpat_logit, a.cpat_bin)
        store = {}
        for oid, prob in scores.items():
            if prob < a.cpat_cut or oid not in by_id:
                continue
            _, tx, aa, ee, codon, klass, _nt, prot = by_id[oid]
            _add(store, prot, f"nuORF|{klass}|{oid}", klass, codon, "cpat", prob,
                 canon_seqs, a.min_aa)
        write_db(out / "db_cpat.fasta", canon, store, out / "db_cpat_class_map.tsv")
        summary["arms"]["cpat"] = {"scored": len(scores), "passing": len(store)}
        print(f"cpat : {len(scores):,} scored -> {len(store):,} novel seqs "
              f"(cut {a.cpat_cut})", file=sys.stderr)

    if "cpc2" in tools:
        res = run_cpc2(nt_fa, out / "cpc2", a.cpc2_python, a.cpc2_script)
        store = {}
        for oid, (prob, label) in res.items():
            if label != "coding" or oid not in by_id:
                continue
            _, tx, aa, ee, codon, klass, _nt, prot = by_id[oid]
            _add(store, prot, f"nuORF|{klass}|{oid}", klass, codon, "cpc2", prob,
                 canon_seqs, a.min_aa)
        write_db(out / "db_cpc2.fasta", canon, store, out / "db_cpc2_class_map.tsv")
        summary["arms"]["cpc2"] = {"scored": len(res), "passing": len(store)}
        print(f"cpc2 : {len(res):,} scored -> {len(store):,} novel seqs (label=coding)",
              file=sys.stderr)

    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"wrote {out}/summary.json", file=sys.stderr)


if __name__ == "__main__":
    main()
