#!/usr/bin/env python3
"""MSFragger parameter templating + run caching (`pgx` step 7).

Two digestion modes, selected by --enzyme, so the same pipeline serves shotgun proteomes and
immunopeptidomes:

  tryptic      trypsin, KR, fully tryptic (num_enzyme_termini = 2), 1 missed cleavage, 6-50 aa
  nonspecific  no enzyme (num_enzyme_termini = 0), 8-14 aa  -- HLA / immunopeptidomics

Everything else (tolerances, modifications, instrument settings) comes from a validated template
in proteogenomics/msfragger/ so the acquisition-specific settings are inherited rather than
reinvented. Override with --template for a different instrument.

Large databases need heap, not a slicing parameter. MSFragger 4.2 sizes its own fragment-index
slices from the JVM heap and REJECTS `num_slices` ("Unknown parameters"); verified on the 2.36M
target near-cognate null, where it built 672,835,480 fragments as 6.27 GB in a single slice. The
only lever is therefore `-Xmx` (the XMX variable in search.sbatch). The target count is reported
here so an undersized heap is easy to spot.

CACHING. The run identity is sha256(database) + sha256(parameters). `pgx_hash` writes that digest
to <out>/.pgx_hash; a search whose digest already matches and whose outputs are present is skipped.
This is what keeps the GENCODE baseline from being re-searched once per arm.

  python -m pgx.search --db db_model.fasta --enzyme tryptic --out-params fragger.params
  python -m pgx.search --db db_model.fasta --enzyme tryptic --print-hash
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

from .refs import tool_path

ENZYME = {
    "tryptic": {
        "search_enzyme_name_1": "trypsin",
        "search_enzyme_cut_1": "KR",
        # `P` = classic trypsin, which does NOT cut before proline. Leaving this empty silently
        # switches the search to trypsin/P and changes the peptide space: measured on BMDM, the
        # GENCODE-only baseline moved 335,548 -> 352,307 PSMs (+5.0%) from this one character,
        # which is enough to make results incomparable with every prior search in this project.
        "search_enzyme_nocut_1": "P",
        "search_enzyme_sense_1": "C",
        "allowed_missed_cleavage_1": "1",
        "num_enzyme_termini": "2",
        "digest_min_length": "6",
        "digest_max_length": "50",
    },
    "nonspecific": {
        "search_enzyme_name_1": "nonspecific",
        "search_enzyme_cut_1": "@",
        "search_enzyme_nocut_1": "",
        "search_enzyme_sense_1": "C",
        "allowed_missed_cleavage_1": "2",
        "num_enzyme_termini": "0",
        "digest_min_length": "8",
        "digest_max_length": "14",
    },
}
# FROZEN templates are the default. The unfrozen originals set calibrate_mass = 2, which re-derives
# search parameters from a first-pass search against WHICHEVER database is given -- so two arms of one
# comparison get scored under different rules. Verified on all three: the macrophage dPSM column was
# inflated by +11,017; on HBL-1 the null arm ran use_topN_peaks=100 / intensity_transform=0 against
# the model arm's 150 / 1 (cutoff 16.742 vs 18.483); on A549 the canonical baseline itself diverged.
# pgx exists to compare databases, so an unfrozen template is never the right default here.
# Pass --template explicitly to reproduce a pre-2026-08-04 run.
DEFAULT_TEMPLATE = {"tryptic": "fragger_macro_lfq_frozen.params",
                    "nonspecific": "fragger_hbl1_hla_frozen.params"}


def count_targets(fasta):
    n = 0
    with open(fasta) as fh:
        for ln in fh:
            if ln.startswith(">") and not ln.startswith(">REV_"):
                n += 1
    return n


def render(template, db, enzyme, extra=None):
    """Rewrite the template's key=value lines. Keys absent from the template are appended."""
    over = {"database_name": str(db), **ENZYME[enzyme]}
    if extra:
        over.update(extra)
    seen, out = set(), []
    for ln in Path(template).read_text().splitlines():
        s = ln.strip()
        if s and not s.startswith("#") and "=" in s:
            key = s.split("=", 1)[0].strip()
            if key in over:
                seen.add(key)
                out.append(f"{key} = {over[key]}")
                continue
        out.append(ln)
    for k, v in over.items():
        if k not in seen:
            out.append(f"{k} = {v}")
    return "\n".join(out) + "\n"


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def pgx_hash(db, params_text):
    """Content identity of a search: sha256(database) + sha256(settings).

    `database_name` is EXCLUDED from the settings digest. It holds the absolute path of the FASTA,
    which is per-run, so including it made two searches over byte-identical databases hash
    differently and defeated --shared-search-root entirely: the model-independent arms
    (gencode, null_atg, null_nc) were re-searched for every checkpoint even though their FASTAs
    had identical md5s. The database's CONTENT is already covered by sha256(db); where the file
    happens to live is not part of what the search is.
    """
    settings = "\n".join(ln for ln in params_text.splitlines()
                         if not ln.strip().startswith("database_name"))
    h = hashlib.sha256()
    h.update(sha256(db).encode())
    h.update(hashlib.sha256(settings.encode()).hexdigest().encode())
    return h.hexdigest()[:32]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True)
    ap.add_argument("--enzyme", required=True, choices=sorted(ENZYME))
    ap.add_argument("--template", default=None, help="MSFragger params template to inherit from")
    ap.add_argument("--out-params", default=None, help="write the rendered params here")
    ap.add_argument("--print-hash", action="store_true", help="print the cache digest and exit")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="additional parameter override, repeatable")
    a = ap.parse_args()

    # --template accepts EITHER a path or a bare basename. The basename form resolves against the
    # project's msfragger dir, the same place DEFAULT_TEMPLATE is looked up. Without this the
    # override resolved against the JOB's cwd while the default resolved against the params dir --
    # two different contracts for one option, which killed all 5 A549 searches at 1 second each with
    # "missing params template: fragger_a549_tmt_frozen.params" for a file that existed.
    pdir = tool_path("msfragger_params_dir")
    if a.template:
        tmpl = Path(a.template)
        if not tmpl.exists() and (pdir / tmpl.name).exists():
            tmpl = pdir / tmpl.name
    else:
        tmpl = pdir / DEFAULT_TEMPLATE[a.enzyme]
    if not tmpl.exists():
        raise SystemExit(f"missing params template: {tmpl}\n"
                         f"(looked for a bare basename under {pdir} as well)")
    extra = dict(kv.split("=", 1) for kv in a.set)
    text = render(tmpl, a.db, a.enzyme, extra)

    if a.print_hash:
        print(pgx_hash(a.db, text))
        return
    print(f"template   : {tmpl}", file=sys.stderr)
    print(f"enzyme     : {a.enzyme} (num_enzyme_termini={ENZYME[a.enzyme]['num_enzyme_termini']}, "
          f"digest {ENZYME[a.enzyme]['digest_min_length']}-{ENZYME[a.enzyme]['digest_max_length']})",
          file=sys.stderr)
    print(f"targets    : {count_targets(a.db):,} (MSFragger sizes its index slices from -Xmx)",
          file=sys.stderr)
    if a.out_params:
        Path(a.out_params).write_text(text)
        print(f"wrote {a.out_params}", file=sys.stderr)
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
