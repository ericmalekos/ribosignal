#!/usr/bin/env python3
"""N-terminal CDS extension calling (`pgx` step 5). Runs in the `ribocode` env.

WHY THIS IS A SEPARATE STEP, AND WHY IT TESTS THE EXTENSION REGION ALONE
-----------------------------------------------------------------------
An N-terminally extended ORF shares its stop codon with the annotated CDS and CONTAINS it, so any
whole-ORF statistic is dominated by the CDS's own translation signal:

* `pred_frame0` (the retired DB criterion) is a RATIO over the ORF interval, so prepending a
  zero-signal upstream region changes neither numerator nor denominator. It is INVARIANT to
  extension, not merely diluted. Measured on BMDM: over 12,719 stop-codon groups containing an
  annotated CDS start plus upstream in-frame starts, 96.2% of 51,646 upstream non-AUG ORFs sat
  within 0.05 of the true CDS's f0; median |delta f0| = 0.0026.
* RiboCode does not help either. `start_check` with only_longest_orf=True simply takes the most 5'
  in-frame start (`start_list[0]`) and tests the WHOLE ORF, so it reports an extension without ever
  asking whether that start is used. And `orf_finder.orf_find` sets alt_flag=0 as soon as any
  in-frame ATG precedes the stop, so a NEAR-COGNATE extension of an annotated CDS is unreachable
  through RiboCode at any parameter setting.

The discriminating question is not "is this ORF in a translated frame" but "is THIS START used".
That is answered by the extension region alone: the segment from the candidate start up to the
annotated CDS start. Real upstream initiation puts periodic P-sites there; a downstream true start
leaves it empty. So this step applies RiboCode's OWN frame statistic (`extract_frame` +
`test_frame`: Wilcoxon f0>f1 and f0>f2, Stouffer-combined) to that segment only, then BH-adjusts
across every candidate tested.

DENSITY. The same construction the calling arm uses (`ribocode_dropin.build_density`), scaled by
the calibrated theta, but with deterministic rounding rather than a Poisson draw: `rint(p * depth *
theta)`. Rounding keeps real zeros (so a smooth low-signal region cannot sneak through RiboCode's
tie logic) while remaining exactly reproducible, with no dependence on the caller's RNG stream.

CANDIDATE GATE. Only transcripts whose annotated CDS was itself called translated in the arm are
extended -- extending an untranslated CDS is meaningless.

OUTPUT also carries `orf_f0`, the whole-ORF frame-0 fraction the retired criterion would have used,
so the contrast between the two tests is visible per candidate rather than asserted.

  ribocode_python -m pgx.extensions --profiles p.npz --calls-collapsed <arm>/..._collapsed.txt \
      --fasta universe.fa --species mouse --theta 0.05 --out <dir>
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from . import rc_io
from .refs import species_refs
from .seqtools import STOPS, bh, norm, parse_starts, read_fasta, translate

try:
    from RiboCode.detectORF import extract_frame, test_frame
except ImportError:  # pragma: no cover
    raise SystemExit("pgx.extensions must run in the `ribocode` env (RiboCode.detectORF missing)")


def upstream_candidates(seq, cds_start0, starts):
    """In-frame candidate starts upstream of the CDS, walking back to the first in-frame stop.

    -> [(start0, codon)] nearest-first. The extension may not cross an in-frame stop codon, so the
    walk terminates there; every returned start shares the CDS reading frame by construction.
    """
    out = []
    i = cds_start0 - 3
    while i >= 0:
        codon = seq[i:i + 3]
        if codon in STOPS:
            break
        if codon in starts:
            out.append((i, codon))
        i -= 3
    return out


def build_density(prof, pred_total, theta):
    """Deterministic integerized density: the arm's rate construction at theta, rounded."""
    return np.rint(np.maximum(prof.astype(np.float64), 0.0) * float(pred_total) * float(theta))


def scan(profiles, collapsed, seqs, starts, theta, qvalue, min_ext_aa, min_nonzero,
         orf_qvalue=0.05):
    """-> (rows, stats). One row per tested candidate, with BH q across all of them."""
    z = np.load(profiles, allow_pickle=False)
    tx_ids = [str(t) for t in z["tx_ids"]]
    lengths = z["lengths"]
    pred_flat = z["pred_flat"]
    pred_total = z["pred_total"]
    off = np.concatenate([[0], np.cumsum(lengths)])
    index = {t: i for i, t in enumerate(tx_ids)}

    # Only extend CDSs the arm called translated.
    targets = [o for o in rc_io.read_collapsed(collapsed)
               if o["orf_type"] == "annotated" and o["qval"] <= orf_qvalue
               and o["cds_start0"] is not None and o["start0"] is not None]
    stats = {"translated_cds": len(targets), "tx_with_seq": 0, "candidates": 0,
             "tested": 0, "skipped_short": 0, "skipped_lowsignal": 0}

    rows = []
    for o in targets:
        tx = o["tx"]
        i = index.get(tx)
        seq = seqs.get(tx)
        if i is None or seq is None:
            continue
        L = int(lengths[i])
        if len(seq) != L:
            continue
        stats["tx_with_seq"] += 1
        seq = norm(seq)
        # RiboCode may already have walked the start 5' of the annotation; anchor on the ANNOTATED
        # CDS start so the extension region is defined against the reference, not against its call.
        cds0 = o["cds_start0"]
        end0 = o["end0"]
        dens = build_density(pred_flat[off[i]:off[i + 1]], pred_total[i], theta)
        whole = dens[cds0:end0]
        orf_f0 = float(whole[0::3].sum() / whole.sum()) if whole.sum() > 0 else 0.0

        for s0, codon in upstream_candidates(seq, cds0, starts):
            stats["candidates"] += 1
            ext_aa = (cds0 - s0) // 3
            if ext_aa < min_ext_aa:
                stats["skipped_short"] += 1
                continue
            region = dens[s0:cds0]
            f0, f1, f2 = extract_frame(region)
            nz = int(np.flatnonzero(f0).size)
            if nz < min_nonzero:
                stats["skipped_lowsignal"] += 1
                continue
            _, _, pv = test_frame(f0, f1, f2, "none")
            pv = float(pv)
            stats["tested"] += 1
            prot = translate(seq[s0:end0 - 3])
            ext_whole = dens[s0:end0]
            rows.append({
                "tx": tx, "gene_id": o["gene_id"], "gene_name": o["gene_name"],
                "ext_start0": s0, "cds_start0": cds0, "end0": end0,
                "start_codon": codon, "ext_len_aa": ext_aa, "aa_len": len(prot),
                "ext_psites_f0": float(f0.sum()), "ext_nonzero_f0": nz,
                "pval": pv,
                # what the retired whole-ORF criterion would have said, for the contrast
                "orf_f0": round(orf_f0, 4),
                "ext_orf_f0": round(float(ext_whole[0::3].sum() / ext_whole.sum()), 4)
                if ext_whole.sum() > 0 else 0.0,
                "aaseq": prot.split("*")[0] if "*" in prot else prot,
            })

    if rows:
        qs = bh([r["pval"] for r in rows])
        for r, q in zip(rows, qs):
            r["qval"] = float(q)
    stats["passing"] = sum(1 for r in rows if r.get("qval", 1.0) <= qvalue)
    return rows, stats


COLS = ("tx", "gene_id", "gene_name", "ext_start0", "cds_start0", "end0", "start_codon",
        "ext_len_aa", "aa_len", "ext_psites_f0", "ext_nonzero_f0", "pval", "qval",
        "orf_f0", "ext_orf_f0", "aaseq")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--calls-collapsed", required=True,
                    help="the arm's RiboCode *_collapsed.txt (supplies the translated CDS set)")
    ap.add_argument("--fasta", required=True, help="universe FASTA (headers = versioned tx id)")
    ap.add_argument("--species", required=True, choices=["human", "mouse"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--theta", type=float, default=None, help="calibrated theta")
    ap.add_argument("--theta-json", default=None, help="theta.json from pgx.calibrate")
    ap.add_argument("--start-codons", default="near_cognate",
                    help="'ATG' | 'near_cognate' | comma list. Extensions are the ONLY route by "
                         "which a near-cognate start reaches a database, since RiboCode's alt-start "
                         "mode is fallback-only per common stop.")
    ap.add_argument("--ext-qvalue", type=float, default=0.05, help="BH q cutoff on the extension test")
    ap.add_argument("--orf-qvalue", type=float, default=0.05, help="q cutoff selecting translated CDSs")
    ap.add_argument("--min-ext-aa", type=int, default=5, help="minimum extension length in codons")
    ap.add_argument("--min-nonzero", type=int, default=5,
                    help="minimum nonzero frame-0 codons in the extension region to run the test")
    a = ap.parse_args()

    theta = a.theta
    if theta is None:
        if not a.theta_json:
            raise SystemExit("need --theta or --theta-json")
        theta = float(json.loads(Path(a.theta_json).read_text())["theta"])

    species_refs(a.species)   # validates species is known
    starts = parse_starts(a.start_codons)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    seqs = read_fasta(a.fasta)
    print(f"start codons ({len(starts)}): {','.join(starts)}   theta={theta:g}", file=sys.stderr)

    rows, stats = scan(a.profiles, a.calls_collapsed, seqs, starts, theta, a.ext_qvalue,
                       a.min_ext_aa, a.min_nonzero, a.orf_qvalue)
    rows.sort(key=lambda r: (r["qval"], r["pval"]))
    keep = [r for r in rows if r["qval"] <= a.ext_qvalue]

    with open(out / "extensions_all.tsv", "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")
    with open(out / "extensions.tsv", "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in keep:
            fh.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")

    by_codon = {}
    for r in keep:
        by_codon[r["start_codon"]] = by_codon.get(r["start_codon"], 0) + 1
    # The contrast that motivates this step: how many candidates the retired whole-ORF f0 cut would
    # have admitted, versus how many the extension-region test actually supports.
    n_f0 = sum(1 for r in rows if r["ext_orf_f0"] >= 0.5)
    stats.update({"theta": theta, "ext_qvalue": a.ext_qvalue, "start_codons": list(starts),
                  "passing_by_codon": dict(sorted(by_codon.items(), key=lambda kv: -kv[1])),
                  "would_pass_whole_orf_f0_0.5": n_f0})
    (out / "extension_summary.json").write_text(json.dumps(stats, indent=2) + "\n")

    print(f"translated CDS targets : {stats['translated_cds']:,}")
    print(f"candidate starts       : {stats['candidates']:,} "
          f"(short {stats['skipped_short']:,}, low-signal {stats['skipped_lowsignal']:,})")
    print(f"tested                 : {stats['tested']:,}")
    print(f"passing q<={a.ext_qvalue}       : {len(keep):,}"
          + (f"  ({len(keep) / stats['tested'] * 100:.1f}% of tested)" if stats["tested"] else ""))
    print(f"  by start codon       : {stats['passing_by_codon']}")
    print(f"whole-ORF f0>=0.5 would have admitted: {n_f0:,} of {len(rows):,} tested "
          "(the retired criterion, shown for contrast)")
    print(f"\nwrote {out}/extensions.tsv (+ _all.tsv, extension_summary.json)")


if __name__ == "__main__":
    main()
