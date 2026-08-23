#!/usr/bin/env python3
"""Build ONE training/eval-ready model-input pack from already-produced per-nt hd5, driven
entirely by CLI paths (no hardcoded project paths, tissue names, or sample-count asserts).

This is the CLI-fied, generalized successor to scripts/heldout/build_heldout_pack.py and the
scripts/pack_target_coverage*.py / pool_rnaseq_coverage*.py family. It pools per-sample Ribo
P-sites (RiboCode `*_psites.hd5`) and RNA coverage (`rnaseq_coverage.py *_coverage.hd5`) onto
a transcript universe and writes the 7-file pack contract, optionally building the per-pack
ORF track (default --kozak none).

Universe: pass EITHER --universe-tx (a fresh tx-id list, one versioned id/line) OR --ref-pack
(reuse an existing pack's tx_order/lengths verbatim so a shared ORF track + embeddings align).

Examples
  # rebuild a Chothani tissue pack on the shared universe (reuse target, new mm1 coverage):
  prepare_pack.py --ribo-psites data/target/HUVEC_psites_pooled.hd5 \
      --rna-coverage data/rnaseq_coverage/per_sample/SRR*.hd5 \
      --ref-pack data/packed --tx2biotype data/tx2biotype.tsv \
      --group HUVEC --species human --out data/packed_HUVEC

  # a future user's fresh dataset (own universe + per-pack no-Kozak ORF track):
  prepare_pack.py --ribo-psites mydata/psites/ --rna-coverage mydata/coverage/ \
      --universe-tx mydata/universe_tx.txt --fasta mydata/universe.fa --build-orf-track \
      --group mycellline --species human --out packs/mycellline
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import packlib  # noqa: E402
from paths import project_root  # noqa: E402


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ribo-psites", nargs="*", default=[],
                    help="RiboCode per-nt *_psites.hd5 file(s) or dir(s); pooled if several. "
                         "Omit with --coverage-only.")
    # NOT required: --reuse-coverage supplies coverage from an existing pack instead. Validated
    # below so the error names the actual missing thing rather than argparse's generic message.
    ap.add_argument("--rna-coverage", nargs="+", default=[],
                    help="rnaseq_coverage.py *_coverage.hd5 file(s) or dir(s); pooled if several. "
                         "Omit only with --reuse-coverage.")
    # universe (+ optionally target) source: exactly one
    uni = ap.add_mutually_exclusive_group(required=True)
    uni.add_argument("--universe-tx", help="fresh tx-id list (one versioned id per line)")
    uni.add_argument("--ref-pack", help="reuse an existing pack's tx_order/lengths verbatim")
    uni.add_argument("--reuse-target",
                     help="reuse an existing pack's universe AND target_counts verbatim; only "
                          "pool + write new coverage (the 'swap RNA-seq' path, e.g. mm20->mm1)")
    uni.add_argument("--reuse-coverage",
                     help="MIRROR of --reuse-target: reuse an existing pack's universe AND "
                          "coverage verbatim, pool + write only a NEW Ribo target. This is the "
                          "'Ribo re-aligned, RNA untouched' path the final recipe requires -- "
                          "PIPELINE_POLICY and rebuild_heldout_canon are explicit that RNA-seq "
                          "stays Local and the universe is never redefined, since the universe came "
                          "from salmon TPM on that RNA. Use when the Ribo BAMs were rebuilt but the "
                          "RNA coverage hd5 no longer exist.")
    ap.add_argument("--group", required=True, help="pack/group name (provenance + coverage_norm)")
    ap.add_argument("--species", default="human", help="written into coverage_norm.json")
    ap.add_argument("--out", required=True, help="output pack directory")
    ap.add_argument("--coverage-only", action="store_true",
                    help="predict-only: no Ribo target (target_counts = zeros)")
    ap.add_argument("--psites-glob", default="*_psites.hd5")
    ap.add_argument("--coverage-glob", default="*_coverage.hd5")
    # pack_meta biotype columns (drop-in for the tissue packs = the eval biotype source)
    ap.add_argument("--tx2biotype", help="build_tx2biotype.py TSV -> 6-col pack_meta")
    ap.add_argument("--gtf", help="GENCODE GTF -> build tx2biotype if --tx2biotype absent")
    # optional per-pack ORF track (FRESH universes; shared-universe packs reuse a shared track)
    ap.add_argument("--build-orf-track", action="store_true")
    ap.add_argument("--fasta", help="universe FASTA (required for --build-orf-track)")
    ap.add_argument("--kozak", default="none", choices=["none", "heuristic", "pwm"],
                    help="ORF-track start context (default none; NEVER set heuristic)")
    ap.add_argument("--orf-mode", default="ext", choices=["atg", "ext", "atgctg"])
    ap.add_argument("--md5", action="store_true", help="also record input md5s (slow)")
    return ap.parse_args()


def main():
    args = parse_args()
    root = project_root()
    out_dir = Path(args.out)

    need_psites = not args.coverage_only and not args.reuse_target
    if need_psites and not args.ribo_psites:
        print("ERROR: --ribo-psites required unless --coverage-only/--reuse-target",
              file=sys.stderr)
        return 2
    if args.build_orf_track and not args.fasta:
        print("ERROR: --build-orf-track requires --fasta", file=sys.stderr)
        return 2

    psites_files = packlib.resolve_hd5(args.ribo_psites, args.psites_glob)
    cov_files = [] if args.reuse_coverage else packlib.resolve_hd5(args.rna_coverage,
                                                                   args.coverage_glob)
    if not cov_files and not args.reuse_coverage:
        print("ERROR: --rna-coverage resolved no coverage hd5, and --reuse-coverage was not given",
              file=sys.stderr)
        return 2
    if need_psites and not psites_files:
        print("ERROR: no psites hd5 resolved from --ribo-psites", file=sys.stderr)
        return 2

    # universe: --reuse-target and --ref-pack both take it (and lengths) from an existing pack
    universe_pack = args.reuse_target or args.reuse_coverage or args.ref_pack
    order, ref_len_of = packlib.load_universe(
        universe_tx=args.universe_tx, ref_pack=universe_pack)
    want = set(order)
    src = (f"reuse-coverage {Path(args.reuse_coverage).name}" if args.reuse_coverage
           else f"reuse-target {Path(args.reuse_target).name}" if args.reuse_target
           else f"ref_pack {Path(args.ref_pack).name}" if args.ref_pack else "fresh universe_tx")
    print(f"group={args.group}  universe={len(order):,} tx ({src})  "
          f"ribo={len(psites_files)} rna={len(cov_files)}", file=sys.stderr)

    print("pooling RNA-seq coverage...", file=sys.stderr)
    if args.reuse_coverage:
        cp = Path(args.reuse_coverage)
        cc = np.load(cp / "coverage.npy", mmap_mode="r")
        coff = np.load(cp / "offsets.npy")
        cov = {t_: np.asarray(cc[coff[i]:coff[i + 1]], dtype=np.int64)
               for i, t_ in enumerate(order)}
        print(f"reused coverage from {cp.name} ({len(order):,} tx) -- RNA is unchanged by the "
              f"final recipe", file=sys.stderr)
    else:
        cov = packlib.pool_per_nt(cov_files, "coverage", want)
    if args.reuse_target:
        tp = Path(args.reuse_target)
        tc = np.load(tp / "target_counts.npy", mmap_mode="r")
        off = np.load(tp / "offsets.npy")
        tgt = {t: np.asarray(tc[off[i]:off[i + 1]], dtype=np.int64) for i, t in enumerate(order)}
        print(f"reused target_counts from {tp.name} ({len(order):,} tx)", file=sys.stderr)
    elif args.coverage_only:
        tgt = {t: np.zeros(cov[t].shape[0], dtype=np.int64) for t in cov}
        print("coverage-only: target = zeros (predict-only)", file=sys.stderr)
    else:
        print("pooling P-site target...", file=sys.stderr)
        tgt = packlib.pool_per_nt(psites_files, "p_sites", want)

    # biotype columns for pack_meta (optional)
    biotype = None
    tx2b_path = args.tx2biotype
    if tx2b_path is None and args.gtf:
        tx2b_path = str(out_dir / "tx2biotype.tsv")
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"building tx2biotype from {args.gtf} -> {tx2b_path}", file=sys.stderr)
        subprocess.run([sys.executable, str(root / "scripts" / "build_tx2biotype.py"),
                        "--gtf", args.gtf, "--out", tx2b_path], check=True)
    if tx2b_path:
        biotype = packlib.load_tx_biotype(tx2b_path)

    norm = packlib.write_pack(
        out_dir, order, tgt, cov,
        group=args.group, species=args.species, biotype=biotype,
        ref_len_of=ref_len_of, n_ribo=len(psites_files), n_rna=len(cov_files))

    # optional per-pack ORF track (fresh universes). Shared-universe packs reuse a shared track
    # via $RIBO_ORF_TRACK and should NOT pass --build-orf-track.
    orf_track = None
    if args.build_orf_track:
        # Encode the Kozak setting in the FILENAME. Project convention is that a bare
        # orf_track_v2.npy is the Kozak-heuristic track and orf_track_v2_nokozak.npy is the
        # --kozak none track; writing no-Kozak content under the bare name is a silent
        # train/inference-mismatch trap (it is exactly how the nokozak Janich dump ended up
        # being fed the Kozak track). Name follows content.
        orf_names = {"atg": "orf_track", "ext": "orf_track_v2", "atgctg": "orf_track_v3"}
        # bare name = the Kozak heuristic (the historical default); matches the existing
        # data/packed*/orf_track_v2_nokozak.npy convention used across the project.
        kozak_suffix = {"heuristic": "", "none": "_nokozak", "pwm": "_pwm"}
        orf_out = out_dir / f"{orf_names[args.orf_mode]}{kozak_suffix[args.kozak]}.npy"
        print(f"building ORF track (--kozak {args.kozak}) -> {orf_out}", file=sys.stderr)
        subprocess.run([sys.executable, str(root / "scripts" / "build_orf_track.py"),
                        "--mode", args.orf_mode, "--kozak", args.kozak,
                        "--pack", str(out_dir), "--fasta", args.fasta,
                        "--out", str(orf_out)], check=True)
        orf_track = str(orf_out)

    # provenance
    def stat(p):
        p = Path(p)
        rec = {"path": str(p), "bytes": p.stat().st_size}
        if args.md5:
            rec["md5"] = packlib.md5(p)
        return rec

    provenance = {
        "group": args.group, "species": args.species,
        "universe": {"ref_pack": args.ref_pack, "universe_tx": args.universe_tx,
                     "n_tx": len(order)},
        "coverage_only": args.coverage_only,
        "ribo_psites_inputs": [stat(f) for f in psites_files],
        "rna_coverage_inputs": [stat(f) for f in cov_files],
        "orf_track": {"built": bool(orf_track), "path": orf_track,
                      "kozak": args.kozak, "mode": args.orf_mode} if args.build_orf_track else None,
        "coverage_norm": norm,
        "argv": sys.argv,
    }
    (out_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    # ready-to-use env hints for the trainer/eval
    onehot_fasta = args.fasta or "<universe.fa>"
    print("\n# pack ready. To train/eval on it, export:", file=sys.stderr)
    print(f"export RIBO_PACK_DIR={out_dir}", file=sys.stderr)
    if orf_track:
        print(f"export RIBO_ORF_TRACK={orf_track}", file=sys.stderr)
    else:
        print("#   (shared-universe pack: set RIBO_ORF_TRACK to the shared track, "
              "e.g. data/packed/orf_track_v2_nokozak.npy)", file=sys.stderr)
    print(f"export RIBO_ONEHOT_FASTA={onehot_fasta}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
