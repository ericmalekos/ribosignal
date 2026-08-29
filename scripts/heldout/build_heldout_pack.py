#!/usr/bin/env python3
"""Pool an EXTERNAL held-out dataset's per-SRR RiboCode P-sites + RNAseq coverage directly
into a model-input pack (posture A), for the cross-study / cross-species held-out eval.

Mirrors pack_target_coverage_tissue.py's output format exactly, but pools straight from the
per-SRR hd5 (no 3.7 GB pooled intermediates -- ceph is tight) and joins strictly by
versioned transcript id (the coverage hd5 @SQ order may differ from the RiboCode order).

Two universe modes:
  ref_pack given (human Ruiz-Orera): reuse the Fibroblast pack's tx_order / offsets / lengths
    verbatim, so the shared orf_track_v2.npy, the one-hot FASTA, and the RiNALMo / Orthrus
    per-token embeddings all align 1:1 -- only target_counts + coverage + coverage_norm +
    pack_meta are dataset-specific, and ALL THREE backends can be evaluated. The build asserts
    every pooled per-nt length equals the reference length (same GENCODE v49 annotation), so a
    structural mismatch fails loudly.
  universe_tx file (mouse Wang): a species-specific tx list (one versioned id per line);
    offsets / lengths are computed fresh from the pooled arrays. One-hot backend only (no
    mouse FM embeddings needed -- the whole point of the one-hot control).

Output (data/packed_heldout_<dataset>/): tx_order.txt, offsets.npy, lengths.npy,
target_counts.npy, coverage.npy, coverage_norm.json, pack_meta.tsv.

Usage: build_heldout_pack.py <dataset>            # dataset in DATASETS below
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")

DATASETS = {
    "human_ruizorera": {
        "psites_dir": NEW / "data" / "heldout_psites" / "human_ruizorera",
        "coverage_dir": NEW / "data" / "heldout_rnaseq_coverage" / "human_ruizorera",
        "ref_pack": NEW / "data" / "packed",          # reuse Fibroblast v49 universe
        "universe_tx": None,
        "species": "human",
    },
    "mouse_wang_liver": {
        "psites_dir": NEW / "data" / "heldout_psites" / "mouse_wang_liver",
        "coverage_dir": NEW / "data" / "heldout_rnaseq_coverage" / "mouse_wang_liver",
        "ref_pack": None,
        "universe_tx": NEW / "data" / "heldout_refs" / "mouse_wang_universe.txt",
        "species": "mouse",
    },
    # GSE120762 BMDM +/- LPS (mouse). Conditions kept SEPARATE: each has its own salmon-derived expressed
    # universe (build_line_universe.py --species mouse), its own mm1 RNA coverage + RiboCode P-sites.
    "mouse_gse120762_nt": {
        "psites_dir": NEW / "data" / "heldout_psites" / "mouse_gse120762_nt",
        "coverage_dir": NEW / "data" / "heldout_rnaseq_coverage" / "mouse_gse120762_nt",
        "ref_pack": None,
        "universe_tx": NEW / "proteogenomics" / "data" / "gse120762" / "gse120762_nt_universe_tx.txt",
        "species": "mouse",
    },
    "mouse_gse120762_lps": {
        "psites_dir": NEW / "data" / "heldout_psites" / "mouse_gse120762_lps",
        "coverage_dir": NEW / "data" / "heldout_rnaseq_coverage" / "mouse_gse120762_lps",
        "ref_pack": None,
        "universe_tx": NEW / "proteogenomics" / "data" / "gse120762" / "gse120762_lps_universe_tx.txt",
        "species": "mouse",
    },
    # GSE155087 mouse CD4+ T cells (WT/control only). Own salmon-derived expressed universe (20,618 tx,
    # mean-TPM>=1 over 3 mRNA reps, pc/lncRNA, <=10K nt, no chrM), own mm1 RNA coverage + pooled RiboCode P-sites.
    "mouse_gse155087_tcell": {
        "psites_dir": NEW / "data" / "heldout_psites" / "mouse_gse155087_tcell",
        "coverage_dir": NEW / "data" / "heldout_rnaseq_coverage" / "mouse_gse155087_tcell",
        "ref_pack": None,
        "universe_tx": NEW / "data" / "heldout_refs" / "mouse_gse155087_tcell_universe_tx.txt",
        "species": "mouse",
    },
    # proteogenomics UC-A: predict-only A549 (cross-project). No Ribo-seq target (zeros); the model
    # gets A549's cross-project RNA-seq coverage on the Fibroblast v49 universe so orf_track/one-hot align.
    "human_a549_encode": {
        "psites_dir": None,
        "coverage_dir": NEW / "proteogenomics" / "data" / "A549_pilot" / "coverage",
        "ref_pack": NEW / "data" / "packed",
        "universe_tx": None,
        "species": "human",
        "coverage_only": True,
    },
    # UC-A FRESH universe: A549's OWN expressed transcriptome (46,307 tx), NOT the Fibroblast universe --
    # cell-type-specific prediction, so A549-specific transcripts (and their novel ORFs) are not skipped.
    "human_a549_fresh": {
        "psites_dir": None,
        "coverage_dir": NEW / "proteogenomics" / "data" / "A549_pilot" / "coverage",
        "ref_pack": None,
        "universe_tx": NEW / "proteogenomics" / "data" / "A549_pilot" / "A549_universe_tx.txt",
        "species": "human",
        "coverage_only": True,
    },
    # UC-B HBL-1 (DLBCL) fresh universe -- its own expressed transcriptome, for the immunopeptidome run
    "human_hbl1": {
        "psites_dir": None,
        "coverage_dir": NEW / "proteogenomics" / "data" / "HBL1_pilot" / "coverage",
        "ref_pack": None,
        "universe_tx": NEW / "proteogenomics" / "data" / "HBL1_pilot" / "HBL1_universe_tx.txt",
        "species": "human",
        "coverage_only": True,
    },
    # MHC-I expansion (Fig 2 pattern): same coverage-only structure, one per cell line.
    "human_dohh2": {
        "psites_dir": None,
        "coverage_dir": NEW / "proteogenomics" / "data" / "DoHH2_pilot" / "coverage",
        "ref_pack": None,
        "universe_tx": NEW / "proteogenomics" / "data" / "DoHH2_pilot" / "DoHH2_universe_tx.txt",
        "species": "human",
        "coverage_only": True,
    },
    "human_sudhl4": {
        "psites_dir": None,
        "coverage_dir": NEW / "proteogenomics" / "data" / "SUDHL4_pilot" / "coverage",
        "ref_pack": None,
        "universe_tx": NEW / "proteogenomics" / "data" / "SUDHL4_pilot" / "SUDHL4_universe_tx.txt",
        "species": "human",
        "coverage_only": True,
    },
}

# Per-population macrophage coverage-only packs (proteogenomics per-population DB build). Auto-generated from
# the 12 built universes + cov_by_pop symlink dirs, so `build_heldout_pack.py macro_<pop>` just works.
_MACRO = NEW / "proteogenomics" / "data" / "macrophage_tissue"
for _u in sorted((_MACRO / "universes").glob("*_universe_tx.txt")):
    _pop = _u.name[: -len("_universe_tx.txt")]
    _cov = _MACRO / "cov_by_pop" / _pop
    if _cov.is_dir():
        DATASETS[f"macro_{_pop}"] = {
            "psites_dir": None, "coverage_dir": _cov, "ref_pack": None,
            "universe_tx": _u, "species": "mouse", "coverage_only": True,
        }


# Cross-species expansion packs, auto-registered from the universes built by
# scripts/xspecies/build_universe_xspecies.py, in the same spirit as the _MACRO block above.
#
# All of these use universe_tx mode with ref_pack=None, i.e. the ONE-HOT backend. That is
# deliberate and it is what makes the expansion tractable: ref_pack mode would require RiNALMo
# and Orthrus per-token embeddings for every transcript of every species, which do not exist
# outside human and mouse. The mouse_wang_liver arm above is the precedent.
#
# The ribo and rna dataset labels differ per species (the two arms often come from different
# BioProjects), so the pairing is explicit rather than inferred from a shared prefix.
_XSP = NEW / "data" / "xspecies_refs"
_XSP_ARMS = {
    # pack name            (ribo dataset,          rna dataset,            species key)
    "xsp_yeast":           ("yeast_gse173654_ribo", "yeast_gse173654_rna",  "yeast"),
    "xsp_celegans":        ("worm_gse52905_ribo",   "worm_gse52861_rna",    "celegans"),
    "xsp_zebrafish":       ("zf_gse46512_ribo",     "zf_gse70549_rna",      "zebrafish"),
    "xsp_gorilla":         ("primate_gg_ribo",      "primate_gg_rna",       "gorilla"),
    "xsp_chimp":           ("primate_pt_ribo",      "primate_pt_rna",       "chimp"),
    "xsp_macaque":         ("primate_rm_ribo",      "primate_rm_rna",       "macaque"),
    "xsp_human":           ("ruizorera_hsCM_ribo",  "ruizorera_hsCM_rna",   "human_refseq"),
}
for _name, (_ribo, _rna, _sp) in _XSP_ARMS.items():
    _uni = _XSP / f"{_ribo}_universe.txt"
    _ps = NEW / "data" / "xspecies_psites" / _ribo
    _cv = NEW / "data" / "xspecies_rna_coverage_mm10" / _rna
    if _uni.exists() and _ps.is_dir() and _cv.is_dir():
        DATASETS[_name] = {
            "psites_dir": _ps, "coverage_dir": _cv, "ref_pack": None,
            "universe_tx": _uni, "species": _sp,
        }


def read_ids(h):
    return h["transcript_ids"].asstr()[:]


def pool_universe(files, ds_name, want):
    """tx -> int64 pooled per-nt array, summed across `files`, restricted to `want` (a set).
    Reads rows in on-disk order per file (sorted picks) so vlen reads stay sequential."""
    acc = {}
    for k, fp in enumerate(files):
        with h5py.File(fp, "r") as h:
            ids = read_ids(h)
            idx = {t: i for i, t in enumerate(ids)}
            picks = sorted((idx[t], t) for t in want if t in idx)
            ds = h[ds_name]
            for n, (i, t) in enumerate(picks):
                a = np.asarray(ds[i], dtype=np.int64)
                if t in acc:
                    if acc[t].shape[0] != a.shape[0]:
                        raise SystemExit(f"per-nt length mismatch {t} in {fp.name}: "
                                         f"{a.shape[0]} vs {acc[t].shape[0]}")
                    acc[t] += a
                else:
                    acc[t] = a
                if (n + 1) % 8000 == 0:
                    print(f"  [{k + 1}/{len(files)}] {fp.name}: {n + 1}/{len(picks)}",
                          file=sys.stderr)
        print(f"pooled {ds_name} from {fp.name} ({len(picks)} universe rows)", file=sys.stderr)
    return acc


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in DATASETS:
        print(f"usage: build_heldout_pack.py <{'|'.join(DATASETS)}>", file=sys.stderr)
        return 2
    name = sys.argv[1]
    cfg = DATASETS[name]
    # RIBO_COVERAGE_DIR / RIBO_PACK_OUT_SUFFIX let a variant coverage posture (e.g. mm10) be packed
    # WITHOUT editing DATASETS and without clobbering the existing pack. Both default to unset, so
    # every existing caller is byte-for-byte unchanged.
    import os as _os
    _cov = _os.environ.get("RIBO_COVERAGE_DIR", "")
    if _cov:
        cfg = dict(cfg)
        cfg["coverage_dir"] = Path(_cov)
        print(f"RIBO_COVERAGE_DIR override -> {_cov}", file=sys.stderr)
        if not Path(_cov).is_dir():
            sys.exit(f"RIBO_COVERAGE_DIR={_cov} is not a directory")
    psites_files = sorted(cfg["psites_dir"].glob("*_psites.hd5")) if cfg.get("psites_dir") else []
    cov_files = sorted(cfg["coverage_dir"].glob("*_coverage.hd5"))
    if not cfg.get("coverage_only"):
        assert psites_files, f"no *_psites.hd5 in {cfg['psites_dir']}"
    assert cov_files, f"no *_coverage.hd5 in {cfg['coverage_dir']}"
    _suf = _os.environ.get("RIBO_PACK_OUT_SUFFIX", "")
    out_dir = NEW / "data" / f"packed_heldout_{name}{_suf}"
    if _suf:
        print(f"RIBO_PACK_OUT_SUFFIX -> {out_dir.name} (existing pack kept)",
              file=sys.stderr)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- universe ----------------------------------------------------------------------
    ref_pack = cfg["ref_pack"]
    if ref_pack is not None:
        order = (ref_pack / "tx_order.txt").read_text().split()
        ref_len = np.load(ref_pack / "lengths.npy")
        ref_len_of = {t: int(ref_len[i]) for i, t in enumerate(order)}
        print(f"dataset={name}  universe={len(order):,} tx (ref_pack {ref_pack.name})  "
              f"ribo={len(psites_files)} rna={len(cov_files)}", file=sys.stderr)
    else:
        order = Path(cfg["universe_tx"]).read_text().split()
        ref_len_of = None
        print(f"dataset={name}  universe={len(order):,} tx (fresh)  "
              f"ribo={len(psites_files)} rna={len(cov_files)}", file=sys.stderr)
    want = set(order)
    assert len(want) == len(order), "duplicate tx in universe"

    # ---- pool ------------------------------------------------------------------------
    print("pooling RNAseq coverage...", file=sys.stderr)
    cov = pool_universe(cov_files, "coverage", want)
    if cfg.get("coverage_only"):
        # predict-only: no Ribo-seq target -> zeros matched to each tx's coverage length
        tgt = {t: np.zeros(cov[t].shape[0], dtype=np.int64) for t in cov}
        print("coverage-only dataset: target = zeros (predict-only)", file=sys.stderr)
    else:
        print("pooling P-site target...", file=sys.stderr)
        tgt = pool_universe(psites_files, "p_sites", want)

    # a universe tx must appear in at least one ribo AND one rna file; else it cannot be packed
    miss_t = [t for t in order if t not in tgt]
    miss_c = [t for t in order if t not in cov]
    if miss_t or miss_c:
        raise SystemExit(f"universe tx absent from inputs: {len(miss_t)} in target "
                         f"(e.g. {miss_t[:3]}), {len(miss_c)} in coverage (e.g. {miss_c[:3]})")

    # ---- lengths / offsets -----------------------------------------------------------
    lengths = np.empty(len(order), dtype=np.int32)
    for k, t in enumerate(order):
        lt, lc = tgt[t].shape[0], cov[t].shape[0]
        if lt != lc:
            raise SystemExit(f"length mismatch {t}: target {lt} vs coverage {lc}")
        if ref_len_of is not None and lt != ref_len_of[t]:
            raise SystemExit(f"held-out length {t}={lt} != ref_pack {ref_len_of[t]} "
                             "(annotation mismatch -> orf_track/embeddings would misalign)")
        lengths[k] = lt
    offsets = np.zeros(len(order) + 1, dtype=np.int64)
    offsets[1:] = np.cumsum(lengths.astype(np.int64))
    total = int(offsets[-1])

    tgt_cat = np.empty(total, dtype=np.int32)
    cov_cat = np.empty(total, dtype=np.int32)
    for k, t in enumerate(order):
        s, e = offsets[k], offsets[k + 1]
        tgt_cat[s:e] = tgt[t].astype(np.int32)
        cov_cat[s:e] = cov[t].astype(np.int32)

    np.save(out_dir / "offsets.npy", offsets)
    np.save(out_dir / "lengths.npy", lengths)
    np.save(out_dir / "target_counts.npy", tgt_cat)
    np.save(out_dir / "coverage.npy", cov_cat)
    (out_dir / "tx_order.txt").write_text("\n".join(order) + "\n")

    global_mean = float(cov_cat.sum(dtype=np.int64)) / int(cov_cat.size)
    (out_dir / "coverage_norm.json").write_text(json.dumps({
        "global_mean_coverage": global_mean,
        "dataset": name,
        "species": cfg["species"],
        "n_ribo_samples": len(psites_files),
        "n_rna_samples": len(cov_files),
        "sum_L": total,
        "coverage_total": int(cov_cat.sum(dtype=np.int64)),
    }, indent=2) + "\n")

    n_scorable = 0
    with (out_dir / "pack_meta.tsv").open("w") as m:
        m.write("tx_id\tlength\ttotal_psites\ttotal_coverage\n")
        for k, t in enumerate(order):
            s, e = offsets[k], offsets[k + 1]
            tp = int(tgt_cat[s:e].sum())
            if tp >= 50:
                n_scorable += 1
            m.write(f"{t}\t{lengths[k]}\t{tp}\t{int(cov_cat[s:e].sum())}\n")

    print(f"wrote {out_dir}/  global_mean_coverage={global_mean:.3f}", file=sys.stderr)
    print(f"  transcripts={len(order):,}  sum_L={total:,}  "
          f"target_total={int(tgt_cat.sum()):,}  coverage_total={int(cov_cat.sum()):,}  "
          f"scorable(>=50 psites)={n_scorable:,}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
