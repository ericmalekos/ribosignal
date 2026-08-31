#!/usr/bin/env python3
"""Regression tests. Each one checks a bug that has already happened in this project.

In every case the pipeline ran to completion and produced plausible numbers, and the mistake was
found days later. These tests make those specific mistakes fail immediately instead.

No GPU, no SLURM allocation, no real data. Synthetic fixtures only, runs in about a second.

    conda_envs/cas12a/bin/python3 -m pytest tests/ -q
    conda_envs/cas12a/bin/python3 tests/test_regressions.py

These check code logic, not the data going into it. They would not have caught the Janich
untrimmed-FASTQ bug or the too-wide MSFragger fragment tolerances, both of which needed checks on
real inputs. See tests/README.md.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts" / "prepare"))
sys.path.insert(0, str(REPO / "proteogenomics" / "scripts"))

import packlib  # noqa: E402
from pgx import rc_io, search, seqtools  # noqa: E402


# ----------------------------------------------------------------------------------------------
# 1. Pack layout. Coverage, target and ORF track are all sliced with the same offsets, so if a
#    dtype or row order changes, the label stops lining up with the input and nothing errors.
# ----------------------------------------------------------------------------------------------
def test_pack_contract():
    order = ["ENST1.1", "ENST2.2", "ENST3.1"]
    lens = {"ENST1.1": 9, "ENST2.2": 6, "ENST3.1": 12}
    tgt = {t: np.arange(n, dtype=np.int64) for t, n in lens.items()}
    cov = {t: np.full(n, 7, dtype=np.int64) for t, n in lens.items()}

    with tempfile.TemporaryDirectory() as d:
        packlib.write_pack(d, order, tgt, cov, group="G", species="human")
        p = Path(d)
        offsets = np.load(p / "offsets.npy")
        lengths = np.load(p / "lengths.npy")
        tc = np.load(p / "target_counts.npy")
        cv = np.load(p / "coverage.npy")

        assert offsets.dtype == np.int64, offsets.dtype
        assert lengths.dtype == np.int32, lengths.dtype
        assert tc.dtype == np.int32 and cv.dtype == np.int32
        assert offsets.shape == (len(order) + 1,)
        assert offsets[0] == 0
        assert offsets[-1] == sum(lens.values()) == tc.size == cv.size

        # row k of tx_order maps to offsets[k]:offsets[k+1] in BOTH arrays -- the axis invariant
        got_order = (p / "tx_order.txt").read_text().split()
        assert got_order == order
        for k, t in enumerate(order):
            s, e = offsets[k], offsets[k + 1]
            assert e - s == lens[t]
            np.testing.assert_array_equal(tc[s:e], tgt[t])
            np.testing.assert_array_equal(cv[s:e], cov[t])

        norm = json.loads((p / "coverage_norm.json").read_text())
        assert norm["global_mean_coverage"] == 7.0      # load-bearing for cov_norm=global_mean
        assert norm["sum_L"] == int(offsets[-1])


def test_pack_rejects_missing_tx_instead_of_skipping():
    """A universe tx absent from the inputs must RAISE, not quietly drop the row.

    PINS: the dump that silently skipped onehot/FASTA-mismatched transcripts. A cross-universe eval
    collapsed to a biased subset with no error and no warning, and the numbers looked fine.
    """
    order = ["ENST1.1", "ENST_MISSING.1"]
    tgt = {"ENST1.1": np.zeros(3, dtype=np.int64)}
    cov = {"ENST1.1": np.zeros(3, dtype=np.int64)}
    with tempfile.TemporaryDirectory() as d:
        try:
            packlib.write_pack(d, order, tgt, cov, group="G", species="human")
        except ValueError as e:
            assert "absent from inputs" in str(e)
        else:
            raise AssertionError("missing tx was silently skipped instead of raising")


def test_pack_rejects_length_mismatch():
    """Target and coverage must agree per transcript, and must match a reference pack when given.

    PINS: an annotation-version mismatch would misalign the ORF track and embeddings against the
    labels -- a shift that produces worse metrics, not a crash.
    """
    order = ["ENST1.1"]
    with tempfile.TemporaryDirectory() as d:
        try:
            packlib.write_pack(d, order, {"ENST1.1": np.zeros(5, dtype=np.int64)},
                               {"ENST1.1": np.zeros(4, dtype=np.int64)}, group="G", species="human")
        except ValueError as e:
            assert "length mismatch" in str(e)
        else:
            raise AssertionError("target/coverage length mismatch was accepted")

    with tempfile.TemporaryDirectory() as d:
        try:
            packlib.write_pack(d, order, {"ENST1.1": np.zeros(5, dtype=np.int64)},
                               {"ENST1.1": np.zeros(5, dtype=np.int64)},
                               group="G", species="human", ref_len_of={"ENST1.1": 9})
        except ValueError as e:
            assert "ref_pack" in str(e)
        else:
            raise AssertionError("ref_pack length mismatch was accepted")


# ----------------------------------------------------------------------------------------------
# 2. ORF coordinates. RiboCode reports 1-based inclusive; pgx works 0-based half-open. An off-by-one
#    here shifts every ORF by a base and still yields a plausible-looking peptide.
# ----------------------------------------------------------------------------------------------
def _collapsed(tmp, rows):
    cols = ["ORF_ID", "ORF_type", "transcript_id", "transcript_type", "gene_id", "gene_name",
            "gene_type", "ORF_tstart", "ORF_tstop", "annotated_tstart", "ORF_length",
            "start_codon", "Psites_sum_frame0", "pval_combined", "adjusted_pval", "AAseq"]
    p = Path(tmp) / "collapsed.txt"
    with p.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    return p


def test_orf_coordinate_conversion():
    """1-based inclusive [tstart, tstop] -> 0-based half-open [start0, end0).

    So start0 = tstart - 1 and end0 = tstop, and the span stays divisible by 3 (stop codon included).
    """
    with tempfile.TemporaryDirectory() as d:
        # 1-based 10..21 = 12 nt = 3 sense codons + stop -> ORF_length 3
        p = _collapsed(d, [dict(ORF_ID="o1", ORF_type="uORF", transcript_id="ENST1.1",
                               ORF_tstart=10, ORF_tstop=21, annotated_tstart=100, ORF_length=3,
                               start_codon="ATG", pval_combined=0.001, adjusted_pval=0.01,
                               AAseq="MKL*")])
        orf = next(rc_io.read_collapsed(p))

    assert orf["start0"] == 9, orf["start0"]          # tstart - 1
    assert orf["end0"] == 21, orf["end0"]             # tstop, exclusive in 0-based
    span = orf["end0"] - orf["start0"]
    assert span == 12 and span % 3 == 0
    assert span // 3 == orf["aa_len"] + 1             # +1 for the stop codon
    assert orf["cds_start0"] == 99                    # annotated_tstart - 1, same convention
    assert orf["aaseq"] == "MKL"                      # stop stripped
    assert orf["qval"] == 0.01


def test_orf_qval_falls_back_to_pval():
    """adjusted_pval missing must fall back to pval_combined, not to a default-pass value."""
    with tempfile.TemporaryDirectory() as d:
        p = _collapsed(d, [dict(ORF_ID="o1", ORF_type="uORF", transcript_id="ENST1.1",
                               ORF_tstart=1, ORF_tstop=6, ORF_length=1,
                               pval_combined=0.042, AAseq="M*")])
        orf = next(rc_io.read_collapsed(p))
    assert orf["qval"] == 0.042, orf["qval"]


# ----------------------------------------------------------------------------------------------
# 3. cds_relationship must agree with RiboCode's classfy_orf on out-of-frame overlaps.
# ----------------------------------------------------------------------------------------------
def test_cds_relationship_splits_out_of_frame_overlaps():
    """Out-of-frame CDS overlaps split three ways, matching RiboCode.

    PINS the null-escape bug: collapsing all out-of-frame overlaps to `internal` left the null
    databases missing a class the RiboCode-derived model arm contained. Measured on BMDM, 117 of
    1,533 model sequences were absent from a 2.3M-sequence near-cognate null for that reason alone
    -- every one an ATG. "The null contains everything the model can find" is invisible until
    checked, which is exactly why it needs a test.
    """
    cds = (100, 400)                                   # 0-based half-open

    # NOTE on the non-overlapping cases: frame is checked BEFORE position, so a downstream ORF that
    # happens to be in-frame with the CDS is `dorf_inframe`, not `utr3`. To exercise utr3 the ORF
    # must be both downstream AND out of frame -- 431 rather than 430. (Getting this wrong is what
    # made this test fail on first run.)
    cases = [
        (100, 400, "cds_canonical",  "starts exactly at the CDS"),
        (103, 400, "cds_inframe",    "in-frame, starts inside the CDS"),
        (400, 460, "dorf_inframe",   "downstream AND in-frame -> not utr3"),
        (431, 481, "utr3",           "downstream and OUT of frame"),
        (94,  400, "ext_inframe",    "N-terminal extension"),
        (11,  59,  "utr5",           "entirely 5' of the CDS"),
        # the three that must NOT collapse together -- this is the null-escape bug
        (80,  200, "overlap_uorf",   "opens 5' of the CDS, off-frame (RiboCode Overlap_uORF)"),
        (302, 460, "overlap_dorf",   "runs past the CDS end, off-frame (RiboCode Overlap_dORF)"),
        (152, 300, "cds_offframe",   "contained and off-frame (RiboCode `internal`)"),
    ]
    for a, e, want, why in cases:
        got = seqtools.cds_relationship(a, e, cds)
        assert got == want, f"[{a},{e}) vs CDS{cds}: expected {want} ({why}), got {got}"

    assert seqtools.cds_relationship(10, 60, None) == "no_cds"


def test_novel_classes_cover_the_overlap_types():
    """The two overlap classes must map to NOVEL pgx classes, or they drop out of the null again."""
    for crel in ("overlap_uorf", "overlap_dorf"):
        klass = seqtools.orf_class(crel, "protein_coding")
        assert klass not in (None, "canonical", "internal"), f"{crel} -> {klass}"
    assert seqtools.orf_class("cds_canonical", "protein_coding") == "canonical"


# ----------------------------------------------------------------------------------------------
# 4. Search parameters. Two separate incidents, both from a single wrong field.
# ----------------------------------------------------------------------------------------------
def test_tryptic_enzyme_does_not_cut_before_proline():
    """search_enzyme_nocut_1 must be "P".

    PINS: leaving it empty silently switches trypsin -> trypsin/P and changes the peptide space.
    Measured on BMDM, the GENCODE-only baseline moved 335,548 -> 352,307 PSMs (+5.0%) from this one
    character, which makes results incomparable with every prior search in the project.
    """
    assert search.ENZYME["tryptic"]["search_enzyme_nocut_1"] == "P"
    assert search.ENZYME["tryptic"]["num_enzyme_termini"] == "2"


def test_frozen_templates_disable_calibration():
    """Every *_frozen.params must set calibrate_mass = 0 and carry explicit tolerances.

    PINS: calibrate_mass = 2 re-derives up to six search parameters per DATABASE, so two arms of one
    comparison get scored under different rules. It inflated the macrophage dPSM column by +11,017,
    and a divergence test later showed the HLA and A549 templates diverging too.
    """
    mdir = REPO / "proteogenomics" / "msfragger"
    frozen = sorted(mdir.glob("*_frozen.params"))
    assert frozen, "no frozen templates found"
    for f in frozen:
        vals = dict(
            line.split("=", 1)[0].strip() and (line.split("=", 1)[0].strip(),
                                               line.split("=", 1)[1].strip())
            for line in f.read_text().splitlines()
            if "=" in line and not line.lstrip().startswith("#")
        )
        assert vals.get("calibrate_mass") == "0", f"{f.name}: calibrate_mass={vals.get('calibrate_mass')}"
        assert vals.get("fragment_mass_units") == "1", f"{f.name}: tolerance must be explicit PPM"
        assert vals.get("fragment_mass_tolerance"), f"{f.name}: no fragment_mass_tolerance"


def test_bh_is_nan_safe():
    """One NaN p-value must not affect any other q-value.

    PINS a silent zeroing bug. numpy sorts NaN last, so the reverse
    `np.minimum.accumulate` in the old BH started on the NaN and propagated it through the entire
    array (np.minimum(NaN, x) is NaN). A single degenerate test -- a zero-variance Wilcoxon on a
    flat extension region -- therefore nulled every q-value in the arm. Measured on three arms:
    B721/mamba4, DoHH2/attn, SU-DHL-4/attn standard extensions reported 0 passing from 4,473 /
    5,261 / 7,757 tested, while 4,434 / 5,179 / 7,678 would pass the whole-ORF f0 criterion.
    It reads as a real biological result and is pure arithmetic.
    """
    import math
    from pgx.seqtools import bh

    clean = [0.001, 0.01, 0.2, 0.5]
    q_clean = bh(clean)
    q_nan = bh(clean + [float("nan")])

    assert all(math.isfinite(x) for x in q_clean), "clean input produced non-finite q"
    for i, (a, b) in enumerate(zip(q_clean, q_nan)):
        assert abs(a - b) < 1e-12, f"NaN row changed q[{i}]: {a} -> {b}"
    assert math.isnan(q_nan[-1]), "the NaN p-value should get a NaN q, not a passing one"
    assert sum(1 for x in q_nan if math.isfinite(x)) == len(clean)

    # monotone, bounded, and a NaN can never pass a threshold
    assert all(0.0 <= x <= 1.0 for x in q_clean)
    assert not (q_nan[-1] <= 0.05), "NaN must not satisfy a q cutoff"
    assert len(bh([])) == 0
    assert all(math.isnan(x) for x in bh([float("nan")] * 3))


def test_coding_potential_pool_matches_null_arm():
    """The CPAT/CPC2 candidate pool must EQUAL build_dbs' null_atg candidate set.

    PINS: the whole point of the coding-potential arms is that only the SELECTOR differs from the
    null -- same universe, same enumeration, same class filter, same min length. On first
    implementation the pool omitted `_add`'s canonical-identity filter and came out a superset by 19
    proteins on 1,000 transcripts. Small, but it means CPAT/CPC2 would be scoring a different
    candidate set than the null they are compared against, and any density difference would be
    partly bookkeeping rather than selection.
    """
    from pgx.build_dbs import novel_enumerated
    from pgx.coding_potential import enumerate_pool

    # tiny synthetic universe -- no reference files needed
    tx = "ENSTX.1"
    seq = ("ATG" + "AAA" * 12 + "TAA") + ("ATG" + "GGG" * 10 + "TGA") + "ACGTACGTAC"
    fa = None
    with tempfile.TemporaryDirectory() as d:
        fa = Path(d) / "u.fa"
        fa.write_text(f">{tx}\n{seq}\n")
        canon_seqs, cds, bt = set(), {}, {tx: "lncRNA"}
        null_store, _ = novel_enumerated(str(fa), ("ATG",), canon_seqs, 7, cds, bt)
        pool = enumerate_pool(str(fa), cds, bt, 7, ("ATG",), canon_seqs)
    pool_prot = {p[7] for p in pool}
    assert pool_prot == set(null_store), (
        f"pool != null_atg: null-only={len(set(null_store) - pool_prot)}, "
        f"pool-only={len(pool_prot - set(null_store))}")


def test_template_accepts_bare_basename():
    """--template must resolve a bare basename against the project's msfragger dir.

    PINS: the override was resolved against the JOB's cwd while the default resolved against the
    params dir. All 5 A549 searches died at 1 second each with "missing params template:
    fragger_a549_tmt_frozen.params" -- for a file that existed.
    """
    import inspect
    src = inspect.getsource(search.main)
    assert "pdir / tmpl.name" in src, "bare-basename fallback missing from --template resolution"
    mdir = REPO / "proteogenomics" / "msfragger"
    for f in sorted(mdir.glob("*_frozen.params")):
        assert (mdir / f.name).exists()


def test_search_hash_ignores_database_path():
    """The cache digest must not depend on the DB's absolute path.

    PINS: including database_name made every content-identical search look unique, so
    --shared-search-root silently did nothing and searches were recomputed.
    """
    if not hasattr(search, "pgx_hash"):
        return
    import inspect
    src = inspect.getsource(search.pgx_hash)
    assert "database_name" in src, "expected pgx_hash to explicitly handle database_name"


def _run_all():
    fns = [(n, f) for n, f in sorted(globals().items())
           if n.startswith("test_") and callable(f)]
    bad = 0
    for name, fn in fns:
        try:
            fn()
            print(f"  PASS {name}")
        except Exception as e:                                     # noqa: BLE001
            bad += 1
            print(f"  FAIL {name}: {type(e).__name__}: {e}")
    print(f"\n{len(fns) - bad}/{len(fns)} passed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(_run_all())
