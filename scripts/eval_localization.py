#!/usr/bin/env python3
"""Localization eval: does the model concentrate predicted P-site signal on TRANSLATED ORFs?

Whole-transcript profile Pearson is CDS-dominated, so it never asks the sharpest question this
project cares about: on a held-out transcript, does the model localize predicted Ribo-seq signal
to the specific ORFs that are genuinely translated -- especially non-canonical ones (uORF / dORF /
novel) -- rather than to the many untranslated candidate AUG ORFs sharing the same transcript?

Truth set = RiboCode ORF calls (<Tissue>_collapsed.txt): 1-based ORF_tstart (the AUG) and ORF_tstop
(1-based, inclusive of the last stop-codon base -- both conventions verified: AUG at [tstart-1,
tstart+2) and stop at [tstop-3, tstop), 100% of calls), typed annotated / uORF / dORF / novel /
internal / Overlap_uORF / Overlap_dORF. RiboCode here was run AUG-only, so every truth ORF starts
at AUG; this is therefore not a non-AUG-start DISCOVERY test but a translation-localization test.

Candidates: every maximal AUG..in-frame-stop ORF on a test transcript (the same construction as the
orf_track occupancy channels), length >= --min_orf_nt. A candidate is POSITIVE if its 0-based start
matches a RiboCode ORF start on that transcript (carrying that ORF's type), else NEGATIVE. Crucially
the orf_track marks candidate positions + reading frames but NOT which candidates translate -- so
separating positives from negatives has to come from the model (FM embedding + Ribo-seq-trained
profile), not from the track it was given.

Per candidate, two length-fair scores from the PREDICTED per-nt profile p (softmax over the tx):
  density = mean_{i in [a,e)} p[i]                                   per-nt predicted P-site density
  frame0  = sum_{i in [a,e), (i-a)%3==0} p[i] / sum_{i in [a,e)} p[i]   periodicity signature (~1/3
            if the model sees no translation there, ->1 if it predicts a clean in-frame ORF)
where [a, e) = [tstart-1, tstop) is the ORF body incl. stop (length divisible by 3).

Metric 1 (discrimination): AUROC of each score separating translated from untranslated candidates,
overall and per ORF_type, with an ORF-length AUROC floor (guards the length confound) and an
observed-profile AUROC ceiling. Metric 2 (fidelity, positives only): predicted-vs-observed frame0
and density correlation, and median predicted frame0 per ORF_type.

Held-out only: within-tissue runs use the fold's test tx + Fibroblast calls; LOTO runs use the
held-out tissue's test tx + that tissue's calls. The mature-mRNA sequence (hence every candidate
ORF) is tissue-independent, so one universe FASTA serves all tissues. Reuses eval_extra's model
load + batched inference. Writes <run>/localization_metrics.json + <run>/localization_orfs.tsv.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dataset import (  # noqa: E402
    BACKENDS,
    PackedStore,
    RiboDataset,
    TokenBudgetSampler,
    collate_pad,
    heldout_pack_dir,
    heldout_test_tx,
    load_split,
    loto_split,
    tissue_pack_dir,
)
from model import RiboSignalModel  # noqa: E402
from train import load_biotype, pearson, spearman  # noqa: E402

NEW = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/riboseq_signal_model")
ECH = Path("/private/groups/carpenterlab/emalekos/RNAZoo_meta/"
           "RNAZoo/experiments/biotype_probe/expression_context_human")
FASTA = NEW / "data" / "fibroblast_universe.fa"
RIBOCODE = ECH / "data" / "ribocode_per_tissue"
STOP = frozenset({"TAA", "TAG", "TGA"})
# ORF_type groups: canonical (annotated CDS) vs the non-canonical classes this project targets.
NONCANON = ("uORF", "Overlap_uORF", "dORF", "Overlap_dORF", "novel", "internal")


def read_fasta(path):
    seqs, h, parts = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if h is not None:
                    seqs[h] = "".join(parts)
                h = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if h is not None:
        seqs[h] = "".join(parts)
    return seqs


def load_ribocode_orfs(tissue, path=None):
    """tx_id -> {start0: orf_type} for RiboCode ORF calls on this tissue.

    Verified conventions: ORF_tstart is 1-based (AUG at [ts-1, ts+2)); ORF_tstop is 1-based
    inclusive of the last stop base; body incl. stop is [ts-1, tstop). start0 = ts-1.

    path overrides the default per-tissue location (external held-out collapsed.txt).
    """
    path = Path(path) if path is not None else RIBOCODE / tissue / f"{tissue}_collapsed.txt"
    out = {}
    with path.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        c = {k: hdr.index(k) for k in ("transcript_id", "ORF_type", "ORF_tstart", "ORF_tstop")}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            tx = f[c["transcript_id"]]
            start0 = int(f[c["ORF_tstart"]]) - 1
            out.setdefault(tx, {})[start0] = f[c["ORF_type"]]
    return out


def load_cds_bounds():
    """tx_id -> (cds_start, cds_end) 0-based half-open, from data/tx2cds.tsv (coding tx only).
    cds_start = utr5_len (0-based position of the CDS's first base = the annotated AUG)."""
    out = {}
    # RIBO_TX2CDS overrides the default human table so a cross-species held-out uses its own CDS.
    p = Path(os.environ.get("RIBO_TX2CDS", str(NEW / "data" / "tx2cds.tsv")))
    if not p.exists():
        return out
    with p.open() as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {k: h.index(k) for k in ("tx_id", "utr5_len", "cds_len")}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            u5, cl = int(f[ix["utr5_len"]]), int(f[ix["cds_len"]])
            if cl > 0:
                out[f[ix["tx_id"]]] = (u5, u5 + cl)
    return out


def cds_relationship(a, e, cds):
    """Where does candidate ORF [a, e) sit relative to the annotated CDS [cs, ce)? The ONLY class
    collapsed (not a distinct ORF) is:

      cds_inframe   INTERNAL in-frame alt start: in the CDS frame, starting strictly inside the CDS
                    (cs < a < ce) -- a downstream AUG re-initiating the same protein truncated.

    Everything else is kept as a genuine candidate:
      cds_canonical the annotated CDS start itself (a == cs, in frame).
      dorf_inframe  in the CDS frame but starting past the stop (a >= ce): an in-frame dORF (KEEP --
                    dORFs are genuine alt-ORFs regardless of frame).
      ext_inframe   in the CDS frame but starting 5' of the CDS (a < cs) and overlapping it: an
                    N-terminal EXTENSION (the key case is a non-AUG CUG extension, MYC / VEGFA).
      cds_offframe  overlaps the CDS in a DIFFERENT frame (out-of-frame overlapping ORF, KEEP).
      utr5 / utr3   entirely 5' / 3' of the CDS (off-frame 3'UTR dORFs land here). no_cds = lncRNA.
    """
    if cds is None:
        return "no_cds"
    cs, ce = cds
    inframe = (a - cs) % 3 == 0
    if inframe and a >= ce:                           # in CDS frame, past the stop = in-frame dORF
        return "dorf_inframe"
    if inframe and a > cs:                            # in CDS frame, strictly inside the CDS
        return "cds_inframe"                          # <- internal in-frame alt start (collapsed)
    if inframe and a == cs:                           # the annotated CDS start
        return "cds_canonical"
    if max(a, cs) < min(e, ce):                       # overlaps the CDS extent
        return "ext_inframe" if inframe else "cds_offframe"   # a < cs, in-frame = N-term extension
    if e <= cs:
        return "utr5"
    if a >= ce:
        return "utr3"
    return "other"


def candidate_orfs(seq, min_nt):
    """All maximal AUG..in-frame-stop ORFs: list of (start0, end0_excl). Every in-frame AUG is a
    candidate start (so nested/internal AUGs are included, matching RiboCode's 'internal' type);
    end0_excl = one past the stop codon, so [start0, end0_excl) has length divisible by 3."""
    s = seq.upper().replace("U", "T")
    L = len(s)
    if L < 6:
        return []
    b = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
    A, G, T = 65, 71, 84
    orfs = []
    for f in range(3):
        n = (L - f) // 3
        if n < 2:
            continue
        idx = f + 3 * np.arange(n)                       # codon-start positions in this frame
        c0, c1, c2 = b[idx], b[idx + 1], b[idx + 2]
        is_atg = (c0 == A) & (c1 == T) & (c2 == G)
        is_stop = ((c0 == T) & (c1 == A) & (c2 == A)) | \
                  ((c0 == T) & (c1 == A) & (c2 == G)) | \
                  ((c0 == T) & (c1 == G) & (c2 == A))
        stop_codon_pos = idx[is_stop]                    # 0-based stop positions (sorted, in-frame)
        for a in idx[is_atg]:
            j = int(np.searchsorted(stop_codon_pos, a, side="right"))  # first stop strictly > a
            if j >= stop_codon_pos.size:
                continue
            e = int(stop_codon_pos[j]) + 3               # one past the stop
            if e - a >= min_nt:
                orfs.append((int(a), e))
    return orfs


def _rankdata(a):
    a = np.asarray(a, dtype=np.float64)
    order = np.argsort(a, kind="mergesort")
    r = np.empty(len(a), dtype=np.float64)
    sa = a[order]
    i, n = 0, len(a)
    while i < n:
        j = i
        while j < n and sa[j] == sa[i]:
            j += 1
        r[order[i:j]] = (i + j - 1) / 2.0 + 1.0          # average rank (1-based) over the tie block
        i = j
    return r


def auroc(scores, labels):
    """AUROC via the Mann-Whitney rank-sum identity (ties handled by average ranks)."""
    scores = np.asarray(scores, dtype=np.float64)
    labels = np.asarray(labels, dtype=bool)
    keep = scores == scores                              # drop nan scores
    scores, labels = scores[keep], labels[keep]
    n_pos = int(labels.sum())
    n_neg = int((~labels).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    r = _rankdata(scores)
    return float((r[labels].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def auroc_stratified(scores, labels, strat, n_bins=5):
    """Length-controlled AUROC: bin candidates by a covariate (ORF length) into quantile strata,
    take the Mann-Whitney AUROC within each stratum (only translated-vs-untranslated pairs of
    SIMILAR length are compared), and average weighted by the pair count per stratum. If this
    stays well above 0.5 the model's score discriminates beyond the length confound."""
    scores = np.asarray(scores, dtype=np.float64)
    labels = np.asarray(labels, dtype=bool)
    strat = np.asarray(strat, dtype=np.float64)
    keep = scores == scores
    scores, labels, strat = scores[keep], labels[keep], strat[keep]
    if labels.sum() == 0 or (~labels).sum() == 0:
        return float("nan")
    qs = np.quantile(strat, np.linspace(0, 1, n_bins + 1))
    qs[-1] += 1e-9
    num = den = 0.0
    for k in range(n_bins):
        m = (strat >= qs[k]) & (strat < qs[k + 1])
        yk = labels[m]
        npos, nneg = int(yk.sum()), int((~yk).sum())
        if npos == 0 or nneg == 0:
            continue
        w = npos * nneg
        num += auroc(scores[m], yk) * w
        den += w
    return float(num / den) if den > 0 else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--min_orf_nt", type=int, default=30, help="min candidate ORF length (nt)")
    ap.add_argument("--min_body_signal", type=int, default=5,
                    help="observed P-sites in the ORF body to keep a POSITIVE (drops called "
                         "ORFs with no signal on this tx's pooled target)")
    ap.add_argument("--no_collapse_inframe_cds", action="store_true",
                    help="keep in-frame AUG ORFs nested in an annotated CDS as separate "
                         "candidates (default collapses them -- 5'-truncations of the same "
                         "protein, not distinct ORFs). Collapse is AUG-specific: in-frame "
                         "non-AUG extensions (CUG MYC/VEGFA) are real alt-ORFs.")
    ap.add_argument("--eval_cap", type=int, default=0, help="0 = all test tx")
    ap.add_argument("--heldout", default=None,
                    help="external held-out dataset name (packed_heldout_<name>); scores the "
                         "trained model's localization against that dataset's RiboCode calls")
    ap.add_argument("--calls", default=None,
                    help="explicit RiboCode collapsed.txt for the held-out (default "
                         "data/heldout_psites/<heldout>/<heldout>_collapsed.txt)")
    ap.add_argument("--min_signal", type=int, default=50,
                    help="held-out scorable floor: min pooled P-sites for a test tx")
    ap.add_argument("--fasta", default=None,
                    help="override the universe FASTA (cross-species held-out)")
    ap.add_argument("--out", default=None, help="output dir for metrics/tsv (default <run>)")
    ap.add_argument("--budget", type=int, default=24000)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    run = Path(args.run)
    if not run.is_absolute():
        run = NEW / run
    cfg = json.loads((run / "args.json").read_text())
    device = args.device if (args.device != "cuda" or torch.cuda.is_available()) else "cpu"
    backend = cfg.get("emb_backend", "rinalmo")
    biotype = load_biotype()

    if args.heldout:
        tissue = args.heldout
        store = PackedStore(pack=heldout_pack_dir(args.heldout), tx_index=BACKENDS[backend])
        test_ids = store.usable(heldout_test_tx(args.heldout, args.min_signal))
        calls_path = args.calls or (NEW / "data" / "heldout_psites" / args.heldout
                                    / f"{args.heldout}_collapsed.txt")
    elif cfg.get("holdout"):
        tissue = cfg["holdout"]
        store = PackedStore(pack=tissue_pack_dir(tissue), tx_index=BACKENDS[backend])
        _, _, _, test_ids = loto_split(
            tissue, train_tissues=cfg.get("resolved_train_tissues"),
            min_signal=cfg.get("min_train_signal", 50), val_fold=cfg.get("val_fold", 0))
        calls_path = None
    else:
        tissue = "Fibroblast"
        store = PackedStore(tx_index=BACKENDS[backend])
        _, _, test_ids = load_split(store, test_fold=cfg.get("test_fold", 0))
        calls_path = None

    seqs = read_fasta(Path(args.fasta) if args.fasta else FASTA)
    ribocode = load_ribocode_orfs(tissue, path=calls_path)
    cds_bounds = load_cds_bounds()
    print(f"run={run.name} tissue={tissue} backend={backend} test_tx={len(test_ids)} "
          f"device={device}", file=sys.stderr)

    use_orf = cfg.get("use_orf_track", False)
    ds = RiboDataset(store, test_ids, input_mode=cfg.get("input_mode", "both"),
                     use_orf=use_orf, cov_norm=cfg.get("cov_norm", "raw"))
    lens = [store.length_of(t) for t in test_ids]
    eval_budget = min(args.budget, cfg.get("budget", args.budget))
    samp = TokenBudgetSampler(lens, budget=eval_budget, shuffle=False)
    d_emb = store.load_emb(test_ids[0]).shape[1]
    extra_in = store.orf_channels if use_orf else 0
    model = RiboSignalModel(d_emb=d_emb, channels=cfg["channels"], n_blocks=cfg["n_blocks"],
                            dropout=0.0, extra_in=extra_in,
                            n_attn_layers=cfg.get("n_attn_layers", 0),
                            n_heads=cfg.get("n_heads", 8),
                            mixer=cfg.get("mixer", "transformer"),
                            mamba_d_state=cfg.get("mamba_d_state", 16),
                            mamba_d_conv=cfg.get("mamba_d_conv", 4),
                            mamba_expand=cfg.get("mamba_expand", 2),
                            learn_start_context=cfg.get("learn_start_context", False),
                            fm_to_mixer=cfg.get("fm_to_mixer", False)).to(device)
    model.load_state_dict(torch.load(run / "best.pt", map_location=device))
    model.eval()

    # per candidate ORF row: (orf_type|"none", is_pos, length, pred_density, pred_frame0,
    #                          obs_density, obs_frame0, body_signal, biotype)
    rows = []
    seen = 0
    with torch.no_grad():
        for batch in torch.utils.data.DataLoader(ds, batch_sampler=samp,
                                                 collate_fn=collate_pad, num_workers=4):
            feats = batch["feats"].to(device)
            mask = batch["mask"].to(device)
            logits, _ = model(feats, mask)
            logp = torch.log_softmax(logits.masked_fill(~mask, -1e9), dim=1)
            P = torch.exp(logp).cpu().numpy()
            counts = batch["counts"].numpy()
            for i, tx in enumerate(batch["tx"]):
                L = int(batch["lengths"][i])
                seq = seqs.get(tx)
                if seq is None or len(seq) != L:
                    continue
                p = P[i, :L]
                c = counts[i, :L].astype(np.float64)
                cands = candidate_orfs(seq, args.min_orf_nt)
                if not cands:
                    continue
                truth = ribocode.get(tx, {})
                bt = biotype.get(tx, "unknown")
                cds = cds_bounds.get(tx)
                for a, e in cands:
                    if e > L:
                        continue
                    crel = cds_relationship(a, e, cds)
                    pm = p[a:e]
                    cm = c[a:e]
                    psum = pm.sum()
                    dens = float(pm.mean())
                    f0 = float(pm[0::3].sum() / psum) if psum > 0 else float("nan")
                    csum = cm.sum()
                    obs_dens = float(cm.mean())
                    obs_f0 = float(cm[0::3].sum() / csum) if csum > 0 else float("nan")
                    otype = truth.get(a)
                    is_pos = otype is not None
                    # a POSITIVE with no observed signal on this tx's pooled target is not a fair
                    # localization target here; drop it (still counts elsewhere in RiboCode's pooled
                    # multi-tx evidence). Negatives are kept regardless (that is the candidate sea).
                    if is_pos and csum < args.min_body_signal:
                        continue
                    rows.append((otype or "none", is_pos, e - a, dens, f0,
                                 obs_dens, obs_f0, int(csum), bt, tx, a, e, crel))
            seen += len(batch["tx"])
            if args.eval_cap and seen >= args.eval_cap:
                break

    # Collapse INTERNAL in-frame AUG alt starts only (default ON): cds_rel == "cds_inframe" is an
    # AUG strictly inside the CDS in its frame -- the same protein re-initiated truncated, not a
    # distinct ORF. Dropped, except an "annotated" positive (kept even if its start looks internal
    # -- a RiboCode / tx2cds coordinate mismatch). Everything else stays: in-frame + off-frame
    # dORFs, out-of-frame overlaps, upstream extensions, uORFs, lncRNA. Full candidate list is still
    # dumped to the tsv. AUG-specific: an in-frame non-AUG start 5' of the CDS (CUG MYC / VEGFA) is
    # a real alt-ORF, so the collapse stays gated on AUG once a non-AUG caller is added.
    collapse = not args.no_collapse_inframe_cds

    def keep_row(r):
        return not (collapse and r[12] == "cds_inframe" and not (r[1] and r[0] == "annotated"))
    metric_rows = [r for r in rows if keep_row(r)]
    n_collapsed = len(rows) - len(metric_rows)

    n_pos = sum(1 for r in metric_rows if r[1])
    n_neg = len(metric_rows) - n_pos
    print(f"candidates: {len(rows):,} enumerated, {n_collapsed:,} in-frame-CDS collapsed -> "
          f"{len(metric_rows):,} scored  pos {n_pos:,}  neg {n_neg:,}", file=sys.stderr)

    lab = np.array([r[1] for r in metric_rows], dtype=bool)
    length = np.array([r[2] for r in metric_rows], dtype=np.float64)
    pdens = np.array([r[3] for r in metric_rows], dtype=np.float64)
    pf0 = np.array([r[4] for r in metric_rows], dtype=np.float64)
    odens = np.array([r[5] for r in metric_rows], dtype=np.float64)
    of0 = np.array([r[6] for r in metric_rows], dtype=np.float64)
    otypes = np.array([r[0] for r in metric_rows], dtype=object)

    def disc(pos_mask):
        """AUROC of each score, pos_mask vs all negatives (untranslated candidate ORFs)."""
        sel = pos_mask | (~lab)
        y = lab[sel]
        return {
            "n_pos": int(pos_mask.sum()), "n_neg": int((~lab).sum()),
            "auroc_pred_density": auroc(pdens[sel], y),
            "auroc_pred_frame0": auroc(pf0[sel], y),
            # length-controlled: model score must beat the length floor WITHIN length strata
            "auroc_pred_frame0_lengthctrl": auroc_stratified(pf0[sel], y, length[sel]),
            "auroc_pred_density_lengthctrl": auroc_stratified(pdens[sel], y, length[sel]),
            "auroc_orf_length_floor": auroc(length[sel], y),
            "auroc_obs_density_ceiling": auroc(odens[sel], y),
            "auroc_obs_frame0_ceiling": auroc(of0[sel], y),
            # length-controlled observed ceiling: the fair reference for the length-controlled
            # predicted scores (both compare translated vs untranslated within length strata)
            "auroc_obs_frame0_ceiling_lengthctrl": auroc_stratified(of0[sel], y, length[sel]),
        }

    discrimination = {"all": disc(lab)}
    discrimination["noncanonical"] = disc(lab & np.isin(otypes, NONCANON))
    for t in ("annotated",) + NONCANON:
        m = lab & (otypes == t)
        if m.sum() >= 20:
            discrimination[t] = disc(m)

    # Metric 2: fidelity on positives -- does predicted periodicity/density track the observed?
    pos = lab
    def corr_block(mask, label):
        pv, ov = pf0[mask], of0[mask]
        k = (pv == pv) & (ov == ov)
        if k.sum() < 8:
            return {"label": label, "n": int(k.sum())}
        return {
            "label": label, "n": int(mask.sum()),
            "median_pred_frame0": float(np.median(pv[pv == pv])),
            "median_obs_frame0": float(np.median(ov[ov == ov])),
            "frame0_pearson": pearson(pv[k], ov[k]),
            "density_pearson": pearson(pdens[mask][k], odens[mask][k]),
            "density_spearman": spearman(pdens[mask][k], odens[mask][k]),
        }

    fidelity = {"all_positives": corr_block(pos, "all positives"),
                "noncanonical": corr_block(pos & np.isin(otypes, NONCANON), "non-canonical"),
                "annotated": corr_block(pos & (otypes == "annotated"), "annotated")}
    for t in NONCANON:
        m = pos & (otypes == t)
        if m.sum() >= 20:
            fidelity[t] = corr_block(m, t)

    out = {
        "run": run.name, "tissue": tissue, "emb_backend": backend,
        "use_orf_track": use_orf, "n_attn_layers": cfg.get("n_attn_layers", 0),
        "min_orf_nt": args.min_orf_nt, "min_body_signal": args.min_body_signal,
        "collapse_inframe_cds": collapse, "n_inframe_cds_collapsed": n_collapsed,
        "n_candidates_enumerated": len(rows),
        "n_candidates": len(metric_rows), "n_positive": n_pos, "n_negative": n_neg,
        "discrimination": discrimination,
        "fidelity": fidelity,
    }
    outdir = Path(args.out) if args.out else run
    outdir.mkdir(parents=True, exist_ok=True)
    suffix = f"_{args.heldout}" if args.heldout else ""
    metrics_path = outdir / f"localization_metrics{suffix}.json"
    orfs_path = outdir / f"localization_orfs{suffix}.tsv"
    metrics_path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))

    with orfs_path.open("w") as fh:
        fh.write("orf_type\tis_translated\tlength_nt\tpred_density\tpred_frame0\t"
                 "obs_density\tobs_frame0\tbody_psites\tbiotype\ttx_id\torf_start\torf_end\tcds_rel\n")
        for r in rows:
            fh.write(f"{r[0]}\t{int(r[1])}\t{r[2]}\t{r[3]:.8f}\t"
                     f"{r[4] if r[4] == r[4] else float('nan'):.6f}\t"
                     f"{r[5]:.6f}\t{r[6] if r[6] == r[6] else float('nan'):.6f}\t{r[7]}\t{r[8]}\t"
                     f"{r[9]}\t{r[10]}\t{r[11]}\t{r[12]}\n")
    print(f"wrote {metrics_path} and {orfs_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
