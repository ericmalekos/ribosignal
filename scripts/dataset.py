#!/usr/bin/env python3
"""Training data for the per-nt Ribo-seq signal model. numpy-only (no h5py), so it runs
inside the CUDA SIF.

Per transcript the model needs three per-nt-aligned arrays of length L:
  - FM embedding (L, d_emb) float16, read from the per-tx .npy indexed by tx_index.tsv
  - RNAseq coverage (L,) int32, sliced from the packed ragged coverage array
  - P-site counts  (L,) int32, sliced from the packed ragged target array (the label)
feats = concat[ emb.float32 , log1p(coverage)[:,None] ] -> (L, d_emb+1).

Splits (gene-disjoint by chromosome, from fibroblast_chrom_kfold.json): fold 0 = test,
fold 1 = val, folds 2/3/4 = train. Batching is by a token budget (sum of lengths per
batch) with length-sorted buckets, so variable-length transcripts (up to 10,000 nt) pack
with little padding and even GPU load.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paths  # noqa: E402  -- every data location resolves through here, never a cluster path

DATA = paths.data_dir()          # $RIBO_DATA_DIR, else <repo>/data -- see scripts/paths.py
PACK = paths.pack_dir()          # $RIBO_PACK_DIR, else <data>/packed
TX_INDEX_RINALMO = paths.emb_dir() / "rinalmo_token_emb" / "tx_index.tsv"
TX_INDEX_ORTHRUS = paths.emb_dir() / "orthrus4t_token_emb" / "tx_index.tsv"
TX_INDEX_HYDRARNA = paths.emb_dir() / "hydrarna_token_emb" / "tx_index.tsv"
# One-hot control backend: no on-disk per-tx emb -- a (L,4) A/C/G/T one-hot computed from the
# universe FASTA on the fly. It measures how much lift the FM embeddings provide over raw sequence
# (same architecture + ORF track + coverage, only the emb input differs). ONEHOT_SENTINEL flags it.
# RIBO_ONEHOT_FASTA overrides the universe FASTA so a cross-species held-out (e.g. mouse vM38) can
# run the one-hot backend on its own sequence, with no FM embeddings.
ONEHOT_UNIVERSE_FASTA = paths.onehot_fasta()
ONEHOT_SENTINEL = "__onehot__"
# Embedding backends: which per-token index file(s) feed the model input. "concat" stacks
# RiNALMo (1280-d) + Orthrus 4-track (512-d) along the feature axis (1792-d). All three
# resolve over the identical universe; usable() intersects across the chosen indices so
# the RiNALMo-vs-Orthrus-vs-concat comparison runs on exactly the same train/val/test tx.
BACKENDS = {
    "rinalmo": [TX_INDEX_RINALMO],
    "orthrus": [TX_INDEX_ORTHRUS],
    "hydrarna": [TX_INDEX_HYDRARNA],
    "concat": [TX_INDEX_RINALMO, TX_INDEX_ORTHRUS],
    "onehot": [ONEHOT_SENTINEL],
}
TX_INDEX = TX_INDEX_RINALMO  # backward-compat default
SPLIT = paths.split_json()
TX2BIOTYPE = paths.tx2biotype()
# Mitochondrial protein-coding genes are translated by the mitoribosome: a different genetic
# code, no 5'UTR (leaderless, so no uORFs), and none of the cytoplasmic 3-nt periodicity the
# ORF-candidate track and the FM embeddings encode. They are not this model's target -- drop the
# 13 chrM mRNAs from every split so they never enter train/val/test. (Mt_rRNA / Mt_tRNA were
# already removed upstream by the ncRNA drop-list; these 13 pc genes slipped through as
# protein_coding.) Impact on the already-reported Fibroblast medians is negligible (pc profile
# median 0.6180 -> 0.6181, uORF unaffected since utr5=0), but this keeps LOTO / any rebuild clean.
EXCLUDE_CHROMS = frozenset({"chrM"})
_EXCLUDED_TX_CACHE = None


def excluded_tx():
    """Versioned tx_ids to drop from all splits (mitochondrial chrom). Cached; empty if the
    tx2biotype table is missing (older packs)."""
    global _EXCLUDED_TX_CACHE
    if _EXCLUDED_TX_CACHE is not None:
        return _EXCLUDED_TX_CACHE
    ex = set()
    if TX2BIOTYPE.exists():
        with TX2BIOTYPE.open() as fh:
            header = fh.readline().rstrip("\n").split("\t")
            ci, ti = header.index("chrom"), header.index("tx_id")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > max(ci, ti) and f[ci] in EXCLUDE_CHROMS:
                    ex.add(f[ti])
    _EXCLUDED_TX_CACHE = ex
    return ex


_REPRESENTATIVE_TX_CACHE = None


def representative_tx():
    """Optional highest-expressed-isoform-per-gene INCLUDE-list for the posture-B multimap
    sensitivity check (results.md Task 21). Env `RIBO_REPRESENTATIVE_TX` = path to a file with one
    versioned tx_id per line (comments '#' ok); when set, every split keeps ONLY those tx (in
    addition to the chrM exclusion). Returns None when unset (normal posture-A behaviour).

    Rationale: for a gene's highest-expressed isoform, posture-A per-nt counts already EQUAL its
    posture-B (one-isoform-commit) counts -- a footprint on that isoform lands in an exon it
    contains, so committing it there does not change the representative's counts; only low-expressed
    siblings lose inherited signal, and those fall below min_signal anyway. So restricting the
    universe to representatives IS the posture-B target, with no target re-derivation."""
    global _REPRESENTATIVE_TX_CACHE
    p = os.environ.get("RIBO_REPRESENTATIVE_TX")
    if not p:
        return None
    if _REPRESENTATIVE_TX_CACHE is not None:
        return _REPRESENTATIVE_TX_CACHE
    keep = set()
    with open(p) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                keep.add(line.split()[0])
    _REPRESENTATIVE_TX_CACHE = frozenset(keep)
    return _REPRESENTATIVE_TX_CACHE


# ACGT -> 0..3, U -> T's column, upper + lower; anything else (e.g. N) -> -1 = all-zero one-hot row.
_ONEHOT_LUT = np.full(256, -1, dtype=np.int8)
for _ch, _j in {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3}.items():
    _ONEHOT_LUT[ord(_ch)] = _j
    _ONEHOT_LUT[ord(_ch.lower())] = _j


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = ln[1:].split()[0]
                buf = []
            else:
                buf.append(ln.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def onehot_encode(seq):
    """(L, 4) float32 one-hot over A,C,G,T (U->T; N / other -> all-zero row). L == len(seq), so it
    aligns 1:1 to the ORF track / coverage / target, exactly like the FM per-token embeddings."""
    codes = _ONEHOT_LUT[np.frombuffer(seq.encode("latin-1", "replace"), dtype=np.uint8)]
    a = np.zeros((codes.shape[0], 4), dtype=np.float32)
    valid = codes >= 0
    a[np.nonzero(valid)[0], codes[valid]] = 1.0
    return a


class PackedStore:
    """Memmapped packed target/coverage + embedding path index. Shared across splits."""

    def __init__(self, pack=PACK, tx_index=TX_INDEX):
        # tx_index may be one index path or a list of paths whose per-token embeddings
        # are concatenated along the feature axis (the "concat" backend).
        self.offsets = np.load(pack / "offsets.npy")
        self.lengths = np.load(pack / "lengths.npy")
        self.counts = np.load(pack / "target_counts.npy", mmap_mode="r")
        self.coverage = np.load(pack / "coverage.npy", mmap_mode="r")
        # Depth-normalizer for the coverage input: dividing per-nt coverage by the dataset's
        # global mean per-nt depth (then log1p) removes absolute sequencing depth (CPM-like), so
        # the input is comparable across datasets of different depth. A new dataset computes its
        # OWN mean over its packed coverage (see data/packed/coverage_norm.json).
        cnp = pack / "coverage_norm.json"
        self.cov_mean = (json.loads(cnp.read_text())["global_mean_coverage"] if cnp.exists()
                         else float(np.asarray(self.coverage, dtype=np.int64).sum())
                         / int(self.coverage.size))
        # (sum_L, C) ORF-candidate track, optional. RIBO_ORF_TRACK selects a variant file
        # (e.g. orf_track_v2.npy with non-AUG starts) without disturbing a running default.
        op = Path(os.environ.get("RIBO_ORF_TRACK", str(pack / "orf_track.npy")))
        self.orf_track = np.load(op, mmap_mode="r") if op.exists() else None
        self.orf_channels = int(self.orf_track.shape[1]) if self.orf_track is not None else 0
        order = (pack / "tx_order.txt").read_text().split()
        self.row = {t: i for i, t in enumerate(order)}
        paths = [tx_index] if isinstance(tx_index, str | Path) else list(tx_index)
        self.emb_index_paths = [str(p) for p in paths]
        # One-hot backend: no on-disk index -- keep every universe-FASTA tx as usable and compute
        # the (L,4) one-hot in load_emb. Value None (never dereferenced) so usable()/emb_path work.
        self.onehot = len(paths) == 1 and str(paths[0]) == ONEHOT_SENTINEL
        if self.onehot:
            self._seqs = read_fasta(ONEHOT_UNIVERSE_FASTA)
            self._emb_maps = [{t: None for t in self._seqs}]
        else:
            self._emb_maps = []
            for ip in paths:
                m = {}
                with open(ip) as fh:
                    fh.readline()
                    for line in fh:
                        t, p = line.rstrip("\n").split("\t")
                        m[t] = p
                self._emb_maps.append(m)
        self.emb_path = self._emb_maps[0]  # backward-compat: first backend's map

    def usable(self, tx_ids):
        """Keep tx present in the pack AND in every embedding index, so a single- or
        multi-backend store (and thus the split) uses one consistent tx set."""
        return [t for t in tx_ids
                if t in self.row and all(t in m for m in self._emb_maps)]

    def load_emb(self, tx):
        """(L, d_emb) float32; feature-axis concat across backends when >1 index. For the onehot
        backend d_emb == 4, computed from the transcript sequence (the raw-sequence, no-FM control)."""
        if self.onehot:
            return onehot_encode(self._seqs[tx])
        arrs = [np.load(m[tx]) for m in self._emb_maps]
        if len(arrs) == 1:
            return arrs[0].astype(np.float32)
        return np.concatenate([a.astype(np.float32) for a in arrs], axis=1)

    def length_of(self, tx):
        return int(self.lengths[self.row[tx]])


_TRAIN_EXCLUDE_CACHE = None


def train_excluded_tx():
    """Optional TRAIN+VAL-only exclusion set. Env `RIBO_TRAIN_EXCLUDE_TX` = path to a file with one
    versioned tx_id per line (comments '#' ok). Returns an empty set when unset.

    Unlike `excluded_tx()` (chrM, dropped everywhere) and `representative_tx()` (an include-list
    applied to every split), these tx are dropped from TRAIN and VAL but KEPT IN TEST. That is the
    point: a transcript whose mm1 Ribo label is missing most of its reads teaches the model an
    artifact (the GTF2I case -- a 1,617 nt CDS window reading 9% of CDS mean at mm1 and 138% at
    mm10), so it must not be a training example, but it is exactly what the resulting model has to
    be SCORED on. Keeping test intact also keeps the test split identical across arms, so an
    exclusion arm and a no-exclusion arm stay directly comparable.

    Val is excluded too, deliberately: checkpoint selection is on val Pearson, so leaving these tx
    in val would make the exclusion arm select checkpoints against a different criterion and add a
    second difference between the arms. Val metrics are therefore NOT comparable across arms; test
    metrics are.
    """
    global _TRAIN_EXCLUDE_CACHE
    if _TRAIN_EXCLUDE_CACHE is not None:
        return _TRAIN_EXCLUDE_CACHE
    path = os.environ.get("RIBO_TRAIN_EXCLUDE_TX", "")
    ex = set()
    if path:
        pth = Path(path)
        if not pth.exists():
            raise FileNotFoundError(f"RIBO_TRAIN_EXCLUDE_TX={path} does not exist")
        for line in pth.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                ex.add(line)
        if not ex:
            raise ValueError(f"RIBO_TRAIN_EXCLUDE_TX={path} contained no tx_ids")
    _TRAIN_EXCLUDE_CACHE = ex
    return ex


def load_split(store: PackedStore, test_fold=0, val_fold=None):
    """Gene-disjoint chromosome-fold split. test = fold[test_fold], val = fold[val_fold]
    (default the next fold, cyclically), train = the remaining folds."""
    d = json.loads(SPLIT.read_text())
    k = d["n_folds"]
    folds = {f["fold"]: f["test"] for f in d["folds"]}
    if val_fold is None:
        val_fold = (test_fold + 1) % k
    ex = excluded_tx()  # mitochondrial (chrM) tx -- dropped from every split
    keep = representative_tx()  # posture-B: None (all) or highest-expressed-isoform include-set
    trex = train_excluded_tx()  # train+val only; test deliberately keeps these
    ok = lambda t: t not in ex and (keep is None or t in keep)  # noqa: E731
    ok_tr = lambda t: ok(t) and t not in trex  # noqa: E731
    test = [t for t in store.usable(folds[test_fold]) if ok(t)]
    val = [t for t in store.usable(folds[val_fold]) if ok_tr(t)]
    train_folds = [j for j in range(k) if j not in (test_fold, val_fold)]
    train = [t for t in store.usable([t for j in train_folds for t in folds[j]]) if ok_tr(t)]
    return train, val, test


# ---- Leave-one-tissue-out (LOTO) ---------------------------------------------------------
# The 9 Chothani tissues share the fixed 36,668-tx universe: FM embeddings + ORF track are
# tissue-independent, so only each tissue's packed target + coverage (data/packed_<Tissue>/,
# Fibroblast = data/packed/) differ. LOTO trains on the held-in tissues and tests on the
# held-out tissue; per-tissue expression enters via the coverage input, not the transcript set.
ALL_TISSUES = ["Fibroblast", "VSMC", "ES", "Fat", "HA_EC", "Brain", "HCAEC",
               "Hepatocytes", "HUVEC"]


def tissue_pack_dir(tissue):
    """Pack dir for a tissue (Fibroblast is the original data/packed/). RIBO_PACK_SUFFIX selects a
    variant pack family: e.g. "union" -> data/packed_union (Fibroblast) + data/packed_union_<T> (others),
    for training on the broader union universe without touching the default Fibroblast-universe packs."""
    suf = os.environ.get("RIBO_PACK_SUFFIX", "")
    if suf:
        base = f"packed_{suf}"
        return DATA / (base if tissue == "Fibroblast" else f"{base}_{tissue}")
    return PACK if tissue == "Fibroblast" else DATA / f"packed_{tissue}"


# ---- External held-out datasets (cross-study / cross-species) --------------------------
# An external dataset (Ruiz-Orera human iPSC-CM, Wang mouse liver) is pooled into a pack by
# scripts/heldout/build_heldout_pack.py. There is no train/val split: a Chothani-trained model
# is applied and every scorable tx (>= min_signal pooled P-sites, chrM excluded) is a test tx.
def heldout_pack_dir(name):
    return DATA / f"packed_heldout_{name}"


def pack_test_tx(pack, min_signal=50):
    """Scorable tx in any pack (total_psites >= min_signal, chrM excluded), sorted.

    An inference-only pack built by build_pack.py without --psites has no P-site column
    worth thresholding -- every total is 0 -- so this returns nothing and the caller is
    expected to pass an explicit transcript list (the pack's expressed_tx.txt) instead.
    """
    ex = excluded_tx()
    keep = representative_tx()  # posture-B include-set (None = all)
    out = []
    with (Path(pack) / "pack_meta.tsv").open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if (int(f[c["total_psites"]]) >= min_signal and f[c["tx_id"]] not in ex
                    and (keep is None or f[c["tx_id"]] in keep)):
                out.append(f[c["tx_id"]])
    return sorted(out)


def heldout_test_tx(name, min_signal=50):
    """Scorable tx of the named held-out dataset (data/packed_heldout_<name>/)."""
    return pack_test_tx(heldout_pack_dir(name), min_signal)


def tissue_signal(tissue):
    """tx_id -> total pooled P-sites in this tissue, from its pack_meta.tsv."""
    d = {}
    with (tissue_pack_dir(tissue) / "pack_meta.tsv").open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[c["tx_id"]]] = int(f[c["total_psites"]])
    return d


def loto_split(holdout, train_tissues=None, min_signal=50, val_fold=0):
    """Tissue-held-out split. Returns (train, val, test_tissue, test_tx):
      train / val : list of (tissue, [tx_ids]) over the held-in tissues, where a tx is kept
                    if it has >= min_signal pooled P-sites in that tissue (learn shape where
                    shape exists) and is not chrM. val = the val_fold chromosomes (gene-disjoint
                    from train for honest early stopping, spanning the training tissues);
                    train = the other chromosomes.
      test_tx     : the held-out tissue's scorable tx (>= min_signal, not chrM), all chromosomes.
    """
    if holdout not in ALL_TISSUES:
        raise ValueError(f"unknown holdout tissue {holdout!r}; expected one of {ALL_TISSUES}")
    if train_tissues is None:
        train_tissues = [t for t in ALL_TISSUES if t != holdout]
    ex = excluded_tx()
    keep = representative_tx()  # posture-B include-set (None = all)
    trex = train_excluded_tx()  # train+val only; test_tx deliberately keeps these
    d = json.loads(SPLIT.read_text())
    folds = {f["fold"]: set(f["test"]) for f in d["folds"]}
    val_chroms = folds[val_fold]

    def scorable(tissue, drop_train_excluded=False):
        return {tx for tx, n in tissue_signal(tissue).items()
                if n >= min_signal and tx not in ex and (keep is None or tx in keep)
                and not (drop_train_excluded and tx in trex)}

    train, val = [], []
    for t in train_tissues:
        sc = scorable(t, drop_train_excluded=True)
        tr = sorted(tx for tx in sc if tx not in val_chroms)
        va = sorted(tx for tx in sc if tx in val_chroms)
        train.append((t, tr))
        val.append((t, va))
    test_tx = sorted(scorable(holdout))
    return train, val, holdout, test_tx


class RiboDataset(Dataset):
    def __init__(self, store: PackedStore, tx_ids, input_mode="both", use_orf=False,
                 cov_norm="raw"):
        # input_mode: "both" | "emb" (zero the coverage channel) | "cov" (zero the
        # embedding channels) -- for input-ablation runs.
        # use_orf: append the sequence-derived ORF-candidate track between the embedding and
        # the coverage channel (coverage stays last). Requires store.orf_track present.
        # cov_norm: "raw" (log1p of raw pooled depth; legacy) or "global_mean" (divide by the
        # dataset global mean depth before log1p -> depth-normalized, transferable input).
        assert input_mode in ("both", "emb", "cov")
        assert cov_norm in ("raw", "global_mean")
        self.store = store
        self.tx_ids = tx_ids
        self.input_mode = input_mode
        self.use_orf = use_orf and store.orf_track is not None
        self.cov_norm = cov_norm

    def __len__(self):
        return len(self.tx_ids)

    def __getitem__(self, k):
        s = self.store
        tx = self.tx_ids[k]
        r = s.row[tx]
        a, b = s.offsets[r], s.offsets[r + 1]
        counts = np.asarray(s.counts[a:b], dtype=np.float32)
        cov = np.asarray(s.coverage[a:b], dtype=np.float32)
        emb = s.load_emb(tx)                                  # (L, d_emb) float32
        L = counts.shape[0]
        if emb.shape[0] != L:
            raise ValueError(f"{tx}: emb L={emb.shape[0]} != target L={L}")
        d = emb.shape[1]
        c = s.orf_channels if self.use_orf else 0
        feats = np.empty((L, d + c + 1), dtype=np.float32)
        feats[:, :d] = 0.0 if self.input_mode == "cov" else emb
        if c:
            orf = np.asarray(s.orf_track[a:b], dtype=np.float32)   # (L, C) 0/1
            # ORF is sequence-derived like the embedding: zero it in the coverage-only ablation
            feats[:, d:d + c] = 0.0 if self.input_mode == "cov" else orf
        if self.input_mode == "emb":
            feats[:, -1] = 0.0
        else:
            cn = cov / s.cov_mean if self.cov_norm == "global_mean" else cov
            feats[:, -1] = np.log1p(cn)
        return {"feats": torch.from_numpy(feats),
                "counts": torch.from_numpy(counts),
                "tx": tx, "length": L}


def collate_pad(items):
    B = len(items)
    Lmax = max(it["length"] for it in items)
    C = items[0]["feats"].shape[1]
    feats = torch.zeros(B, Lmax, C, dtype=torch.float32)
    counts = torch.zeros(B, Lmax, dtype=torch.float32)
    mask = torch.zeros(B, Lmax, dtype=torch.bool)
    txs, lens = [], []
    for i, it in enumerate(items):
        L = it["length"]
        feats[i, :L] = it["feats"]
        counts[i, :L] = it["counts"]
        mask[i, :L] = True
        txs.append(it["tx"])
        lens.append(L)
    return {"feats": feats, "counts": counts, "mask": mask, "tx": txs, "lengths": lens}


class TokenBudgetSampler(torch.utils.data.Sampler):
    """Yield batches (index lists) whose summed length <= budget, length-bucketed to
    minimize padding; batch order shuffled each epoch. Caps per-batch B so a pile of very
    short transcripts doesn't make a huge batch."""

    def __init__(self, lengths, budget=24000, max_batch=64, shuffle=True, seed=0):
        self.lengths = list(lengths)
        self.budget = budget
        self.max_batch = max_batch
        self.shuffle = shuffle
        self.epoch = 0
        self.seed = seed
        self._batches = self._build()

    def _build(self):
        order = sorted(range(len(self.lengths)), key=lambda i: self.lengths[i])
        batches, cur, cur_max = [], [], 0
        for i in order:
            li = self.lengths[i]
            newmax = max(cur_max, li)
            if cur and (newmax * (len(cur) + 1) > self.budget
                        or len(cur) >= self.max_batch):
                batches.append(cur)
                cur, cur_max = [i], li
            else:
                cur.append(i)
                cur_max = newmax
        if cur:
            batches.append(cur)
        return batches

    def set_epoch(self, e):
        self.epoch = e

    def __iter__(self):
        batches = self._batches
        if self.shuffle:
            g = torch.Generator().manual_seed(self.seed + self.epoch)
            perm = torch.randperm(len(batches), generator=g).tolist()
            batches = [batches[i] for i in perm]
        yield from batches

    def __len__(self):
        return len(self._batches)


if __name__ == "__main__":
    # smoke against whatever pack exists (dry-run pack acceptable)
    store = PackedStore()
    train, val, test = load_split(store)
    print(f"train {len(train):,}  val {len(val):,}  test {len(test):,}")
    ds = RiboDataset(store, train[:50])
    lens = [store.length_of(t) for t in train[:50]]
    samp = TokenBudgetSampler(lens, budget=24000)
    print(f"50 tx -> {len(samp)} token-budget batches")
    batch = collate_pad([ds[i] for i in next(iter(samp))])
    print("batch feats", tuple(batch["feats"].shape), "mask true",
          int(batch["mask"].sum()), "counts total", float(batch["counts"].sum()))
