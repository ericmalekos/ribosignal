# Training data

The model is trained on the **Chothani et al. human Ribo-seq + matched RNA-seq atlas** (GENCODE v49), across
eight tissues / primary-cell types. Each tissue contributes one **pack**: aligned, pooled, per-nucleotide
P-site counts (the label) and RNA coverage (an input) on a shared transcript universe.

## Tissues

```{list-table}
:header-rows: 1
:widths: 40 20 40

* - Held-in (train + validation)
  - Held-out (test)
  - Dropped
* - Fibroblast, VSMC, ES, Fat, HA_EC, HCAEC, HUVEC
  - **Hepatocytes**
  - Brain (noisy / leaderless)
```

Training holds out **Hepatocytes** entirely (leave-one-tissue-out, LOTO). Within the held-in tissues, the
train/validation split is **gene-disjoint by chromosome** so validation genes are never seen in training.
Mitochondrial (chrM) transcripts are excluded everywhere -- the mitoribosome uses a different code.

## The pack (what the model reads)

`scripts/prepare/` turns aligned BAMs into a pack: everything is a flat array indexed by `tx_order.txt`.

```{list-table}
:header-rows: 1
:widths: 30 30 40

* - File
  - Shape
  - Content
* - `tx_order.txt`
  - n_tx
  - versioned transcript IDs, defining row order
* - `offsets.npy` / `lengths.npy`
  - n_tx (+1)
  - where each transcript starts / its length in the flat arrays
* - `target_counts.npy`
  - sum_L
  - per-nt **P-sites** (the label; pooled over the tissue's Ribo-seq reps)
* - `coverage.npy`
  - sum_L
  - per-nt **RNA depth** (the input; unique-mapper STAR, "mm1")
* - `pack_meta.tsv`
  - per-tx
  - biotype, total P-sites (a transcript is *scorable* at >= 50 P-sites)
* - ORF track
  - sum_L x 5
  - the reading-frame prior (shared across tissues; `--kozak none`)
```

```{admonition} mm1 RNA coverage
:class: note
Coverage uses **unique mappers only** (`--outFilterMultimapNmax 1`, "mm1"). Multi-mapper coverage inflates
isoform signal ~18x and misleads the count head; per-transcript normalisation and mm1 both guard against it.
```

## Transcript universe: fibroblast vs union

Two universes define *which* transcripts the model trains and is scored on. Both are protein-coding + lncRNA,
<= 10 kb, no chrM.

```{list-table}
:header-rows: 1
:widths: 22 18 22 38

* - Universe
  - Transcripts
  - Hepatocytes test set
  - Definition
* - **Fibroblast** (default)
  - 36,668
  - 33,918
  - salmon TPM >= 1 **in Fibroblast**
* - **Union** (broad)
  - 84,472
  - 70,883
  - salmon TPM >= 1 in **any** of the 8 tissues
```

The union universe adds ~48k transcripts -- including ~2,100 expressed only in Hepatocytes -- so it can
**discover ~80% more novel ORFs** on the held-out tissue, at a small canonical-precision cost. A universe is
**never reused across datasets**: each new dataset gets its own salmon quantification, or it imports that
dataset's expression bias. See {doc}`test-sets`.
