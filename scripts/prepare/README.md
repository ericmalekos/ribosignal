# `prepare/` -- CLI data-prep for riboseq_signal_model packs

Build a training / eval-ready **pack** from your own data by passing paths on the command
line. No source edits, no hardcoded project root, tissue names, or sample-count asserts.

A "pack" is what the model trains and predicts on: per-nt Ribo-seq P-site **target** +
per-nt RNA-seq **coverage**, pooled onto a transcript **universe**, plus a sequence-derived
**ORF track**. The trainer/eval read a pack via `$RIBO_PACK_DIR` (+ `$RIBO_ORF_TRACK`,
`$RIBO_ONEHOT_FASTA`).

## The three tools

| Tool | Use when | Entry |
|---|---|---|
| `prepare_from_bams.py` | you have **STAR transcriptome BAMs** | `--ribo-bam`, `--rna-bam` |
| `prepare_pack.py` | you already have RiboCode `*_psites.hd5` + `*_coverage.hd5` | `--ribo-psites`, `--rna-coverage` |
| `prepare_packs.py` | **many groups** at once | `--samplesheet groups.csv` |

`paths.py` (root resolution) and `packlib.py` (shared pool+pack primitives) are libraries.

## Universe: `--ref-pack` vs `--universe-tx`

- `--ref-pack <dir>`: reuse an existing pack's `tx_order`/`lengths` verbatim. The shared ORF
  track + embeddings + one-hot FASTA then align 1:1, so you do **not** build a per-pack track.
  Use for same-annotation data (e.g. another human GENCODE v49 sample on the Fibroblast universe).
- `--universe-tx <file>`: a fresh transcript-id list (one versioned id per line, e.g. from
  `build_line_universe.py --salmon quant.sf`). Offsets are computed fresh; pass
  `--build-orf-track --fasta universe.fa` to ship a per-pack no-Kozak ORF track.

The BAM `@SQ` reference names **must** be versioned transcript ids from the **same GENCODE
annotation** as the universe, or `prepare_pack` fails loudly (annotation mismatch).

## Examples

Build a pack from STAR transcriptome BAMs (fresh universe + no-Kozak ORF track):
```bash
export RIBOCODE_BIN=/path/to/ribocode/bin
export RIBO_PYSAM_PYTHON=/path/to/riboseq_env/bin/python   # has pysam
python scripts/prepare/prepare_from_bams.py \
  --ribo-bam A.toTranscriptome.bam B.toTranscriptome.bam \
  --rna-bam  rna1.toTranscriptome.bam rna2.toTranscriptome.bam \
  --annot /path/ribocode_annot --ncrna-tx ncrna_tx.txt --tx-to-gene tx2gene.tsv \
  --universe-tx universe_tx.txt --fasta universe.fa --build-orf-track \
  --gtf gencode.v49.annotation.gtf \
  --group mycellline --species human --out packs/mycellline
# add --dry-run first to print the command plan without running anything
```

From already-pooled hd5 (skip RiboCode / coverage), reusing an existing universe:
```bash
python scripts/prepare/prepare_pack.py \
  --ribo-psites psites/ --rna-coverage coverage/ \
  --ref-pack data/packed --tx2biotype data/tx2biotype.tsv \
  --group mysample --species human --out packs/mysample
```

Swap only the RNA-seq coverage against a fixed target/universe (e.g. mm20 -> mm1):
```bash
python scripts/prepare/prepare_pack.py \
  --reuse-target data/packed_HUVEC --rna-coverage newcov/SRR*.hd5 \
  --tx2biotype data/tx2biotype.tsv --group HUVEC --species human --out /tmp/huvec_new
```

Many groups from a samplesheet (`group,assay,path`; assay in {ribo,rna}; path = BAM or hd5):
```bash
python scripts/prepare/prepare_packs.py --samplesheet groups.csv \
  --ref-pack data/packed --gtf gencode.v49.annotation.gtf --species human --outdir packs/ \
  --annot /path/ribocode_annot --ncrna-tx nc.txt --tx-to-gene t2g.tsv --ribocode-bin $RIBOCODE_BIN
```

## Then train / predict on the pack

`prepare_pack` prints the exact exports. Typical:
```bash
export RIBO_PACK_DIR=packs/mycellline
export RIBO_ORF_TRACK=packs/mycellline/orf_track_v2.npy    # per-pack track (fresh universe)
#   or the shared track for a --ref-pack build, e.g. data/packed/orf_track_v2_nokozak.npy
export RIBO_ONEHOT_FASTA=universe.fa
```
Then run `train.py` / `eval_localization.py` / `dump_pred_profiles.py` as usual.

## Environment / external tools (resolved from flags or env, never hardcoded)

- `--ribocode-bin` / `$RIBOCODE_BIN` -- dir with RiboCode `metaplots` + `RiboCode` (BAM path only).
- `--pysam-python` / `$RIBO_PYSAM_PYTHON` -- a python with `pysam` for `filter_tx_heldout.py`
  and `rnaseq_coverage.py` (BAM path only).
- `--samtools` / `$RIBO_SAMTOOLS` -- `samtools` (BAM path only; default: on `PATH`).
- Project root: `$RIBOSEQ_SIGNAL_MODEL_ROOT`, else auto-detected (the dir holding `scripts/`
  and `data/`).

## Pack contract (what gets written)

`tx_order.txt`, `offsets.npy` (int64), `lengths.npy` (int32), `target_counts.npy` (int32),
`coverage.npy` (int32), `coverage_norm.json` (`global_mean_coverage`), `pack_meta.tsv`,
`provenance.json`, and -- with `--build-orf-track` -- `orf_track_v2.npy` (float16, `--kozak
none`). All per-nt arrays are row-aligned to `tx_order`.

**Kozak:** the ORF track defaults to `--kozak none` and you should keep it there;
the hand-picked heuristic is redundant with what the sequence backbone learns (see
`results.md` Task 20).
