# S_metagene -- positional P-site metagene, start and stop anchors, 8 training tissues

**What it shows:** the plot every Ribo-seq audience expects and that this project did not have. Density
jumps at the start codon, runs in 3-nt periodicity through the CDS, pauses at the stop codon and
falls off in the 3'UTR. `S_riboseq_qc` had per-library QC *scalars* (read-length distribution,
frame-0 fraction, P-site offset per read length) but no positional profile: RiboCode's `metaplots`
computes these quantities to choose offsets and does not persist them.

## THE POSITION IS THE P-SITE, NOT THE READ 5' END

State this in any caption. `target_counts.npy` in each pack is RiboCode `process_bam` output with the
per-read-length P-site offset **already applied**, so:

- **start anchor, offset 0** = the first nucleotide of the start codon (the A of AUG)
- **stop anchor, offset 0** = the first nucleotide of the stop codon

A 5'-end metagene would sit roughly 12 nt to the left of this one. That offset 0 is also where the
peak actually lands is an independent check that the offsets were applied correctly.

Coordinates are **transcript**, not genomic, so no splicing correction is needed or applied.

## Panel -> data

| element | source |
|---|---|
| all four panels | `figures/S_metagene/metagene_pooled.tsv`, `metagene_by_tissue.tsv` |
| produced by | `scripts/build_metagene.py` (~7 min over 8 packs) |
| P-sites | `data/packed_canon*/target_counts.npy` (8 final-recipe training packs) |
| CDS coordinates | `data/tx2cds.tsv`, restricted to `has_start_codon == 1` |
| semantics + frame stats | `metagene_meta.json` |

Windows: **-50 to +150 nt** around the start, **-150 to +50 nt** around the stop.

**Brain is excluded.** It exists in the source study but was dropped from training for low
periodicity, so including it in a "training data" metagene would misrepresent what was trained on.
The 8 packs are Fibroblast (`packed_canon`), ES, Fat, HA_EC, HCAEC, Hepatocytes, HUVEC, VSMC.

## `n_tx` VARIES WITH OFFSET -- divide by the per-row value

A transcript contributes to an offset only if that offset exists within it, so a transcript with a
20 nt 5'UTR contributes nothing at offset -50. `n_tx` is therefore reported **per row**. Dividing the
pooled `psites` by one global transcript count would manufacture a decay at the window edges purely
from missing UTR, which would look like biology and is not.

Both panels A and B plot `psites / n_tx` for this reason.

## Two normalisations, and which to use

| column | meaning | use when |
|---|---|---|
| `psites` | raw pooled sum | you want the conventional metagene, dominated by highly translated genes |
| `psites_norm` | each transcript scaled to sum 1 inside its own window before pooling | you want shape independent of a few very deep genes |

Panels C/D instead normalise each **tissue** to sum 1, so 8 libraries spanning 28.8M to 330.8M
windowed P-sites are comparable in shape.

## Numbers

| quantity | value |
|---|--:|
| transcripts contributing | 71,116 (of 220,781 with an annotated CDS + start codon) |
| start-codon peak | 83.9 P-sites per transcript |
| mean 5'UTR (-50..-4) | 1.71 |
| mean CDS (0..+149) | 7.83 |
| **start peak / 5'UTR** | **49x** |
| CDS / 5'UTR | 4.6x |
| **CDS frame-0, first 150 nt** | **82.3%** (frame1 9.1%, frame2 8.6%) |

Frame-0 at 82.3% sits inside the 74.6-87.7% per-library range reported by `S_riboseq_qc`, computed
by a completely separate path, which is a useful consistency check on both.

## Reading it

- **The start-codon peak is 49x the 5'UTR baseline**, and the stop-codon peak is comparable. Both are
  the expected initiation and termination pauses.
- **All 8 tissues overlay almost exactly** in panels C/D despite an 11x depth spread. Shape is a
  property of the assay here, not of the library.
- **The 3'UTR falls to near zero** within ~10 nt of the stop, which is the cleanest single indicator
  that the ncRNA/cross-gene filtering and offset assignment are behaving.

## Caveats a caption must carry

- Transcripts **without** an annotated start codon are excluded entirely. A CDS inferred without one
  would blur offset 0, which is the position the figure exists to show.
- A CDS annotation that overruns the packed transcript is skipped rather than clipped, so the pack
  and the annotation can never be silently forced into the wrong frame.
- This is **observed** Ribo-seq only. No model output appears in this figure.
- Depth is NOT matched across tissues; use panels C/D for any cross-tissue statement.

## Regenerate

```bash
PY=/private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3
$PY scripts/build_metagene.py --out results/metagene
cd figures/S_metagene && $PY make_metagene.py
```
