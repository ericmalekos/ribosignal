# SUPPLEMENTAL S_riboseq_qc -- Ribo-seq library quality, every dataset used

**What it shows:** all 15 Ribo-seq library groups the study depends on, scored on the same three
axes, so a reviewer can check in one place that no dataset is bad enough to poison the labels or the
validation. Closes FIGURES_PLAN decision **D5**.

- **(a) read-length distribution** -- monosome footprints peak at 28-31 nt. An off-peak or flat
  distribution means the gel cut or the adapter trimming was wrong.
- **(b) 3-nt periodicity (f0 share at the annotated CDS)** -- the signal the model is trained to
  reproduce. 33% is the no-periodicity floor.
- **(c) P-site offset vs read length** -- should be near-constant within a library. A jumping offset
  means the per-nucleotide labels are smeared across codons.

---

## `make_riboseq_qc.py` -> `riboseq_qc.{png,pdf}` + `riboseq_qc_values.json`

### Data inputs

Every number comes from RiboCode `metaplots` output (`*_pre_config.txt`), which is produced as a
required step before any `RiboCode` call in this project. The commented block above the sample table
in each file holds, per read length: share of mapped reads, the fitted P-site offset, the frame
0/1/2 P-site sums at annotated start codons, and the periodicity p-values. No new alignment or
computation is performed by this figure.

**Training (Chothani, human, 9 tissues)** -- from the upstream expression-context project, which is
where the training P-sites were called:

    ../../results/chothani_regeneration/_canonical_<T>/canon_<T>_pre_config.txt
    T in {ES, Fat, Fibroblast, HA_EC, HCAEC, Hepatocytes, HUVEC, VSMC}

**Held-out validation (5)**:

    data/heldout_psites/human_ruizorera/human_ruizorera_pre_config.txt
    data/heldout_psites/mouse_wang_liver/mouse_wang_liver_pre_config.txt
    data/heldout_psites/mouse_gse120762_nt/mouse_gse120762_nt_pre_config.txt
    data/heldout_psites/mouse_gse120762_lps/mouse_gse120762_lps_pre_config.txt
    data/heldout_psites/mouse_gse155087_tcell/mouse_gse155087_tcell_pre_config.txt

**Proteogenomics cell lines (2)**:

    proteogenomics/data/HBL1_pilot/ribocode/HBL1_DMSO_pre_config.txt
    proteogenomics/data/B721_pilot/ribocode/B721_pre_config.txt

### Why RiboCode metaplots and not ribotish

FIGURES_PLAN D5 named ribotish. ribotish needs a **genome** BAM plus a GTF; this project keeps
transcriptome BAMs for everything except the training set, so honouring the letter of D5 would mean
re-aligning a dozen studies to regenerate a plot whose content already exists on disk. RiboCode
`metaplots` computes the same quantities (per-read-length offset, frame sums at start codons,
periodicity test) and has already been run on every dataset. Reading those files yields a strictly
larger panel at zero compute. The axes shown are the ones ribotish would have shown.

### Aggregation rules (so the numbers are reproducible)

- A `pre_config.txt` can describe **many libraries**; every one is parsed. `n_libraries` counts the
  distinct BAMs.
- `f0_frac` is computed from **summed** frame counts across all read lengths and libraries, not as a
  mean of per-length percentages, so a rare read length cannot swing it.
- `mode_read_length` and `mean_read_length` are weighted by each read length's share of mapped reads.
- Panel (a) normalises each row to sum to 1, so rows are comparable regardless of depth.
- Datasets whose `pre_config.txt` is absent are **printed as MISSING** and recorded in the JSON's
  `missing` key rather than silently dropped.

### Row order and colour

Rows are grouped by role (training -> held-out -> proteogenomics) and sorted by `f0_frac` within a
group. Panel (a) is an `imshow`, which places row 0 at the top; panels (b) and (c) set
`set_ylim(n-0.4, -0.6)` to match. **Do not replace this with `invert_yaxis()`** -- an earlier version
did, which flipped only (b) and (c) and silently mislabelled every row against the shared y-axis
labels in (a).

---

## What the panel actually shows (2026-08-07)

Every library clears the bar. f0 ranges 74.6% to 87.7% against a 33% no-periodicity floor, and the
fitted P-site offset is 12 nt in the majority of libraries. There is no dataset in this study whose
quality would justify excluding it.

Two observations worth stating in the caption rather than hiding:

- **HBL-1 is the weakest library in the project on two axes at once.** Its footprints peak at 34 nt,
  outside the 28-31 nt monosome window every other dataset sits in, and it contributes only 128,641
  P-sites at annotated CDS -- roughly 150x fewer than B721.221's 18.8 M. It also has the lowest f0
  (74.8%). This is consistent with its known ~2% unique mapping (~98% rRNA) and is a reason to weight
  the HBL-1 proteogenomics result below the other immunopeptidomes, not to drop it.
- **B721.221 is by far the deepest** (18.8 M P-sites at CDS across 7 libraries, f0 78.7%, a single
  consistent 12 nt offset). That depth is what makes it usable as the measured-translation reference
  in the drop-in ground-truth check.

### Rebuild

    /private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3 \
        figures/S_riboseq_qc/make_riboseq_qc.py

Deterministic: no sampling, no fitting, no randomness. Re-running on unchanged inputs reproduces the
figure and the JSON byte-for-byte.


## Repointed to the canonical configs (2026-08-14)

The training tissues previously read the source study's own
`experiments/biotype_probe/expression_context_human/data/ribocode_per_tissue/`, which no longer
exists: the final-recipe alignment redo replaced it. They now read
`results/chothani_regeneration/_canonical_<T>/canon_<T>_pre_config.txt`, produced by
`scripts/chothani_canonical_diverge.sbatch`.

**The failure mode was silent and worth recording.** The nine training entries dropped into the
generator's `missing` list and the panel rendered anyway -- 7 held-out/pgx arms and NOT ONE training
tissue -- while its own README row still claimed "16 dataset arms (9 training + ...)". A QC panel
whose entire job is to establish training-library quality had stopped showing any training library,
and nothing failed. The `missing` list is printed, but nothing gates on it.

**Brain is correctly absent**, so 8 training arms is the right count, not 9: it was dropped from the
source study's 9-tissue panel for low periodicity (period_obs 0.044) before any pack was built, so
there is no canonical config to read and nothing to draw.

Per-tissue library counts now reconcile exactly with P3 (Fibroblast 32, VSMC 11, ES/Fat/HA_EC 6,
HCAEC/Hepatocytes 5, HUVEC 3 = 69 training + 5 held out).

## Expanded to 18 arms 2026-08-15 (was 15)

Added at the poster session's request, filling every frame-0 hole the poster's dataset census had:

| arm | role | pooled frame-0 | libs |
|---|---|--:|--:|
| THP-1 GSE208041 | heldout | 84.0% | 5 |
| CAR-T GSE304796 | heldout | 78.5% | 3 |
| Brain | **dropped** (new role) | 76.2% | 5 |

**`dropped` is a new role**, grey, holding the project's two EXCLUDED Ribo-seq datasets. Brain is
`dropped` and not `train`: it appears in no training pack, so colouring it as training would
misstate what the model saw.

**THP-1 GSE39561 cannot be added.** Its `pre_config.txt` is 957 bytes of `#/path` and bare
`#read_length` headers with ZERO data rows -- RiboCode wrote no per-read-length table because it
selected no periodic read lengths. It lands in `missing` and the script says so. Its frame-0 comes
from prose (40.7%, results.md:2813), which is the SAME quantity as `f0_frac`: pooled frame
concentration at annotated CDS starts.

**Brain is not a marginal library on these axes**, which matters for how its exclusion is captioned:
frame-0 **76.2%** is ABOVE HUVEC's 74.6%, and its pack holds **194M** P-sites against HUVEC's 94M --
and HUVEC is trained on. Brain's exclusion rests entirely on `period_obs` 0.044, a transcriptome-wide
statistic this panel does not plot. It was excluded on signal strength, not because the library is
broken.

The old comment claiming Brain was "dropped before any pack was built, so 8 is the correct count,
not 9" was stale twice over: `data/packed_Brain` exists (2026-07-11), and Brain was realigned and
metaplotted on the final recipe 2026-08-15 (5/5 libraries selected offsets). **8 training tissues and
9 Chothani tissues are both correct**, for different counts.

## `psites_at_cds` renamed to `psites_at_start_window` (2026-08-15)

The value never changed; the NAME was wrong and cost a downstream figure. It is the sum of the
frame-0/1/2 counts in RiboCode's metaplots table -- P-sites in the window around ANNOTATED START
CODONS, not across the CDS and not the library. It runs ~0.5% of real pack depth (Fibroblast 5.03M
here vs 1,069.8M in the pack).

**It is not usable as a depth proxy either.** Its ratio to real pack depth swings 0.33-0.65% across
the eight Chothani tissues, enough to invert the ES/Fat ordering: this panel says ES 2.10M > Fat
1.86M, while real depth says ES 494M < Fat 565M. For depth, read the pack.
