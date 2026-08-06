# Test sets

The model is judged on data it never trained on, at three widening levels of generalisation: a held-out
**tissue**, a held-out **study** (same species), and held-out **species**. Each external dataset is aligned
and RiboCode-called the same way as training, giving an *observed* ORF-call set to score predictions against.

```{list-table}
:header-rows: 1
:widths: 22 26 16 36

* - Dataset
  - Sample
  - Species
  - Generalisation tested
* - **Hepatocytes** (LOTO)
  - Chothani hepatocytes
  - human
  - held-out **tissue**, same atlas (in-distribution)
* - **Ruiz-Orera**
  - iPSC-cardiomyocytes
  - human
  - held-out **study** (different lab, protocol)
* - **Wang** / **Janich**
  - liver
  - mouse
  - held-out **species**
* - **GSE120762**
  - BMDM +/- LPS
  - mouse
  - held-out species + a **stimulus** contrast
* - **GSE155087**
  - CD4+ T cells (WT)
  - mouse
  - held-out species + a distinct **cell type**
```

## How a held-out is scored (the drop-in)

The deployed model predicts the P-site profile for the held-out transcripts; RiboCode then calls ORFs on
three densities (`scripts/ribocode_dropin.py`), and we compare the call sets:

```{list-table}
:header-rows: 1
:widths: 26 74

* - Density
  - What it isolates
* - `real`
  - RiboCode on the **observed** P-sites -- the truth/reference for this dataset
* - `pred_obsdepth`
  - predicted **shape** scaled to the real per-tx depth -- did the model get the ORF distribution right,
    independent of the count head?
* - `pred_preddepth`
  - fully **standalone**: predicted shape x predicted depth -- the true no-experiment prediction
```

## What generalisation looks like

The deployed human model transfers cleanly to held-out human tissue, a held-out human study, and -- most
tellingly -- to mouse cell types it never saw:

```{list-table}
:header-rows: 1
:widths: 30 22 24 24

* - Held-out
  - CDS F1 (standalone)
  - novel precision
  - non-canon F1
* - Ruiz-Orera (human iPSC-CM)
  - ~0.93
  - n/a
  - drop-in F1 0.931
* - Wang (mouse liver)
  - ~0.93
  - n/a
  - drop-in F1 0.929
* - GSE155087 (mouse T-cell)
  - **0.966**
  - 0.412 ($\theta$=1) / 0.592 (Poisson)
  - 0.495
```

```{admonition} Data provenance
:class: note
Every external dataset is logged before use (accession, genome build, access tier, raw-vs-BAM, matching,
md5) and processed by the same rules as training (adapter trim, unique mappers, ncRNA/chrM drop). The mouse
sets align to GRCm39/vM38; the human sets to GRCh38/v49.
```
