# Environments & references

Three conda environments and two Singularity images cover the whole pipeline; alignment references are built
once and reused.

## Compute environments

```{list-table}
:header-rows: 1
:widths: 26 20 54

* - Environment
  - Used for
  - Key tools
* - `conda_envs/riboseq`
  - alignment
  - `cutadapt`, `STAR`, `salmon`, `samtools`
* - `conda_envs/ribocode`
  - ORF calling
  - RiboCode `metaplots` + `detectORF` (the P-site + ORF-call reference, and the drop-in)
* - `conda_envs/cas12a`
  - analysis / plots
  - `numpy`, `h5py`, `matplotlib` (pack building, metrics, figures)
* - `singularity_cache/…rnazoo-rinalmo…`
  - train / predict
  - PyTorch; the model + one-hot and RiNALMo/Orthrus/HydraRNA embedding backends
* - `singularity_cache/…rnazoo-orthrus…`
  - Mamba variants
  - PyTorch **+ `mamba_ssm`** (required for the `--mixer mamba` runs)
```

```{note}
The model is ~5M parameters and **runs on CPU**. GPU is faster for training and the profile dump, but a jammed
GPU queue is not a blocker -- the dump and drop-in run on idle CPUs (`--device cpu`).
```

## Alignment references (built once, reused)

```{list-table}
:header-rows: 1
:widths: 30 20 50

* - Reference
  - Species
  - Notes
* - `STAR_indexes/star_index_grch38_v49`
  - human
  - GENCODE v49 primary assembly
* - `STAR_indexes/star_index_grcm39_vM38`
  - mouse
  - GENCODE vM38 (held-out mouse datasets)
* - `Salmon_indexes/salmon_index_decoy_v49` / `…_vM38`
  - human / mouse
  - **decoy-aware** (genome as decoys) -- required, transcriptome-only misattributes intronic reads
* - `ribocode_annot` / `ribocode_annot_mouse`
  - human / mouse
  - RiboCode annotation for `metaplots` + `detectORF`
```

## Ribo-seq alignment rules

Every Ribo-seq footprint alignment follows three non-negotiable rules (they change the P-site counts the model
learns from):

- **Adapter-trim first.** Many SRA/ENA submissions ship reads with the 3' TruSeq adapter attached; without
  `cutadapt -a AGATCGGAAGAGC --minimum-length 20 --maximum-length 40` you get ~0% unique mapping.
- **Unique mappers only.** `STAR --outFilterMultimapNmax 1`. RiboCode's caller does not drop multimappers, so
  they would inflate P-site counts and break periodicity.
- **Drop ncRNA + mitochondrial loci.** rRNA/tRNA/miRNA and the 13 chrM ORFs are removed before the caller.
