# Environments

Captured 2026-08-04 from the live envs on `prism`. Nothing here was reconstructed by hand -- the
`.spec.txt` / `.pip.txt` files are dumps of what actually ran.

## The important part: training happens in CONTAINERS, not conda

The conda envs do CPU work (analysis, plotting, RiboCode, alignment, MS). **Model training runs in
Singularity**, and the two released models use **different images**, because only one of them has the
Mamba CUDA kernels:

| model | image | torch | mamba_ssm | causal_conv1d |
|---|---|---|---|---|
| `orf_v2_attn` (transformer) | `ghcr.io-ericmalekos-rnazoo-rinalmo-latest.img` | 2.2.0+cu118 | **absent** | **absent** |
| `orf_v2_mamba4` (primary) | `ghcr.io-ericmalekos-rnazoo-orthrus-latest.img` | 2.2.0+cu118 | 1.2.0.post1 | 1.2.0.post2 |

Both Python 3.10.20, CUDA 11.8. Images live in
`/private/groups/carpenterlab/emalekos/singularity_cache/` (5.6 GB and 6.8 GB).
`scripts/train_loto_union.sbatch` picks the image from `MIXER`; see its `case` block.

Trying to train or run mamba4 in the rinalmo image fails with `ModuleNotFoundError: mamba_ssm`.
That is the image being wrong, not the model.

Invocation pattern (the flags are load-bearing -- see below):

```
singularity exec --nv --no-home --cleanenv --env PYTHONNOUSERSITE=1 <image> python3 scripts/train_loto.py ...
```

`--no-home --cleanenv --env PYTHONNOUSERSITE=1` are REQUIRED. Without them the host's
`~/.local` site-packages leak into the container and shadow its pinned packages.

## Conda environments

| file | env | role | pkgs |
|---|---|---|---|
| `cas12a.spec.txt` | `cas12a` | main analysis / plotting / pgx pipeline (62 refs) | 156 + 70 pip |
| `ribocode.spec.txt` | `ribocode` | **RiboCode ORF calling** (28 refs) | 110 + 35 pip |
| `riboseq.spec.txt` | `riboseq` | alignment: STAR, salmon, cutadapt, samtools | 136 + 31 pip |
| `msfragger.spec.txt` | `msfragger` | MSFragger 4.2 + openjdk 25.0.2 | 82 + 41 pip |
| `ms2rescore.spec.txt` | `ms2rescore` | MS2Rescore / mokapot rescoring | 304 + 132 pip |
| `proteomics.spec.txt` | `proteomics` | Comet | 18 |

Recreate:

```
micromamba create -p <prefix> --file env/<name>.spec.txt
<prefix>/bin/pip install -r env/<name>.pip.txt     # where a .pip.txt exists
```

Key pins: python 3.11.15 / numpy 2.4.6 (cas12a), STAR 2.7.11b, salmon 2.1.2, cutadapt 5.2,
samtools 1.23.1, MSFragger 4.2.

## Two inconsistencies

**1. `cas12a` carries `torch==2.12.0+cpu`.** It is a CPU-only torch and is NOT the training
environment -- training is the containers above. Use `cas12a` for CPU inference and analysis only.

**2. Two RiboCode versions are installed:**

| env | RiboCode |
|---|---|
| `ribocode` | **1.2.15** |
| `riboseq` | 1.2.13 |

The pipeline calls RiboCode through `conda_envs/ribocode/bin/python` (28 references) and
`ribocode/bin/prepare_transcripts`, so **1.2.15 is the version behind the ORF calls**. `riboseq`'s
1.2.13 is incidental -- that env is used for STAR/salmon/cutadapt, not calling. Worth removing or
pinning to match, since two versions of the ORF caller in one project is a trap: a future script that
reaches for `riboseq/bin/python` would silently call ORFs with a different version.

## Not captured

- **Exact image digests.** The two `.img` files are `latest` tags pulled into the local cache, so the
  filename does not pin a version. For external release the images need immutable digests
  (`ghcr.io/...@sha256:...`) recorded here.
- `gedi` (PRICE) and `perturb` envs, referenced 2x and 1x respectively -- peripheral to the main
  pipeline.
