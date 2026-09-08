# Example SLURM scripts

Adapt these, do not run them as-is. They take every path from an environment variable and have no
site-specific defaults, but the `--partition` and `--account` lines are placeholders and the
container path assumes Singularity.

```bash
export RIBO_IMAGE=/path/to/riboseq-model.sif   # singularity build from ghcr.io/ericmalekos/riboseq-model
export RIBO_WORK=/path/to/workdir              # the demo working directory
export RIBO_REF=$RIBO_WORK/ref                 # annotation, genome, STAR index
sbatch --partition=<yours> --account=<yours> examples/slurm/01_align.sbatch
```

The tutorial covers chr22 on one machine. These exist for a whole-genome run, where the STAR index
alone needs about 32 GB.
