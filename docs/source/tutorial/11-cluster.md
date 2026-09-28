# Running it on a cluster

The tutorial runs on chr22 on one machine. For a whole-genome run on SLURM, adapt the example
scripts in `examples/slurm/`.

```bash
ls examples/slurm/
```

They are examples to adapt, not a supported entry point. Each takes its paths from environment
variables with no cluster-specific defaults, so set them before submitting.

```bash
export RIBO_REF=/path/to/reference
export RIBO_WORK=/path/to/workdir
export RIBO_IMAGE=/path/to/riboseq-model.sif
sbatch --partition=<yours> --account=<yours> examples/slurm/02_predict.sbatch
```

The scripts carry no `--partition` or `--account` line, since both are site-specific; pass them on
the `sbatch` command line as above. They run the image through Singularity, not Docker.
