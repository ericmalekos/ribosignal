# Predict the profile

The weights are already in the image at `$RIBO_WEIGHTS`. Verify them, then predict with both
architectures.

```bash
(cd $RIBO_WEIGHTS && sha256sum -c SHA256SUMS)
```

Both architectures run on CPU and on GPU. `mamba_ssm` has no CPU kernel, so on CPU `mamba4` falls
back to `scripts/mamba_ref.py`, a transcription of upstream's own pure-PyTorch reference; it
reproduces the CUDA kernel to 1.19e-07 max absolute difference, so the CPU and GPU outputs are the
same profile.

```bash
python $RIBO_SCRIPTS/dump_pred_profiles.py --run $RIBO_WEIGHTS --arch attn \
       --pack pack/demo --out pred/attn --device cpu

python $RIBO_SCRIPTS/dump_pred_profiles.py --run $RIBO_WEIGHTS --arch mamba4 \
       --pack pack/demo --out pred/mamba4 --device cpu     # or --device cuda
```

On a GPU box set `RIBO_MAMBA_IMPL=cuda` so a missing `mamba_ssm` is an error rather than a silent
fall back to the reference path, which is far slower.

## How long it takes

Measured on this chr22 pack, 6,583 scored transcripts. CPU is 16 threads on one compute node; GPU
is a single A100.

| architecture | CPU, 16 threads | GPU (A100) | GPU speedup |
|---|--:|--:|--:|
| `attn` | 30 min (1,818 s) | 59.6 s | 30x |
| `mamba4` | 3 h 10 min (11,431 s) | 35.4 s | 323x |

**Which architecture is faster depends on the device, and the ordering flips.** On a GPU `mamba4`
is the quicker of the two. On CPU it is 6.3x slower than `attn`, because the reference path
replaces one fused kernel with an explicit scan. On CPU prefer `attn`; on GPU either is cheap.

Pin the thread count before starting, or torch takes every core it can see and the wall time is
not reproducible. The demo script exports these itself.

```bash
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16 NUMEXPR_NUM_THREADS=16
```

**More threads is not faster here.** Left unpinned on a 160-core machine, torch took 67 cores and
`mamba4` ran 5 h 36 min, against 3 h 10 min on 16. `attn` was unchanged. The reference scan is
sequential along the transcript, so the extra threads contend rather than help. If you have many
cores, prefer a GPU or run `attn`; do not raise this number expecting a speedup.

`pred_total` differs slightly between CPU and GPU: about 0.2% in counts, median relative
difference 7e-04, correlation 0.9999998. The count head emits `log1p(total)`, so `expm1` amplifies
ordinary float32 accumulation-order differences. It is far below the Poisson dial, where theta
0.05 is a 20x change, and does not move the ORF calls.
