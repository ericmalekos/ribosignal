# Predict the profile

The weights are already in the image at `$RIBO_WEIGHTS`. Verify them, then predict with both
architectures.

```bash
(cd $RIBO_WEIGHTS && sha256sum -c SHA256SUMS)
```

`attn` runs on CPU. `mamba4` needs a GPU, because `mamba_ssm` has no CPU kernel.

```bash
python $RIBO_SCRIPTS/dump_pred_profiles.py --run $RIBO_WEIGHTS --arch attn \
       --pack pack/demo --out pred/attn --device cpu

python $RIBO_SCRIPTS/dump_pred_profiles.py --run $RIBO_WEIGHTS --arch mamba4 \
       --pack pack/demo --out pred/mamba4 --device cuda
```

On a GPU box set `RIBO_MAMBA_IMPL=cuda` so a missing `mamba_ssm` is an error rather than a silent
fall back to a pure-PyTorch path that is roughly 100x slower.
