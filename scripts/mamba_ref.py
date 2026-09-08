#!/usr/bin/env python3
"""A dependency-free Mamba mixer, so the mamba4 checkpoint is not GPU-only.

`mamba_ssm` cannot be installed without CUDA: `causal-conv1d`'s setup reads
`torch.version.cuda` to pick a prebuilt wheel and dies on `None`, and `mamba-ssm`
hard-requires `triton`, which publishes no macOS wheel at all. So on any machine without
a CUDA toolchain the primary released checkpoint could not be loaded, let alone run --
`scripts/dump_pred_profiles.py`'s `RIBO_MAMBA_CPU=1` path does not help, because it
patches kernels *inside* an already-imported `mamba_ssm`.

This module is a transcription of upstream's own pure-PyTorch reference path, so a
checkpoint trained against the CUDA kernels loads and runs here unchanged:

  mamba_ssm/modules/mamba_simple.py   Mamba.__init__ and Mamba.forward   (v1.2.0.post1)
  mamba_ssm/ops/selective_scan_interface.py   selective_scan_ref, mamba_inner_ref
  causal_conv1d/causal_conv1d_interface.py    causal_conv1d_ref          (v1.2.0.post2)

Copyright (c) 2023-2024 Tri Dao, Albert Gu; Apache-2.0, as published with those files.

Parameter names and shapes are byte-for-byte what `mamba_ssm.Mamba` registers -- in_proj,
conv1d, x_proj, dt_proj, A_log, D, out_proj -- so `load_state_dict(strict=True)` is itself
the check that this is the same module. `selftest()` at the bottom asserts that, and
compares against the real `mamba_ssm` whenever one is importable.

Deliberately NOT supported: inference_params / step() incremental decoding, complex A, and
the fused fast path. This model runs one full transcript at a time, so none apply.

Cost: the selective scan is a sequential recurrence over transcript positions, and here it
is a Python loop rather than a fused kernel. Expect a few seconds per kilobase per block.
Use the attn checkpoint when that matters; it is the CPU release for this reason.
"""
from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F

# Positions per iteration of the scan. Only bounds peak memory: the recurrence is carried
# across chunk boundaries, so the result does not depend on this value. The intermediates
# are (batch, d_inner, chunk, d_state) float32, hence ~134 MB at batch 8 / d_inner 512.
SCAN_CHUNK = 512


def selective_scan(u, delta, A, B, C, D=None, z=None, delta_bias=None,
                   delta_softplus=False, chunk=SCAN_CHUNK):
    """Upstream `selective_scan_ref`, chunked over length to bound peak memory.

      u, delta : (batch, d_inner, L)
      A        : (d_inner, d_state), real, negative
      B, C     : (batch, d_state, L)   -- input-dependent
      D        : (d_inner,)
      z        : (batch, d_inner, L)   -- the SiLU gate

    Everything is computed in float32 and cast back at the end, exactly as upstream does.
    """
    dtype_in = u.dtype
    u = u.float()
    delta = delta.float()
    if delta_bias is not None:
        delta = delta + delta_bias[..., None].float()
    if delta_softplus:
        delta = F.softplus(delta)
    B = B.float()
    C = C.float()

    batch, dim, L = u.shape
    state = A.new_zeros((batch, dim, A.shape[1]))            # (batch, d_inner, d_state)
    ys = []
    for lo in range(0, L, chunk):
        hi = min(lo + chunk, L)
        d = delta[:, :, lo:hi]                               # (b, d, l)
        # deltaA[b,d,l,n] = exp(delta[b,d,l] * A[d,n]); A broadcast as (1, d_inner, 1, d_state)
        deltaA = torch.exp(d.unsqueeze(-1) * A.unsqueeze(0).unsqueeze(2))
        # deltaB_u[b,d,l,n] = delta[b,d,l] * B[b,n,l] * u[b,d,l]
        deltaB_u = (d * u[:, :, lo:hi]).unsqueeze(-1) * B[:, :, lo:hi].permute(0, 2, 1).unsqueeze(1)
        Cc = C[:, :, lo:hi]                                  # (b, n, l)
        for i in range(hi - lo):
            state = deltaA[:, :, i] * state + deltaB_u[:, :, i]
            ys.append((state * Cc[:, :, i].unsqueeze(1)).sum(-1))     # (b, d)
    y = torch.stack(ys, dim=2)                               # (batch, d_inner, L)
    out = y if D is None else y + u * D.unsqueeze(-1)
    if z is not None:
        out = out * F.silu(z.float())
    return out.to(dtype=dtype_in)


class Mamba(nn.Module):
    """`mamba_ssm.Mamba` with the CUDA kernels replaced by the reference path.

    The constructor signature accepts, and ignores, the fused-path arguments so a caller
    can swap this in for the real class without a conditional.
    """

    def __init__(self, d_model, d_state=16, d_conv=4, expand=2, dt_rank="auto",
                 conv_bias=True, bias=False, use_fast_path=False, layer_idx=None,
                 device=None, dtype=None, **_ignored):
        factory_kwargs = {"device": device, "dtype": dtype}
        super().__init__()
        self.d_model = d_model
        self.d_state = d_state
        self.d_conv = d_conv
        self.expand = expand
        self.d_inner = int(expand * d_model)
        self.dt_rank = math.ceil(d_model / 16) if dt_rank == "auto" else dt_rank
        self.layer_idx = layer_idx

        self.in_proj = nn.Linear(self.d_model, self.d_inner * 2, bias=bias, **factory_kwargs)
        self.conv1d = nn.Conv1d(self.d_inner, self.d_inner, kernel_size=d_conv, bias=conv_bias,
                                groups=self.d_inner, padding=d_conv - 1, **factory_kwargs)
        self.activation = "silu"
        self.act = nn.SiLU()
        self.x_proj = nn.Linear(self.d_inner, self.dt_rank + self.d_state * 2, bias=False,
                                **factory_kwargs)
        self.dt_proj = nn.Linear(self.dt_rank, self.d_inner, bias=True, **factory_kwargs)
        # S4D-real init. Only reached when this module is trained from scratch; loading a
        # released checkpoint overwrites all of it.
        A = torch.arange(1, self.d_state + 1, dtype=torch.float32,
                         device=device).repeat(self.d_inner, 1).contiguous()
        self.A_log = nn.Parameter(torch.log(A))
        self.A_log._no_weight_decay = True
        self.D = nn.Parameter(torch.ones(self.d_inner, device=device))
        self.D._no_weight_decay = True
        self.out_proj = nn.Linear(self.d_inner, self.d_model, bias=bias, **factory_kwargs)

    def forward(self, hidden_states, inference_params=None):
        """(B, L, d_model) -> (B, L, d_model)."""
        if inference_params is not None:
            raise NotImplementedError(
                "mamba_ref.Mamba has no incremental-decoding cache; it scores whole "
                "transcripts in one pass. Install mamba_ssm for step() decoding.")
        batch, seqlen, _ = hidden_states.shape

        # (B, L, d) -> (B, 2*d_inner, L). Upstream fuses the matmul with the transpose to
        # save one copy; nn.Linear then transpose is the same arithmetic, bias included.
        xz = self.in_proj(hidden_states).transpose(1, 2)
        x, z = xz.chunk(2, dim=1)

        # Causal depthwise conv + SiLU. conv1d pads by d_conv-1 on both sides, so trimming
        # to seqlen leaves exactly the causal window -- identical to causal_conv1d_ref.
        x = self.act(self.conv1d(x)[..., :seqlen])

        x_dbl = self.x_proj(x.transpose(1, 2).reshape(batch * seqlen, self.d_inner))
        dt, B, C = torch.split(x_dbl, [self.dt_rank, self.d_state, self.d_state], dim=-1)
        dt = (self.dt_proj.weight @ dt.t()).reshape(self.d_inner, batch, seqlen).transpose(0, 1)
        B = B.reshape(batch, seqlen, self.d_state).transpose(1, 2).contiguous()
        C = C.reshape(batch, seqlen, self.d_state).transpose(1, 2).contiguous()

        A = -torch.exp(self.A_log.float())
        y = selective_scan(x, dt, A, B, C, self.D.float(), z=z,
                           delta_bias=self.dt_proj.bias.float(), delta_softplus=True)
        return self.out_proj(y.transpose(1, 2))


def selftest(seed=0, tol=2e-4):
    """Shapes, finiteness, chunk-invariance, and -- where mamba_ssm exists -- equality.

    The chunk-invariance check is the one that catches a broken scan without a GPU in the
    room: the recurrence must give the same answer whatever the chunking, which it does
    only if the carried state is right.
    """
    torch.manual_seed(seed)
    m = Mamba(d_model=64, d_state=16, d_conv=4, expand=2).eval()
    x = torch.randn(2, 97, 64)
    with torch.no_grad():
        y = m(x)
    assert y.shape == x.shape, y.shape
    assert torch.isfinite(y).all(), "non-finite output"

    global SCAN_CHUNK
    keep, SCAN_CHUNK = SCAN_CHUNK, 7
    try:
        with torch.no_grad():
            y_small = m(x)
    finally:
        SCAN_CHUNK = keep
    d = (y - y_small).abs().max().item()
    assert d < 1e-5, f"chunking changed the result by {d:.3g}"
    print(f"mamba_ref: shapes OK, finite, chunk-invariant (max diff {d:.2e})")

    try:
        from mamba_ssm import Mamba as CudaMamba
    except Exception as e:                                    # noqa: BLE001
        print(f"mamba_ssm not importable ({type(e).__name__}); skipped the equality check")
        return
    ref = CudaMamba(d_model=64, d_state=16, d_conv=4, expand=2, use_fast_path=False).eval()
    ref.load_state_dict(m.state_dict(), strict=True)          # names/shapes must match exactly
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    with torch.no_grad():
        y_ref = ref.to(dev)(x.to(dev)).cpu()
    d = (y - y_ref).abs().max().item()
    assert d < tol, f"mamba_ssm disagrees by {d:.3g}"
    print(f"mamba_ref: matches mamba_ssm on {dev} (max abs diff {d:.2e})")


if __name__ == "__main__":
    selftest()
