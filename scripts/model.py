#!/usr/bin/env python3
"""Dual-head dilated-CNN for per-nucleotide Ribo-seq P-site prediction.

Design (see methods.md section 5B): the per-nt FM embedding already carries long-range
sequence context, so the network body only needs an O(L) local operator to mix in the
RNAseq coverage and sharpen to a base-resolution, frame-aware profile. A stack of residual
dilated 1-D convolutions (BPNet-style, Avsec et al. 2021) does this. Two heads:

  - profile head: per-nt logits -> (at loss time) log-softmax over the transcript ->
    multinomial NLL against per-nt counts. Scale-free, captures 3-nt periodicity and
    CDS localization, and is the transferable target across datasets.
  - count head: a scalar log-total from masked-mean body features + log1p(total coverage),
    regressed to log1p(total P-sites). Magnitude; dataset-depth-dependent.

Input per nt is [FM emb (D_emb) ; log1p(coverage) (1)] -> D_emb+1 channels.
Convolutions operate on (B, C, L); the dataset provides (B, L, C) so we transpose.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class ResidualDilatedBlock(nn.Module):
    """Conv -> norm -> GELU -> Conv -> norm, added back to the input (same length)."""

    def __init__(self, channels: int, kernel_size: int, dilation: int, dropout: float):
        super().__init__()
        pad = dilation * (kernel_size - 1) // 2
        self.conv1 = nn.Conv1d(channels, channels, kernel_size,
                               padding=pad, dilation=dilation)
        self.conv2 = nn.Conv1d(channels, channels, kernel_size,
                               padding=pad, dilation=dilation)
        # GroupNorm is batch-size independent (we use token-budget batches of varying B).
        self.norm1 = nn.GroupNorm(8, channels)
        self.norm2 = nn.GroupNorm(8, channels)
        self.drop = nn.Dropout(dropout)

    def forward(self, x):
        h = self.drop(F.gelu(self.norm1(self.conv1(x))))
        h = self.norm2(self.conv2(h))
        return F.gelu(x + h)


class BiMambaBlock(nn.Module):
    """Pre-norm residual block with a BIDIRECTIONAL Mamba mixer (forward + reversed scan).

    seq2ribo's polisher (Kaynar & Kingsford 2026) refines an sTASEP-simulated profile with a
    single CAUSAL Mamba, which suffices there because the simulator already injects downstream
    (exclusion-zone) context. This model has no simulator, so a causal scan alone would lose the
    downstream context that the bidirectional transformer provides; we run Mamba both directions
    and sum, matching the transformer's full-transcript receptive field at O(L) instead of
    O(L^2). mamba_ssm is imported lazily so the default transformer path needs no CUDA Mamba build.
    """

    def __init__(self, d_model: int, d_state: int = 16, d_conv: int = 4,
                 expand: int = 2, dropout: float = 0.1):
        super().__init__()
        from mamba_ssm import Mamba  # lazy: only imported when a mamba mixer is actually built
        self.norm = nn.LayerNorm(d_model)
        self.fwd = Mamba(d_model=d_model, d_state=d_state, d_conv=d_conv, expand=expand)
        self.bwd = Mamba(d_model=d_model, d_state=d_state, d_conv=d_conv, expand=expand)
        self.drop = nn.Dropout(dropout)

    def forward(self, x, mask):
        # x (B, L, d); mask (B, L) bool. Zero padded positions before the scan so pads do not
        # seed the SSM state (they would otherwise start the reversed scan). Both scans read the
        # same normed input; the assignment to h happens only after the RHS is fully evaluated.
        residual = x
        h = self.norm(x) * mask.unsqueeze(-1).to(x.dtype)
        h = self.fwd(h) + self.bwd(h.flip(1)).flip(1)
        return residual + self.drop(h)


class MambaBody(nn.Module):
    """Stack of bidirectional Mamba blocks; a drop-in alternative to the TransformerEncoder."""

    def __init__(self, d_model: int, n_layers: int, d_state: int = 16, d_conv: int = 4,
                 expand: int = 2, dropout: float = 0.1):
        super().__init__()
        self.blocks = nn.ModuleList([
            BiMambaBlock(d_model, d_state, d_conv, expand, dropout) for _ in range(n_layers)
        ])

    def forward(self, x, mask):
        for blk in self.blocks:
            x = blk(x, mask)
        return x


class RiboSignalModel(nn.Module):
    def __init__(self, d_emb: int = 1280, channels: int = 256, n_blocks: int = 10,
                 kernel_size: int = 3, dropout: float = 0.1,
                 dilations: tuple[int, ...] | None = None,
                 extra_in: int = 0, n_attn_layers: int = 0, n_heads: int = 8,
                 mixer: str = "transformer", mamba_d_state: int = 16,
                 mamba_d_conv: int = 4, mamba_expand: int = 2,
                 learn_start_context: bool = False, fm_to_mixer: bool = False):
        super().__init__()
        self.d_emb = d_emb
        self.extra_in = extra_in            # e.g. sequence-derived ORF-candidate channels
        self.n_attn_layers = n_attn_layers
        self.mixer = mixer                  # global-context mixer: "transformer" or "mamba"
        # Q3 (Kozak ablation): replace the hand-picked start-context heuristic with a LEARNED gate.
        # A conv over the one-hot window [-6,+4] around each nt yields a per-nt gate in (0,1); the
        # (no-Kozak) codon-weight start channel is remodulated w*(0.5+0.5*gate) before in_proj, so the
        # 4x10 kernel is a learned nucleotide x position Kozak matrix (extracted post-training). The
        # start channel is ORF-track channel 3 -> feats index d_emb+3. One-hot backend only.
        self.learn_start_context = learn_start_context
        self.start_ch_idx = d_emb + 3
        if learn_start_context:
            assert d_emb == 4 and extra_in >= 4, \
                "learn_start_context requires the onehot backend (d_emb=4) + the ORF track"
            self.start_ctx_conv = nn.Conv1d(4, 1, kernel_size=10, bias=True)
            self.start_ctx_pad = (6, 3)     # transcript window [i-6, i+3] = Kozak -6..+4
        if dilations is None:
            dilations = tuple(2 ** (i % 10) for i in range(n_blocks))  # 1,2,...,512,1,...
        # input per nt is [FM emb (d_emb) ; ORF track (extra_in) ; log1p(coverage) (1)];
        # coverage stays the LAST channel so the count head can read feats[..., -1].
        self.in_proj = nn.Conv1d(d_emb + extra_in + 1, channels, 1)
        self.blocks = nn.ModuleList([
            ResidualDilatedBlock(channels, kernel_size, d, dropout) for d in dilations
        ])
        # Optional global-context stack: the dilated body has a ~4 kb receptive field, so a
        # position deep in a long CDS cannot see the start codon that sets its frame. Self-
        # attention gives full-transcript context (torch 2.x uses a memory-efficient kernel).
        if n_attn_layers > 0:
            if mixer == "mamba":
                # bidirectional Mamba stack (O(L)); same layer count as the transformer for a
                # controlled swap. Built only here, so mamba_ssm is required just for this path.
                self.mamba_body = MambaBody(channels, n_attn_layers, mamba_d_state,
                                            mamba_d_conv, mamba_expand, dropout)
            else:
                layer = nn.TransformerEncoderLayer(
                    d_model=channels, nhead=n_heads, dim_feedforward=channels * 2,
                    dropout=dropout, activation="gelu", batch_first=True, norm_first=True)
                self.attn = nn.TransformerEncoder(layer, num_layers=n_attn_layers)
        # Direct FM -> mixer path. By default the FM embedding reaches the mixer only after being
        # (a) linearly squeezed to `channels` by the SHARED in_proj, jointly with the ORF track and
        # coverage, and (b) rewritten by all n_blocks dilated convs. For RiNALMo that is 1280 -> 256
        # in a single 1x1 conv before any spatial processing -- a 5x bottleneck the 4-d one-hot
        # backend never pays, which confounds "does the FM help" with "does the FM survive fusion".
        #
        # This adds a dedicated projection of the embedding straight onto the mixer input:
        #     x_mixer = cnn_out + gate * GroupNorm(fm_proj(emb))
        # `gate` is a scalar initialised to ZERO (LayerScale / ReZero style), so the network starts
        # numerically identical to the baseline and must actively open the path. The learned |gate|
        # is therefore a direct readout of how much the model wants the un-convolved embedding --
        # a gate that stays ~0 is itself an interpretable negative result.
        self.fm_to_mixer = fm_to_mixer
        if fm_to_mixer:
            assert n_attn_layers > 0, "fm_to_mixer requires a mixer stack (n_attn_layers > 0)"
            self.fm_proj = nn.Conv1d(d_emb, channels, 1)
            self.fm_norm = nn.GroupNorm(8, channels)
            self.fm_gate = nn.Parameter(torch.zeros(1))
        self.profile_head = nn.Conv1d(channels, 1, 1)
        # count head: masked-mean body features (channels) + log1p(total coverage) (1)
        self.count_head = nn.Sequential(
            nn.Linear(channels + 1, channels), nn.GELU(), nn.Dropout(dropout),
            nn.Linear(channels, 1),
        )

    def forward(self, feats, mask):
        """feats (B, L, d_emb+extra_in+1) float; mask (B, L) bool (True = real nt).

        Returns profile_logits (B, L) with -1e9 at padded positions, and pred_logcount
        (B,) the predicted log1p(total P-sites).
        """
        if self.learn_start_context:
            # learned Kozak: gate the codon-weight start channel by a conv over the one-hot context
            oh = feats[:, :, :4].transpose(1, 2)                  # (B, 4, L) one-hot (onehot backend)
            gate = torch.sigmoid(
                self.start_ctx_conv(F.pad(oh, self.start_ctx_pad))).squeeze(1)   # (B, L) in (0,1)
            idx = self.start_ch_idx
            w = feats[:, :, idx]                                  # (B, L) codon-weight start (no Kozak)
            new_start = (w * (0.5 + 0.5 * gate)).unsqueeze(-1)    # (B, L, 1)
            feats = torch.cat([feats[:, :, :idx], new_start, feats[:, :, idx + 1:]], dim=2)
        x = feats.transpose(1, 2)                 # (B, C_in, L)
        x = self.in_proj(x)
        m = mask.unsqueeze(1).to(x.dtype)         # (B, 1, L)
        for blk in self.blocks:
            x = blk(x * m)                        # zero padded positions before each conv
        x = x * m
        if self.n_attn_layers > 0:
            if self.fm_to_mixer:
                # un-convolved embedding, own projection, zero-init gate (see __init__)
                fm = feats[:, :, :self.d_emb].transpose(1, 2)      # (B, d_emb, L)
                x = x + self.fm_gate * self.fm_norm(self.fm_proj(fm)) * m
            xt = x.transpose(1, 2)                # (B, L, C)
            if self.mixer == "mamba":
                xt = self.mamba_body(xt, mask)
            else:
                xt = self.attn(xt, src_key_padding_mask=~mask)
            x = xt.transpose(1, 2) * m           # (B, C, L), re-zero padded positions
        logits = self.profile_head(x).squeeze(1)  # (B, L)
        logits = logits.masked_fill(~mask, -1e9)

        # count head: masked mean over valid positions
        denom = mask.sum(dim=1, keepdim=True).clamp(min=1).to(x.dtype)   # (B,1)
        feat_mean = (x * m).sum(dim=2) / denom                          # (B, C)
        logcov = torch.log1p(feats[:, :, -1].clamp(min=0).sum(dim=1, keepdim=True))
        pred_logcount = self.count_head(torch.cat([feat_mean, logcov], dim=1)).squeeze(1)
        return logits, pred_logcount


def profile_multinomial_nll(logits, counts, mask):
    """Per-transcript multinomial NLL, equal-weighted across transcripts.

    logits, counts, mask: (B, L). Returns scalar mean over the batch of
    -sum_i (c_i / N) * log softmax(logits)_i  (cross-entropy of the empirical profile),
    which is scale-free (independent of a transcript's total depth N).
    """
    log_p = F.log_softmax(logits, dim=1)               # -inf at padded (logit -1e9)
    log_p = torch.where(mask, log_p, torch.zeros_like(log_p))
    counts = counts * mask                             # zero padded counts
    n = counts.sum(dim=1).clamp(min=1.0)               # (B,)
    per_tx = -(counts * log_p).sum(dim=1) / n          # (B,)
    return per_tx.mean()


def count_mse(pred_logcount, counts, mask):
    """MSE between predicted log-count and log1p(observed total P-sites)."""
    target = torch.log1p((counts * mask).sum(dim=1))
    return F.mse_loss(pred_logcount, target)


def profile_entropy_gap(logits, counts, mask):
    """Anti-smoothing penalty (O3): how much SMOOTHER the predicted profile is than the
    observed one, per transcript, averaged over the batch.

    The multinomial NLL is minimized by matching the empirical profile, but a smooth
    prediction still scores well because cross-entropy tolerates a hedged distribution.
    That smoothness is why RiboCode over-calls short non-canonical ORFs on the predicted
    density (the shape features that separate real from spurious ORFs are washed out).

    This term penalizes the Shannon-entropy excess H(p) - H(q), where p = softmax(logits)
    is the predicted profile and q = counts / N the empirical one, both over the masked
    transcript positions. relu() keeps only over-smoothing (H(p) > H(q)); a prediction
    peakier than the target is never pushed back (avoids chasing shot noise). Because the
    pooled Chothani targets are deep (~6.6 P-sites/nt), H(q) reflects real 3-nt periodicity
    and CDS localization, not sparse-sampling spikes. Returns a scalar >= 0 (nats).
    """
    log_p = F.log_softmax(logits, dim=1)                  # padded logits (-1e9) -> ~ -inf
    zero = torch.zeros_like(log_p)
    p = torch.where(mask, log_p.exp(), zero)              # zero the padded positions
    ent_p = -(p * torch.where(mask, log_p, zero)).sum(dim=1)          # (B,) H(p)

    counts = counts * mask
    n = counts.sum(dim=1).clamp(min=1.0)                  # (B,)
    q = counts / n.unsqueeze(1)
    log_q = torch.where(q > 0, torch.log(q.clamp(min=1e-12)), zero)
    ent_q = -(q * log_q).sum(dim=1)                       # (B,) H(q)

    return F.relu(ent_p - ent_q).mean()


if __name__ == "__main__":
    # smoke: shapes + finite loss on random data
    torch.manual_seed(0)
    B, L, D = 3, 512, 1280
    feats = torch.randn(B, L, D + 1)
    feats[:, :, -1] = feats[:, :, -1].abs()            # coverage >= 0
    mask = torch.ones(B, L, dtype=torch.bool)
    mask[0, 300:] = False                              # a shorter transcript
    counts = (torch.rand(B, L) * 5).floor() * mask
    m = RiboSignalModel(d_emb=D, channels=64, n_blocks=6)
    logits, pc = m(feats, mask)
    loss = profile_multinomial_nll(logits, counts, mask) + 0.1 * count_mse(pc, counts, mask)
    print("logits", tuple(logits.shape), "pred_logcount", tuple(pc.shape),
          "loss", float(loss))
    loss.backward()
    print("backward OK; params:", sum(p.numel() for p in m.parameters()))
