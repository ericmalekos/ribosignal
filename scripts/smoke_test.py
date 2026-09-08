#!/usr/bin/env python3
"""End-to-end check: synthetic pack -> ORF track -> both checkpoints -> density.

This is the test that would have caught the release bug the day it appeared. The README
quickstart could not run for anyone -- the loader resolved a relative --run against a path
that existed only on the author's cluster, then looked for `args.json` and `best.pt`, which
Hugging Face has never published -- and nothing in the repository executed the chain, so
nothing said so.

It runs on CPU with no real data and no GPU: `make_demo_pack.py` invents a handful of
transcripts, `build_orf_track.py` builds their track, and each released checkpoint predicts
a profile for them. mamba4 included -- scripts/mamba_ref.py means it no longer needs CUDA.

    python scripts/smoke_test.py --download        # fetch the weights, then run everything
    python scripts/smoke_test.py --weights weights/

What is asserted is PLUMBING, not accuracy. The sequences are random, so no number here
says anything about how well the model predicts real translation; held-out metrics live in
release/*/test_metrics.json. What the frame check does establish is that the model is
reading its inputs rather than emitting noise: the ORF-candidate track states the reading
frame, so a working model concentrates its density on it, and a wiring mistake -- a
transposed axis, a checkpoint loaded into the wrong architecture, a track that is off by
one -- lands the density back at chance.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
CHANCE = 1.0 / 3.0          # a profile ignoring frame puts a third of its mass in any frame
FRAME_FLOOR = 0.40          # observed 0.59 for both checkpoints on the default demo pack


def run(cmd, env=None):
    print(f"\n$ {' '.join(str(c) for c in cmd)}", flush=True)
    subprocess.run([str(c) for c in cmd], check=True, env=env)


def frame_share(pred, obs):
    """Share of predicted density in the observed reading frame.

    The frame is read off the observed counts' start peak rather than assumed, so this
    stays honest if make_demo_pack.py's layout ever changes.
    """
    f = int(np.argmax(obs)) % 3
    return float(pred[np.arange(len(pred)) % 3 == f].sum())


def check_profiles(npz_path, n_expected):
    z = np.load(npz_path, allow_pickle=False)
    lengths = z["lengths"]
    off = np.concatenate([[0], np.cumsum(lengths)])
    assert len(lengths) == n_expected, f"{npz_path}: {len(lengths)} tx, expected {n_expected}"
    assert z["pred_flat"].shape[0] == off[-1], "pred_flat does not match the summed lengths"
    assert z["obs_flat"].shape[0] == off[-1], "obs_flat does not match the summed lengths"

    shares, corrs = [], []
    for j in range(len(lengths)):
        p = z["pred_flat"][off[j]:off[j + 1]].astype(np.float64)
        o = z["obs_flat"][off[j]:off[j + 1]].astype(np.float64)
        assert np.isfinite(p).all(), f"tx {j}: non-finite profile"
        assert (p >= 0).all(), f"tx {j}: negative probability"
        assert abs(p.sum() - 1.0) < 1e-4, f"tx {j}: profile sums to {p.sum():.6f}, not 1"
        shares.append(frame_share(p, o))
        corrs.append(float(np.corrcoef(p, o)[0, 1]))
    tot = z["pred_total"]
    assert np.isfinite(tot).all() and (tot > 0).all(), f"pred_total not positive-finite: {tot}"

    mean_share = float(np.mean(shares))
    print(f"  {len(lengths)} tx | profiles sum to 1 | pred_total finite and > 0")
    print(f"  in-frame density share  mean {mean_share:.3f}  min {min(shares):.3f}  "
          f"(chance {CHANCE:.3f}, floor {FRAME_FLOOR})")
    print(f"  Pearson r vs observed   mean {np.mean(corrs):+.3f}  min {min(corrs):+.3f}")
    assert mean_share > FRAME_FLOOR, (
        f"mean in-frame share {mean_share:.3f} is at chance; the model is not using the "
        f"ORF track -- suspect the checkpoint/config pairing or the feature layout")
    return z


def check_density_flags():
    """build_density must refuse a flag it would silently discard (see validate_density_flags)."""
    sys.path.insert(0, str(SCRIPTS))
    from ribocode_dropin import build_density

    prof = np.array([0.1, 0.2, 0.7], np.float32)
    obs = np.array([1, 2, 7], np.int32)
    for kwargs in ({"variant": "pred_obsdepth", "pred_scale": 0.05},
                   {"variant": "real", "pred_scale": 0.05},
                   {"variant": "real", "poisson": True}):
        try:
            build_density(prof=prof, obs=obs, ptot=10.0, **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"build_density silently accepted {kwargs}")
    assert build_density("real", prof, obs, 10.0).tolist() == [1.0, 2.0, 7.0]
    assert build_density("pred_obsdepth", prof, obs, 10.0).tolist() == [1.0, 2.0, 7.0]
    print("  ignored-flag combinations rejected; real / pred_obsdepth densities unchanged")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--weights", default=None,
                    help="directory holding <arch>_best.pt and <arch>_config.json")
    ap.add_argument("--download", action="store_true",
                    help="fetch emalek/RiboSignal into <workdir>/weights first")
    ap.add_argument("--workdir", default=None, help="default: a temporary directory")
    ap.add_argument("--arch", nargs="+", default=["attn", "mamba4"],
                    help="checkpoints to run (default both)")
    ap.add_argument("--n-tx", type=int, default=6)
    args = ap.parse_args()

    tmp = None
    if args.workdir:
        work = Path(args.workdir)
        work.mkdir(parents=True, exist_ok=True)
    else:
        tmp = tempfile.TemporaryDirectory()
        work = Path(tmp.name)

    weights = Path(args.weights) if args.weights else work / "weights"
    if args.download:
        from huggingface_hub import snapshot_download
        print(f"downloading emalek/RiboSignal -> {weights}", flush=True)
        snapshot_download("emalek/RiboSignal", local_dir=str(weights))
    if not weights.is_dir():
        sys.exit(f"no weights at {weights}; pass --weights or --download")

    data = work / "data"
    pack = data / "packed_heldout_demo"
    fasta = data / "universe.fa"
    track = pack / "orf_track_v2_nokozak.npy"
    env = {**__import__("os").environ,
           "RIBO_DATA_DIR": str(data),
           "RIBO_ONEHOT_FASTA": str(fasta),
           "RIBO_ORF_TRACK": str(track)}

    run([sys.executable, SCRIPTS / "make_demo_pack.py", "--out", pack,
         "--fasta", fasta, "--n-tx", args.n_tx], env=env)
    run([sys.executable, SCRIPTS / "build_orf_track.py", "--mode", "ext", "--kozak", "none",
         "--pack", pack, "--fasta", fasta], env=env)

    for arch in args.arch:
        out = work / f"out_{arch}"
        run([sys.executable, SCRIPTS / "dump_pred_profiles.py", "--run", weights,
             "--arch", arch, "--heldout", "demo", "--device", "cpu", "--out", out], env=env)
        print(f"\n--- {arch} ---")
        check_profiles(out / "pred_profiles.npz", args.n_tx)

    print("\n--- ribocode_dropin.build_density ---")
    check_density_flags()

    print(f"\nSMOKE TEST PASSED ({', '.join(args.arch)})")
    if tmp:
        tmp.cleanup()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
