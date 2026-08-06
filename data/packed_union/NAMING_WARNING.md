# WARNING: `orf_track_v2.npy` here holds NO-KOZAK content despite the bare name

The convention (`prepare_pack.py`, since 2026-07-31) is:

| filename | content |
|---|---|
| `orf_track_v2.npy` (bare) | heuristic Kozak start-context factor |
| `orf_track_v2_nokozak.npy` | `--kozak none` |
| `orf_track_v2_pwm.npy` | learned PWM |

This pack was built 2026-07-25, BEFORE that naming fix, so its bare-named file actually contains
`--kozak none` content. `orf_track_v2_meta.json` is authoritative and says `"kozak": "none"`.

**ALWAYS read the meta.json rather than trusting the filename in this directory.** Reading the bare
name as "heuristic" is exactly the mistake that previously fed a nokozak-trained model the Kozak
track (see the ROOT FIX note in the project memory).

`orf_track_v2_nokozak.npy` is a symlink to the same file, added so the correct name resolves. Both
paths give identical bytes. Existing scripts that reference the bare name keep working.

The released union models (`orf_v2_{mamba4,attn}_onehot_union_noBrain_nokozak_mm1_holdout_Hepatocytes`)
were trained on THIS track, so `--kozak none` is the correct inference setting for them.

Recommended cleanup (not done, to avoid touching 2.1 GB referenced by live scripts): make
`orf_track_v2_nokozak.npy` the real file and the bare name the symlink, so the canonical name holds
the bytes.
