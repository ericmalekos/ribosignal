#!/usr/bin/env python3
"""Flag scripts that are not registered in docs/CANONICAL_ENTRYPOINTS.md.

A rule that lives only in prose does not get followed -- docs/PIPELINE_POLICY.md already said "one
pipeline" on the day a redundant training script was written anyway. This is the mechanical version:
it lists every `*.sbatch` (and the `scripts/prepare/*.py` tools) that the registry does not mention,
so an unregistered script is visible rather than merely discouraged.

It deliberately does NOT try to detect duplication automatically. Deciding whether two scripts do the
same job needs judgement; what a checker can do reliably is refuse to let a new entry point be
invisible.

Exit status is 1 when unregistered scripts exist, so it can gate a commit hook. Scripts carrying a
`SUPERSEDED` header are reported in their own section and do NOT fail the check -- they are retained
for provenance and are meant to be absent from the registry.

cas12a env. Run: `python3 scripts/check_entrypoints.py`
"""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
REGISTRY = ROOT / "docs" / "CANONICAL_ENTRYPOINTS.md"
# Scripts that predate the registry. Flagging all 159 of them would make the checker noise, and a
# checker that always fails is a checker nobody runs. The baseline freezes the back catalogue so the
# check has exactly one job: catch a NEW entry point that skipped the registry.
BASELINE = ROOT / "docs" / "entrypoints_baseline.txt"


def main():
    if not REGISTRY.exists():
        sys.exit(f"missing registry {REGISTRY}")
    reg = REGISTRY.read_text()
    # A script counts as registered if its path appears anywhere in the registry prose.
    named = set(re.findall(r"scripts/[A-Za-z0-9_/]+\.(?:sbatch|py|sh)", reg))

    candidates = sorted(
        [p for p in (ROOT / "scripts").rglob("*.sbatch")]
        + [p for p in (ROOT / "scripts" / "prepare").glob("*.py")],
        key=lambda p: str(p))

    base = set()
    if BASELINE.exists():
        base = {l.strip() for l in BASELINE.read_text().splitlines()
                if l.strip() and not l.startswith("#")}
    unregistered, superseded, grandfathered = [], [], []
    for p in candidates:
        rel = str(p.relative_to(ROOT))
        if rel in named:
            continue
        head = "\n".join(p.read_text(errors="replace").splitlines()[:15])
        if "SUPERSEDED" in head:
            superseded.append(rel)
        elif rel in base:
            grandfathered.append(rel)
        else:
            unregistered.append(rel)

    print(f"  registered: {len(candidates)-len(superseded)-len(grandfathered)-len(unregistered)}"
          f" | grandfathered: {len(grandfathered)} | superseded: {len(superseded)}")
    if superseded:
        print(f"  SUPERSEDED, retained for provenance ({len(superseded)}) -- not a failure")
    if unregistered:
        print(f"\n  UNREGISTERED ({len(unregistered)}) -- add to docs/CANONICAL_ENTRYPOINTS.md, or")
        print("  extend an existing entry instead of keeping a new script:")
        for r in unregistered:
            print(f"    {r}")
        print(f"\n  {len(unregistered)} NEW unregistered script(s). Either extend an existing entry,")
        print("  or add it to the registry with a one-line reason no entry fit.")
        return 1
    print("\n  OK: no NEW unregistered entry points.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
