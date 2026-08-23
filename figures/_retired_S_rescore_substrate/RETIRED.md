# S_rescore_substrate -- RETIRED 2026-08-15

This figure cannot exist without A549, which was removed from the project on 2026-08-15.

**Its entire claim was a CONTRAST between two substrates**: MS2Rescore gave opposite answers on the
A549 tryptic whole proteome (no gain -- the novel class is untrainable) and the HBL-1 HLA-I
immunopeptidome (a large gain -- mokapot learns the novel class). A549 was the only tryptic dataset
in the project. With it gone there is one substrate left, so there is no contrast to draw and the
panel has nothing to say.

It is retired rather than rebuilt because no rebuild is possible: the missing half is a whole
substrate class, not a missing number.

**What is lost.** The panel was the justification for a pipeline-wide decision -- reporting RAW
hyperscore class-specific FDR everywhere instead of rescored scores. That decision is still correct
and is still documented in `docs/PIPELINE_POLICY.md` and `proteogenomics/methods.md`; it simply no
longer has a figure demonstrating why.

Directories prefixed `_` are skipped by `scripts/build_poster_package.py`, so this no longer ships.
