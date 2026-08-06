# Minimal Makefile for the Ribo-seq signal-model tutorial.
# Uses the sphinx-docs conda env by default; override SPHINXBUILD to point elsewhere.
SPHINXBUILD ?= /private/groups/carpenterlab/emalekos/conda_envs/sphinx-docs/bin/sphinx-build
AUTOBUILD   ?= /private/groups/carpenterlab/emalekos/conda_envs/sphinx-docs/bin/sphinx-autobuild
SOURCEDIR    = source
BUILDDIR     = build

PLOTPY ?= /private/groups/carpenterlab/emalekos/conda_envs/cas12a/bin/python3

html:                         ## build the static site (regenerates the benchmark figure first)
	$(PLOTPY) ../scripts/tutorial/plot_benchmarks.py || true
	rm -rf "$(BUILDDIR)/doctrees" "$(BUILDDIR)/html/_images"   # force figures to always re-copy fresh
	$(SPHINXBUILD) -b html -d "$(BUILDDIR)/doctrees" "$(SOURCEDIR)" "$(BUILDDIR)/html"

strict:                       ## fail the build on broken xrefs / missing images
	$(SPHINXBUILD) -W -b html "$(SOURCEDIR)" "$(BUILDDIR)/html"

serve:                        ## live-reloading local preview
	$(AUTOBUILD) "$(SOURCEDIR)" "$(BUILDDIR)/html"

clean:
	rm -rf "$(BUILDDIR)"

.PHONY: html strict serve clean
