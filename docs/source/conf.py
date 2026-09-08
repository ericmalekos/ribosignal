"""Sphinx configuration for the RiboSignal documentation."""

project = "RiboSignal"
author = "Eric Malekos"
copyright = "2026, Eric Malekos"

extensions = ["myst_parser"]
myst_enable_extensions = ["colon_fence", "deflist"]

source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
master_doc = "index"
exclude_patterns = ["_build"]

html_theme = "furo"
html_title = "RiboSignal"
html_static_path = []

# Warnings are errors in CI (.readthedocs.yaml fail_on_warning), so a broken cross-reference
# fails the build instead of shipping a dead link.
nitpicky = False
