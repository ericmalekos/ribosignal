# Sphinx configuration for the Ribo-seq signal-model tutorial.
import os
import re
import shutil

project = "Ribo-seq signal model"
author = "Eric Malekos"
copyright = "2026, Eric Malekos"
release = "0.1"

extensions = [
    "myst_parser",            # author pages in Markdown
    "sphinx_copybutton",      # copy button on code blocks
    "sphinx_design",          # cards / grids / tabs
    "sphinxcontrib.mermaid",  # architecture / pipeline diagrams
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "attrs_inline",
    "tasklist",
    "fieldlist",
    "dollarmath",             # inline/blocks of LaTeX for the metric definitions
]
myst_heading_anchors = 3

source_suffix = {".md": "markdown"}
exclude_patterns = ["build", "_build", "_code", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_title = "Ribo-seq signal model"
html_static_path = ["_static"]
html_theme_options = {
    "navigation_depth": 3,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "titles_only": False,
}

mermaid_init_js = (
    "mermaid.initialize({startOnLoad:true, theme:'neutral', "
    "flowchart:{curve:'basis', nodeSpacing:45, rankSpacing:45, useMaxWidth:true}});"
)

# --- Build hooks ----------------------------------------------------------------------------------
# 1. Mirror scripts/ -> source/_code/scripts/ with machine-specific paths scrubbed, so {literalinclude}
#    can show real code without leaking the working-directory layout.
# 2. Mirror figures/<name>/*.png|pdf -> source/img/<name>/ so pages embed figures with
#    {figure} /img/<name>/<file>.png. Both regenerate every build (source/img is gitignored).
_CONF_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_CONF_DIR, "..", ".."))
_PATH_SUBS = [
    ("/private/home/emalekos/.local/bin/micromamba", "/path/to/micromamba"),
    ("/private/groups/carpenterlab/emalekos/conda_envs", "/path/to/conda_envs"),
    ("/private/groups/carpenterlab/emalekos/singularity_cache", "/path/to/singularity_cache"),
    ("/private/groups/carpenterlab/emalekos/STAR_indexes", "/path/to/STAR_indexes"),
    ("/private/groups/carpenterlab/emalekos/Salmon_indexes", "/path/to/Salmon_indexes"),
    ("/private/groups/carpenterlab/emalekos/genomes", "/path/to/genomes"),
    ("/private/groups/carpenterlab/emalekos/RNAZoo_meta", "/path/to/project"),
    ("/private/groups/carpenterlab/emalekos", "/path/to/data"),
    ("/private/home/emalekos", "/path/to/home"),
    ("/data/tmp/emalekos", "/path/to/tmp"),
]
_SUBS = [(re.compile(re.escape(a)), b) for a, b in _PATH_SUBS]


def _sanitize_code(app=None):
    src = os.path.join(_REPO, "scripts")
    dst = os.path.join(_CONF_DIR, "_code", "scripts")
    if not os.path.isdir(src):
        return
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    for root, _dirs, files in os.walk(src):
        if "__pycache__" in root:
            continue
        out = os.path.join(dst, os.path.relpath(root, src))
        os.makedirs(out, exist_ok=True)
        for f in files:
            try:
                with open(os.path.join(root, f), "r", encoding="utf-8") as fh:
                    text = fh.read()
            except (UnicodeDecodeError, IsADirectoryError, OSError):
                continue
            for rx, repl in _SUBS:
                text = rx.sub(repl, text)
            with open(os.path.join(out, f), "w", encoding="utf-8") as fh:
                fh.write(text)


def _mirror_figures(app=None):
    src = os.path.join(_REPO, "figures")
    dst = os.path.join(_CONF_DIR, "img")
    if not os.path.isdir(src):
        return
    os.makedirs(dst, exist_ok=True)
    for root, _dirs, files in os.walk(src):
        rel = os.path.relpath(root, src)
        for f in files:
            if not f.lower().endswith((".png", ".svg")):
                continue
            outdir = os.path.join(dst, rel) if rel != "." else dst
            os.makedirs(outdir, exist_ok=True)
            try:
                shutil.copy2(os.path.join(root, f), os.path.join(outdir, f))
            except OSError:
                continue


_sanitize_code()
_mirror_figures()


def setup(app):
    app.connect("builder-inited", _sanitize_code)
    app.connect("builder-inited", _mirror_figures)
