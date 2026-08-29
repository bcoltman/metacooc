from pathlib import Path
import os
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/metacooc-matplotlib")

from metacooc.version import __version__  # noqa: E402


project = "MetaCoOc"
copyright = "2026, MetaCoOc contributors"
author = "MetaCoOc contributors"
release = __version__

extensions = [
    "myst_parser",
    "sphinxarg.ext",
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
]
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "sphinx_rtd_theme"
source_suffix = {".md": "markdown"}
myst_heading_anchors = 3
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
autodoc_typehints = "description"
