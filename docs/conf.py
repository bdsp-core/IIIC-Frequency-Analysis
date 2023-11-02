# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "pdDetector"
copyright = "2023, Dr. Chris McGraw, Shashank Manjunath, Samvrit Rao"
author = "Dr. Chris McGraw, Shashank Manjunath, Samvrit Rao"

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join("..")))
sys.path.insert(0, os.path.abspath(os.path.join("..", "pdDetector")))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_rtd_theme",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Napoleon settings

napoleon_custom_sections = [("Returns", "params_style")]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
#  html_static_path = ["_static"]
