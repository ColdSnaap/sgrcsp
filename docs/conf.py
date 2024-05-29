# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('../src'))  # Adjust the path to include your project module if necessary

# -- Project information -----------------------------------------------------

project = 'SGRCSP'
author = 'ColdSnaap'
release = '0.0.1'
version = '0.0.1'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',    # Include documentation from docstrings
    'sphinx.ext.napoleon',   # Support for Google and NumPy style docstrings
    'sphinx.ext.viewcode',   # Add links to highlighted source code
    'sphinx.ext.githubpages' # Create .nojekyll file to publish the docs on GitHub Pages
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

html_theme = 'alabaster'
html_static_path = ['_static']

# -- Options for autodoc -----------------------------------------------------

autoclass_content = 'both'  # Include both class docstring and __init__ docstring