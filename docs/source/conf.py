# Configuration file for the Sphinx documentation builder.
# URI BEM Analysis Pipeline

import os
import sys

# -- Project information -----------------------------------------------------
project = 'URI BEM Analysis Pipeline'
copyright = '2025, University of Rhode Island'
author = 'URI Ocean Engineering'
version = '1.0'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinxcontrib.bibtex',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.mathjax',
    'sphinx_copybutton',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['css/custom.css']
html_title = 'URI BEM Analysis Pipeline'

html_theme_options = {
    'navigation_depth': 4,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'titles_only': False,
}

# -- Figure and equation numbering ------------------------------------------
numfig = True
numfig_format = {
    'figure': 'Figure %s',
    'table': 'Table %s',
    'code-block': 'Listing %s',
}

# -- Bibliography (sphinxcontrib-bibtex) -------------------------------------
bibtex_bibfiles = ['bibliography.bib']
bibtex_default_style = 'unsrt'

# -- Autosection labels ------------------------------------------------------
autosectionlabel_prefix_document = True
autosectionlabel_maxdepth = 4

# -- MathJax -----------------------------------------------------------------
mathjax3_config = {
    'tex': {
        'tags': 'ams',
        'macros': {
            'bm': [r'\boldsymbol{#1}', 1],
        },
    },
}
