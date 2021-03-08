# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

import os, sys
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(0, os.path.abspath('../../../tmo'))
sys.path.insert(0, os.path.abspath('../../../bst'))
del os, sys

# -- Project information -----------------------------------------------------

import time, qsdsan

project = 'QSDsan'
author = 'Quantitative Sustainable Design Group'
copyright = f'2020-{time.gmtime().tm_year}, Quantitative Sustainable Design Group'
# version = qsdsan.__version__
# The full version, including alpha/beta/rc tags
release = '0.0.1' if not qsdsan.__version__ else qsdsan.__version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
'nbsphinx',
'sphinx.ext.autodoc',
'sphinx.ext.intersphinx',
'sphinx.ext.mathjax',
'sphinx.ext.napoleon',
]

# Allow exceptions to occur in notebooks
nbsphinx_allow_errors = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'manni'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['css/qsdsan.css']

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

# -- Extension settings -------------------------------------------------------
# napoleon_custom_sections = [
# 'Reference documents']

# -- External mapping -------------------------------------------------------
intersphinx_mapping = {
	'BioSTEAM': ('https://biosteam.readthedocs.io/en/latest', None),
	'Thermosteam': ('https://thermosteam.readthedocs.io/en/latest', None),
	'chemicals': ('https://chemicals.readthedocs.io/en/latest', None),
    'SALib': ('https://salib.readthedocs.io/en/latest', None),
}