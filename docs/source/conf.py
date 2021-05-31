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

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = '0.2.10'
# The full version, including alpha/beta/rc tags.
release = version

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
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
# html_theme_options = {}

# Add any paths that contain custom themes here, relative to this directory.
# html_theme_path = []

# The name for this set of Sphinx documents.
# "<project> v<release> documentation" by default.
html_title = f'QSDsan {version}'

# A shorter title for the navigation bar.  Default is the same as html_title.
html_short_title = 'QSDsan'

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = '_static/logo_dark_small.png'

# The name of an image file (relative to this directory) to use as a favicon of
# the docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = '_static/favicon.png'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_css_files = [
	'css/qsdsan.css',
	'css/theme_overrides.css',
	]

# -- Extension settings -------------------------------------------------------
# napoleon_custom_sections = ['Tips']

# -- External mapping -------------------------------------------------------
intersphinx_mapping = {
	'biosteam': ('https://biosteam.readthedocs.io/en/latest', None),
	'thermosteam': ('https://thermosteam.readthedocs.io/en/latest', None),
	'BioSTEAM': ('https://biosteam.readthedocs.io/en/latest', None),
	'Thermosteam': ('https://thermosteam.readthedocs.io/en/latest', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference/', None),
    'SALib': ('https://salib.readthedocs.io/en/latest', None),
}
