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
version = qsdsan.__version__
# The full version, including alpha/beta/rc tags.
release = version

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
'nbsphinx',
'sphinx_copybutton',
'sphinx_design',
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
# pygments_style = 'manni'

# Webpage dark mode
# default_dark_mode = False

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'furo'

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
 	'css/copybutton.css',
# 	'css/theme_overrides.css',
 	]
html_js_files = [
    'js/unit-operation-filters.js',
]

# Docs chatbot widget (internal-first): register the assets only on the dedicated
# 'chatbot' Read the Docs version, and on local builds for development. Set
# CHATBOT_WIDGET_VERSIONS (comma-separated) to override which RTD versions show it.
import os as _os
_chatbot_versions = _os.environ.get("CHATBOT_WIDGET_VERSIONS", "chatbot").split(",")
_rtd_version = _os.environ.get("READTHEDOCS_VERSION")
if _rtd_version is None or _rtd_version in _chatbot_versions:
    html_css_files.append('css/chatbot.css')
    html_js_files.append('js/chatbot.js')
del _os, _chatbot_versions, _rtd_version

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
# html_logo = 'logo_dark_small.png'
html_theme_options = {
    'dark_logo': 'logo_dark_mode_m.png',
    'light_logo': 'logo_light_mode_m.png',
	'sidebar_hide_name': True,
	'top_of_page_button': 'edit', # only edit or None is supported
    'announcement': (
        '📣 Stay in the loop: join the '
        '<a href="https://groups.google.com/g/qsdsan" target="_blank" rel="noopener">'
        'QSDsan Google Group</a> for updates.'
    ),
    'footer_icons': [
        {
            'name': 'GitHub',
            'url': 'https://github.com/QSD-Group/QSDsan',
            'html': """<svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16"><title>GitHub</title><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path></svg>""",
            'class': '',
        },
        {
            'name': 'PyPI',
            'url': 'https://pypi.org/project/qsdsan/',
            'html': """<svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 24 24"><title>PyPI</title><path d="M12 2 2 7v10l10 5 10-5V7L12 2zm0 2.18L18.5 7.5 12 10.82 5.5 7.5 12 4.18zM4 9.18l7 3.5v7.14l-7-3.5V9.18zm9 10.64v-7.14l7-3.5v7.14l-7 3.5z"></path></svg>""",
            'class': '',
        },
        {
            'name': 'YouTube',
            'url': 'https://www.youtube.com/channel/UC8fyVeo9xf10KeuZ_4vC_GA',
            'html': """<svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 24 24"><title>YouTube</title><path d="M23.498 6.186a3.016 3.016 0 0 0-2.122-2.136C19.505 3.545 12 3.545 12 3.545s-7.505 0-9.377.505A3.017 3.017 0 0 0 .502 6.186C0 8.07 0 12 0 12s0 3.93.502 5.814a3.016 3.016 0 0 0 2.122 2.136c1.871.505 9.376.505 9.376.505s7.505 0 9.377-.505a3.015 3.015 0 0 0 2.122-2.136C24 15.93 24 12 24 12s0-3.93-.502-5.814zM9.545 15.568V8.432L15.818 12l-6.273 3.568z"></path></svg>""",
            'class': '',
        },
        {
            'name': 'Email',
            'url': 'mailto:quantitative.sustainable.design@gmail.com',
            'html': """<svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 24 24"><title>Email</title><path d="M20 4H4c-1.1 0-1.99.9-1.99 2L2 18c0 1.1.9 2 2 2h16c1.1 0 2-.9 2-2V6c0-1.1-.9-2-2-2zm0 4-8 5-8-5V6l8 5 8-5v2z"></path></svg>""",
            'class': '',
        },
        {
            'name': 'Google Group',
            'url': 'https://groups.google.com/g/qsdsan',
            'html': """<svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 24 24"><title>Google Group</title><path d="M16 11c1.66 0 2.99-1.34 2.99-3S17.66 5 16 5c-1.66 0-3 1.34-3 3s1.34 3 3 3zm-8 0c1.66 0 2.99-1.34 2.99-3S9.66 5 8 5C6.34 5 5 6.34 5 8s1.34 3 3 3zm0 2c-2.33 0-7 1.17-7 3.5V19h14v-2.5c0-2.33-4.67-3.5-7-3.5zm8 0c-.29 0-.62.02-.97.05 1.16.84 1.97 1.97 1.97 3.45V19h6v-2.5c0-2.33-4.67-3.5-7-3.5z"></path></svg>""",
            'class': '',
        },
    ],
}


# Canonical site URL (docs are served at qsdsan.com); sets <link rel="canonical">.
html_baseurl = 'https://qsdsan.com'

# Output file base name for HTML help builder.
htmlhelp_basename = 'QSDsan_documentation'


# -- Extension settings -------------------------------------------------------
# napoleon_custom_sections = ['Tips']
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True


# -- External mapping -------------------------------------------------------
intersphinx_mapping = {
	'biosteam': ('https://biosteam.readthedocs.io/en/latest/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy', None),
    'SALib': ('https://salib.readthedocs.io/en/latest/', None),
}
