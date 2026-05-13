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
import datetime as dt
import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'Pepper'
copyright = f'{dt.datetime.now().year}, Pepper Maintainers'
author = 'Pepper Maintainers'
mattermost_channel_url = "https://mattermost.web.cern.ch/cms-exp/channels"
gitlab_repository_url = "https://gitlab.cern.ch/cms-analysis/general/pepper/pepper"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.extlinks",
    "sphinx_autodoc_typehints",
    "sphinx_design",
    "sphinx_copybutton",
]


autosummary_generate = True
napoleon_google_docstring = True
napoleon_numpy_docstring = True
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'special-members': '__init__',
    'show-inheritance': True,
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**/.venv/**', '**/site-packages/**']

# COPY code
copybutton_exclude = '.linenos, .gp'
copybutton_prompt_text = "Copy to clipboard"
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

html_title = f"{project} Documentation"
html_favicon = "_static/icons/pepper_logo.svg"
html_theme = 'pydata_sphinx_theme'
html_theme_options = {
    "logo": {
        "text": f"{project}",      # Show name next to logo (important!)
        "image_light": "_static/icons/pepper_logo_light.svg",
        "image_dark": "_static/icons/pepper_logo_dark.svg",
    },
    "logo_link": "/",             # clicking logo -> homepage
    "show_toc_level": 2,
    "navigation_depth": 2,
    "collapse_navigation": True,
    "icon_links": [
        {
            "name": "GitLab",
            "url": gitlab_repository_url,
            "icon": "fab fa-gitlab",      # FontAwesome icon
            "type": "fontawesome",
        },
        {
            "name": "CERN Mattermost",
            "url": f"{mattermost_channel_url}/pepper-users",
            "icon": "fa-custom fa-mattermost",
            "type": "fontawesome",
        }
    ],
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


html_css_files = [
    "css/custom.css",
]

html_js_files = [
   "js/mattermost-icon.js",
]

extlinks = {
    "repo": (f"{gitlab_repository_url}/-/blob/master/%s", "%s"),
    "mattermost": (f"{mattermost_channel_url}/%s", "%s")
}