# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ONETEP Tutorials'
copyright = '2024, ONETEP Collective'
author = 'ONETEP Collective'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx_rtd_dark_mode","sphinx_favicon"]
default_dark_mode = False

templates_path = ['_templates']
exclude_patterns = []
latex_engine = 'lualatex'
latex_elements = {
'papersize': 'a4paper',
'pointsize': '12pt',
'figure_align': '!htbp',
}

# jd: For proper referencing of figures.
numfig = True


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
html_logo = '_static/onetep_logo.svg'
html_theme_options = {
    'logo_only': True,
    'display_version': True,
}

# --- Favicon ---
favicons = [
    "favicon.ico",
]
