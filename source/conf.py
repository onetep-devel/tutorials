# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Onetep Tutorials'
copyright = '2023, Onetep Collective'
author = 'Onetep Collective'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx_rtd_dark_mode"]
default_dark_mode = False

templates_path = ['_templates']
exclude_patterns = []
latex_engine = 'lualatex'
latex_elements = {
'papersize': 'a4paper',
'pointsize': '12pt',
'figure_align': '!htbp',
}




# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
