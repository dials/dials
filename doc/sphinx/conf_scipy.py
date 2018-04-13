from __future__ import absolute_import, division, print_function

import datetime

# needs_sphinx = '1.1'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.autosummary',
    'matplotlib.sphinxext.plot_directive',
    'sphinx.ext.viewcode',
    'libtbx.sphinx.phil',
    'libtbx.sphinx.python_string'
]

# Add CDN path for mathjax script, converting Latex to readable text on the fly.
# (The Sphinx builtin path is deprecated.)
mathjax_path = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'

project = u'DIALS'
copyright = u'%d, Diamond Light Source, Lawrence Berkeley National Laboratory and STFC' % datetime.datetime.now().year

version = ''
release = ''
exclude_patterns = ['_build', 'scipy-sphinx-theme']
pygments_style = 'sphinx'

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
  'http://docs.python.org/': None,
  'http://cci.lbl.gov/cctbx_docs/': None
}

# -- Options for HTML output ---------------------------------------------------

html_theme = 'scipy'
html_theme_path = ['scipy-sphinx-theme/_theme']
#html_logo = '_static/scipyshiny_small.png'
html_static_path = ['_static']
html_theme_options = {
    "edit_link": "true",
    "sidebar": "right",
    "scipy_org_logo": "false",
    "rootlinks": [("dials://diamond.ac.uk/doc/", "DIALS"),]
}

imgmath_latex_preamble = r"""
\usepackage{color}
\definecolor{textgray}{RGB}{51,51,51}
\color{textgray}
"""
imgmath_use_preview = True
imgmath_dvipng_args = ['-gamma 1.5', '-D 96', '-bg Transparent']

#------------------------------------------------------------------------------
# Plot style
#------------------------------------------------------------------------------

plot_pre_code = """
import numpy as np
import scipy as sp
np.random.seed(123)
"""
plot_include_source = True
plot_formats = [('png', 96), 'pdf']
plot_html_show_formats = False

import math
phi = (math.sqrt(5) + 1)/2

font_size = 13*72/96.0  # 13 px

plot_rcparams = {
    'font.size': font_size,
    'axes.titlesize': font_size,
    'axes.labelsize': font_size,
    'xtick.labelsize': font_size,
    'ytick.labelsize': font_size,
    'legend.fontsize': font_size,
    'figure.figsize': (3*phi, 3),
    'figure.subplot.bottom': 0.2,
    'figure.subplot.left': 0.2,
    'figure.subplot.right': 0.9,
    'figure.subplot.top': 0.85,
    'figure.subplot.wspace': 0.4,
    'text.usetex': False,
}
