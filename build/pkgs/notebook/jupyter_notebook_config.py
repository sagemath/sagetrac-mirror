# Configuration file for Sage's builtin Jupyter notebook server

# Note for distributors: Sage uses mathjax, so the notebook server
# needs to have the mathjax_url set to wherever your distribution
# installs mathjax.

c.NotebookApp.mathjax_url = 'nbextensions/mathjax/MathJax.js'

# Disable XSRF checking to fix Thebe. See
# * https://trac.sagemath.org/ticket/22458
# * https://github.com/oreillymedia/thebe/issues/93

c.NotebookApp.disable_check_xsrf = True
