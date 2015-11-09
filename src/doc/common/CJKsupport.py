"""
Hack around Sphinx's current lack of support for CJK
"""

import os
from subprocess import check_output


CJK_LANGUAGES = ['zh_CN']

def enable_if_necessary(lang, latex_directory):
    """
    Check whether the CJK hack is necessary, and perform it optionally.
    """
    if lang not in CJK_LANGUAGES:
        return False
    if not have_xelatex():
        print "Disabling pdf documentation for {0} because it requires xelatex,".format(lang)
        print "which you do not have installed"
        disable(latex_directory, lang)
        return False
    fixup_makefile(os.path.join(latex_directory, 'Makefile'), lang)
    return True


HAVE_XELATEX = None

def have_xelatex():
    """
    Return whether XeTeX is installed

    EXAMPLES::

        sage: import os, sys; sys.path.append(os.environ['SAGE_DOC']+'/common/'); import CJKsupport
        sage: CJKsupport.have_xelatex() in [True, False]
        True
    """
    global HAVE_XELATEX
    if HAVE_XELATEX is None:
        try:
            HAVE_XELATEX = 'XeTeX' in check_output(['xelatex', '-v'])
        except OSError:
            HAVE_XELATEX = False
    return HAVE_XELATEX


def fixup_makefile(makefile, lang):
    """
    Switch to XeLaTeX
    
    Sphinx hardcodes pdflatex, but for CJK you pretty much have to use
    xelatex.
    """
    with open(makefile, 'r') as f:
        src = f.read()
    src = src.replace('pdflatex', 'xelatex', -1)
    with open(makefile, 'w') as f:
        f.write(src)


DISABLED_LATEX_FILE = r"""
\documentclass{article}
\begin{document}

Running LaTeX for documents containing CJK fonts requires XeTeX, which
you do not have installed. Consult your latex installation manual to
find out how you can add it.
\begin{itemize}
\item 
  If you are using TeXlive, you can use 
  \texttt{tlmgr install common-xetex}
\end{itemize}

\end{document}
"""

def disable(latex_dir, lang):
    """
    Disable pdf build
    
    Sage still expects a pdf output, so we write a dummy file
    explaining how to enable it.
    """
    for name in os.listdir(latex_dir):
        if not name.endswith('.tex'):
            continue
        os.rename(
            os.path.join(latex_dir, name),
            os.path.join(latex_dir, name+'-disabled'))
    with open(os.path.join(latex_dir, 'disabled.tex'), 'w') as f:
        f.write(DISABLED_LATEX_FILE)
