# -*- encoding: utf-8 -*-
"""
Latex Compilers
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import os
from sage.interfaces.cmdline import Tool, ToolAdapter, CompilerTool


SAMPLE_LATEX = r"""
\documentclass{article}
\begin{document}
Hello World!
\begin{equation}
  x_{1,2} =
  \frac{
    -b \pm \sqrt{b^2-4ac}
  }{
    2a
  }
\end{equation}
\end{document}
"""


def test_compiles_to_pdf(tool):
    """
    Test that a tool compiles TeX to PDF

    INPUT:

    - ``tool`` -- :class:`sage.interfaces.cmdline.tool.Tool`
      instance. The tool to test.

    EXAMPLES::

        sage: from sage.interfaces.cmdline.posix import cat
        sage: from sage.interfaces.cmdline.latex import test_compiles_to_pdf
        sage: test_compiles_to_pdf(cat)
        Traceback (most recent call last):
        ...
        AssertionError: cat did not return valid pdf
    """
    pdf = tool(SAMPLE_LATEX)
    assert len(pdf) > 0, \
        '{0} ran without error but did not produce output'.format(tool)
    assert pdf.startswith('%PDF-'), \
        '{0} did not return valid pdf'.format(tool)


def test_compiles_to_dvi(tool):
    """
    Test that a tool compiles TeX to DVI

    INPUT:

    - ``tool`` -- :class:`sage.interfaces.cmdline.tool.Tool`
      instance. The tool to test.

    EXAMPLES::

        sage: from sage.interfaces.cmdline.posix import cat
        sage: from sage.interfaces.cmdline.latex import test_compiles_to_dvi
        sage: test_compiles_to_dvi(cat)
        Traceback (most recent call last):
        ...
        AssertionError: cat did not return valid dvi
    """
    dvi = tool(SAMPLE_LATEX)
    assert len(dvi) > 0, \
        '{0} ran without error but did not produce output'.format(tool)
    assert 'TeX output' in dvi, \
        '{0} did not return valid dvi'.format(tool)


def test_convert_to_pdf(tool):
    """
    Test that a tool converts DVI to PDF

    INPUT:

    - ``tool`` -- :class:`sage.interfaces.cmdline.tool.Tool`
      instance. The tool to test.

    EXAMPLES::

        sage: from sage.interfaces.cmdline.posix import cat
        sage: from sage.interfaces.cmdline.latex import test_convert_to_pdf
        sage: test_convert_to_pdf(cat)
        Traceback (most recent call last):
        ...
        AssertionError: cat did not return valid pdf
    """
    from sage.env import SAGE_EXTCODE
    filename = os.path.join(SAGE_EXTCODE, 'doctest', 'rich_output', 'example.dvi')
    with open(filename, 'r') as f:
        dvi = f.read()
    pdf = tool(dvi)
    assert len(pdf) > 0, \
        '{0} ran without error but did not produce output'.format(tool)
    assert pdf.startswith('%PDF-'), \
        '{0} did not return valid pdf'.format(tool)


pdflatex = CompilerTool(
    'pdflatex',
    ['pdflatex', '-interaction', 'nonstopmode', 'document.tex'],
    'Provided by TeXlive (https://www.tug.org/texlive/)',
    test_compiles_to_pdf,
    filename=('document.tex', 'document.pdf'),
    rpm='texlive-latex',
)


latex = CompilerTool(
    'latex',
    ['latex', '-interaction', 'nonstopmode', 'document.tex'],
    'Provided by TeXlive (https://www.tug.org/texlive/)',
    test_compiles_to_dvi,
    filename=('document.tex', 'document.dvi'),
    rpm='texlive-latex',
)


xelatex = CompilerTool(
    'xelatex',
    ['xelatex', '-interaction', 'nonstopmode', 'document.tex'],
    'Provided by TeXlive (https://www.tug.org/texlive/)',
    test_compiles_to_pdf,
    filename=('document.tex', 'document.pdf'),
    rpm='texlive-xelatex',
)


dvipdf = ToolAdapter(
    'dvipdf',
    ['dvipdf', '{input}', '{output}'],
    'Part of Ghostscript (http://www.ghostscript.com)',
    test_convert_to_pdf,
    rpm='ghostscript',
)


convert_pdf = ToolAdapter(
    'convert dvi->pdf',
    ['convert', 'dvi:{input}', 'pdf:{output}'],
    'Part of the ImageMagick suite of tools',
    test_convert_to_pdf,
    rpm='ImageMagick',
)



def compile_latex(source, implementation=None):
    """
    Compile LaTeX source to PDF

    INPUT:

    - ``source`` -- string. The LaTeX source code.

    - ``implementation`` -- string. The LaTeX implementation to
      use. Must be one of:

        - ``None`` -- default. Choose automatically.

        - ``latex`` -- use the obsolete latex -> dvi -> pdf conversion
          chain.

        - ``pdflatex`` -- use the new pdflatex. Preferred implementation.

        - ``xelatex`` -- use xelatex.

    OUTPUT:

    A PDF Document wrapped in a
    :class:`~sage.repl.portable_document_format.PortableDocumentFormat`
    instance.

    EXAMPLES::

        sage: src = '\n'.join([
        ....:     r'\documentclass{article}',
        ....:     r'\begin{document}',
        ....:     r'  This is a test',
        ....:     r'\end{document}',
        ....:  ])
        sage: from sage.interfaces.cmdline.latex import compile_latex
        sage: compile_latex(src)
        PDF document (... bytes)
    """
    if implementation is None:
        compiler = Tool.require(pdflatex, xelatex, latex)
    elif implementation == 'pdflatex':
        compiler = Tool.require(pdflatex)
    elif implementation == 'xelatex':
        compiler = Tool.require(xelatex)
    elif implementation == 'latex':
        compiler = Tool.require(latex)
    else:
        raise ValueError('unknown implementation: {0}'.format(implementation))
    result = compiler(source)
    if implementation == latex:
        converter = Tool.require(dvipdf)
        pdf = converter(result, convert_pdf)
    else:
        pdf = result
    from sage.repl.portable_document_format import PortableDocumentFormat
    return PortableDocumentFormat(pdf)

    
