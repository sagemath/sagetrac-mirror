# -*- encoding: utf-8 -*-
"""
Commandline Tools for PDF Documents
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
from cStringIO import StringIO
from sage.interfaces.cmdline import Tool, ToolAdapter


def test_convert_to_svg(tool):
    """
    Test that a tool converts to SVG

    INPUT:

    - ``tool`` -- :class:`sage.interfaces.cmdline.tool.Tool`
      instance. The tool to test.

    EXAMPLES::

        sage: from sage.interfaces.cmdline.posix import cat
        sage: from sage.interfaces.cmdline.pdf import test_convert_to_svg
        sage: test_convert_to_svg(cat)
        Traceback (most recent call last):
        ...
        Error: Syntax error at line 1: illegal data at start of file
    """
    from sage.env import SAGE_EXTCODE
    filename = os.path.join(SAGE_EXTCODE, 'doctest', 'rich_output', 'example.pdf')
    with open(filename, 'r') as f:
        pdf = f.read()
    svg = tool(pdf)
    assert len(svg) > 0, \
        '{0} ran without error but did not produce output'.format(tool)
    from xmllib import XMLParser
    XMLParser().feed(svg)   # will raise error if svg is invalid
    

def test_convert_to_png(tool):
    """
    Test that a tool converts to PNG

    INPUT:

    - ``tool`` -- :class:`sage.interfaces.cmdline.tool.Tool`
      instance. The tool to test.

    EXAMPLES::

        sage: from sage.interfaces.cmdline.posix import cat
        sage: from sage.interfaces.cmdline.pdf import test_convert_to_png
        sage: test_convert_to_png(cat)
        Traceback (most recent call last):
        ...
        IOError: cannot identify image file ...
    """
    from sage.env import SAGE_EXTCODE
    filename = os.path.join(SAGE_EXTCODE, 'doctest', 'rich_output', 'example.pdf')
    with open(filename, 'r') as f:
        pdf = f.read()
    # We use an odd dpi value to check that it is used correctly
    png = tool(pdf, dpi=123)
    from PIL import Image
    img = Image.open(StringIO(png))
    img.verify()
    assert img.size == (112, 54), \
        'image has the wrong resolution, {0} did not respect dpi'.format(tool)


def test_output_is_pdf(tool):
    """
    Test that a tool returns PDF

    INPUT:

    - ``tool`` -- :class:`sage.interfaces.cmdline.tool.Tool`
      instance. The tool to test.

    EXAMPLES::

        sage: from sage.interfaces.cmdline.posix import cat
        sage: from sage.interfaces.cmdline.pdf import test_output_is_pdf
        sage: test_output_is_pdf(cat)
    """
    from sage.env import SAGE_EXTCODE
    filename = os.path.join(SAGE_EXTCODE, 'doctest', 'rich_output', 'example.pdf')
    with open(filename, 'r') as f:
        pdf = f.read()
    pdf = tool(pdf)
    assert len(pdf) > 0, \
        '{0} ran without error but did not produce output'.format(tool)
    assert pdf.startswith('%PDF-'), \
        '{0} did not return valid pdf'.format(tool)
    

pdfjam = ToolAdapter(
    'pdfjam',
    ['pdfjam', '--outfile', '{output}', '--', '{input}'],
    'The homepage of PDFjam is at http://go.warwick.ac.uk/pdfjam',
    test_output_is_pdf,
    rpm='texlive-pdfjam',
    file_ext=['.pdf', '.pdf'],
)


# Bug: "pdftocairo -png -singlefile - -" should work but does not
# pdftocairo_png = Tool(
#     'pdftocairo pdf->png',
#     ['pdftocairo', '-singlefile', '-png', '-', '-'],
#     'Part of the Poppler PDF library utilities',
#     test_convert_to_png,
#     rpm='poppler-utils',
# )


pdftocairo_svg = Tool(
    'pdftocairo pdf->svg',
    ['pdftocairo', '-svg', '-', '-'],
    'Part of the Poppler PDF library utilities',
    test_convert_to_svg,
    rpm='poppler-utils',
)


pdf2svg = ToolAdapter(
    'pdf2svg',
    ['pdf2svg', '{input}', '{output}'],
    'Homepage at http://www.cityinthesky.co.uk/opensource/pdf2svg/',
    test_convert_to_svg,
    rpm='pdf2svg',
)
    

convert_svg = Tool(
    'convert pdf->svg',
    ['convert', 'pdf:', 'svg:-'],
    'Part of the ImageMagick suite of tools',
    test_convert_to_svg,
    rpm='ImageMagick',
)


convert_png = Tool(
    'convert pdf->png',
    ['convert', '-units', 'PixelsPerInch', '-density', '{dpi}', 'pdf:', 'png:-'],
    'Part of the ImageMagick suite of tools',
    test_convert_to_png,
    rpm='ImageMagick',
)


ghostscript_png = ToolAdapter(
    'ghostscript pdf->png',
    ['gs', '-sDEVICE=png16m', '-dNOPAUSE', '-dBATCH', '-dSAFER',
     '-r{dpi}', '-sOutputFile={output}', '{input}'],
    'Ghostscript (http://www.ghostscript.com)',
    test_convert_to_png,
    rpm='ghostscript',
)


pdfcrop = ToolAdapter(
    'pdfcrop',
    ['pdfcrop', '{input}', '{output}'],
    'PDFCrop by by Heiko Oberdiek',
    test_output_is_pdf,
    rpm='texlive-pdfcrop',
    deb='texlive-extra-utils',
)


