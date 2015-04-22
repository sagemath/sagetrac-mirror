# -*- encoding: utf-8 -*-
r"""
Guess the rich output type given a file.
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
import imghdr

import sage.repl.rich_output.output_catalog as catalog



def guess_output(filename):
    """
    Guess a rich output container for a file
    
    Uses a combination of file magics and extension to guess a
    suitable rich output type. Clearly this can only work for rich
    output file formats that consist of a single file.

    INPUT:

    - ``filename`` -- string. The name of a file. The file must exist.

    OUTPUT:

    A subclass of
    :class:`~sage.repl.rich_output.output_base.OutputBase` or ``None``
    if no suitable output type can be identified.

    EXAMPLES::

        sage: from sage.repl.rich_output.guess_output import *
        sage: pdf_file = catalog.OutputImagePdf.example().pdf.filename()
        sage: guess_output(pdf_file)
    """
    # Use file magics
    file_type = imghdr.what(filename)
    if file_type == 'gif':
        return catalog.OutputImageGif
    elif file_type == 'png':
        return catalog.OutputImagePng
    elif file_type == 'jpeg':
        return catalog.OutputImageJpg
    # Otherwise try file extension
    base, ext = os.path.splitext(filename.lower())
    if ext == 'gif':
        return catalog.OutputImageGif
    elif ext == 'png':
        return catalog.OutputImagePng
    elif ext in ['jpg', 'jpeg']:
        return catalog.OutputImageJpg
    elif ext == 'svg':
        return catalog.OutputImageSvg
    elif ext == 'pdf':
        return catalog.OutputImagePdf
    elif ext == 'dvi':
        return catalog.OutputImageDvi
    # Give up
    return None
