# -*- encoding: utf-8 -*-
r"""
Graphics File Formats
"""
#*****************************************************************************
#       Copyright (C) 2016 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function

import os
import collections


GraphicsFormat = collections.namedtuple('GraphicsFormat', ['ext', 'matplotlib'])

EPS  = GraphicsFormat('.eps', matplotlib='eps')
PDF  = GraphicsFormat('.pdf', matplotlib='pdf')
PGF  = GraphicsFormat('.pgf', matplotlib=None)
PNG  = GraphicsFormat('.png', matplotlib='png')
PS   = GraphicsFormat('.ps',  matplotlib='ps')
SOBJ = GraphicsFormat('.sobj', matplotlib=None)
SVG  = GraphicsFormat('.svg', matplotlib='svg')

ALL_FORMATS = [
    EPS, PDF, PGF, PNG, PS, SOBJ, SVG
]

BY_EXTENSION = dict([fmt.ext, fmt] for fmt in ALL_FORMATS)

ALLOWED_EXTENSIONS = sorted([fmt.ext for fmt in ALL_FORMATS])



def graphics_format(filename):
    """
    Figure out graphics format from filename

    INPUT:

    - ``filename`` -- string or file-like object. The filename to use.

    OUTPUT: 

    :class:`GraphicsFormat` instance.

    EXAMPLES::

        sage: from sage.plot.graphics_format import graphics_format
        sage: graphics_format('/tmp/foo.svg')
        GraphicsFormat(ext='.svg', matplotlib='svg')
        sage: graphics_format('/tmp/foo.bar')
        Traceback (most recent call last):
        ...
        ValueError: allowed file extensions for images are ...

    Fallback for missing extension in string or no ``name`` attribute
    on file-like objects is sobj::

        sage: graphics_format('/tmp/foo')
        GraphicsFormat(ext='.sobj', matplotlib=None)
        sage: import io
        sage: graphics_format(io.BytesIO())
        GraphicsFormat(ext='.sobj', matplotlib=None)
    """
    if hasattr(filename, 'write'):
        # is a file-like object
        name = getattr(filename, 'name', '')
    else:
        # is a string
        name = filename
    ext = os.path.splitext(name)[1].lower()
    if not ext:
        ext = '.sobj'
    try:
        return BY_EXTENSION[ext]
    except KeyError:
        raise ValueError("allowed file extensions for images are '"
                         + "', '".join(ALLOWED_EXTENSIONS) + "'!")
    
