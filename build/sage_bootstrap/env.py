# -*- coding: utf-8 -*-
"""
Environment Variables

This module defines the following subset of the Sage environment
variables:

* ``SAGE_ROOT``
* ``SAGE_VPATH``
* ``SAGE_SRC``
* ``SAGE_DISTFILES``
"""


#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os


try:
    SAGE_ROOT = os.environ['SAGE_ROOT']
except KeyError:
    SAGE_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))

try:
    SAGE_VPATH = os.environ['SAGE_VPATH']
except KeyError:
    SAGE_VPATH = SAGE_ROOT

SAGE_SRC = os.environ.get('SAGE_SRC',
    os.path.join(SAGE_VPATH, 'src'))
SAGE_DISTFILES = os.environ.get('SAGE_DISTFILES',
    os.path.join(SAGE_ROOT, 'upstream'))

assert os.path.isfile(os.path.join(SAGE_VPATH, 'configure.ac')), SAGE_VPATH
# Check that SAGE_ROOT is the root of either a configured build directory
# or the root of an unconfigured non-VPATH source directory.
assert (os.path.isfile(os.path.join(SAGE_ROOT, 'config.status'))
        or os.path.isfile(os.path.join(SAGE_ROOT, 'configure.ac'))), SAGE_ROOT

try:
    # SAGE_DISTFILES does not exist in a fresh git clone
    os.mkdir(SAGE_DISTFILES)
except OSError:
    pass

assert os.path.isdir(SAGE_DISTFILES)
