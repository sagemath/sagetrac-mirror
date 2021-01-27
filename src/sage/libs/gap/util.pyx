"""Utility functions for GAP interface."""

#*****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import os

import sage.env


# For backwards-compatibility
from gappy.exceptions import GAPError


def gap_root():
    """
    Find the location of the GAP root install which is stored in the gap
    startup script.

    EXAMPLES::

        sage: from sage.libs.gap.util import gap_root
        sage: gap_root()   # random output
        '/home/vbraun/opt/sage-5.3.rc0/local/gap/latest'
    """
    if os.path.exists(sage.env.GAP_ROOT_DIR):
        return sage.env.GAP_ROOT_DIR

    # Attempt to figure out the appropriate GAP_ROOT by reading the
    # local/bin/gap shell script; this is an ugly hack that exists for
    # historical reasons; the best approach to setting where Sage looks for
    # the appropriate GAP_ROOT is to set the GAP_ROOT_DIR variable
    SAGE_LOCAL = sage.env.SAGE_LOCAL
    with open(os.path.join(SAGE_LOCAL, 'bin', 'gap')) as f:
        gap_sh = f.read().splitlines()
    gapdir = next(x for x in gap_sh if x.strip().startswith('GAP_ROOT'))
    gapdir = gapdir.split('"')[1]
    gapdir = gapdir.replace('$SAGE_LOCAL', SAGE_LOCAL)
    return gapdir
