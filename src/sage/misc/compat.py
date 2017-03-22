"""Cross-platform compatibility routines and wrappers."""

#*****************************************************************************
#       Copyright (C) 2017 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import subprocess
import sys

from sage.env import SAGE_LOCAL
from sage.misc.decorators import sage_wraps


#################################################################
# Replacements (as needed) for Python stdlib functions to provide
# better platform compatibility
#################################################################
from ctypes.util import find_library as _find_library
if sys.platform == 'cygwin':
    # find_library that works in cygwin adapted from
    # http://cygwin-ports.svn.sourceforge.net/viewvc/cygwin-ports/ports/trunk/lang/python/2.5.2-ctypes-util-find_library.patch?revision=8245&view=markup
    @sage_wraps(_find_library)
    def _find_library(name):
        for libdir in [os.path.join(SAGE_LOCAL, 'lib'),
                       '/usr/local/lib', '/usr/lib']:
            for libext in ['dll.a', 'a']:
                implib = os.path.join(libdir,
                                      'lib{0}.{1}'.format(name, libext))
                if not os.path.exists(implib):
                    continue

                cmd = ['dlltool', '-I', implib]

                p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE,
                                          universal_newlines=True)

                stdout, stderr = p.communicate()

                if p.returncode == 0:
                    return stdout.strip()


@sage_wraps(_find_library)
def find_library(name):
    """
    Returns the shared library filename for a given library.

    The library name is given without any prefixes or suffixes--(e.g.
    just "Singular", not "libSingular", as shared library naming is
    platform-specific.

    This does ''not'' currently return the absolute path of the file on most
    platforms; see https://bugs.python.org/issue21042

    EXAMPLES::

        sage: from sage.misc.compat import find_library
        sage: find_library('Singular')
        ...Singular...

    """

    return _find_library(name)
