"""
Management of the .sage directory
"""
# ****************************************************************************
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import os
import platform
from pathlib import Path
from functools import lru_cache


@lru_cache
def dot_sage():
    """
    Return the ".sage" directory (and create it if it does not exist yet).

    This is done with restrictive permissions, since otherwise
    possibly just anybody can easily see every command you type.

    EXAMPLES::

        sage: from sage.misc.dot_sage import dot_sage
        sage: d = dot_sage()
        sage: d.is_dir()
        True
    """
    d = os.environ.get('DOT_SAGE')
    if d is not None:
        path = Path(d)
    else:
        path = Path(os.environ.get("HOME"))
        if ' ' in str(path):
            if platform.system().lower().startswith('cygwin'):
                # on windows/cygwin it is typical for the home directory
                # to have a space in it. Fortunately, users also have
                # write privileges to c:\cygwin\home, so we just put
                # .sage there.
                path = Path("/home")
            else:
                print("Your home directory name has a space in it."
                      "This will probably break some functionality of Sage."
                      "A workaround is to set the environment "
                      "variable DOT_SAGE to a directory with no spaces "
                      "that you have write permissions to before you "
                      "start sage."
        path = path / ".sage"

    path.mkdir(mode=700, exist_ok=True)
    return path
