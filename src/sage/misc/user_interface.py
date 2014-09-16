r"""
User interface

AUTHORS:

- Thierry Monteil

The function :func:`user_interface` tells which user interface is currently
used.


"""

#*****************************************************************************
#       Copyright (C) 2014 Thierry Monteil <sage!lma.metelu.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.server.support import EMBEDDED_MODE
import __builtin__

def user_interface():
    r"""
    Tell which user interface is currently used.

    OUTPUT:

    - a string, among the following:

        - ``sage_command_line``
        - ``ipython_command_line``
        - ``sage_notebook``
        - ``ipython_notebook``
        - ``other``

    Sage can currently be used through various user interfaces, which do not
    handle output formatting the same way (e.g. HTML, LaTeX).

    This function is mainly a transitional tool. While this function could be
    useful to let Sage behave differently depending in which user interface is
    used, you should refrain using it for ipython notebook specific rendering.
    Instead, you should try the ipython way and define ``_repr_*_()`` methods,
    as explained in http://ipython.org/ipython-doc/dev/config/integrating.html.

    EXAMPLES::

        sage: from sage.misc.user_interface import user_interface
        sage: user_interface()      # depends on the current user interface
        'other'
    """
    ipython = hasattr(__builtin__, 'get_ipython')
    if ipython:
        get_ipython = str(__builtin__.get_ipython())
        if 'ZMQInteractiveShell' in get_ipython:
            return 'ipython_notebook'
        elif 'SageInteractiveShell' in get_ipython:
            return 'sage_command_line'
        elif 'TerminalInteractiveShell' in get_ipython:
            return 'ipython_command_line'
    elif EMBEDDED_MODE:
        return 'sage_notebook'
    return 'other'


