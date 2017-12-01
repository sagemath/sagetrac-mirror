# -*- encoding: utf-8 -*-

"""A home for miscellaneous context managers."""

#*****************************************************************************
#       Copyright (C) 2017 Erik M. Bray <erik.bray@lri.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import


import atexit

import six


__all__ = ['restore_atexit']


if six.PY2:
    # These are helper functions for implementing restore_atexit.  The only
    # reason they're broken out from the class itself is because these require
    # very different implements on Python 2 and 3
    def _get_exithandlers():
        """Return list of exit handlers registered with the atexit module."""
        return atexit._exithandlers[:]

    def _set_exithandlers(exithandlers):
        """
        Replace the list of exit handlers registered with the atexit module
        with a new list.
        """
        atexit._exithandlers[:] = exithandlers

    def _clear_exithandlers():
        """Clear the atexit module of all registered exit handlers."""
        # It generally shouldn't matter, but this keeps the same list
        # object in place rather than replacing it with a new one in
        # case any other code is accessing this list directly.
        del atexit._exithandlers[:]
else:
    # Python 3 implementations of the same helper functions
    from ._atexit_py3 import (_get_exithandlers, _set_exithandlers,
                              _clear_exithandlers)


# Wrapper class for the two restore_atexit implementations so that they
# can share a docstring
class restore_atexit(object):
    """
    Context manager that restores the state of the atexit module to its
    previous state when exiting the context.

    INPUT:

    - ``clear`` (`bool`, default `False`) -- if `True`, clear already
      registered atexit handlers upon entering the context.

    EXAMPLES::

    For this example we will wrap the entire example with
    ``restore_atexit(True)`` so as to start with a fresh atexit module state
    for the sake of the example.

    Note that the function ``atexit._run_exitfuncs()`` runs all registered
    handlers, and then clears the list of handlers, so we can use it to test
    manipulation of the ``atexit`` state.

    sage: import atexit
    sage: from sage.cpython.atexit import restore_atexit
    sage: def handler(*args, **kwargs):
    ....:     print((args, kwargs))
    sage: with restore_atexit(clear=True):
    ....:     atexit._run_exitfuncs()  # Should be none registered
    ....:     atexit.register(handler, 1, 2, c=3)
    ....:     with restore_atexit():
    ....:         atexit._run_exitfuncs()  # Run just registered handler
    ....:     atexit._run_exitfuncs()  # Handler should be run again
    <function handler at 0x...>
    ((1, 2), {'c': 3})
    ((1, 2), {'c': 3})
    """

    def __init__(self, clear=False):
        self._clear = clear
        self._exithandlers = None

    def __enter__(self):
        self._exithandlers = _get_exithandlers()
        if self._clear:
            _clear_exithandlers()

        return self

    def __exit__(self, *exc):
        _set_exithandlers(self._exithandlers)
