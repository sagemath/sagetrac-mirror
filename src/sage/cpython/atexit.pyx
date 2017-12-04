# -*- encoding: utf-8 -*-

"""Utilities for interfacing with the standard library's atexit module."""

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


__all__ = ['restore_atexit']


cdef class restore_atexit(object):
    """
    Context manager that restores the state of the atexit module to its
    previous state when exiting the context.

    INPUT:

    - ``clear`` (`bool`, default `False`) -- if `True`, clear already
      registered atexit handlers upon entering the context.

    EXAMPLES:

    For this example we will wrap the entire example with
    ``restore_atexit(True)`` so as to start with a fresh atexit module state
    for the sake of the example.

    Note that the function ``atexit._run_exitfuncs()`` runs all registered
    handlers, and then clears the list of handlers, so we can use it to test
    manipulation of the ``atexit`` state::

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

    cdef list _exithandlers
    cdef int _clear

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


# These are helper functions for implementing restore_atexit.  The only reason
# they're broken out from the class itself is because these require very
# different implementations on Python 2 and 3
IF PY_MAJOR_VERSION == 2:
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
ELSE:
    from cpython.ref cimport PyObject

    # Internal structures defined in the CPython source in
    # Modules/atexitmodule.c and subject to (but unlikely to) change.  Watch
    # https://bugs.python.org/issue32082 for a request to (eventually)
    # re-expose more of the atexit module's internals to Python
    ctypedef struct atexit_callback:
        PyObject* func
        PyObject* args
        PyObject* kwargs


    ctypedef struct atexitmodule_state:
        atexit_callback** atexit_callbacks
        int ncallbacks
        int callback_len


    cdef extern from "Python.h":
        void* PyModule_GetState(object module)


    def _get_exithandlers():
        """Return list of exit handlers registered with the atexit module."""
        cdef atexitmodule_state* state
        cdef atexit_callback callback
        cdef list exithandlers
        cdef int idx
        cdef object kwargs

        state = <atexitmodule_state*>PyModule_GetState(atexit)

        if not state:
            raise RuntimeError("atexit module state missing or corrupt")

        exithandlers = []

        for idx in range(state.ncallbacks):
            callback = state.atexit_callbacks[idx][0]
            if callback.kwargs:
                kwargs = <object>callback.kwargs
            else:
                kwargs = {}
            exithandlers.append((<object>callback.func,
                                 <object>callback.args,
                                 kwargs))
        return exithandlers


    def _set_exithandlers(list exithandlers):
        """
        Replace the list of exit handlers registered with the atexit module
        with a new list.
        """

        # We could do this more efficiently by directly rebuilding the array
        # of atexit_callbacks, but this is much simpler
        for callback in exithandlers:
            atexit.register(callback[0], *callback[1], **callback[2])


    def _clear_exithandlers():
        """Clear the atexit module of all registered exit handlers."""
        atexit._clear()
