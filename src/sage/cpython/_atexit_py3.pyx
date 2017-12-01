# -*- encoding: utf-8 -*-

"""Internal module with C implementations of some context managers."""

from __future__ import absolute_import

import atexit

from cpython.ref cimport PyObject


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


def _set_exithandlers(exithandlers):
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
