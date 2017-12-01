# -*- encoding: utf-8 -*-

"""Internal module with C implementations of some context managers."""

from __future__ import absolute_import

import atexit

from cpython.ref cimport PyObject


__all__ = ['restore_atexit']


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


cdef class restore_atexit:
    cdef int _clear
    cdef list _exithandlers

    def __init__(self, clear=False):
        self._clear = clear
        self._exithandlers = []

    def __enter__(self):
        self.copy_state()
        if self._clear:
            atexit._clear()

        return self

    def __exit__(self, *exc):
        """
        Restore the original atexit module state by first clearing it, then
        re-registering all the saved callbacks.
        """

        atexit._clear()

        for callback in self._exithandlers:
            atexit.register(callback[0], *callback[1], **callback[2])

        del self._exithandlers[:]

    cdef copy_state(self):
        cdef atexitmodule_state* state
        cdef atexit_callback callback
        cdef int idx
        cdef object kwargs

        state = <atexitmodule_state*>PyModule_GetState(atexit)

        if not state:
            raise RuntimeError("atexit module state missing or corrupt")

        for idx in range(state.ncallbacks):
            callback = state.atexit_callbacks[idx][0]
            if callback.kwargs:
                kwargs = <object>callback.kwargs
            else:
                kwargs = {}
            self._exithandlers.append((<object>callback.func,
                                       <object>callback.args,
                                       kwargs))
