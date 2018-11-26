from cpython.exc cimport PyErr_Clear
from cpython.sequence cimport PySequence_Check


cdef extern from "Python.h":
    # Note: The declaration of PySequence_Length from cpython.sequence raises
    # an exception, which we don't want
    Py_ssize_t PySequence_Length(object o)


cpdef inline bint issequence(obj):
    """
    Check whether the given object is sequence-like.

    Here "sequence-like" is defined more closely to the Sequence
    ABC (abstract base class) than the C-level ``PySequence`` API
    in that the object's type must define ``__getitem__`` and
    ``__len__``, and should be possible to iterate over
    monotonically from ``0`` to ``len(obj)`` like
    ``for idx in range(len(obj)): obj[idx]``.

    Note: For custom types the ``__getitem__`` implementation *must* perform
    its own bounds checking, and raise `IndexError` on out of bounds indices.
    Otherwise, Python may go into an infinite loop or run out of memory, as
    it does not use ``__len__`` to perform bounds checking on sequences.

    EXAMPLES::

        sage: from sage.cpython.sequence import issequence
        sage: import sys

        sage: issequence([])
        True
        sage: issequence(())
        True
        sage: issequence('abc')
        True
        sage: issequence(range(10))
        True
        sage: issequence(1)
        False
        sage: issequence(vector([1, 2, 3]))
        True

    Just implementing ``__getitem__`` does *not*  make a sequence for the
    purposes of this function::

        sage: class MyNonSequence(object):
        ....:     def __getitem__(self, x): return x
        sage: issequence(MyNonSequence())
        False

    A sequence must implement ``__len__`` and must implement proper bounds
    checking in its ``__getitem__``::

        sage: class MySequence(object):
        ....:     def __len__(self): return 10
        ....:     def __getitem__(self, x):
        ....:         # You may also implement negative indices here but that
        ....:         # is omitted for simplicity
        ....:         if x >= 0 and x < len(self):
        ....:             return x
        ....:         raise IndexError(x)
        ....:
        sage: s = MySequence()
        sage: issequence(s)
        True
        sage: list(s)
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    """

    cdef bint ret = PySequence_Check(obj) and PySequence_Length(obj) >= 0
    PyErr_Clear()
    return ret
