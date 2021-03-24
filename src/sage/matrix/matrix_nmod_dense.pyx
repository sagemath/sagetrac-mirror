# distutils: extra_compile_args = -D_XPG6
"""
FLINT nmod_mat class wrapper

AUTHORS:

- Edgar Costa (2021) Initial version.
"""

from cpython.sequence cimport *

from cysignals.signals cimport sig_str, sig_off

from sage.structure.sage_object cimport SageObject
from sage.libs.flint.nmod_mat cimport *

from .args cimport SparseEntry, MatrixArgs_init

import sage.matrix.matrix_space as matrix_space


cdef class Matrix_nmod_dense(Matrix_dense):
    def __cinit__(self):
        sig_str("FLINT exception")
        nmod_mat_init(self._matrix, self._nrows, self._ncols, self._modulus.int64)
        sig_off()


    def __dealloc__(self):
        nmod_mat_clear(self._matrix)

    def __init__(self, parent, entries=None, bint coerce=True):
        self._parent = parent # MatrixSpace over IntegerMod_int or IntegerMod_int64
        self._modulus = parent._base._pyx_order
        ma = MatrixArgs_init(parent, entries)

        cdef long z
        for t in ma.iter(coerce, True): #????
            se = <SparseEntry>t
            z = <long>se.entry
            nmod_mat_set_entry(self._matrix, se.i, se.j, z.ivalue)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x):
        self.set_unsafe_si(i, j, x.ivalue)

    cdef void set_unsafe_si(self, Py_ssize_t i, Py_ssize_t j, long value):
        nmod_mat_set_entry(self._matrix, i, j, value)

    cdef long get_unsafe_si(self, Py_ssize_t i, Py_ssize_t j):
        return nmod_mat_get_entry(self._matrix, i, j)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        return self._parent._base._new_c(self.get_unsafe_si(i, j))

    ########################################################################
    # LEVEL 1 helpers:
    #   These function support the implementation of the level 1 functionality.
    ########################################################################
    cdef Matrix_nmod_dense _new(self, Py_ssize_t nrows, Py_ssize_t ncols):
        """
        Return a new matrix over the parent from given parent
        All memory is allocated for this matrix, but its
        entries have not yet been filled in.
        """
        if nrows == self._nrows and ncols == self._ncols:
            P = self._parent
        else:
            P = matrix_space.MatrixSpace(self._parent._base, nrows, ncols, sparse=False)
        cdef Matrix_nmod_dense ans = Matrix_nmod_dense.__new__(Matrix_nmod_dense, P, None, None, None)
        return ans

    cdef void strong_echelon_form(self):
        nmod_mat_strong_echelon_form(self._matrix)

    cdef void howell_form(self):
        nmod_mat_howell_form(self._matrix)


