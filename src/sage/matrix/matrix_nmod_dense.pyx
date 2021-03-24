# distutils: extra_compile_args = -D_XPG6
"""
FLINT nmod_mat class wrapper

AUTHORS:

- Edgar Costa (2021) Initial version.
"""

from cpython.sequence cimport *

from cysignals.signals cimport sig_str, sig_off

from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport Element, Matrix
from sage.libs.flint.nmod_mat cimport *

from .args cimport SparseEntry, MatrixArgs_init

import sage.matrix.matrix_space as matrix_space
from sage.rings.finite_rings.integer_mod cimport IntegerMod_int, IntegerMod_int64

cdef class Matrix_nmod_dense(Matrix_dense):
    def __cinit__(self, parent, *args, **kwds):
        self._modulus = parent._base._pyx_order
        sig_str("FLINT exception")
        nmod_mat_init(self._matrix, self._nrows, self._ncols, self._modulus.int64)
        sig_off()


    def __dealloc__(self):
        nmod_mat_clear(self._matrix)

    def __init__(self, parent, entries=None, bint coerce=True):
        self._parent = parent # MatrixSpace over IntegerMod_int or IntegerMod_int64
        ma = MatrixArgs_init(parent, entries)

        cdef long z
        for t in ma.iter(coerce, True): #????
            se = <SparseEntry>t
            z = <long>se.entry
            nmod_mat_set_entry(self._matrix, se.i, se.j, z)
        nmod_mat_print_pretty(self._matrix)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x):
        if self._modulus.element_class() is IntegerMod_int:
            self.set_unsafe_si(i, j, (<IntegerMod_int>x).ivalue)
        else:
            self.set_unsafe_si(i, j, (<IntegerMod_int64>x).ivalue)

    cdef void set_unsafe_si(self, Py_ssize_t i, Py_ssize_t j, long value):
        nmod_mat_set_entry(self._matrix, i, j, value)

    cdef long get_unsafe_si(self, Py_ssize_t i, Py_ssize_t j):
        return nmod_mat_get_entry(self._matrix, i, j)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        # FIXME is this correct?
        return self._parent._base(self.get_unsafe_si(i, j))

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
        cdef Matrix_nmod_dense ans = Matrix_nmod_dense.__new__(Matrix_nmod_dense, P)
        return ans

    cpdef _add_(self, _right):
        cdef Matrix_nmod_dense right = _right
        cdef Matrix_nmod_dense M = self._new(self._nrows, self._ncols)
        nmod_mat_add(M._matrix, self._matrix, right._matrix)
        return M

    cpdef _sub_(self, _right):
        cdef Matrix_nmod_dense right = _right
        cdef Matrix_nmod_dense M = self._new(self._nrows, self._ncols)
        nmod_mat_sub(M._matrix, self._matrix, right._matrix)
        return M

    cdef Matrix _matrix_times_matrix_(left, Matrix _right):
        if left._ncols != _right._nrows:
            raise IndexError("Number of columns of self must equal number of rows of right.")
        cdef Matrix_nmod_dense right = _right
        cdef Matrix_nmod_dense M = left._new(left._nrows, right._ncols)
        nmod_mat_mul(M._matrix, left._matrix, right._matrix)
        return M

    cpdef _lmul_(self, Element right):
        cdef Matrix_nmod_dense M = self._new(self._nrows, self._ncols)
        print("Hello")
        if self._modulus.element_class() is IntegerMod_int:
            nmod_mat_scalar_mul(M._matrix, self._matrix, (<IntegerMod_int?>right).ivalue)
        else:
            nmod_mat_scalar_mul(M._matrix, self._matrix, (<IntegerMod_int64?>right).ivalue)
        return M

    def strong_echelon_form(self):
        if self._nrows >= self._ncols:
            nmod_mat_strong_echelon_form(self._matrix)
        else:
            raise ValueError("Matrix must have at least as many rows as columns.")


    def howell_form(self):
        if self._nrows >= self._ncols:
            nmod_mat_howell_form(self._matrix)
        else:
            raise ValueError("Matrix must have at least as many rows as columns.")


