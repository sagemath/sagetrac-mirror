# distutils: extra_compile_args = -D_XPG6
"""
FLINT nmod_mat class wrapper

AUTHORS:

- Edgar Costa (2021) Initial version.
"""

from cpython.sequence cimport *

from cysignals.signals cimport sig_str, sig_on, sig_off

from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport Element, Matrix
from sage.structure.richcmp cimport rich_to_bool
from sage.libs.flint.nmod_mat cimport *
from sage.libs.flint.nmod_poly cimport nmod_poly_set
from sage.libs.flint.ulong_extras cimport n_CRT
from sage.rings.polynomial.polynomial_zmod_flint cimport Polynomial_zmod_flint

from .args cimport SparseEntry, MatrixArgs_init

import sage.matrix.matrix_space as matrix_space
from sage.rings.finite_rings.integer_mod cimport IntegerMod_abstract, IntegerMod_int, IntegerMod_int64

cdef class Matrix_nmod_dense(Matrix_dense):
    ########################################################################
    # LEVEL 1 helpers:
    #   These function support the implementation of the level 1 functionality.
    ########################################################################
    def __cinit__(self, parent, *args, **kwds):
        self._modulus = parent._base._pyx_order
        sig_str("FLINT exception")
        nmod_mat_init(self._matrix, self._nrows, self._ncols, self._modulus.int64)
        sig_off()



    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        self._parent = parent # MatrixSpace over IntegerMod_int or IntegerMod_int64
        ma = MatrixArgs_init(parent, entries)

        cdef long z
        for t in ma.iter(coerce, True): #????
            se = <SparseEntry>t
            z = <long>se.entry
            nmod_mat_set_entry(self._matrix, se.i, se.j, z)
        nmod_mat_print_pretty(self._matrix)

    def __dealloc__(self):
        nmod_mat_clear(self._matrix)

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
        cdef type t = self._modulus.element_class()
        cdef IntegerMod_abstract x = t.__new__(t)
        x._parent = self._parent._base
        x.__modulus = self._modulus
        x.set_from_ulong_fast(nmod_mat_get_entry(self._matrix, i, j))
        return x

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

    ########################################################################
    # LEVEL 2 helpers:
    #   These function support the implementation of the level 2 functionality.
    ########################################################################

    cpdef _add_(self, _right):
        cdef Matrix_nmod_dense right = _right
        cdef Matrix_nmod_dense M = self._new(self._nrows, self._ncols)
        nmod_mat_add(M._matrix, self._matrix, right._matrix)
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
        if self._modulus.element_class() is IntegerMod_int:
            nmod_mat_scalar_mul(M._matrix, self._matrix, (<IntegerMod_int?>right).ivalue)
        else:
            nmod_mat_scalar_mul(M._matrix, self._matrix, (<IntegerMod_int64?>right).ivalue)
        return M


    def __copy__(self):
        """
        Return a copy of this matrix. Changing the entries of the copy will
        not change the entries of this matrix.
        """
        cdef Matrix_nmod_dense M = self._new(self._nrows, self._ncols)
        sig_on()
        nmod_mat_set(M._matrix, self._matrix)
        sig_off()
        return M

    def __neg__(self):
        r"""
        Return the negative of this matrix.

        TESTS::


        """
        cdef Matrix_nmod_dense M = self._new(self._nrows, self._ncols)
        sig_on()
        nmod_mat_neg(M._matrix, self._matrix)
        sig_off()
        return M

    cpdef _richcmp_(self, right, int op):
        r"""
        Compare ``self`` with ``right``, examining entries in
        lexicographic (row major) ordering.

        EXAMPLES::


        """
        cdef Py_ssize_t i, j
        cdef int k

        sig_on()
        for i in range(self._nrows):
            for j in range(self._ncols):
                if nmod_mat_get_entry(self._matrix,i,j) != nmod_mat_get_entry((<Matrix_nmod_dense>right)._matrix,i,j):
                    sig_off()
                    return rich_to_bool(op, 1)
        sig_off()
        return rich_to_bool(op, 0)


    ########################################################################
    # LEVEL 3 helpers:
    #   These function support the implementation of the level 2 functionality.
    ########################################################################

    def __nonzero__(self):
        return not nmod_mat_is_zero(self._matrix)

    cpdef _sub_(self, _right):
        cdef Matrix_nmod_dense right = _right
        cdef Matrix_nmod_dense M = self._new(self._nrows, self._ncols)
        nmod_mat_sub(M._matrix, self._matrix, right._matrix)
        return M

    def __invert__(self):
        r"""
        Return the inverse of this matrix.

        Raises a ``ZeroDivisionError`` if the determinant is not a unit,
        and raises an ``ArithmeticError`` if the
        inverse doesn't exist because the matrix is nonsquare.

        EXAMPLES::

            sage: 
        """
        return self.inverse_of_unit()

    def inverse_of_unit(self, algorithm="crt"):
        if not self.is_square():
            raise ArithmeticError("inverse only defined for square matrix")
        if not self.nrows():
            return self

        cdef Matrix_nmod_dense M
        R = self._parent._base
        cdef Py_ssize_t i, j, n = self._nrows
        cdef long k, e, b
        cdef mp_limb_t p, q, N = 1
        cdef nmod_mat_t accum, A, B, combo
        cdef bint lift_required, crt_required, ok
        if R.is_field():
            M = self._new(n, n)
            ok = nmod_mat_inv(M._matrix, self._matrix)
            if not ok:
                raise ZeroDivisionError("input matrix must be nonsingular")
            return M
        else:
            if algorithm == "lift":
                return (~self.lift_centered()).change_ring(R)
            elif algorithm == "crt":
                # The modulus is small, so factoring is feasible.
                # We find inverses modulo each prime dividing the modulus, lift p-adically, then CRT the results together
                M = self._new(n, n)
                F = R.factored_order()
                lift_required = any(ez > 1 for pz, ez in F)
                crt_required = len(F) > 1
                try:
                    nmod_mat_init(A, n, n, 1)
                    nmod_mat_init(B, n, n, 1)
                    if lift_required:
                        nmod_mat_init(combo, n, n, 1)
                    if crt_required:
                        nmod_mat_init(accum, n, n, 1)
                    for pz, ez in F:
                        q = p = pz
                        e = ez
                        _nmod_mat_set_mod(A, p)
                        _nmod_mat_set_mod(B, p)
                        for i in range(n):
                            for j in range(n):
                                nmod_mat_set_entry(A, i, j, nmod_mat_get_entry(self._matrix, i, j) % p)
                        ok = nmod_mat_inv(B, A)
                        if not ok:
                            raise ZeroDivisionError("input matrix must be nonsingular")
                        # Since the modulus fits in a word, we don't worry about optimizing intermediate precision to keep it low, since the cost of operations is constant
                        # So, for example, if e=10 we do 1, 2, 4, 8, 10 rather than 1, 2, 3, 5, 10.
                        k = 1
                        while k < e:
                            k = min(2*k, e)
                            q = p**k
                            _nmod_mat_set_mod(combo, q)
                            _nmod_mat_set_mod(B, q)
                            _nmod_mat_set_mod(A, q)
                            for i in range(n):
                                for j in range(n):
                                    nmod_mat_set_entry(A, i, j, nmod_mat_get_entry(self._matrix, i, j) % q)
                            # Suppose B is an approximation to the inverse of A.
                            # I - AB = p^k E
                            # (I - AB)(I - AB) = I - A(2B - BAB) = p^(2k)E^2
                            # So 2B - BAB is a better approximation
                            nmod_mat_mul(combo, A, B)
                            nmod_mat_mul(combo, B, combo)
                            nmod_mat_scalar_mul(B, B, 2)
                            nmod_mat_sub(B, B, combo)
                        if not crt_required:
                            nmod_mat_set(M._matrix, B)
                            return M
                        _nmod_mat_set_mod(accum, N*q)
                        for i in range(n):
                            for j in range(n):
                                nmod_mat_set_entry(accum, i, j, n_CRT(nmod_mat_get_entry(accum, i, j), N, nmod_mat_get_entry(B, i, j), q))
                        N *= q
                    nmod_mat_set(M._matrix, accum)
                    return M
                finally:
                    nmod_mat_clear(A)
                    nmod_mat_clear(B)
                    if lift_required:
                        nmod_mat_clear(combo)
                    if crt_required:
                        nmod_mat_clear(accum)

    cdef swap_columns_c(self, Py_ssize_t c1, Py_ssize_t c2):
        nmod_mat_swap_cols(self._matrix, NULL, c1, c2)

    cdef swap_rows_c(self, Py_ssize_t r1, Py_ssize_t r2):
        nmod_mat_swap_rows(self._matrix, NULL, r1, r2)

    # Extra

    def strong_echelon_form(self):
        # Edgar: fix so not in place
        if self._nrows >= self._ncols:
            nmod_mat_strong_echelon_form(self._matrix)
        else:
            raise ValueError("Matrix must have at least as many rows as columns.")


    def howell_form(self):
        # Edgar: fix so not in place
        if self._nrows >= self._ncols:
            nmod_mat_howell_form(self._matrix)
        else:
            raise ValueError("Matrix must have at least as many rows as columns.")

    # random matrix generation (David)
    # charpoly and minpoly (David)
    # solve (nmod_mat_can_solve) (David)



    # transpose (Edgar)
    # nmod_mat_pow (Edgar)
    # nmod_mat_trace (Edgar)
    # rank and det (only primes) (Egar)
    # right_kernel_matrix (nmod_mat_nullspace) (Edgar)
    # row reduction (nmod_mat_rref) (Edgar)
