r"""
Arbitrary precision complex ball matrices using Arb

AUTHORS:

- Clemens Heuberger (2014-10-25): Initial version.

This is a rudimentary binding to the `Arb library
<http://arblib.org>`_; it may be useful to refer to its
documentation for more details.

TESTS::

    sage: mat = matrix(CBF, 2, 2, range(4))
    sage: x = polygen(QQ)
    sage: pol = x^3 + 2
    sage: pol(mat)
    [8.000000000000000 11.00000000000000]
    [22.00000000000000 41.00000000000000]
"""

#*****************************************************************************
#       Copyright (C) 2014 Clemens Heuberger <clemens.heuberger@aau.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.object cimport Py_EQ, Py_NE
from cysignals.signals cimport sig_on, sig_str, sig_off

from sage.arith.power cimport generic_power_pos
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_mat cimport *
from sage.libs.gmp.mpz cimport mpz_fits_ulong_p, mpz_get_ui
from sage.matrix.constructor import matrix
from sage.matrix.matrix_generic_sparse cimport Matrix_generic_sparse
from .args cimport SparseEntry, MatrixArgs_init
from sage.rings.complex_interval_field import ComplexIntervalField_class, ComplexIntervalField
from sage.rings.complex_interval cimport ComplexIntervalFieldElement
from sage.rings.complex_arb cimport (
    ComplexBall,
    ComplexIntervalFieldElement_to_acb,
    acb_to_ComplexIntervalFieldElement)
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_complex_arb cimport Polynomial_complex_arb
from sage.structure.element cimport Element, RingElement, Matrix
from sage.structure.parent cimport Parent

from sage.rings.integer_ring import ZZ
from sage.rings.polynomial import polynomial_ring_constructor


cdef void matrix_to_acb_mat(acb_mat_t target, source):
    """
    Convert a matrix containing :class:`ComplexIntervalFieldElement` to an ``acb_mat_t``.

    INPUT:

    - ``target`` -- an ``acb_mat_t``

    - ``source`` -- a matrix consisting of :class:`ComplexIntervalFieldElement`

    OUTPUT:

    None.
    """
    cdef unsigned long nrows, ncols, r, c, precision

    nrows = acb_mat_nrows(target)
    ncols = acb_mat_ncols(target)

    for r in range(nrows):
        for c in range(ncols):
            ComplexIntervalFieldElement_to_acb(acb_mat_entry(target, r, c),
                                               source[r][c])

cdef ComplexIntervalFieldElement _to_CIF(acb_t source, ComplexIntervalFieldElement template):
    cdef ComplexIntervalFieldElement result
    result = template._new()
    acb_to_ComplexIntervalFieldElement(
        result, source)
    return result

cdef Matrix_generic_dense acb_mat_to_matrix(
    acb_mat_t source, Parent CIF):
    """
    Convert an ``acb_mat_t`` to a matrix containing :class:`ComplexIntervalFieldElement`.

    INPUT:

    - ``source`` -- an ``acb_mat_t``

    - ``precision`` -- a positive integer.

    OUTPUT:

    A :class:`~sage.matrix.matrix_generic_dense.Matrix_generic_dense`
    containing :class:`ComplexIntervalFieldElement`.
    """
    cdef unsigned long nrows, ncols, r, c
    cdef ComplexIntervalFieldElement template

    nrows = acb_mat_nrows(source)
    ncols = acb_mat_ncols(source)
    template = CIF(0)

    return matrix(
                  [[_to_CIF(acb_mat_entry(source, r, c), template)
                    for c in range(ncols)]
                   for r in range(nrows)])

cdef inline long prec(Matrix_complex_ball_dense mat):
    return mat._base_ring._prec

cdef class Matrix_complex_ball_dense(Matrix_dense):
    """
    Matrix over a complex ball field. Implemented using the
    ``acb_mat`` type of the Arb library.

    EXAMPLES::

        sage: MatrixSpace(CBF, 3)(2)
        [2.000000000000000                 0                 0]
        [                0 2.000000000000000                 0]
        [                0                 0 2.000000000000000]
        sage: matrix(CBF, 1, 3, [1, 2, -3])
        [ 1.000000000000000  2.000000000000000 -3.000000000000000]
    """
    def __cinit__(self,
                  parent,
                  entries,
                  coerce,
                  copy):
        """
        Create and allocate memory for the matrix.

        INPUT:

        -  ``parent, entries, coerce, copy`` - as for
           ``__init__``.

        EXAMPLES::

            sage: from sage.matrix.matrix_complex_ball_dense import Matrix_complex_ball_dense
            sage: a = Matrix_complex_ball_dense.__new__( # indirect doctest
            ....:     Matrix_complex_ball_dense, Mat(CBF, 2), 0, 0, 0)
            sage: type(a)
            <type 'sage.matrix.matrix_complex_ball_dense.Matrix_complex_ball_dense'>
        """
        self._parent = parent
        self._base_ring = parent.base_ring()
        self._nrows = parent.nrows()
        self._ncols = parent.ncols()
        sig_str("Arb exception")
        acb_mat_init(self.value, self._nrows, self._ncols)
        sig_off()

    def __dealloc__(self):
        """
        Free all the memory allocated for this matrix.

        EXAMPLES::

            sage: a = Matrix(CBF, 2, [1, 2, 3, 4]) # indirect doctest
            sage: del a
        """
        acb_mat_clear(self.value)

    cdef Matrix_complex_ball_dense _new(self, Py_ssize_t nrows, Py_ssize_t ncols):
        r"""
        Return a new matrix over the same base ring.
        """
        cdef Parent P
        if nrows == self._nrows and ncols == self._ncols:
            P = self._parent
        else:
            P = self.matrix_space(nrows, ncols)
        return Matrix_complex_ball_dense.__new__(Matrix_complex_ball_dense, P, None, None, None)

    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        Initialize a dense matrix over the complex ball field.

        INPUT:

        - ``parent`` -- a matrix space over a complex ball field

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if False, assume without checking that the
          entries lie in the base ring

        EXAMPLES:

        The ``__init__`` function is called implicitly in each of the
        examples below to actually fill in the values of the matrix.

        We create a `2 \times 2` and a `1\times 4` matrix::

            sage: matrix(CBF, 2, 2, range(4))
            [                0 1.000000000000000]
            [2.000000000000000 3.000000000000000]
            sage: Matrix(CBF, 1, 4, range(4))
            [                0 1.000000000000000 2.000000000000000 3.000000000000000]

        If the number of columns isn't given, it is determined from the
        number of elements in the list. ::

            sage: matrix(CBF, 2, range(4))
            [                0 1.000000000000000]
            [2.000000000000000 3.000000000000000]
            sage: matrix(CBF, 2, range(6))
            [                0 1.000000000000000 2.000000000000000]
            [3.000000000000000 4.000000000000000 5.000000000000000]

        Another way to make a matrix is to create the space of matrices and
        convert lists into it. ::

            sage: A = Mat(CBF, 2); A
            Full MatrixSpace of 2 by 2 dense matrices over
            Complex ball field with 53 bits of precision
            sage: A(range(4))
            [                0 1.000000000000000]
            [2.000000000000000 3.000000000000000]

        Actually it is only necessary that the input can be converted to a
        list, so the following also works::

            sage: v = reversed(range(4)); type(v)
            <...iterator'>
            sage: A(v)
            [3.000000000000000 2.000000000000000]
            [1.000000000000000                 0]

        Matrices can have many rows or columns (in fact, on a 64-bit
        machine they could have up to `2^{63}-1` rows or columns)::

            sage: v = matrix(CBF, 1, 10^5, range(10^5))
            sage: v.parent()
            Full MatrixSpace of 1 by 100000 dense matrices over
            Complex ball field with 53 bits of precision

        TESTS::

            sage: MatrixSpace(CBF, 0, 0).one()
            []
            sage: Matrix(CBF, 0, 100)
            0 x 100 dense matrix over Complex ball field with 53 bits
            of precision (use the '.str()' method to see the entries)
            sage: Matrix(CBF, 100, 0)
            100 x 0 dense matrix over Complex ball field with 53 bits
            of precision (use the '.str()' method to see the entries)
        """
        ma = MatrixArgs_init(parent, entries)
        cdef ComplexBall z
        for t in ma.iter(coerce, True):
            se = <SparseEntry>t
            z = <ComplexBall>se.entry
            acb_set(acb_mat_entry(self.value, se.i, se.j), z.value)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x):
        """
        Set position ``i``, ``j`` of this matrix to ``x``.

        The object ``x`` must be of type ``ComplexBall``.

        INPUT:

        - ``i`` -- row

        - ``j`` -- column

        - ``x`` -- must be ComplexBall! The value to set self[i,j] to.

        EXAMPLES::

            sage: a = matrix(CBF, 2, 3, range(6)); a
            [                0 1.000000000000000 2.000000000000000]
            [3.000000000000000 4.000000000000000 5.000000000000000]
            sage: a[0, 0] = 10
            sage: a
            [10.00000000000000  1.000000000000000  2.000000000000000]
            [3.000000000000000  4.000000000000000  5.000000000000000]
        """
        acb_set(acb_mat_entry(self.value, i, j), (<ComplexBall> x).value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Return ``(i, j)`` entry of this matrix as a new ComplexBall.

        .. warning::

           This is very unsafe; it assumes ``i`` and ``j`` are in the right
           range.

        EXAMPLES::

            sage: a = MatrixSpace(CBF, 3)(range(9)); a
            [                0 1.000000000000000 2.000000000000000]
            [3.000000000000000 4.000000000000000 5.000000000000000]
            [6.000000000000000 7.000000000000000 8.000000000000000]
            sage: a[1, 2]
            5.000000000000000
            sage: a[4, 7]
            Traceback (most recent call last):
            ...
            IndexError: matrix index out of range
            sage: a[-1, 0]
            6.000000000000000
        """
        cdef ComplexBall z = ComplexBall.__new__(ComplexBall)
        z._parent = self._base_ring
        acb_set(z.value, acb_mat_entry(self.value, i, j))
        return z

    cpdef _richcmp_(left, right, int op):
        r"""
        EXAMPLES::

            sage: a = matrix(CBF, [[1,2],[3,4]])
            sage: b = matrix(CBF, [[1,2],[3,4]])
            sage: a == b
            True
            sage: a + 1/3 == b + 1/3
            False
            sage: a < b
            Traceback (most recent call last):
            ...
            TypeError: no order is defined on complex ball matrices

        TESTS::

            sage: a = matrix(CBF, [1/3])
            sage: b = matrix(CBF, [1/3])
            sage: a == a or b == b or a[0,0] == a[0,0] or a[0,0] == b[0,0]
            False
        """
        cdef Matrix_complex_ball_dense lt = <Matrix_complex_ball_dense> left
        cdef Matrix_complex_ball_dense rt = <Matrix_complex_ball_dense> right
        if op == Py_EQ:
            return acb_mat_eq(lt.value, rt.value)
        elif op == Py_NE:
            return acb_mat_ne(lt.value, rt.value)
        else:
            raise TypeError("no order is defined on complex ball matrices")

    def identical(self, Matrix_complex_ball_dense other):
        r"""
        Test if the corresponding entries of two complex ball matrices
        represent the same balls.

        EXAMPLES::

            sage: a = matrix(CBF, [[1/3,2],[3,4]])
            sage: b = matrix(CBF, [[1/3,2],[3,4]])
            sage: a == b
            False
            sage: a.identical(b)
            True
        """
        return acb_mat_equal(self.value, other.value)

    def overlaps(self, Matrix_complex_ball_dense other):
        r"""
        Test if two matrices with complex ball entries represent overlapping
        sets of complex matrices.

        EXAMPLES::

            sage: b = CBF(0, RBF(0, rad=0.1r)); b
            [+/- 0.101]*I
            sage: matrix(CBF, [0, b]).overlaps(matrix(CBF, [b, 0]))
            True
            sage: matrix(CBF, [1, 0]).overlaps(matrix(CBF, [b, 0]))
            False
        """
        return acb_mat_overlaps(self.value, other.value)

    def contains(self, Matrix_complex_ball_dense other):
        r"""
        Test if the set of complex matrices represented by ``self`` is
        contained in that represented by ``other``.

        EXAMPLES::

            sage: b = CBF(0, RBF(0, rad=.1r)); b
            [+/- 0.101]*I
            sage: matrix(CBF, [0, b]).contains(matrix(CBF, [0, 0]))
            True
            sage: matrix(CBF, [0, b]).contains(matrix(CBF, [b, 0]))
            False
            sage: matrix(CBF, [b, b]).contains(matrix(CBF, [b, 0]))
            True
        """
        return acb_mat_contains(self.value, other.value)

    def __neg__(self):
        r"""
        TESTS::

            sage: -matrix(CBF, [[1,2]])
            [-1.000000000000000 -2.000000000000000]
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, self._ncols)
        sig_on()
        acb_mat_neg(res.value, self.value)
        sig_off()
        return res

    cpdef _add_(self, other):
        r"""
        TESTS::

            sage: matrix(CBF, [[1,2]])._add_(matrix(CBF, [3,4]))
            [4.000000000000000 6.000000000000000]
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, self._ncols)
        sig_on()
        acb_mat_add(res.value, self.value, (<Matrix_complex_ball_dense> other).value, prec(self))
        sig_off()
        return res

    cpdef _sub_(self, other):
        r"""
        TESTS::

            sage: matrix(CBF, [[1,2]])._sub_(matrix(CBF, [3,4]))
            [-2.000000000000000 -2.000000000000000]
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, self._ncols)
        sig_on()
        acb_mat_sub(res.value, self.value, (<Matrix_complex_ball_dense> other).value, prec(self))
        sig_off()
        return res

    cpdef _lmul_(self, Element a):
        r"""
        TESTS::

            sage: matrix(CBF, [[1,2]])._lmul_(CBF(I))
            [1.000000000000000*I 2.000000000000000*I]
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, self._ncols)
        sig_on()
        acb_mat_scalar_mul_acb(res.value, self.value, (<ComplexBall> a).value, prec(self))
        sig_off()
        return res

    cpdef _rmul_(self, Element a):
        r"""
        TESTS::

            sage: matrix(CBF, [[1,2]])._rmul_(CBF(I))
            [1.000000000000000*I 2.000000000000000*I]
        """
        return self._lmul_(a)

    cdef _matrix_times_matrix_(self, Matrix other):
        r"""
        TESTS::

            sage: matrix(CBF, [[1,2]])*matrix([[3], [4]]) # indirect doctest
            [11.00000000000000]
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, other._ncols)
        sig_on()
        acb_mat_mul(res.value, self.value, (<Matrix_complex_ball_dense> other).value, prec(self))
        sig_off()
        return res

    cpdef _pow_int(self, n):
        r"""
        Return the ``n``-th power of this matrix.

        EXAMPLES::

            sage: mat = matrix(CBF, [[1/2, 1/3], [1, 1]])
            sage: mat**2
            [[0.5833333333333...] [0.500000000000000 +/- 1.95e-16]]
            [               1.500000000000000 [1.333333333333333 +/- 5.37e-16]]
            sage: mat**(-2)
            [ [48.00000000000...] [-18.00000000000...]]
            [[-54.0000000000...]  [21.000000000000...]]

        TESTS::

            sage: mat**(0r)
            [1.000000000000000                 0]
            [                0 1.000000000000000]

            sage: mat**(1/2)
            Traceback (most recent call last):
                ...
            NotImplementedError: non-integral exponents not supported

            sage: (-(matrix(CBF, [2])**(-2**100))[0,0].log(2)).log(2)
            [100.000000000000 +/- 7.13e-14]
            sage: (-(matrix(CBF, [2])**(-2**64+1))[0,0].log(2)).log(2)
            [64.0000000000000 +/- 1.34e-14]
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, self._ncols)
        cdef Matrix_complex_ball_dense tmp
        cdef unsigned long expo
        n = Integer(n)
        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix")

        neg = (n < 0)
        if neg:
            n = -n
        if mpz_fits_ulong_p((<Integer>n).value):
            expo = mpz_get_ui((<Integer>n).value)
            sig_on()
            acb_mat_pow_ui(res.value, self.value, expo, prec(self))
            sig_off()
        else:
            tmp = generic_power_pos(self, n)
            acb_mat_set(res.value, tmp.value)
        if neg:
            sig_on()
            acb_mat_inv(res.value, res.value, prec(self))
            sig_off()

        return res

    def __invert__(self):
        r"""
        TESTS::

            sage: ~matrix(CBF, [[1/2, 1/3], [1, 1]])
            [ [6.00000000000000 +/- 3.78e-15] [-2.00000000000000 +/- 1.89e-15]]
            [[-6.00000000000000 +/- 3.78e-15]  [3.00000000000000 +/- 1.89e-15]]
            sage: ~matrix(CBF, [[1/2, 1/3]])
            Traceback (most recent call last):
            ...
            ArithmeticError: self must be a square matrix
        """
        if not self.is_square():
            raise ArithmeticError("self must be a square matrix")
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, self._ncols)
        sig_on()
        acb_mat_inv(res.value, self.value, prec(self))
        sig_off()
        return res

    def _solve_right_nonsingular_square(self, Matrix_complex_ball_dense rhs, check_rank=None):
        r"""
        TESTS::

            sage: matrix(CBF, [[1/2, 1/3], [1, 1]]) \ vector([-1, 1])
            ([-8.00000000000000 +/- ...], [9.00000000000000 +/- ...])
            sage: matrix(CBF, 2, 2, 0) \ vector([-1, 1])
            Traceback (most recent call last):
            ...
            ValueError: matrix equation has no solutions
            sage: b = CBF(0, RBF(0, rad=.1r))
            sage: matrix(CBF, [[1, 1], [0, b]]) \ vector([-1, 1])
            Traceback (most recent call last):
            ...
            ValueError: unable to invert this matrix
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, rhs._ncols)
        sig_on()
        success = acb_mat_solve(res.value, self.value, rhs.value, min(prec(self), prec(rhs)))
        sig_off()
        if success:
            return res
        else:
            raise ValueError("unable to invert this matrix")

    def determinant(self):
        r"""
        Compute the determinant of this matrix.

        EXAMPLES::

            sage: matrix(CBF, [[1/2, 1/3], [1, 1]]).determinant()
            [0.1666666666666667 +/- 7.04e-17]
            sage: matrix(CBF, [[1/2, 1/3], [1, 1]]).det()
            [0.1666666666666667 +/- 7.04e-17]
            sage: matrix(CBF, [[1/2, 1/3]]).determinant()
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix
        """
        cdef ComplexBall res = ComplexBall.__new__(ComplexBall)
        res._parent = self._base_ring
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        sig_on()
        acb_mat_det(res.value, self.value, prec(self))
        sig_off()
        return res

    def trace(self):
        r"""
        Compute the trace of this matrix.

        EXAMPLES::

            sage: matrix(CBF, [[1/3, 1/3], [1, 1]]).trace()
            [1.333333333333333 +/- 5.37e-16]
            sage: matrix(CBF, [[1/2, 1/3]]).trace()
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix
        """
        cdef ComplexBall res = ComplexBall.__new__(ComplexBall)
        res._parent = self._base_ring
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        sig_on()
        acb_mat_trace(res.value, self.value, prec(self))
        sig_off()
        return res

    def charpoly(self, var='x', algorithm=None):
        r"""
        Compute the characteristic polynomial of this matrix.

        EXAMPLES::

            sage: from sage.matrix.benchmark import hilbert_matrix
            sage: mat = hilbert_matrix(5).change_ring(ComplexBallField(10))
            sage: mat.charpoly()
            x^5 + ([-1.8 +/- 0.0258])*x^4 + ([0.3 +/- 0.0567])*x^3 +
            ([+/- 0.0212])*x^2 + ([+/- 0.0266])*x + [+/- 0.0285]

        TESTS::

            sage: mat.charpoly(algorithm="hessenberg")
            x^5 + ([-1.8 +/- 0.0445])*x^4 + ([0.3 +/- 0.0828])*x^3
            + ([+/- 0.0163])*x^2 + ([+/- 5.95e-4])*x + [+/- 6.83e-6]
            sage: mat.charpoly('y')
            y^5 + ([-1.8 +/- 0.0258])*y^4 + ([0.3 +/- 0.0567])*y^3 +
            ([+/- 0.0212])*y^2 + ([+/- 0.0266])*y + [+/- 0.0285]
        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        if algorithm is not None:
            return super(Matrix_dense, self).charpoly(var=var, algorithm=algorithm)
        Pol = polynomial_ring_constructor._single_variate(self.base_ring(), var)
        cdef Polynomial_complex_arb res = Polynomial_complex_arb(Pol)
        sig_on()
        acb_mat_charpoly(res.__poly, self.value, prec(self))
        sig_off()
        return res

    def exp(self):
        r"""
        Compute the exponential of this matrix.

        EXAMPLES::

            sage: matrix(CBF, [[i*pi, 1], [0, i*pi]]).exp()
            [[-1.00000000000000 +/- 8.04e-16] + [+/- 7.05e-16]*I [-1.00000000000000 +/- 8.04e-16] + [+/- 7.05e-16]*I]
            [                                                  0 [-1.00000000000000 +/- 8.04e-16] + [+/- 7.05e-16]*I]
            sage: matrix(CBF, [[1/2, 1/3]]).exp()
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix
        """
        cdef Matrix_complex_ball_dense res = self._new(self._nrows, self._ncols)
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        sig_on()
        acb_mat_exp(res.value, self.value, prec(self))
        sig_off()
        return res
