r"""
Arbitrary precision :class:`RealBallField` matrices.

This module contains :class:`Matrix_real_ball_dense`, a minimal
subclass of :class:`Matrix_dense` that can solving linear interval
systems.
"""
from cysignals.signals cimport sig_on, sig_str, sig_off
from sage.libs.arb.arb cimport *
from sage.libs.arb.arb_mat cimport *
from .args cimport SparseEntry, MatrixArgs_init
from sage.structure.parent cimport Parent
from sage.rings.real_arb cimport (RealBall)

cdef class Matrix_real_ball_dense(Matrix_dense):
    r"""
    A matrix whose entries live in a `RealBallField`. At the moment,
    the only nontrivial method it implements
    :meth:`_solve_right_nonsingular_square`, to avoid the issues that
    arise when the naive superclass implementation is used to perform
    interval arithmetic. This is part of :trac:`29729`.

    The underlying matrix type is ``arb_mat`` of the Arb library.

    TESTS:

    Ensure that we can create a matrix over a field with a random
    number of bits of precision, to catch corner cases::

        sage: set_random_seed()
        sage: bits = ZZ.random_element(2,1000)
        sage: R = RealBallField(bits)
        sage: nrows = ZZ.random_element(100)
        sage: ncols = ZZ.random_element(100)
        sage: matrix.random(R, nrows, ncols).base_ring()
        Real ball field with...

    """
    def __cinit__(self):
        r"""
        Allocate memory for and initialize my underlying Arb matrix C
        structure.
        """
        sig_str("Arb exception")
        arb_mat_init(self.value, self._nrows, self._ncols)
        sig_off()


    def __dealloc__(self):
        """
        Free the memory allocated for my underlying Arb matrix C structure.

        TESTS:

        Ensure that nothing outrageous happens if we call ``del`` on
        an instance of this class::

            sage: set_random_seed()
            sage: bits = ZZ.random_element(2,100)
            sage: R = RealBallField(bits)
            sage: nrows = ZZ.random_element(10)
            sage: ncols = ZZ.random_element(10)
            sage: A = matrix.random(R, nrows, ncols).base_ring()
            sage: del(A)
            sage: A
            Traceback (most recent call last):
            ...
            NameError: name 'A' is not defined

        """
        arb_mat_clear(self.value)


    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        Initialize this dense matrix over a real ball field.

        INPUT:

        - ``parent`` -- a matrix space over a real ball field

        - ``entries`` -- the entries of this matrix, see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if ``False``, assume without checking that the
          entries lie in the base ring (this can crash SageMath if
          you're wrong)

        EXAMPLES:

        We create a two-by-two matrix with the default precision::

            sage: from sage.matrix.matrix_real_ball_dense import *
            sage: MS = MatrixSpace(RBF, 2, 2)
            sage: Matrix_real_ball_dense(MS, range(4))
            [                0 1.000000000000000]
            [2.000000000000000 3.000000000000000]

        With ``coerce=False``, nonsense (or worse...) can result::

            sage: from sage.matrix.matrix_real_ball_dense import *
            sage: MS = MatrixSpace(RBF, 1, 1)
            sage: Matrix_real_ball_dense(MS, [pi], coerce=False)
            [[+/- inf]]

        """
        ma = MatrixArgs_init(parent, entries)
        cdef RealBall z
        for t in ma.iter(coerce, True):
            se = <SparseEntry>t
            z = <RealBall>se.entry
            arb_set(arb_mat_entry(self.value, se.i, se.j), z.value)


    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x):
        r"""
        Set entry ``i,j1` of this matrix to ``x``, unsafely.

        The new value ``x`` must be of type ``RealBall`` and
        ``i,j`` must be within bounds.

        INPUT:

        - ``i`` -- row index

        - ``j`` -- column index

        - ``x`` -- a ``RealBall``; the value to set self[i,j] to.

        TESTS:

        The unsafe setter is indirectly used by the safe one:

            sage: A = matrix(RBF, 1, 1, 0); A
            [0]
            sage: A[0,0] = e; A
            [[2.718281828459045 +/- 5.35e-16]]

        """
        arb_set(arb_mat_entry(self.value, i, j), (<RealBall> x).value)


    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Get entry ``i,j`` of this matrix, unsafely.

        The indices ``i,j`` must be within bounds.

        INPUT:

        - ``i`` -- row index

        - ``j`` -- column index

        .. NOTE::

           The unsafe getter is indirectly tested by the safe getter that
           is used everywhere (in particular, to display a matrix).

        """
        cdef RealBall z = RealBall.__new__(RealBall)
        z._parent = self._base_ring
        arb_set(z.value, arb_mat_entry(self.value, i, j))
        return z


    def _solve_right_nonsingular_square(self,
                                        Matrix_real_ball_dense rhs,
                                        check_rank=None):
        r"""
        Override the superclass method with a square/nonsingular
        algorithm (from Arb) tailored to interval systems.

        EXAMPLES:

        Numerical noise should not trick SageMath into thinking that
        this system is nonsingular (:trac:`29729`)::

            sage: A = matrix(RBF, [[2/3, 1], [2/5, 3/5]])
            sage: b = vector(RBF, [1, 1])
            sage: A.solve_right(b)
            Traceback (most recent call last):
            ...
            ValueError: unable to invert this matrix

        """
        result_space = self.matrix_space(self._nrows, rhs._ncols)
        cdef Matrix_real_ball_dense res
        res = Matrix_real_ball_dense(result_space, coerce=False)
        soln_prec = min(self._base_ring._prec, rhs._base_ring._prec)

        sig_on()
        success = arb_mat_solve(res.value, self.value, rhs.value, soln_prec)
        sig_off()
        if success:
            return res
        else:
            raise ValueError("unable to invert this matrix")
