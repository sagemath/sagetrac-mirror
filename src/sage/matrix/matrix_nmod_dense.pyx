# distutils: extra_compile_args = -D_XPG6
# flags chosen from libs/flint/fmpz_poly.pyx
r"""
FLINT nmod_mat class wrapper

This file implements matrices over `\ZZ/N\ZZ` for `N < 2^{63}`.
It adds some capabilities for composite `N` that are not present in FLINT.

AUTHORS:

- Edgar Costa, David Roe (2021) Initial version.
"""

from cpython.sequence cimport *
from cysignals.signals cimport sig_on, sig_str, sig_off

from sage.structure.element cimport Element, Matrix
from sage.structure.element import is_Vector
from sage.categories.cartesian_product import cartesian_product
from sage.structure.richcmp cimport rich_to_bool, Py_EQ, Py_NE
from sage.structure.factorization import Factorization
from sage.structure.proof.all import linear_algebra as linalg_proof
from sage.misc.prandom import randrange
from sage.arith.all import gcd
from sage.arith.power cimport generic_power
from sage.arith.long cimport integer_check_long_py
from sage.rings.polynomial.polynomial_zmod_flint cimport Polynomial_zmod_flint
from sage.rings.integer cimport Integer
from sage.rings.finite_rings.integer_mod cimport (
    IntegerMod_abstract,
    IntegerMod_int,
    IntegerMod_int64,
    IntegerMod_gmp
)
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.special import identity_matrix
from .args cimport SparseEntry, MatrixArgs_init

from sage.libs.flint.nmod_mat cimport *
from sage.libs.flint.nmod_poly cimport (
    nmod_poly_set,
    nmod_poly_set_coeff_ui,
    nmod_poly_get_coeff_ui,
    nmod_poly_fit_length
)
from sage.libs.flint.ulong_extras cimport (
    n_precompute_inverse,
    n_preinvert_limb,
    n_pow,
    n_addmod,
    n_submod,
    n_invmod,
    n_mod2_preinv,
    n_mulmod2_preinv,
    n_remove2_precomp,
    n_CRT,
    n_div2_preinv,
    n_divrem2_precomp,
)
from sage.libs.gmp.mpz cimport mpz_sgn,  mpz_fits_ulong_p, mpz_get_ui, mpz_get_si


def poly_crt(S, polys, moduli):
    r"""
    Implements the Chinese remainder theorem for polynomials over `\ZZ/m\ZZ`.

    INPUT:

    - ``S`` -- a polynomial ring over `\ZZ/N\ZZ` with `N < 2^{63}`.

    - ``polys`` -- a list of polynomials modulo smaller integers `m`.

    - ``moduli`` -- the moduli `m` for the polynomials.  The list must have the same length, and the product must be `N`.

    OUTPUT:

    The polynomial in S reducing to each polynomial in ``polys`` modulo the integers given in ``moduli``.

    EXAMPLES::

        sage: from sage.matrix.matrix_nmod_dense import poly_crt
        sage: moduli = [4, 9, 25, 49]
        sage: N = prod(moduli)
        sage: S.<x> = Zmod(N)[]
        sage: polys = [Zmod(m)['x']([sqrt(m)] + [0]*(sqrt(m)-1) + [1]) for m in moduli]
        sage: f = poly_crt(S, polys, moduli); f
        27000*x^7 + 15876*x^5 + 34300*x^3 + 11025*x^2 + 40530
        sage: all(g == f.change_ring(Zmod(m)) for (g, m) in zip(polys, moduli))
        True

    TESTS:

    The algorithm involves sorting the polynomials; we check that the result is preserved by reordering::

        sage: polys.reverse()
        sage: moduli.reverse()
        sage: f == poly_crt(S, polys, moduli)
        True
    """
    if len(polys) == 1:
        return polys[0]
    cdef Py_ssize_t i, j, d
    cdef mp_limb_t N
    # sort polys and moduli in parallel
    polys, moduli = zip(*sorted(zip(polys, moduli), key=lambda pair: pair[0].degree()))
    cdef Polynomial_zmod_flint f = S()
    nmod_poly_fit_length(&f.x, polys[-1].degree()+1)
    N = 1
    d = 0
    for i in range(len(polys)):
        d = polys[i].degree()
        for j in range(d+1):
            nmod_poly_set_coeff_ui(&f.x, j, n_CRT(nmod_poly_get_coeff_ui(&f.x, j), N, nmod_poly_get_coeff_ui(&(<Polynomial_zmod_flint?>polys[i]).x, j), <unsigned long>moduli[i]))
        N *= moduli[i]
    return f


cdef class Matrix_nmod_dense(Matrix_dense):
    r"""
    Matrices modulo `N` for `N < 2^{63}`

    EXAMPLES::

        sage: A = matrix(Zmod(36), 3, 3, range(9))
        sage: type(A)
        <class 'sage.matrix.matrix_nmod_dense.Matrix_nmod_dense'>
    """
    ########################################################################
    # LEVEL 1 helpers:
    #   These function support the implementation of the level 1 functionality.
    ########################################################################
    def __cinit__(self, parent, *args, **kwds):
        """
        Memory initialization

        EXAMPLES::

            sage: M = MatrixSpace(Zmod(36), 2, 2)
            sage: M()
            [0 0]
            [0 0]
        """
        self._modulus = parent._base._pyx_order
        sig_str("FLINT exception")
        nmod_mat_init(self._matrix, self._nrows, self._ncols, mpz_get_ui(self._modulus.sageInteger.value))
        sig_off()

    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        """
        Initialization

        EXAMPLES::

            sage: A = matrix(Zmod(36), 3, 3, range(9))
            sage: TestSuite(A).run()
        """
        self._parent = parent
        ma = MatrixArgs_init(parent, entries)
        cdef SparseEntry se
        for se in ma.iter(coerce, sparse=True):
            nmod_mat_set_entry(self._matrix, se.i, se.j, se.entry)

    def __dealloc__(self):
        """
        Memory deallocation

        EXAMPLES::

            sage: A = matrix(Zmod(36), 3, 3, range(9))
            sage: del A
        """
        sig_on()
        nmod_mat_clear(self._matrix)
        sig_off()

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x):
        """
        Low level interface for setting entries, used by generic matrix code.
        """
        e = self._modulus.element_class()
        cdef mp_limb_t ivalue
        if e is IntegerMod_int:
            ivalue = (<IntegerMod_int>x).ivalue
        elif e is IntegerMod_int64:
            ivalue = (<IntegerMod_int64>x).ivalue
        else:
            ivalue = mpz_get_ui((<IntegerMod_gmp>x).value)
        nmod_mat_set_entry(self._matrix, i, j, ivalue)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Low level interface for getting entries as an element of the base ring,
        used by generic matrix code.
        """
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
            P = MatrixSpace(self._parent._base, nrows, ncols, sparse=False)
        cdef Matrix_nmod_dense ans = Matrix_nmod_dense.__new__(Matrix_nmod_dense, P)
        ans._parent = P
        return ans

    ########################################################################
    # LEVEL 2 helpers:
    #   These function support the implementation of the level 2 functionality.
    ########################################################################

    cpdef _add_(self, _right):
        """
        Addition

        EXAMPLES::

            sage: A = matrix(Zmod(36), 3, 3, range(9))
            sage: B = matrix(Zmod(36), 3, 3, range(2, 20, 2))
            sage: A + B
            [ 2  5  8]
            [11 14 17]
            [20 23 26]
        """
        cdef Matrix_nmod_dense right = _right
        cdef Matrix_nmod_dense M = self._new(self._nrows, self._ncols)
        sig_on()
        nmod_mat_add(M._matrix, self._matrix, right._matrix)
        sig_off()
        return M

    cdef Matrix _matrix_times_matrix_(left, Matrix _right):
        """
        Multiplication

        EXAMPLES::

            sage: A = matrix(Zmod(36), 2, 3, range(6))
            sage: B = matrix(Zmod(36), 3, 2, range(6))
            sage: A * B
            [10 13]
            [28  4]
        """
        if left._ncols != _right._nrows:
            raise IndexError("Number of columns of self must equal number of rows of right.")
        cdef Matrix_nmod_dense right = _right
        cdef Matrix_nmod_dense M = left._new(left._nrows, right._ncols)
        sig_on()
        nmod_mat_mul(M._matrix, left._matrix, right._matrix)
        sig_off()
        return M

    cpdef _lmul_(self, Element right):
        """
        Scalar multiplication

        EXAMPLES::

            sage: A = matrix(Zmod(36), 3, 3, range(0, 18, 2))
            sage: 9*A
            [ 0 18  0]
            [18  0 18]
            [ 0 18  0]
        """
        cdef Matrix_nmod_dense M = self._new(self._nrows, self._ncols)
        e = self._modulus.element_class()
        cdef mp_limb_t ivalue
        if e is IntegerMod_int:
            ivalue = (<IntegerMod_int>right).ivalue
        elif e is IntegerMod_int64:
            ivalue = (<IntegerMod_int64>right).ivalue
        else:
            ivalue = mpz_get_ui((<IntegerMod_gmp>right).value)
        sig_on()
        nmod_mat_scalar_mul(M._matrix, self._matrix, ivalue)
        sig_off()
        return M


    def __copy__(self):
        """
        Return a copy of this matrix. Changing the entries of the copy will
        not change the entries of this matrix.

        EXAMPLES::

            sage: A = matrix(Zmod(36), 2, 2, range(4))
            sage: B = copy(A)
            sage: B[1,1] = 35
            sage: A
            [0 1]
            [2 3]
        """
        cdef Matrix_nmod_dense M = self._new(self._nrows, self._ncols)
        sig_on()
        nmod_mat_set(M._matrix, self._matrix)
        sig_off()
        return M

    def __neg__(self):
        r"""
        Return the negative of this matrix.

        EXAMPLES::

            sage: A = matrix(Zmod(36), 2, 2, range(4))
            sage: -A
            [ 0 35]
            [34 33]
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

            sage: A = matrix(Zmod(36), 2, 2, range(4))
            sage: B = matrix(Zmod(36), 2, 2, [0, 1, 3, 2])
            sage: A == A
            True
            sage: A != A
            False
            sage: A == B
            False
            sage: A != B
            True
            sage: A < B
            True
        """
        cdef Py_ssize_t i, j
        cdef int k

        if op == Py_EQ:
            return bool(nmod_mat_equal(self._matrix, (<Matrix_nmod_dense>right)._matrix))
        elif op == Py_NE:
            return not bool(nmod_mat_equal(self._matrix, (<Matrix_nmod_dense>right)._matrix))
        else:
            sig_on()
            for i in range(self._nrows):
                for j in range(self._ncols):
                    k = nmod_mat_get_entry(self._matrix,i,j) - nmod_mat_get_entry((<Matrix_nmod_dense>right)._matrix,i,j)
                    if k:
                        sig_off()
                        if k < 0:
                            return rich_to_bool(op, -1)
                        else:
                            return rich_to_bool(op, 1)
            sig_off()
            return rich_to_bool(op, 0)

    ########################################################################
    # LEVEL 3 helpers:
    #   These function support the implementation of the level 2 functionality.
    ########################################################################

    def __nonzero__(self):
        """
        EXAMPLES::

            sage: bool(matrix(Zmod(36), 2, 2, 1))
            True
            sage: bool(matrix(Zmod(36), 0, 0))
            False
        """
        return not nmod_mat_is_zero(self._matrix)

    cpdef _sub_(self, _right):
        """
        Subtraction

        EXAMPLES::

            sage: A = matrix(Zmod(36), 3, 3, range(9))
            sage: B = matrix(Zmod(36), 3, 3, range(2, 20, 2))
            sage: A - B
            [34 33 32]
            [31 30 29]
            [28 27 26]
        """
        cdef Matrix_nmod_dense right = _right
        cdef Matrix_nmod_dense M = self._new(self._nrows, self._ncols)
        nmod_mat_sub(M._matrix, self._matrix, right._matrix)
        return M

    def __invert__(self):
        r"""
        Return the inverse of this matrix.

        EXAMPLES::

            sage: A = matrix(Zmod(36), [[28, 32, 19], [25, 24, 2], [15, 11, 30]])
            sage: B = ~A
            sage: A * B == B * A == identity_matrix(Zmod(36), 3, 3)
            True
        """
        return self.inverse_of_unit()

    def inverse_of_unit(self, algorithm="crt"):
        """
        Return the inverse of this matrix.

        Raises a ``ZeroDivisionError`` if the determinant is not a unit,
        and raises an ``ArithmeticError`` if the
        inverse doesn't exist because the matrix is nonsquare.

        INPUT:

        - ``algorithm`` -- either ``"crt"`` or ``"lift"``.
            In the first case, CRT is used to assemble results from prime powers
            dividing the modulus, with quadratic Hensel lifting for prime powers.
            In the second case, the inverse is computed over the integers and
            then reduced.

        EXAMPLES::

            sage: A = matrix(Zmod(36), [[28, 32, 19], [25, 24, 2], [15, 11, 30]])
            sage: A.inverse_of_unit()
            [14  5  4]
            [ 0 15 23]
            [23 28 16]

            sage: for N in [5, 625, 36, 2^24, 2^6*3^9, 2^63-1]:
            ....:     for n in [3, 6]:
            ....:         A = random_matrix(Zmod(N), n, n)
            ....:         while not A.det().is_unit():
            ....:             A = random_matrix(Zmod(N), n, n)
            ....:         assert A.inverse_of_unit('crt') == A.inverse_of_unit('lift')

            sage: matrix(Zmod(36), [[2, 0], [0, 2]]).inverse_of_unit()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: input matrix must be nonsingular
            sage: matrix(Zmod(36), [[2], [2]]).inverse_of_unit()
            Traceback (most recent call last):
            ...
            ArithmeticError: inverse only defined for square matrix
        """
        if not self.is_square():
            raise ArithmeticError("inverse only defined for square matrix")
        if not self.nrows():
            return self

        cdef Matrix_nmod_dense M
        R = self._parent._base
        cdef Py_ssize_t i, j, n = self._nrows
        cdef long k, e, b, nlifts
        cdef mp_limb_t p, N = 1
        cdef nmod_mat_t inv, A, B, tmp
        cdef bint lift_required, ok
        if R.is_field():
            M = self._new(n, n)
            ok = nmod_mat_inv(M._matrix, self._matrix)
            if not ok:
                raise ZeroDivisionError("input matrix must be nonsingular")
            return M
        if algorithm == "lift":
            return (~self.lift_centered()).change_ring(R)
        if algorithm == "crt":
            # The modulus is small, so factoring is feasible.
            # We find inverses modulo each prime dividing the modulus, lift p-adically, then CRT the results together
            M = self._new(n, n)
            F = R.factored_order()
            maxe = max(fac[1] for fac in F)
            lift_required = (maxe > 1)
            if lift_required:
                nlifts = Integer(maxe - 1).nbits()
                lift_mods = [1 for _ in range(nlifts)]
            try:
                nmod_mat_init(A, n, n, 1)
                nmod_mat_init(inv, n, n, 1)
                if lift_required:
                    nmod_mat_init(tmp, n, n, 1)
                for pz, ez in F:
                    p = pz
                    _nmod_mat_set_mod(A, p)
                    for i in range(n):
                        for j in range(n):
                            nmod_mat_set_entry(A, i, j, nmod_mat_get_entry(self._matrix, i, j) % p)
                    ok = nmod_mat_inv(A, A)
                    if not ok:
                        raise ZeroDivisionError("input matrix must be nonsingular")
                    _nmod_mat_set_mod(inv, N*p)
                    for i in range(n):
                        for j in range(n):
                            nmod_mat_set_entry(inv, i, j, n_CRT(nmod_mat_get_entry(inv, i, j), N, nmod_mat_get_entry(A, i, j), p))
                    N *= p
                    # Now inv is accurate modulo N
                    if lift_required:
                        e = ez
                        # Update the moduli for lifting so that the exponent of p at most doubles
                        for k in reversed(range(nlifts)):
                            lift_mods[k] *= p**e
                            e = (e + 1) // 2
                if lift_required:
                    for k in range(nlifts):
                        N = lift_mods[k]
                        _nmod_mat_set_mod(tmp, N)
                        _nmod_mat_set_mod(inv, N)
                        _nmod_mat_set_mod(A, N)
                        for i in range(n):
                            for j in range(n):
                                nmod_mat_set_entry(A, i, j, nmod_mat_get_entry(self._matrix, i, j) % N)
                        # Suppose B is an approximation to the inverse of A.
                        # I - AB = p^k E
                        # (I - AB)(I - AB) = I - A(2B - BAB) = p^(2k)E^2
                        # So 2B - BAB is a better approximation
                        nmod_mat_mul(tmp, A, inv)
                        nmod_mat_mul(tmp, inv, tmp)
                        nmod_mat_scalar_mul(inv, inv, 2)
                        nmod_mat_sub(inv, inv, tmp)
                        # Now inv is accurate mod N
                nmod_mat_set(M._matrix, inv)
                return M
            finally:
                nmod_mat_clear(A)
                nmod_mat_clear(inv)
                if lift_required:
                    nmod_mat_clear(tmp)
        raise ValueError("Unrecognized algorithm '%s'" % algorithm)

    def hessenberg_form(self):
        """
        The Hessenberg form of a matrix is almost upper triangular,
        and allows for efficient computation of the characteristic polynomial.

        Requires that the base ring have prime power modulus.

        EXAMPLES::

            sage: A = matrix(Zmod(125), 4, [2,1,1,-2,2,2,-1,-1,-1,1,2,3,4,5,6,7])
            sage: H = A.hessenberg_form(); H
            [  2  59  88   1]
            [  2  63  34 124]
            [  0  15 104   8]
            [  0   0  90  94]
            sage: H.charpoly() == A.charpoly()
            True
        """
        B = self.__copy__()
        B.hessenbergize()
        return B

    def hessenbergize(self):
        """
        Transform this matrix into Hessenberg form.

        The hessenberg form of a matrix `A` is a matrix that is
        similar to `A`, so has the same characteristic polynomial
        as `A`, and is upper triangular except possible for entries
        right below the diagonal.

        Requires that the base ring have prime power modulus.

        The algorithm is a modification of the standard one over fields,
        where rather than swapping in an arbitrary nonzero element in each column
        one instead picks an element of minimal valuation.

        EXAMPLES::

            sage: A = matrix(Zmod(125), 4, [2,1,1,-2,2,2,-1,-1,-1,1,2,3,4,5,6,7])
            sage: A.hessenbergize(); A
            [  2  59  88   1]
            [  2  63  34 124]
            [  0  15 104   8]
            [  0   0  90  94]

        You can't Hessenbergize an immutable matrix::

            sage: A = matrix(Zmod(125), 3, range(9))
            sage: A.set_immutable()
            sage: A.hessenbergize()
            Traceback (most recent call last):
            ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).
        """
        R = self._parent._base
        F = R.factored_order()
        if len(F) != 1:
            raise ValueError("Hessenbergize only supported for prime power modulus")
        if not self.is_square():
            raise TypeError("self must be square")
        self.check_mutability()

        pz, ez = F[0]
        cdef Py_ssize_t i, j, k, m, n, r
        n = self._nrows

        self.check_mutability()

        cdef bint found_nonzero
        cdef mp_limb_t u, minu, inv, q, x, y, prepe, prepe1, preq, pe, p = pz
        cdef int v, minv, e = ez
        cdef double pinv = n_precompute_inverse(p)
        pe = n_pow(p, e)
        prepe = n_preinvert_limb(pe)
        prepe1 = n_preinvert_limb(n_pow(p, e-1))

        for m in range(1, n-1):
            i = -1
            minv = e
            for r in range(m, n):
                u = nmod_mat_get_entry(self._matrix, r, m-1)
                if u != 0:
                    v = n_remove2_precomp(&u, p, pinv)
                    if v < minv:
                        i = r
                        minv = v
                        minu = u
            if i != -1:
                # i is the first row below m with an entry in colunmn m-1 with minimal valuation
                q = n_pow(p, e-minv)
                if minv == 0:
                    preq = prepe
                elif minv == 1:
                    preq = prepe1
                else:
                    preq = n_preinvert_limb(q)
                inv = n_invmod(minu, q)
                if i > m:
                    self.swap_rows_c(i, m)
                    # We must do the corresponding column swap to
                    # maintain the characteristic polynomial (which is
                    # an invariant of Hessenberg form)
                    self.swap_columns_c(i, m)
                # Now the nonzero entry in position (m,m-1) is minu * p^minv.
                # Use it to clear the entries in column m-1 below m.
                for j in range(m+1, n):
                    u = nmod_mat_get_entry(self._matrix, j, m-1)
                    if u != 0:
                        # multiply by the inverse of the (m,m-1) entry
                        if minv == 0:
                            # Simpler case if minv=0
                            u = n_mulmod2_preinv(u, inv, q, preq)
                        else:
                            # extract powers of p from u, multiply by the unit
                            # and put the right number of powers of p back.
                            v = n_remove2_precomp(&u, p, pinv)
                            u = n_pow(p, v-minv) * n_mulmod2_preinv(u, inv, q, preq)
                        # self.add_multiple_of_row(j, m, -u, 0)
                        for k in range(n):
                            x = nmod_mat_get_entry(self._matrix, j, k)
                            y = nmod_mat_get_entry(self._matrix, m, k)
                            y = n_mulmod2_preinv(y, u, pe, prepe)
                            nmod_mat_set_entry(self._matrix, j, k, n_submod(x, y, pe))
                        # To maintain charpoly, do the corresponding column operation,
                        # which doesn't mess up the matrix, since it only changes
                        # column m, and we're only worried about column m-1 right now.
                        # Add u*column_j to column_m.
                        for k in range(n):
                            x = nmod_mat_get_entry(self._matrix, k, m)
                            y = nmod_mat_get_entry(self._matrix, k, j)
                            y = n_mulmod2_preinv(y, u, pe, prepe)
                            nmod_mat_set_entry(self._matrix, k, m, n_addmod(x, y, pe))

    def charpoly(self, var='x', algorithm=None):
        """
        Return the characteristic polynomial of this matrix, as a polynomial over the base ring.

        INPUT:

        - ``var`` -- a string, the variable name for the polynomial returned.

        - ``algorithm`` -- either ``"flint"``, ``"lift"`` or ``"crt"``.
            In the first case, use FLINT's characteristic polynomial function,
            in the second case, lift to the integers and compute there,
            in the third form, compute the hessenberg form for each prime power
            dividing the modulus.
            If not given, defaults to ``"flint"`` over fields and ``"crt"`` otherwise.

        EXAMPLES::

            sage: A = matrix(Zmod(36), [[28, 32, 19], [25, 24, 2], [15, 11, 30]])
            sage: A.charpoly()
            x^3 + 8*x^2 + 8
            sage: A = matrix(Zmod(43^10), 6, [0, 5664461354126771, 12212357361910300, 15947020959157478, 0, 16792952041597449, 14690359073749623, 11237259451999884, 5117434014428142, 15157488677243483, 9004103062307752, 20761679499270441, 4620722392655416, 5445142895231681, 6605357538252496, 7608812697273777, 18542817615638637, 18194689690271501, 0, 20341333098836812, 12117922812876054, 1270149214447437, 0, 10999401748338075, 4620722392655416, 10891113386038365, 956055025271903, 2162842206467093, 18542817615638637, 1143972982339214, 13128267973348003, 15817056104759912, 20531311511260484, 13598045280630823, 7585782589268305, 14053895308766769])
            sage: A.charpoly()
            x^6 + 13124967810747524*x^5 + 20067912494391006*x^4 + 11204731077775359*x^3

            sage: for N in [5, 625, 36, 2^24, 2^6*3^9, 2^63-1]:
            ....:     for n in [3, 6]:
            ....:         A = random_matrix(Zmod(N), n, n)
            ....:         assert A.charpoly(algorithm='crt') == A.charpoly(algorithm='lift')
        """
        if not self.is_square():
            raise TypeError("self must be square")
        R = self._parent._base
        cdef Polynomial_zmod_flint f
        cdef Py_ssize_t i, j, m, jstart, n = self._nrows
        cdef mp_limb_t pe, prepe, preN, x, y, scalar, t, N = R.order()
        preN = n_preinvert_limb(N)
        cdef Matrix_nmod_dense A
        cdef nmod_mat_t H, c
        F = R.factored_order()
        if algorithm is None:
            if len(F) == 1 and F[0][1] == 1: # N prime
                algorithm = "flint"
            else:
                algorithm = "crt"
        S = PolynomialRing(R, var, implementation="FLINT")
        if algorithm == "flint":
            f = S()
            sig_on()
            nmod_mat_charpoly(&f.x, self._matrix)
            sig_off()
            return f
        if algorithm == "lift":
            return self.lift_centered().charpoly(var).change_ring(R)
        if algorithm == "crt":
            # We hessenbergize at each prime, CRT, then compute the charpoly from the Hessenberg form
            try:
                nmod_mat_init(H, n, n, 1)
                nmod_mat_init(c, n+1, n+1, N)
                # Now we use N cummulatively as we work through the factorization
                N = 1
                for pz, ez in F:
                    pe = n_pow(pz, ez)
                    prepe = n_preinvert_limb(pe)
                    R = Zmod(pe)
                    R.factored_order.set_cache(Factorization([(pz, ez)]))
                    A = MatrixSpace(R, self._nrows, self._ncols, sparse=False)()
                    for i in range(n):
                        for j in range(n):
                            x = nmod_mat_get_entry(self._matrix, i, j)
                            x = n_mod2_preinv(x, pe, prepe)
                            nmod_mat_set_entry(A._matrix, i, j, x)
                    A.hessenbergize()
                    _nmod_mat_set_mod(H, N*pe)
                    for i in range(n):
                        jstart = i - 1
                        if i == 0:
                            jstart = 0
                        for j in range(jstart, n):
                            nmod_mat_set_entry(H, i, j, n_CRT(nmod_mat_get_entry(H, i, j), N, nmod_mat_get_entry(A._matrix, i, j), pe))
                    N *= pe
                # We represent the intermediate polynomials that come up in
                # the calculations as rows of an (n+1)x(n+1) matrix, since
                # we've implemented basic arithmetic with such a matrix.
                # Also see Cohen's first GTM, Algorithm 2.2.9.
                nmod_mat_set_entry(c, 0, 0, 1)
                for m in range(1, n+1):
                    # Set the m-th row of c to (x - H[m-1,m-1])*c[m-1] = x*c[m-1] - H[m-1,m-1]*c[m-1]
                    # We do this by hand by setting the m-th row to c[m-1]
                    # shifted to the right by one.  We then add
                    # -H[m-1,m-1]*c[m-1] to the resulting m-th row.
                    for i in range(1, n+1):
                        nmod_mat_set_entry(c, m, i, nmod_mat_get_entry(c, m-1, i-1))
                    # self.sub_multiple_of_row(m, m-1, -H[m-1,m-1], 0)
                    scalar = nmod_mat_get_entry(H, m-1, m-1)
                    for j in range(n+1):
                        x = nmod_mat_get_entry(c, m, j)
                        y = nmod_mat_get_entry(c, m-1, j)
                        y = n_mulmod2_preinv(y, scalar, N, preN)
                        nmod_mat_set_entry(c, m, j, n_submod(x, y, N))
                    t = 1
                    for i in range(1, m):
                        t = n_mulmod2_preinv(t, nmod_mat_get_entry(H, m-i, m-i-1), N, preN)
                        scalar = n_mulmod2_preinv(t, nmod_mat_get_entry(H, m-i-1, m-1), N, preN)
                        for j in range(n+1):
                            x = nmod_mat_get_entry(c, m, j)
                            y = nmod_mat_get_entry(c, m-i-1, j)
                            y = n_mulmod2_preinv(y, scalar, N, preN)
                            nmod_mat_set_entry(c, m, j, n_submod(x, y, N))
                # the answer is now the n-th row of c.
                f = S()
                for j in range(n+1):
                    nmod_poly_set_coeff_ui(&f.x, j, nmod_mat_get_entry(c, n, j))
                return f
            finally:
                nmod_mat_clear(H)
                nmod_mat_clear(c)
        raise ValueError("Unknown algorithm '%s'" % algorithm)

    cpdef _shift_mod(self, mp_limb_t modulus, mp_limb_t shift=1, bint mul=True, bint domod=True):
        """
        A fast method for returning a copy of this matrix with
        different modulus or scaled by an integer.

        It is the caller's responsibility to ensure that the combination of
        shift and modulus yields a result that is mathematically meaningful.

        INPUT:

        - ``modulus`` -- an integer, the modulus for the result.

        - ``shift`` -- an integer to multiply or divide by (default 1)

        - ``mul`` -- boolean.  If true, multiply by ``shift``, otherwise divide.
            In the case of division, floor division is used; the intended purpose
            is for when all entries are divisible by ``shift`` and the modulus
            will be divided as well.

        - ``domod`` -- boolean.  If false, assumes that the entries
            don't need to be reduced after scaling.

        EXAMPLES::

            sage: A = matrix(Zmod(12), 2, [2,4,6,8])
            sage: A._shift_mod(6, 2, mul=False, domod=False)
            [1 2]
            [3 4]
            sage: A._shift_mod(24, 2, mul=True, domod=False)
            [ 4  8]
            [12 16]
            sage: A._shift_mod(12, 3, mul=True)
            [6 0]
            [6 0]

        If not all entries are divisible by shift, the result may not be sensible::

            sage: A._shift_mod(4, 3, mul=False)
            [0 1]
            [2 2]

        It's possible to get unreduced matrices::

            sage: A._shift_mod(12, 2, domod=False)
            [ 4  8]
            [12 16]

        Note that you also need to be careful about not reducing if there's a
        possibility of overflow::

            sage: N = 2^63 - 5
            sage: A = matrix(Zmod(N), 1, [-1])
            sage: A._shift_mod(N, 3, domod=False)
            [9223372036854775790]
            sage: (3*(N-1)) % (2^64)
            9223372036854775790
            sage: A.lift() * 3
            [27670116110564327406]
            sage: (A.lift() * 3).change_ring(Zmod(N))
            [9223372036854775800]

        This works correctly if you use reduction::

            sage: A._shift_mod(N, 3, domod=True)
            [9223372036854775800]
        """
        cdef Py_ssize_t i, j, m = self._nrows, n = self._ncols
        cdef mp_limb_t x, preN, preshift, N = modulus
        preN = n_preinvert_limb(N)
        if not mul:
            preshift = n_preinvert_limb(shift)
        R = Zmod(N)
        P = MatrixSpace(R, m, n, sparse=False)
        cdef Matrix_nmod_dense ans = Matrix_nmod_dense.__new__(Matrix_nmod_dense, P)
        ans._parent = P
        for i in range(m):
            for j in range(n):
                x = nmod_mat_get_entry(self._matrix, i, j)
                if shift == 1:
                    if domod:
                        x = n_mod2_preinv(x, N, preN)
                elif mul:
                    if domod:
                        x = n_mulmod2_preinv(x, shift, N, preN)
                    else:
                        x *= shift
                else:
                    x = n_div2_preinv(x, shift, preshift)
                    if domod:
                        x = n_mod2_preinv(x, N, preN)
                nmod_mat_set_entry(ans._matrix, i, j, x)
        return ans

    def minpoly_ideal(self, var='x', proof=None, **kwds):
        """
        The ideal of polynomials over the base ring that vanish on this matrix.

        When the base ring is not a field, this ideal is not necessarily principal.

        INPUT:

        - ``var`` -- the variable name for the polynomial ring

        - ``proof`` -- boolean, whether to check that the answer computed by using a minimal polynomial is correct.  If not specificied, uses the linear alegbra default proof state.  Note that if the modulus is composite and divisible by small primes the probability of an incorrect result is substantial.  The result is not cached when proof is ``False``, so this function can be called multiple times to get a desired level of certainty.

        EXAMPLES::

            sage: A = matrix(Zmod(36), 2, [6, 0, 0, -6])
            sage: A.minpoly_ideal()
            Ideal (x^2, 3*x + 18) of Univariate Polynomial Ring in x over Ring of integers modulo 36

            sage: A = matrix(Zmod(43^10), 6, [0, 5664461354126771, 12212357361910300, 15947020959157478, 0, 16792952041597449, 14690359073749623, 11237259451999884, 5117434014428142, 15157488677243483, 9004103062307752, 20761679499270441, 4620722392655416, 5445142895231681, 6605357538252496, 7608812697273777, 18542817615638637, 18194689690271501, 0, 20341333098836812, 12117922812876054, 1270149214447437, 0, 10999401748338075, 4620722392655416, 10891113386038365, 956055025271903, 2162842206467093, 18542817615638637, 1143972982339214, 13128267973348003, 15817056104759912, 20531311511260484, 13598045280630823, 7585782589268305, 14053895308766769])
            sage: I = A.minpoly_ideal(); I
            Ideal (x^4 + 3889828461*x^3 + 5524445805550*x^2 + 182758215997616*x, 6321363049*x^3 + 1630911666642*x^2 + 140258403331212*x, 11688200277601*x^2, 502592611936843*x) of Univariate Polynomial Ring in x over Ring of integers modulo 21611482313284249

        We check that all generators vanish on A::

            sage: all(f(A) == 0 for f in I.gens())
            True

        We check that the ideal is preserved by conjugation::

            sage: B = random_matrix(Zmod(43^10), 6)
            sage: while not B.is_unit():
            ....:     B = random_matrix(Zmod(43^10), 6)
            sage: J = (~B * A * B).minpoly_ideal()
            sage: J == I
            True

        We check that the polynomials for random matrices vanish for various moduli::

            sage: for N in [5, 625, 36, 2^24, 2^6*3^9, 2^63-1]:
            ....:     for n in [3, 6]:
            ....:         A = random_matrix(Zmod(N), n, n)
            ....:         I = A.minpoly_ideal()
            ....:         assert all(f(A) == 0 for f in I.gens())
            ....:         assert I.gens()[0].is_monic()
        """
        if not self.is_square():
            raise ValueError("Minimal polynomial ideal not defined for non-square matrices")
        if proof is None:
            proof = linalg_proof()
        R = self.base_ring()
        Rx = PolynomialRing(R, var, implementation="FLINT")
        key = "minpoly_ideal"
        ans = self.fetch(key)
        if ans is not None:
            if ans.ring().variable_name() != var:
                ans = Rx.ideal([Rx(g) for g in ans.gens()])
            return ans
        n = self.nrows()
        F = R.factored_order()
        P = self.parent()
        cdef mp_limb_t piv, c, N = self.base_ring().order()
        cdef Matrix_nmod_dense C
        cdef Polynomial_zmod_flint f, mpoly
        cdef Py_ssize_t i, j, jj, k, d
        if len(F) == 1:
            p, e = F[0]
            # We start with the minimal polynomial mod p
            S = Zmod(p)
            Sx = PolynomialRing(S, var, implementation="FLINT")
            C = self.change_ring(S)
            mpoly = Sx()
            sig_on()
            nmod_mat_minpoly(&mpoly.x, C._matrix)
            sig_off()
            if e == 1:
                ans = Rx.ideal(mpoly)
                self.cache(key, ans)
                return ans
            d = mpoly.degree()
            if d == n:
                # The minimal polynomial is the same as the characteristic polynomial
                ans = Rx.ideal(self.charpoly(var=var))
                self.cache(key, ans)
                return ans
        # We pick a random vector and iteratively multiply by the matrix n times
        cdef Matrix_nmod_dense v = self._new(n, 1), B = self._new(n+1, n+1)
        pows = [self.parent().identity_matrix()]
        def check(polys, pows):
            d = polys[0].degree()
            while d >= len(pows):
                pows.append(pows[-1] * self)
            for f in polys:
                C = self._new(n, n)
                for k in range(f.degree() + 1):
                    C += f[k] * pows[k]
                if C:
                    return False
            return True
        while True:
            for i in range(n):
                c = randrange(N)
                nmod_mat_set_entry(v._matrix, i, 0, c)
                nmod_mat_set_entry(B._matrix, i, n, c)
            # It would be nice to do the kernel computation in parallel with computing self*v so that we know if we can stop early
            for j in range(n-1, -1, -1):
                v = self * v
                for i in range(n):
                    nmod_mat_set_entry(B._matrix, i, j, nmod_mat_get_entry(v._matrix, i, 0))
            # _right_kernel_matrix caches echelon form, and we're changing B
            B.clear_cache()
            C = B._right_kernel_matrix()
            polys = []
            # scan for last nonzero row
            for i in range(n, -1, -1):
                for j in range(n):
                    piv = nmod_mat_get_entry(C._matrix, i, j)
                    if piv:
                        f = Rx()
                        nmod_poly_fit_length(&f.x, n+1-j)
                        for jj in range(n+1-j):
                            nmod_poly_set_coeff_ui(&f.x, jj, nmod_mat_get_entry(C._matrix, i, n-jj))
                        polys.append(f)
                        break
                if piv == 1:
                    # reached the monic generator, so try to return an answer
                    polys.reverse()
                    # We only know that the polynomials vanish on a random vector.
                    # We need to check that they actually vanish on the matrix
                    if proof and check(polys, pows) or not proof:
                        ans = Rx.ideal(polys)
                        if proof: self.cache(key, ans)
                        return ans
                    break

    def minpoly(self, var='x', proof=None, **kwds):
        """
        Returns the minimal polynomial of this matrix.

        Note that, in general, the minimal polynomial of a matrix over a ring with composite modulus is not well defined.

        INPUT:

        - ``var`` -- the variable name for the polynomial ring

        - ``proof`` -- boolean, whether to check that the answer computed by using a minimal polynomial is correct.  If not specificied, uses the linear alegbra default proof state.  Note that if the modulus is composite and divisible by small primes the probability of an incorrect result is substantial.  The result is not cached when proof is ``False``, so this function can be called multiple times to get a desired level of certainty.

        EXAMPLES::

            sage: A = matrix(Zmod(36), 6, 6, 1)
            sage: A.minpoly()
            x + 35

            sage: A = matrix(Zmod(43^4), 4, 4, [1547380, 1211593, 1900136, 2872531, 3010350, 3130824, 2590437, 420249, 1236104, 1628956, 1300526, 145927, 1817994, 1412525, 1305313, 2775020])
            sage: f = A.minpoly(); f
            x^4 + 1502653*x^3 + 2556099*x^2 + 595635*x + 865202
            sage: f(A) == 0
            True

            sage: A = matrix(Zmod(36), 2, [6, 0, 0, -6])
            sage: A.minpoly()
            Traceback (most recent call last)
            ...
            ValueError: Matrix does not have a minimal polynomial; try minpoly_ideal

        We build matrices with a given minimal polynomial, then check that the correct result is returned::

           sage: R.<x> = Zmod(2^63-1)[]
           sage: f = x^3 + R.random_element(2)
           sage: g = x^2 + R.random_element(1)
           sage: A = block_diagonal_matrix([companion_matrix(f), companion_matrix(g), companion_matrix(f*g)])
           sage: B = random_matrix(Zmod(2^63-1), 10, 10)
           sage: while not B.is_unit():
           ....:     B = random_matrix(Zmod(2^63-1), 10, 10)
           sage: C = ~B * A * B
           sage: C.minpoly() == f*g
        """
        I = self.minpoly_ideal(var, proof=proof)
        gens = I.gens()
        if len(gens) != 1:
            raise ValueError("Matrix does not have a minimal polynomial; try minpoly_ideal")
        return gens[0]

    cdef swap_columns_c(self, Py_ssize_t c1, Py_ssize_t c2):
        """
        Swap two columns using FLINT.
        """
        nmod_mat_swap_cols(self._matrix, NULL, c1, c2)

    cdef swap_rows_c(self, Py_ssize_t r1, Py_ssize_t r2):
        """
        Swap two rows using FLINT.
        """
        nmod_mat_swap_rows(self._matrix, NULL, r1, r2)

    def _solve_right_modn(self, B, check=True):
        """
        Solve the matrix equation `A X = B`.

        sage: for N in [5, 625, 36, 2^24, 2^6*3^9, 2^63-1]: # indirect doctest
        ....:     for n in [3, 6]:
        ....:         for m in [4, 5, 6]:
        ....:             A = random_matrix(Zmod(N), m, n)
        ....:             for k in [1, 2, 4]:
        ....:                 X0 = random_matrix(Zmod(N), n, k)
        ....:                 B = A * X0
        ....:                 X = A.solve_right(B)
        ....:                 assert A * X == B
        """
        b_is_vec = is_Vector(B)
        cdef Py_ssize_t i, j, ii, jj, n, Cm, An = self.ncols()
        cdef mp_limb_t piv, q, r, tmp, entry, N = self.base_ring().order(), Ninv = n_preinvert_limb(N)
        cdef double pivinv
        if b_is_vec:
            n = 1
        else:
            n = B.ncols()
        cdef Matrix_nmod_dense C, X
        R = self.base_ring()
        X = self._new(An, n)
        if R.is_field():
            C = B.column() if b_is_vec else B
            sig_on()
            nmod_mat_can_solve(X._matrix, self._matrix, C._matrix)
            sig_off()
        else:
            # We use Howell form to solve the system
            C = self.augment(B)
            if C.nrows() < C.ncols():
                C = C.stack(C.new_matrix(nrows=C.ncols() - C.nrows()))
            C.echelonize()
            # save the echelon form if not already cached
            if self.fetch('echelon_form') is None:
                E = C.submatrix(max(self.ncols(), self.nrows()), self.ncols())
                self.cache('echelon_form', E)
            Cm = C.nrows()
            for i in range(Cm - 1, -1, -1):
                for j in range(An):
                    piv = nmod_mat_get_entry(C._matrix, i, j)
                    if not piv:
                        continue
                    elif piv != 1:
                        pivinv = n_precompute_inverse(piv)
                    for jj in range(n):
                        if piv == 1:
                            q = nmod_mat_get_entry(C._matrix, i, An + jj)
                        else:
                            entry = nmod_mat_get_entry(C._matrix, i, An + jj)
                            r = n_divrem2_precomp(&q, entry, piv, pivinv)
                            if r:
                                raise ValueError("matrix equation has no solutions")
                        nmod_mat_set_entry(X._matrix, j, jj, q)
                        # We now update the right augmented block based on the entries above the pivot.
                        for ii in range(i - 1, -1, -1):
                            entry = n_mulmod2_preinv(q, nmod_mat_get_entry(C._matrix, ii, j), N, Ninv)
                            entry = n_submod(nmod_mat_get_entry(C._matrix, ii, An + jj), entry, N)
                            nmod_mat_set_entry(C._matrix, ii, An + jj, entry)
                    break
                else:
                    # full zero row on the left augmented block; check that it's also zero on the right
                    for jj in range(n):
                        if nmod_mat_get_entry(C._matrix, i, An + jj):
                            raise ValueError("matrix equation has no solution")
        if b_is_vec:
            # Convert back to a vector
            return X.column(0)
        else:
            return X

    def tranpose(self):
        """
        Return the transpose of this matrix, without changing this matrix.

        EXAMPLES::

            sage: matrix(Zmod(36), 2, 2, [0, 1, 0, 0]).transpose()
            [0 0]
            [1 0]
        """
        cdef Matrix_nmod_dense M = self._new(self._nrows, self._ncols)
        sig_on()
        nmod_mat_transpose(M._matrix, self._matrix)
        sig_off()
        if self._subdivisions is not None:
            row_divs, col_divs = self.subdivisions()
            M.subdivide(col_divs, row_divs)
        return M


    def echelonize(self, algorithm="default", **kwds):
        """
        Transform this matrix into echelon form in place, over the same base ring.

        .. NOTE::

            Over a field, transforms to standard reduced row echelon form.
            Otherwise, uses Howell form, which may increase the number of nonzero rows.  In this case, we require that the number of rows is at least the number of columns.

        INPUT:

        - ``algorithm`` -- ignored (always uses FLINT)

        - ``transformation`` -- boolean. Whether to return a matrix `T` so that left multiplication by `T` transforms the original matrix into echelon form.

        OUTPUT:

        This matrix is put into echelon form.  Nothing is returned unless
        the keyword option ``transformation=True`` is specified, in
        which case the transformation matrix is returned.

        EXAMPLES::

            sage: A = matrix(Zmod(625), 4, 3, [[404, 355, 133], [375, 482, 448], [506, 115,  77], [370, 384, 66]])
            sage: A
            [404 355 133]
            [375 482 448]
            [506 115  77]
            [370 384  66]
            sage: A.echelonize()
            sage: A
            [1 0 2]
            [0 1 4]
            [0 0 5]
            [0 0 0]
        """
        self.check_mutability()
        self.clear_cache()
        cdef bint transformation = 'transformation' in kwds and kwds['transformation']
        cdef Matrix_nmod_dense aug, trans
        cdef Py_ssize_t i, j, m, n
        if self._parent._base.is_field():
            if transformation:
                m = self.nrows()
                n = self.ncols()
                aug = self.augment(identity_matrix(self.base_ring(), m))
                trans = self._new(m, m)
                sig_on()
                nmod_mat_rref(aug._matrix)
                sig_off()
                for i in range(m):
                    for j in range(n):
                        nmod_mat_set_entry(self._matrix, i, j, nmod_mat_get_entry(aug._matrix, i, j))
                    for j in range(m):
                        nmod_mat_set_entry(trans._matrix, i, j, nmod_mat_get_entry(aug._matrix, i, n+j))
            else:
                sig_on()
                rank = nmod_mat_rref(self._matrix)
                sig_off()
                self.cache('rank', Integer(rank))
        else:
            trans = self._howell_form(transformation)
        self.cache('in_echelon_form', True)
        if transformation:
            return trans

    def echelon_form(self, algorithm="default", **kwds):
        """
        Return the echelon form of this matrix.

        Note that when not over a field, the echelon form can have more rows than the input.

        INPUT:

        - ``algorithm`` -- ignored (always uses FLINT)

        - ``transformation`` -- boolean. Whether to also return the
            transformation matrix.

        OUTPUT:

        If the base ring is a field, the reduced row echelon form
        of this matrix, with the same dimensions.  If the base ring
        is not a field, the Howell form of this matrix, which will
        have extra rows added so that the number of rows is at least
        the number of columns.

        This matrix is not changed by this function, and the returned
        echelon form is immutable.  Use :meth:`echelonize` to transform
        this matrix to echelon form in place.

        If the optional parameter ``transformation=True`` is
        specified, the output consists of a pair `(E,T)` of matrices
        where `E` is the echelon form of this matrix and `T` is the
        transformation matrix.

        EXAMPlES::

            sage: A = matrix(Zmod(625), 4, 3, [[404, 355, 133], [375, 482, 448], [506, 115,  77], [370, 384, 66]])
            sage: A.echelon_form()
            [1 0 2]
            [0 1 4]
            [0 0 5]
            [0 0 0]
            sage: E, T = A.echelon_form(transformation=True); T
            [  2  17  23 564]
            [  4  22 429 488]
            [  3   4 188 543]
            [  5   8 510 316]
            sage: E == T * A
            True
        """
        transformation = kwds.get('transformation')
        key='echelon_form'
        E = self.fetch(key)
        if E is not None:
            if not transformation:
                return E
            T = self.fetch('echelon_transformation')
            if T is not None:
                return E, T
        if self._nrows < self._ncols and not self.base_ring().is_field():
            M = self._new(self._ncols - self._nrows, self._ncols)
            E = M.stack(self)
        else:
            E = self.__copy__()
        T = E.echelonize(algorithm, **kwds)
        if T is not None:
            self.cache('echelon_transformation', T)
        E.set_immutable()
        self.cache(key, E)
        if transformation:
            return E, T
        else:
            return E

    def _pivots(self):
        """
        When not over a field, it is possible to have leading entries
        that are zero divisors.  This function distinguishes
        between such pivots and pivots that are 1.

        OUTPUT:

        - a list of columns containing a leading 1

        - a list of columns containing a leading nonzero zero divisor

        EXAMPLES::

            sage: R = Zmod(625)
            sage: A = matrix(Zmod(625), 3, 4, [1, 2, 3, 4, 0, 5, 5, 6, 0, 0, 0, 25])
            sage: A
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0 25]
            sage: A._pivots()
            ((0,), (1, 3))
        """
        key = 'pivots'
        ans = self.fetch(key)
        if ans is not None:
            return ans
        cdef Matrix_nmod_dense E
        cdef Py_ssize_t i, j, k = 0
        # howell form has all the zero rows at the bottom
        E = self.echelon_form()
        p = []
        zdp = []
        for i in range(E._nrows):
            for j in range(k, E._ncols):
                if nmod_mat_get_entry(E._matrix, i, j) != 0:  # nonzero position
                    if nmod_mat_get_entry(E._matrix, i, j) == 1:
                        p.append(j)
                    else:
                        zdp.append(j) # is a zero divisor
                    k = j+1  # so start at next position next time
                    break
        ans = (tuple(p), tuple(zdp))
        self.cache(key, ans)
        return ans

    def pivots(self):
        """
        Return the columns containing a leading 1 in the echelon form of this matrix.

        When the base ring is not a field, there may be other rows with leading entries
        a zero divisor.  These columns are available using the :meth:`_pivots` method.

        EXAMPLES::

            sage: R = Zmod(625)
            sage: A = matrix(Zmod(625), 3, 4, [1, 2, 3, 4, 0, 5, 5, 6, 0, 0, 0, 25])
            sage: A
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0 25]
            sage: A.pivots()
            (0,)
        """
        return self._pivots()[0]

    def rank(self):
        """
        Return the rank of this matrix.

        When the base ring is not a field, this is defined as the number of leading 1 pivots.
        This choice has the benefit that a square matrix will be invertible exactly when
        the rank is the same as the number of rows.

        .. SEEALSO:

        - :meth:`pivots`

        - :meth:`_pivots`

        EXAMPlES::

            sage: entries = [1, 2, 3, 4, 0, 5, 5, 6, 0, 0, 0, 25]
            sage: A = matrix(GF(17), 3, 4, entries)
            sage: A
            [1 2 3 4]
            [0 5 5 6]
            [0 0 0 8]
            sage: A.rank()
            3
            sage: A = matrix(Zmod(625), 3, 4, entries)
            sage: A
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0 25]
            sage: A.rank()
            1
        """
        key = 'rank'
        ans = self.fetch(key)
        if ans is not None:
            return ans
        if not self._parent._base.is_field():
            ans = len(self.pivots())
        else:
            p = self.fetch('pivots')
            if ans is not None:
                ans = len(p[0])
            else:
                sig_on()
                ans = nmod_mat_rank(self._matrix)
                sig_off()
        self.cache(key, ans)
        return ans


    def determinant(self, algorithm=None):
        """
        Return the determinant of this matrix.

        EXAMPLES::

            sage: A = matrix(ZZ, 3, [1, 6, 11, -2, 3, 1, 3, 1, 0])
            sage: all(A.change_ring(Zmod(N)).det() == A.det() for N in range(2, 100))

            sage: A = random_matrix(Zmod(36), 6)
            sage: all(A.det()^n == (A^n).det() for n in range(2, 10))
            True
            sage: A.det("charpoly") == A.det("lift")
            True
            sage: B = random_matrix(Zmod(36), 6)
            sage: while not B.is_unit():
            ....:     B = random_matrix(Zmod(36), 6)
            sage: A.det() == (~B * A * B).det()
            True
        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        key = 'det'
        d = self.fetch(key)
        if d is not None:
            return d
        if algorithm is None:
            algorithm = "flint" if self.base_ring().is_field() else "charpoly"

        # If charpoly known, then det is easy.
        f = self.fetch('charpoly')
        if f is not None:
            d = f[0]
            if self._ncols % 2:
                d = -d
        elif algorithm == "flint":
            sig_on()
            d = nmod_mat_det(self._matrix)
            sig_off()
        elif algorithm == "charpoly":
            f = self.charpoly()
            d = f[0]
            if self._ncols % 2:
                d = -d
        elif algorithm == "lift":
            d = self.lift_centered().det()
        else:
            raise ValueError("Unknown algorithm '%s'" % algorithm)

        d = self._parent._base(d)
        self.cache(key, d)
        return d

    def trace(self):
        """
        Return the trace of this matrix, which is the sum of the diagonal entries.

        The input must be square.

        EXAMPLES::

            sage: A = matrix(ZZ, 3, [1, 6, 11, -2, 3, 1, 3, 1, 0])
            sage: all(A.change_ring(Zmod(N)).trace() == A.trace() for N in range(2, 100))
            True

        If the characteristic polynomial is known, negating the second coefficient
        is used to recover the trace::

            sage: A = matrix(Zmod(36), 3, [12, 12, 10, 6, 23, 24, 35, 24, 30])
            sage: A.charpoly()
            x^3 + 7*x^2 + 4*x + 22
            sage: A.trace()
            29
        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        key = 'trace'
        ans = self.fetch(key)
        if ans is not None:
            return ans
        f = self.fetch('charpoly')
        if f is not None:
            ans = -f[self._nrows - 1]
        else:
            sig_on()
            ans = nmod_mat_trace(self._matrix)
            sig_off()
            ans = self.base_ring()(ans)
        self.cache(key, ans)
        return ans

    def strong_echelonize(self):
        """
        Transform this matrix in place into its strong echelon form.

        The strong echelon form obtained from the Howell form by permuting rows:
        rather than having all zero rows at the bottom,
        they are placed so that the pivots occur on the diagonal.

        The input must have at least as many rows as columns.

        EXAMPLES::

            sage: entries = [1, 2, 3, 4, 0, 5, 5, 6, 0, 0, 0, 25, 0, 0, 0, 0]
            sage: A = matrix(Zmod(625), 4, entries); A
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0 25]
            [ 0  0  0  0]
            sage: A.strong_echelonize(); A
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0  0]
            [ 0  0  0 25]
        """
        if self._nrows < self._ncols:
            raise ValueError("self must must have at least as many rows as columns.")
        self.check_mutability()
        self.clear_cache()
        sig_on()
        nmod_mat_strong_echelon_form(self._matrix)
        sig_off()

    def strong_echelon_form(self):
        """
        Strong echelon form of this matrix.

        This is obtained from the Howell form (the form used
        as echelon form when not over a field) by permuting rows:
        rather than having all zero rows at the bottom,
        they are placed so that the pivots occur on the diagonal.

        The result will always have at least as many rows as columns, with zero
        rows added if necessary to accomplish this.

        EXAMPLES::

            sage: entries = [624, 353, 352, 497, 182, 374, 556, 365, 271, 557, 203, 327]
            sage: A = matrix(Zmod(625), 3, 4, entries)
            sage: A.strong_echelon_form()
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0  0]
            [ 0  0  0 25]
        """
        key = 'strong_echelon_form'
        ans = self.fetch(key)
        if ans is not None:
            return ans
        if self._nrows < self._ncols:
            M = self._new(self._ncols - self._nrows, self._ncols)
            ans = self.stack(M)
        else:
            ans = self.__copy__()
        ans.strong_echelonize()
        self.cache(key, ans)
        return ans

    def howellize(self, transformation=False):
        r"""
        Transform this matrix in place into its Howell form.

        See :meth:`howell_form` for the definition.

        Since the Howell form may have more rows than `A`, to transform `A` in place
        we require that it has at least as many rows as columns.

        EXAMPLES::

            sage: entries = [624, 353, 352, 497, 182, 374, 556, 365, 271, 557, 203, 327, 181, 102, 283, 237]
            sage: A = matrix(Zmod(625), 4, entries)
            sage: A.howellize(); A
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0 25]
            [ 0  0  0  0]
        """
        if self._nrows < self._ncols:
            raise ValueError("matrix must have at least as many rows as columns.")
        self.check_mutability()
        cdef Matrix_nmod_dense aug, trans
        cdef Py_ssize_t i, j, m, n
        self.clear_cache()
        if transformation:
            m = self.nrows()
            n = self.ncols()
            aug = self.augment(identity_matrix(self.base_ring(), m)).stack(self._new(n, m + n))
            trans = self._new(m, m)
            sig_on()
            nmod_mat_howell_form(aug._matrix)
            sig_off()
            for i in range(m):
                for j in range(n):
                    nmod_mat_set_entry(self._matrix, i, j, nmod_mat_get_entry(aug._matrix, i, j))
                for j in range(m):
                    nmod_mat_set_entry(trans._matrix, i, j, nmod_mat_get_entry(aug._matrix, i, n+j))
            return trans
        else:
            sig_on()
            nmod_mat_howell_form(self._matrix)
            sig_off()

    def howell_form(self):
        """
        Return the Howell form of this matrix.

        The Howell form is an echelon form with a few additional properties
        that make it unique and useful for solving equations modulo `N`.

        For a matrix `A`, let `S(A)` be the module spanned by the rows, and
        `S_j(A)` the submodule consisting of vectors with the first `j` entries
        zero.  Recall that elements `a, b` of a ring `R` are associates if there
        is a unit `u \in R` with `a = ub`; this is an equivalence relation.
        We can choose a set `X` of representatives in `\ZZ/N\ZZ` by taking
        all products of the primes dividing `N`, raised to arbitrary powers.

        The Howell form of `A` is a matrix `H` with the same number of columns
        as `A` so that

        - it is an echelon form: the zero rows are at the bottom, and the leading
          entries (pivots) in each row each occur down and to the right of the previous.

        - Each pivot lies in `X`, and the entries above a pivot are smaller
          (as in Hermite normal form).

        - If `(i, j_i)` is the location of a pivot, then rows `i+1`, ... generate
          `S_{j_i}(A)`.

        if `r` is the number of nonzero rows of `H`, then the first `r` rows
          of `H` are nonzero

        .. NOTE::

            The Howell form will always have at least as many rows as columns.

        EXAMPLES::

            sage: entries = [624, 353, 352, 497, 182, 374, 556, 365, 271, 557, 203, 327]
            sage: A = matrix(Zmod(625), 3, 4, entries)
            sage: A.howell_form()
            [ 1  2  3  4]
            [ 0  5  5  6]
            [ 0  0  0 25]
            [ 0  0  0  0]

        The following example illustrates why the number of rows needs to be increased::

            sage: A = matrix(Zmod(625), [125, 25, 5, 1])
            sage: A.howell_form()
            [125  25   5   1]
            [  0 125  25   5]
            [  0   0 125  25]
            [  0   0   0 125]
        """
        key='howell_form'
        ans = self.fetch(key)
        if ans is not None:
            return ans

        if self._nrows < self._ncols:
            M = self._new(self._ncols - self._nrows, self._ncols)
            ans = self.stack(M)
        else:
            ans = self.__copy__()
        ans.howellize()
        self.cache(key, ans)
        return ans

    def __pow__(sself, n, dummy):
        """
        Exponentiation.

        EXAMPLES::

        
        """
        cdef Matrix_nmod_dense self = <Matrix_nmod_dense?>sself

        if dummy is not None:
            raise ValueError
        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix")

        cdef unsigned long e
        cdef long e_sgn
        cdef int err

        if integer_check_long_py(n, &e_sgn, &err) and not err:
            if e_sgn < 0:
                return (~self) ** (-n)
            e = <unsigned long>e_sgn
        else:
            if not isinstance(n, Integer):
                try:
                    n = Integer(n)
                except TypeError:
                    from sage.symbolic.expression import Expression
                    if isinstance(n, Expression):
                        from sage.matrix.matrix2 import _matrix_power_symbolic
                        return _matrix_power_symbolic(self, n)
                    else:
                        raise NotImplementedError("the given exponent is not supported")
            if mpz_sgn((<Integer>n).value) < 0:
                return (~self) ** (-n)

            if mpz_fits_ulong_p((<Integer>n).value):
                e = mpz_get_ui((<Integer>n).value)
            else:
                return generic_power(self, n)

        if e == 0:
            return self._parent.identity_matrix()
        if e == 1:
            return self

        cdef Matrix_nmod_dense M = self._new(self._nrows, self._ncols)
        sig_on()
        nmod_mat_pow(M._matrix, self._matrix, e)
        sig_off()
        return M

    def _right_kernel_matrix(self, zero_divisors_are_pivots=False):
        cdef Matrix_nmod_dense X, ans, E
        cdef Py_ssize_t i, j, k, l
        cdef long x, y, s
        E = self.echelon_form()
        if self._parent._base.is_field():
            # nmod_mut_nullspace will do this regardless
            # so we are better off to start in echelon form to have the rank
            X = self._new(self.ncols, self._nrows - self.rank())
            ans = self._new(self._nrows - self.rank(), self._ncols)
            sig_on()
            nmod_mat_nullspace(X._matrix, E._matrix) # columns of X form a basis
            nmod_mat_transpose(ans._matrix, X._matrix)
            sig_off()
            return ans
        else:
            zero = self._parent._base.zero()
            one = self._parent._base.one()
            p, zdp = map(set, self._pivots())
            set_pivots = p.union(zdp)
            pivot = sorted(set_pivots)
            N = mpz_get_si(self._modulus.sageInteger.value)
            basis = []
            k = 0
            for j in range(self._ncols):
                if j in p:
                    k += 1
                    continue
                v = [zero] * self._ncols
                i = k
                if j in zdp:
                    k += 1
                    if zero_divisors_are_pivots:
                        continue
                    v[j] = self._parent._base(N//nmod_mat_get_entry(E._matrix, i, j))
                else:
                    v[j] = one

                # figure out the remaining coefficients of v
                # note that v might need to be rescaled several times
                for l in reversed(range(i)):
                    x = nmod_mat_get_entry(E._matrix, l, pivot[l])
                    # solve
                    # v[pivot[l]] E[l, pivot[l]]  + y*sum(E[l,m] * v[m] for m in range(pivot[l] + 1, j)) = 0
                    s = sum(v[m]*nmod_mat_get_entry(E._matrix, l, m) for m in range(pivot[l] + 1, j + 1)) % N
                    if s % x != 0: # make sure we can work mod N/x
                        y = x // gcd(s, x)
                        s *= y # now s is divisible by x
                        for m in range(pivot[l] + 1, j + 1):
                            v[m] *= y
                        assert v[j] % N != 0
                    # QUESTION: this is correct modulo N/x, does one need to consider the various lifts?
                    # FIXME, this feels wrong
                    v[pivot[l]] = self._parent._base(-s//x)
                basis.append(v)
            # FIXME, this feels wrong
            ans = self._new(len(basis), self._ncols)
            ma = MatrixArgs_init(ans._parent, basis)
            for t in ma.iter(False, True): #????
                se = <SparseEntry>t
                x = <long>se.entry
                nmod_mat_set_entry(ans._matrix, se.i, se.j, x)
            if zero_divisors_are_pivots:
                return ans
            else:
                return ans.howell_form()


    # transpose (Edgar) x
    # nmod_mat_pow (Edgar) x
    # nmod_mat_trace (Edgar) X
    # rank and det (only primes) (Edgar) X
    # rank generic (Edgar)
    # right_kernel_matrix (nmod_mat_nullspace) (Edgar) x
    # row reduction (nmod_mat_rref) (Edgar) ~
    # richcmp x
