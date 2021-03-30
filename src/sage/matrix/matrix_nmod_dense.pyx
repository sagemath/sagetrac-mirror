# distutils: extra_compile_args = -D_XPG6
# flags chosen from libs/flint/fmpz_poly.pyx
"""
FLINT nmod_mat class wrapper

AUTHORS:

- Edgar Costa (2021) Initial version.
"""

from cpython.sequence cimport *

from cysignals.signals cimport sig_on, sig_str, sig_off

from sage.arith.all import gcd
from sage.arith.power cimport generic_power
from sage.arith.long cimport integer_check_long_py
from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport Element, Matrix
from sage.libs.flint.nmod_mat cimport *
from sage.libs.flint.nmod_poly cimport nmod_poly_set, nmod_poly_set_coeff_ui
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
)
from sage.rings.polynomial.polynomial_zmod_flint cimport Polynomial_zmod_flint
from sage.libs.gmp.mpz cimport mpz_sgn,  mpz_fits_ulong_p, mpz_get_ui, mpz_get_si
from sage.rings.integer cimport Integer
from sage.structure.factorization import Factorization

from .args cimport SparseEntry, MatrixArgs_init

import sage.matrix.matrix_space as matrix_space
from sage.rings.finite_rings.integer_mod cimport IntegerMod_abstract, IntegerMod_int, IntegerMod_int64, IntegerMod_gmp
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.structure.richcmp cimport rich_to_bool, Py_EQ, Py_NE

cdef class Matrix_nmod_dense(Matrix_dense):
    ########################################################################
    # LEVEL 1 helpers:
    #   These function support the implementation of the level 1 functionality.
    ########################################################################
    def __cinit__(self, parent, *args, **kwds):
        self._modulus = parent._base._pyx_order
        sig_str("FLINT exception")
        nmod_mat_init(self._matrix, self._nrows, self._ncols, mpz_get_ui(self._modulus.sageInteger.value))
        sig_off()



    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        self._parent = parent
        ma = MatrixArgs_init(parent, entries)
        cdef long z
        for t in ma.iter(coerce, True): #????
            se = <SparseEntry>t
            z = <long>se.entry
            nmod_mat_set_entry(self._matrix, se.i, se.j, z)

    def _print(self):
        #FIXME: For debugging; remove when ready
        nmod_mat_print_pretty(self._matrix)

    def __dealloc__(self):
        nmod_mat_clear(self._matrix)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x):
        e = self._modulus.element_class()
        cdef mp_limb_t ivalue
        if e is IntegerMod_int:
            ivalue = (<IntegerMod_int>x).ivalue
        elif e is IntegerMod_int64:
            ivalue = (<IntegerMod_int64>x).ivalue
        else:
            ivalue = mpz_get_ui((<IntegerMod_gmp>x).value)
        self.set_unsafe_si(i, j, ivalue)

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
        ans._parent = P
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
        e = self._modulus.element_class()
        cdef mp_limb_t ivalue
        if e is IntegerMod_int:
            ivalue = (<IntegerMod_int>right).ivalue
        elif e is IntegerMod_int64:
            ivalue = (<IntegerMod_int64>right).ivalue
        else:
            ivalue = mpz_get_ui((<IntegerMod_gmp>right).value)
        nmod_mat_scalar_mul(M._matrix, self._matrix, ivalue)
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
        else:
            if algorithm == "lift":
                return (~self.lift_centered()).change_ring(R)
            elif algorithm == "crt":
                # The modulus is small, so factoring is feasible.
                # We find inverses modulo each prime dividing the modulus, lift p-adically, then CRT the results together
                M = self._new(n, n)
                F = R.factored_order()
                maxe = max(fac[1] for fac in F)
                lift_required = (maxe > 1)
                if lift_required:
                    nlifts = (maxe - 1).nbits()
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
                            for k in reversed(range(maxe)):
                                lift_mods[k] *= p**e
                                e = (e + 1) // 2
                    if lift_required:
                        for k in range(maxe):
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

    def hessenberg_form(self):
        B = self.__copy__()
        B.hessenbergize()
        return B

    def hessenbergize(self):
        """
        """
        R = self._parent._base
        F = R.factored_order()
        if len(F) != 1:
            raise ValueError("Hessenbergize only supported for prime power modulus")
        if not self.is_square():
            raise TypeError("self must be square")

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
                algorithm = "FLINT"
            else:
                algorithm = "crt"
        if algorithm == "FLINT":
            f = R[var]()
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
                    A = matrix_space.MatrixSpace(R, self._nrows, self._ncols, sparse=False)()
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
                f = R[var]()
                for j in range(n+1):
                    nmod_poly_set_coeff_ui(&f.x, j, nmod_mat_get_entry(c, n, j))
                return f
            finally:
                nmod_mat_clear(H)
                nmod_mat_clear(c)
        raise ValueError("Unknown algorithm '%s'" % algorithm)

    cdef swap_columns_c(self, Py_ssize_t c1, Py_ssize_t c2):
        nmod_mat_swap_cols(self._matrix, NULL, c1, c2)

    cdef swap_rows_c(self, Py_ssize_t r1, Py_ssize_t r2):
        nmod_mat_swap_rows(self._matrix, NULL, r1, r2)

    # Extra

    def tranpose(self):
        cdef Matrix_nmod_dense M = self._new(self._nrows, self._ncols)
        sig_on()
        nmod_mat_transpose(M._matrix, self._matrix)
        sig_off()
        return M


    def echelonize(self):
        #FIXME: add a warning regarding the output potentially having more rows than input
        """
        Echelon form in place
        """
        self.check_mutability()
        self.clear_cache()
        if self._parent._base.is_field():
            sig_on()
            rank = nmod_mat_rref(self._matrix)
            sig_off()
            self.cache('rank', rank)
        else:
            self._howell_form()

    def echelon_form(self):
        #FIXME: add a warning regarding the output potentially having more rows than input
        key='echelon_form'
        ans = self.fetch(key)
        if ans is not None:
            return ans
        if self._nrows < self._ncols:
            M = self._new(self._ncols - self._nrows, self._ncols)
            ans = M.stack(self)
        else:
            ans = self.__copy__()
        ans.echelonize()
        self.cache(key, ans)
        return ans


    def _pivots(self):
        key = 'pivots'
        ans = self.fetch(key)
        if ans is not None:
            return ans
        cdef Matrix_nmod_dense E
        cdef Py_ssize_t i, j, k
        # howell form has all the zero rows at the bottom
        E = self.echelon_form()
        p = []
        zdp = []
        k = 0
        for i in range(E._nrows):
            for j in range(k, E._ncols):
                if nmod_mat_get_entry(E._matrix, i, j) != 0:  # nonzero position
                    if nmod_mat_get_entry(E._matrix, i, j) == 1:
                        p.append(j)
                    else:
                        zdp.append(j) # is a zero divisor
                    k = j+1  # so start at next position next time
                    break
        ans = (p, zdp)
        self.cache(key, ans)
        return ans

    def pivots(self):
        return self._pivots()[0]


    def rank(self):
        key = 'rank'
        ans = self.fetch(key)
        if ans is not None:
            return ans
        if not self._parent._base.is_field():
            p = self.pivots()
            ans = len(p[0])
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


    def determinant(self):
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        key = 'det'
        ans = self.fetch(key)
        if ans is not None:
            return ans

        # If charpoly known, then det is easy.
        f = self.fetch('charpoly')
        if f is not None:
            c = f[0]
            if self._ncols % 2:
                c = -c
            d = self._coerce_element(c)
            self.cache('det', d)
            return d

        if self._parent._base.is_field():
            sig_on()
            ans = nmod_mat_det(self._matrix)
            sig_off()
        else:
            ans = self.apply_map(lambda x : x.lift_centered()).det()

        ans = self._parent._base(ans)
        self.cache(key, ans)
        return ans

    def trace(self):
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        key = 'trace'
        ans = self.fetch(key)
        if ans is not None:
            return ans
        sig_on()
        ans = nmod_mat_trace(self._matrix)
        sig_off()
        self.cache(key, ans)
        return ans

    def _strong_echelon_form(self):
        """
        In place strong echelon form of self
        """
        if self._nrows < self._ncols:
            raise ValueError("self must must have at least as many rows as columns.")
        self.check_mutability()
        self.clear_cache()
        sig_on()
        nmod_mat_strong_echelon_form(self._matrix)
        sig_off()

    def strong_echelon_form(self):
        #FXIME add warning
        """
        Strong echelon form of self
        """
        key='strong_echelon_form'
        ans = self.fetch(key)
        if ans is not None:
            return ans
        if self._nrows < self._ncols:
            M = self._new(self._ncols - self._nrows, self._ncols)
            ans = self.stack(M)
        else:
            ans = self.__copy__()
        ans._strong_echelon_form()
        self.cache(key, ans)
        return ans

    def _howell_form(self):
        """
        In place Howell form of self
        """
        if self._nrows < self._ncols:
            raise ValueError("self must have at least as many rows as columns.")
        self.check_mutability()
        cdef Matrix_nmod_dense M
        self.clear_cache()
        sig_on()
        nmod_mat_howell_form(self._matrix)
        sig_off()

    def howell_form(self):
        #FXIME add warning
        """
        Howell form of self
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
        ans._howell_form()
        self.cache(key, ans)
        return ans

    def __pow__(sself, n, dummy):
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
                # it is very likely that the following will never finish except
                # if self has only eigenvalues 0, 1 or -1.
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
            pivot = sorted(list(set_pivots))
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
                        y = x//gcd(s, x)
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


    # random matrix generation (David)
    # charpoly and minpoly (David)
    # solve (nmod_mat_can_solve) (David)



    # transpose (Edgar) x
    # nmod_mat_pow (Edgar) x
    # nmod_mat_trace (Edgar) X
    # rank and det (only primes) (Edgar) X
    # rank generic (Edgar)
    # right_kernel_matrix (nmod_mat_nullspace) (Edgar) x
    # row reduction (nmod_mat_rref) (Edgar) ~
    # richcmp x
