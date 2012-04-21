from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.polynomial.all import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.arith import binomial, bernoulli
from sage.modules.free_module_element import vector, zero_vector
from sage.matrix.matrix cimport Matrix
from sage.matrix.all import matrix
from sage.misc.prandom import random
from sage.functions.other import floor
from sage.structure.element cimport RingElement, Element
import operator
#from sage.modular.overconvergent.pollack.S0p import S0
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.padics.padic_capped_absolute_element cimport pAdicCappedAbsoluteElement
from sage.rings.padics.padic_capped_relative_element cimport pAdicCappedRelativeElement
from sage.rings.padics.padic_fixed_mod_element cimport pAdicFixedModElement
from sage.rings.integer cimport Integer
from sage.misc.misc import verbose, cputime

cdef extern from "zn_poly/zn_poly.h":
    pass
from sage.libs.flint.zmod_poly cimport *, zmod_poly_t
from sage.libs.flint.long_extras cimport *

from fund_domain import M2ZSpace
cdef long overflow = 1 << (4*sizeof(long)-1)
cdef long underflow = -overflow

include "stdsage.pxi"
include "cdefs.pxi"

def get_dist_classes(p, prec_cap, base, symk):
    if symk or p is None or base.is_field() or (isinstance(base, pAdicGeneric) and base.degree() > 1):
        return Dist_vector, WeightKAction_vector
    if 7*p**(prec_cap) < ZZ(2)**(4*sizeof(long)-1):
        return Dist_long, WeightKAction_long
    else:
        return Dist_vector, WeightKAction_vector

cdef class Dist(ModuleElement):
    cpdef normalize(self):
        raise NotImplementedError

    def scale(self,left):
        if isinstance(self, Dist_long) and isinstance(left, (Integer, pAdicCappedRelativeElement, pAdicCappedAbsoluteElement, pAdicFixedModElement)):
            return self._lmul_(left)
        R = left.parent()
        base = self.parent().base_ring()
        if base is R:
            return self._lmul_(left)
        elif base.has_coerce_map_from(R):
            return self._lmul_(base(left))
        else:
            from sage.categories.pushout import pushout
            new_base = pushout(base, R)
            V = self.parent().change_ring(new_base)
            scalar = new_base(left)
            return V([scalar * new_base(self.moment(i)) for i in range(self.precision_absolute())])

    def find_scalar(self, other, p, M, check=True):
        """
        Returns an alpha with other = self * alpha, or raises a ValueError.

        self should be nonzero
        """
        p = self.parent().prime()
        i = 0
        n = self.precision_absolute()
        a = self.moment(i)
        if p is None:
            while a == 0:
                if other.moment(i) != 0:
                    raise ValueError("not a scalar multiple")
                i += 1
                a = self.moment(i)
            alpha = other.moment(i) / a
            if check:
                i += 1
                while i < n:
                    if self.moment(i) != alpha * other.moment(i):
                        raise ValueError("not a scalar multiple")
                    i += 1
        else:
            v = a.valuation(p)
            while v >= n - i:
                i += 1
                a = self.moment(i)
                v = a.valuation(p)
            relprec = n - i - v
            alpha = other.moment(i) / a % p**(n-i)
            i += 1
            while i < n:
                a = self.moment(i)
                if check and a % p**(n-i) != alpha * other.moment(i) % p**(n-i):
                    raise ValueError("not a scalar multiple")
                v = a.valuation(p)
                if n - i - v > relprec:
                    relprec = n - i - v
                    alpha = other.moment(i) / a % p**(n-i)
        try:
            return self.parent().base_ring()(alpha)
        except ValueError:
            return alpha

    cpdef ModuleElement _rmul_(self, RingElement _left):
        return self._lmul_(_left)

    def valuation(self, p=None):
        r"""
        Returns the highest power of `p` which divides all moments of the distribution
        """
        if p is None:
            p = self.parent()._p
        return min([self.moment(a).valuation(p) for a in
            range(self.precision_absolute())])

    def specialize(self, new_base_ring=None):
        r"""
        Specializes to weight `k` -- i.e. projects to `Sym^k`
        """
        k=self.parent()._k
        if k < 0:
            raise ValueError("negative weight")
        if self.precision_absolute() < k+1:
            raise ValueError("not enough moments")
        V = self.parent().specialize(new_base_ring)
        new_base_ring = V.base_ring()
        return V([binomial(k, j) * (-1)**j * new_base_ring(self.moment(j)) for j in range(k+1)])

    def lift(self, p=None, M=None, new_base_ring=None):
        V = self.parent().lift(p, M, new_base_ring)
        k = V._k
        p = V.prime()
        M = V.precision_cap()
        R = V.base_ring()
        moments = [self.moment(j) * (-1)**j / binomial(k, j) for j in range(k+1)]
        zero = R(0)
        moments.extend([zero] * (M - k - 1))
        mu = V(moments)
        val = mu.valuation()
        if val < 0:
            # This seems unnatural
            print "scaling by %s^%s to keep things integral"%(p, -val)
            mu *= p**(-val)
        return mu

    def act_right(self,gam):
        r"""
        Return self|gam
        """
        return self.parent()._act(self, gam)

cdef class Dist_vector(Dist):
    def __init__(self,moments,parent,check=True):
        r"""
        A distribution is stored as a vector whose `j`-th entry is the `j`-th moment of the distribution.

        The `j`-th entry is stored modulo `p^(N-j)` where `N` is the total number of moments.
        (This is the accuracy that is maintained after acting by `\Gamma_0(p)`.)

        INPUTS:

        - ``parent`` -- a :class:`distributions.Distributions_Zp` instance
        - ``moments`` -- the list of moments given as a vector
        - ``check`` -- boolean, whether to coerce the vector into the appropriate module
        """
        Dist.__init__(self,parent)
        if check:
            base = parent.base_ring()
            try:
                M = len(moments)
            except TypeError:
                M = 1
                moments = [moments]
            moments = parent.approx_module(M)(moments)
        self.moments = moments

    cdef Dist_vector _new_c(self):
        cdef Dist_vector ans = PY_NEW(Dist_vector)
        ans._parent = self._parent
        return ans

    def _repr_(self):
        r"""
        Displays the moments of the distribution
        """
        self.normalize()
        if len(self.moments) == 1:
            return repr(self.moments[0])
        else:
            return repr(self.moments)

    def _rational_(self):
        """
        Convert to a rational number.
        
        EXAMPLES::

            sage: D = Distributions(0); d = D(4/3); d
            4/3
            sage: QQ(d)
            4/3

        We get a TypeError if there is more than 1 moment::

            sage: D = Distributions(1); d = D([1,2]); d
            (1, 2)
            sage: QQ(d)
            Traceback (most recent call last):
            ...
            TypeError: k must be 0
        """
        if len(self.moments) == 1:
            return QQ(self.moments[0])
        raise TypeError, "k must be 0"

    def moment(self,n):
        r"""
        Returns the `n`-th moment
        """
        return self.moments[n]

    cpdef ModuleElement _add_(self, ModuleElement _right):
        cdef Dist_vector ans = self._new_c()
        cdef Dist_vector right = _right
        if len(self.moments) == len(right.moments):
            ans.moments = self.moments + right.moments
        elif len(self.moments) < len(right.moments):
            ans.moments = self.moments + right.moments[:len(self.moments)]
        else:
            ans.moments = self.moments[:len(right.moments)] + right.moments
        return ans

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        cdef Dist_vector ans = self._new_c()
        cdef Dist_vector right = _right
        if len(self.moments) == len(right.moments):
            ans.moments = self.moments - right.moments
        elif len(self.moments) < len(right.moments):
            ans.moments = self.moments - right.moments[:len(self.moments)]
        else:
            ans.moments = self.moments[:len(right.moments)] - right.moments
        return ans

    cpdef ModuleElement _lmul_(self, RingElement right):
        cdef Dist_vector ans = self._new_c()
        ans.moments = self.moments * right
        return ans

    def precision_absolute(self):
        r"""
        Returns the number of moments of the distribution
        """
        return len(self.moments)

    cdef int _cmp_c_impl(left, Element right) except -2:
        return cmp(left.moments, right.moments)

    def zero(self):
        cdef Dist_vector ans = self._new_c()
        ans.moments = 0 * self.moments # could make this faster
        return ans

    cpdef normalize(self):
        r"""
        Reduce modulo `Fil^N` -- that is the `i`-th moment is reduced modulo `p^(N-i)`
        """
        p = self.parent()._p
        if not self.parent().is_symk(): # non-classical
            assert self.valuation() >= 0, "moments not integral in normalization"
            V = self.moments.parent()
            R = V.base_ring()
            n = self.precision_absolute()
            if isinstance(R, pAdicGeneric):
                self.moments = V([self.moment(i).add_bigoh(n-i) for i in range(n)])
            else:
                self.moments = V([self.moment(i)%(p**(n-i)) for i in range(n)])
        return self

    def reduce_precision(self, M):
        r"""
        Only holds on to `M` moments.
        """
        assert M<=self.precision_absolute(),"not enough moments"

        cdef Dist_vector ans = self._new_c()
        ans.moments = self.moments[:M]
        return ans

    def solve_diff_eqn(self):
        r"""
        Solves the difference equation.
        """
        # assert self.moments[0][0]==0, "not total measure zero"
        # print "result accurate modulo p^",self.moment(0).valuation(self.p)
        #v=[0 for j in range(0,i)]+[binomial(j,i)*bernoulli(j-i) for j in range(i,M)]
        M = self.precision_absolute()
        K = self.parent().base_ring().fraction_field()
        V = self.moments.parent()
        v = [K(0) for i in range(M)]
        for m in range(1,M):
            scalar = K(self.moment(m) / m) # division can take us out of the base ring.
            for j in range(m-1,M):
                v[j] += binomial(j,m-1) * bernoulli(j-m+1) * scalar
        cdef Dist_vector ans = self._new_c()
        ans.moments = V(v)
        return ans

    #def lift(self):
    #    r"""
    #    Increases the number of moments by `1`.
    #    """
    #    n = len(self.moments)
    #    if n >= self.parent()._prec_cap:
    #        raise ValueError("Cannot lift above precision cap")
    #    cdef Dist_vector ans = self._new_c()
    #    ans.moments = self.parent().approx_module(n+1)(list(self.moments) + [0])
    #    return ans

cdef class Dist_long(Dist):
    def __init__(self, moments, parent, check=True):
        Dist.__init__(self, parent)
        p = parent._p
        cdef int i
        if check:
            try:
                M = len(moments)
            except TypeError:
                M = 1
                moments = [moments]
            if M > 100 or 7*p**M > ZZ(2)**(4*sizeof(long) - 1): # 6 is so that we don't overflow on gathers
                raise ValueError("moments too long")
            
        for i in range(M):
            # TODO: shouldn't be doing the conversion to ZZ when check=False?
            self.moments[i] = ZZ(moments[i])
        self.prec = M
        self.prime_pow = <PowComputer_long?>parent.prime_pow
        #gather = 2**(4*sizeof(long)-1) // p**len(moments)
        #if gather >= len(moments):
        #    gather = 0
        #self._gather = gather

    cdef Dist_long _new_c(self):
        cdef Dist_long ans = PY_NEW(Dist_long)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        return ans

    def _repr_(self):
        self.normalize()
        if self.prec == 1:
            return repr(self.moments[0])
        else:
            return "(" + ", ".join([repr(self.moments[i]) for i in range(self.prec)]) + ")"

    cdef int quasi_normalize(self) except -1:
        cdef int i
        for i in range(self.prec):
            if self.moments[i] > overflow:
                self.moments[i] = self.moments[i] % self.prime_pow.small_powers[self.prec-i]
            elif self.moments[i] < underflow:
                self.moments[i] = self.moments[i] % self.prime_pow.small_powers[self.prec-i]
                self.moments[i] += self.prime_pow.small_powers[self.prec-i]

    cpdef normalize(self):
        cdef int i
        for i in range(self.prec):
            if self.moments[i] < 0:
                self.moments[i] = self.moments[i] % self.prime_pow.small_powers[self.prec-i]
                self.moments[i] += self.prime_pow.small_powers[self.prec-i]
            elif self.moments[i] >= self.prime_pow.small_powers[self.prec-i]:
                self.moments[i] = self.moments[i] % self.prime_pow.small_powers[self.prec-i]
        return self

    def moment(self, _n):
        cdef int n = _n
        if n < 0:
            n += self.prec
        if n < 0 or n >= self.prec:
            raise IndexError("list index out of range")
        return self.moments[n]

    cpdef ModuleElement _add_(self, ModuleElement _right):
        cdef Dist_long ans = self._new_c()
        cdef Dist_long right = _right
        ans.prec = self.prec if self.prec < right.prec else right.prec
        cdef int i
        # The following COULD overflow, but it would require 2^32
        # additions (on a 64-bit machine), since we restrict p^k to be
        # less than 2^31/7.
        for i in range(ans.prec):
            ans.moments[i] = self.moments[i] + right.moments[i]
        return ans

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        cdef Dist_long ans = self._new_c()
        cdef Dist_long right = _right
        ans.prec = self.prec if self.prec < right.prec else right.prec
        cdef int i
        # The following COULD overflow, but it would require 2^32
        # additions (on a 64-bit machine), since we restrict p^k to be
        # less than 2^31/7.
        for i in range(ans.prec):
            ans.moments[i] = self.moments[i] - right.moments[i]
        return ans

    cpdef ModuleElement _lmul_(self, RingElement _right):
        cdef Dist_long ans = self._new_c()
        ans.prec = self.prec
        self.quasi_normalize()
        cdef long scalar, absprec
        cdef Integer iright
        cdef pAdicCappedAbsoluteElement pcaright
        cdef pAdicCappedRelativeElement pcrright
        cdef pAdicFixedModElement pfmright
        if PY_TYPE_CHECK(_right, Integer):
            iright = <Integer>_right
            if mpz_fits_slong_p(iright.value):
                scalar = mpz_get_si(iright.value) % self.prime_pow.small_powers[self.prec]
            else:
                scalar = mpz_fdiv_ui(iright.value, self.prime_pow.small_powers[self.prec])
        elif PY_TYPE_CHECK(_right, pAdicCappedAbsoluteElement):
            pcaright = <pAdicCappedAbsoluteElement>_right
            if pcaright.absprec <= self.prec:
                scalar = mpz_get_si(pcaright.value)
            else:
                scalar = mpz_fdiv_ui(pcaright.value, self.prime_pow.small_powers[self.prec])
        elif PY_TYPE_CHECK(_right, pAdicCappedRelativeElement):
            pcrright = <pAdicCappedRelativeElement>_right
            absprec = pcrright.ordp + pcrright.relprec
            if pcrright.ordp < 0:
                raise NotImplementedError
            if absprec <= self.prec:
                scalar = mpz_get_si(pcrright.unit) * self.prime_pow.small_powers[pcrright.ordp]
            else:
                scalar = mpz_fdiv_ui(pcrright.unit, self.prime_pow.small_powers[self.prec - pcrright.ordp]) * self.prime_pow.small_powers[pcrright.ordp]
        elif PY_TYPE_CHECK(_right, pAdicFixedModElement):
            pfmright = <pAdicFixedModElement>_right
            scalar = mpz_get_si(pfmright.value)
        cdef int i
        for i in range(self.prec):
            ans.moments[i] = self.moments[i] * scalar
        ans.quasi_normalize()
        return ans

    def precision_absolute(self):
        return self.prec

    cdef int _cmp_c_impl(left, Element _right) except -2:
        cdef int i
        cdef Dist_long right = _right
        for i in range(left.prec):
            if left.moments[i] < right.moments[i]:
                return -1
            if left.moments[i] > right.moments[i]:
                return 1
        return 0

    def zero(self):
        cdef Dist_long ans = self._new_c()
        ans.prec = self.prec
        cdef int i
        for i in range(self.prec):
            ans.moments[i] = 0
        return ans

    def change_precision(self, M):
        if M > self.prec: raise ValueError("not enough moments")
        if M < 0: raise ValueError("precision must be non-negative")
        cdef Dist_long ans = self._new_c()
        ans.prec = M
        cdef int i
        for i in range(ans.prec):
            ans.moments[i] = self.moments[i]
        return ans

    def solve_diff_eqn(self):
        raise NotImplementedError

    #def lift(self):
    #    if self.prec >= self.parent()._prec_cap:
    #        raise ValueError("Cannot lift above precision cap")
    #    cdef Dist_long ans = self._new_c()
    #    ans.prec = self.prec + 1
    #    cdef int i
    #    for i in range(self.prec):
    #        ans.moments[i] = self.moments[i]
    #    return ans

cdef class WeightKAction(Action):
    def __init__(self, Dk, character, tuplegen, on_left):
        self._k = Dk._k
        if self._k < 0: raise ValueError("k must not be negative")
        if tuplegen is None:
            tuplegen = lambda g: (g[0,0], g[0,1], g[1,0], g[1,1])
        self._tuplegen = tuplegen
        if isinstance(character, tuple):
            if len(character) != 2:
                raise ValueError("character must have length 2")
            if character[0] is None and character[1] is None:
                self._character = None
                self._Np = Dk._p
            else:
                chr = character[0]
                dettwist = character[1]
                if chr is None:
                    self._character = lambda a, b, c, d: (a*d - b*c)**dettwist
                elif dettwist is None:
                    self._character = lambda a, b, c, d: chr(a)
                else:
                    self._character = lambda a, b, c, d: chr(a) * (a*d - b*c)**dettwist
                if chr is None:
                    self._Np = Dk._p
                else:
                    self._Np = Dk._p * chr.conductor()
        else:
            self._character = character
            self._Np = Dk._p # need to get conductor somehow in the case character = lambda g: ...
        self._p = Dk._p
        self._symk = Dk._symk
        self._actmat = {}
        self._maxprecs = {}
        Action.__init__(self, M2ZSpace, Dk, on_left, operator.mul)

    def clear_cache(self):
        self._actmat = {}
        self._maxprecs = {}

    cpdef acting_matrix(self, g, M):
        if not self._maxprecs.has_key(g):
            self._maxprecs[g] = M
            A = self._compute_acting_matrix(g, M)
            self._actmat[g] = {M:A}
            return A
        else:
            mats = self._actmat[g]
            if mats.has_key(M):
                return mats[M]
            maxprec = self._maxprecs[g]
            if M < maxprec:
                A = mats[maxprec][:M,:M] # submatrix; might want to reduce precisions
                mats[M] = A
                return A
            if M < 2*maxprec:
                maxprec = 2*maxprec
            else:
                maxprec = M
            self._maxprecs[g] = maxprec
            mats[maxprec] = self._compute_acting_matrix(g, maxprec) # could lift from current maxprec
            if M == maxprec:
                return mats[maxprec]
            A = mats[maxprec][:M,:M] # submatrix; might want to reduce precisions
            mats[M] = A
            return A

    cpdef _check_mat(self, a, b, c, d):
        if a*d == b*c:
            raise ValueError("zero determinant")
        if not self._symk:
            if self._p.divides(a):
                raise ValueError("p divides a")
            if not self._Np.divides(c):
                raise ValueError("Np does not divide c")

    cpdef _compute_acting_matrix(self, g, M):
        """
        Forms a large M x M matrix say G such that if v is the vector of
        moments of a distribution mu, then v*G is the vector of moments of
        mu|[a,b;c,d]
        """
        raise NotImplementedError

cdef class WeightKAction_vector(WeightKAction):
    cpdef _compute_acting_matrix(self, g, M):
        tim = verbose("Starting")
        a, b, c, d = self._tuplegen(g)
        self._check_mat(a, b, c, d)
        k = self._k
        if self._symk:
            base_ring = QQ
        else:
            base_ring = Zmod(self._p**M)
        R = PowerSeriesRing(base_ring, 'y', default_prec = M)
        y = R.gen()
        tim = verbose("Checked, made R",tim)
        # special case for small precision, large weight
        scale = (b+d*y)/(a+c*y)
        t = (a+c*y)**k # will already have precision M
        if self._character is not None:
            t *= self._character(a, b, c, d)
        cdef Matrix B = matrix(base_ring,M,M)
        cdef long row, col
        tim = verbose("Made matrix",tim)
        for col in range(M):
            for row in range(M):
                B.set_unsafe(row, col, t[row])
            t *= scale
        verbose("Finished loop",tim)
        # the changering here is annoying, but otherwise we have to change ring each time we multiply
        return B.change_ring(self.codomain().base_ring())

    cpdef _call_(self, _v, g):
        # if g is a matrix it needs to be immutable
        # hashing on arithmetic_subgroup_elements is by str
        cdef Dist_vector v = <Dist_vector?>_v
        cdef Dist_vector ans = v._new_c()
        verbose(str(v.moments))
        verbose(str(self.acting_matrix(g, len(v.moments))))
        ans.moments = v.moments * self.acting_matrix(g, len(v.moments))
        return ans

cdef inline long mymod(long a, unsigned long pM):
    a = a % pM
    if a < 0:
        a += pM
    return a

cdef class SimpleMat(SageObject):
    def __cinit__(self, unsigned long M):
        self._inited = False
        self.M = M
        self._mat = <long*>sage_malloc(M*M*sizeof(long))
        if self._mat == NULL:
            raise MemoryError
        self._inited = True

    def __getitem__(self, i):
        cdef Py_ssize_t r, c, Mnew, Morig = self.M
        cdef SimpleMat ans
        if PyTuple_Check(i) and PyTuple_Size(i) == 2:
            a, b = i
            if PySlice_Check(a) and PySlice_Check(b):
                r0, r1, rs = a.indices(Morig)
                c0, c1, cs = b.indices(Morig)
                if r0 != 0 or c0 != 0 or rs != 1 or cs != 1: raise NotImplementedError
                Mr = r1
                Mc = c1
                if Mr != Mc: raise ValueError("result not square")
                Mnew = Mr
                if Mnew > Morig: raise IndexError("index out of range")
                ans = SimpleMat(Mnew)
                for r in range(Mnew):
                    for c in range(Mnew):
                        ans._mat[Mnew*c + r] = self._mat[Morig*c + r]
                return ans
        raise NotImplementedError

    def __dealloc__(self):
        sage_free(self._mat)

cdef class WeightKAction_long(WeightKAction):
    cpdef _compute_acting_matrix(self, g, _M):
        _a, _b, _c, _d = self._tuplegen(g)
        if self._character is not None: raise NotImplementedError
        self._check_mat(_a, _b, _c, _d)
        cdef long k = self._k
        cdef Py_ssize_t row, col, M = _M
        cdef zmod_poly_t t, scale, xM, bdy
        cdef unsigned long pM = self._p**M
        cdef long a, b, c, d
        a = mymod(_a, pM)
        b = mymod(_b, pM)
        c = mymod(_c, pM)
        d = mymod(_d, pM)
        cdef double pMinv = pM
        pMinv = 1.0 / pMinv
        zmod_poly_init2_precomp(t, pM, pMinv, M)
        zmod_poly_init2_precomp(scale, pM, pMinv, M)
        zmod_poly_init2_precomp(xM, pM, pMinv, M+1)
        zmod_poly_init2_precomp(bdy, pM, pMinv, 2)
        zmod_poly_set_coeff_ui(xM, M, 1)
        zmod_poly_set_coeff_ui(t, 0, a)
        zmod_poly_set_coeff_ui(t, 1, c)
        zmod_poly_newton_invert(scale, t, M)
        zmod_poly_set_coeff_ui(bdy, 0, b)
        zmod_poly_set_coeff_ui(bdy, 1, d)
        zmod_poly_mul_trunc_n(scale, scale, bdy, M) # scale = (b+dy)/(a+cy)
        zmod_poly_powmod(t, t, k, xM) # t = (a+cy)^k
        cdef SimpleMat B = SimpleMat(M)
        for col in range(M):
            for row in range(M):
                B._mat[M*col + row] = zmod_poly_get_coeff_ui(t, row)
            if col < M - 1:
                zmod_poly_mul_trunc_n(t, t, scale, M)
        return B

    cpdef _call_(self, _v, g):
        cdef Dist_long v = <Dist_long?>_v
        cdef Dist_long ans = v._new_c()
        cdef long M = v.prec
        cdef long pM = self._p**M
        cdef SimpleMat B = <SimpleMat>self.acting_matrix(g, M)
        cdef long row, col, entry = 0
        for col in range(M):
            ans.moments[col] = 0
            for row in range(M):
                ans.moments[col] += mymod(B._mat[entry] * v.moments[row], pM)
                entry += 1
        ans.prec = M
        return ans

cdef class iScale(Action):
    def __init__(self, Dk, on_left):
        Action.__init__(self, ZZ, Dk, on_left, operator.mul)

    cpdef _call_(self, a, b):
        if PY_TYPE_CHECK(a, Dist):
            return (<Dist>a)._lmul_(b)
        else:
            return (<Dist?>b)._lmul_(a)

#@cached_function
#def eta(Dk, i, M):
#    """
#    Helper function in solving the difference equation -- see Lemma 4.4 of [PS]
#    """
#    v=[0 for j in range(0,i)]+[binomial(j,i)*bernoulli(j-i) for j in range(i,M)]
#    return dist(p,k,vector(v))
#
#def random_dist(p,k,M):
#    """
#    Returns a random distribution with prime p, weight k, and M moments
#    """
#    moments=vector([ZZ(floor(random()*(p**M))) for i in range(1,M+1)])
#
#    return dist(p,k,moments)


