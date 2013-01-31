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

from fund_domain import M2ZSpace,M2Z
cdef long overflow = 1 << (4*sizeof(long)-1)
cdef long underflow = -overflow

include "stdsage.pxi"
include "cdefs.pxi"

def get_dist_classes(p, prec_cap, base, symk):
    r"""
    Determines the element and action classes to be used for given inputs.

    INPUT:

    - ``p``        -- prime

    - ``prec_cap`` -- The p-adic precision cap

    - ``base``     -- The base ring

    - ``symk``     -- An element of Symk

    OUTPUT:

    - Either a Dist_vector and WeightKAction_vector, or a Dist_vector_long
       and WeightKAction_vector_long

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.dist import get_dist_classes
        sage: pass
    """
    if symk or p is None or base.is_field() or (isinstance(base, pAdicGeneric) and base.degree() > 1):
        return Dist_vector, WeightKAction_vector
    if 7*p**(prec_cap) < ZZ(2)**(4*sizeof(long)-1):
        return Dist_long, WeightKAction_long
    else:
        return Dist_vector, WeightKAction_vector

cdef class Dist(ModuleElement):
    r"""
        The main p-adic distribution class, implemented as per the paper
        'Overconvergent Modular Symbols and p-adic L-functions' by Pollack
        & Stevens
    """
    cpdef normalize(self):
        r"""
        Normalize so that the precision of the `i`-th moment is `n-i`,
        where `n` is the number of moments.

        OUTPUT:

        - Normalized entries of the distribution

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(5, 7, 15)
            sage: D
            Space of 7-adic distributions with k=5 action and precision cap 15
            sage: v = D([1,2,3,4,5]); v
            (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
            sage: v.normalize()
            (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
        """
        raise NotImplementedError

    def scale(self,left):
        r"""
        Scales the moments of the distribution by `left`

        INPUT:

        - ``left`` -- scalar

        OUTPUT:

        - Scales the moments by `left`

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(5, 7, 15)
            sage: v = D([1,2,3,4,5]); v
            (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
            sage: v.scale(2)
            (2 + O(7^5), 4 + O(7^4), 6 + O(7^3), 1 + 7 + O(7^2), 3 + O(7))
        """
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

    def is_zero(self, p=None, M=None):
        r"""
        Returns True if all of the moments are either identically zero or zero modulo p^M

        INPUT:

        - ``p`` -- prime

        - ``M`` -- precision

        OUTPUT:

        - True/False

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(5, 7, 15)
            sage: v = D([1,2,3,4,5]); v
            (1 + O(7^5), 2 + O(7^4), 3 + O(7^3), 4 + O(7^2), 5 + O(7))
            sage: v.is_zero()
            False
            sage: v = D(5*[0])
            sage: v.is_zero()
            True
        """
        n = self.precision_absolute()
        if M is None:
            return all([self.moment(a).is_zero() for a in range(n)])
        else:
            return all([self.moment(a).valuation(p) >= M for a in range(n)])

    def find_scalar(self, other, p, M = None, check=True):
        r"""
        Returns an ``alpha`` with ``other = self * alpha``, or raises a ValueError.

        It will also raise a ValueError if this distribution is zero.

        INPUT:

        - ``other`` -- another distribution

        - ``p`` -- an integral prime (only used if the parent is not a Symk)

        - ``M`` -- (default: None) an integer, the relative precision
          to which the scalar must be determined

        - ``check`` -- (default: True) boolean, whether to validate
          that ``other`` is actually a multiple of this element.

        OUTPUT:

        - A scalar ``alpha`` with ``other = self * alpha``.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(5, 7, 15)
            sage: v = D([1,2,3,4,5])
            sage: w = D([3,6,9,12,15])
            sage: v.find_scalar(w,p=7)
            3 + O(7^5)

            sage: u = D([1,4,9,16,25])
            sage: v.find_scalar(u,p=7)
            Traceback (most recent call last):
            ...
            ValueError: not a scalar multiple

        """
        i = 0
        n = self.precision_absolute()
        if n != other.precision_absolute():
            raise ValueError("other should have the same number of moments")
        verbose("n = %s"%n)
        verbose("moment 0")
        a = self.moment(i)
        verbose("a = %s"%(a))
        padic = isinstance(a.parent(), pAdicGeneric)
        if self.parent().is_symk():
            while a == 0:
                if other.moment(i) != 0:
                    raise ValueError("not a scalar multiple")
                i += 1
                verbose("moment %s"%i)
                try:
                    a = self.moment(i)
                except IndexError:
                    raise ValueError("self is zero")
            alpha = other.moment(i) / a
            if check:
                i += 1
                while i < n:
                    verbose("comparing moment %s"%i)
                    if self.moment(i) != alpha * other.moment(i):
                        raise ValueError("not a scalar multiple")
                    i += 1
        else:
            p = self.parent().prime()
            v = a.valuation(p)
            while v >= n - i:
                i += 1
                verbose("p moment %s"%i)
                try:
                    a = self.moment(i)
                except IndexError:
                    raise ValueError("self is zero")
                v = a.valuation(p)
            relprec = n - i - v
            verbose("p=%s, n-i=%s\nself.moment=%s, other.moment=%s"%(p, n-i, a, other.moment(i)),level=2)
            if padic:
                alpha = (other.moment(i) / a).add_bigoh(n-i)
            else:
                alpha = (other.moment(i) / a) % p**(n-i)
            verbose("alpha = %s"%(alpha))
            while i < n-1:
                i += 1
                verbose("comparing p moment %s"%i)
                a = self.moment(i)
                if check:
                    verbose("self.moment=%s, other.moment=%s"%(a, other.moment(i)))
                    if (padic and other.moment(i) != alpha * a) or \
                       (not padic and other.moment() % p**(n-i) != alpha * a % p**(n-i)):
                        raise ValueError("not a scalar multiple")
                v = a.valuation(p)
                if n - i - v > relprec:
                    verbose("Reseting alpha: relprec=%s, n-i=%s, v=%s"%(relprec, n-i, v))
                    relprec = n - i - v
                    if padic:
                        alpha = (other.moment(i) / a).add_bigoh(n-i)
                    else:
                        alpha = (other.moment(i) / a) % p**(n-i)
                    verbose("alpha=%s"%(alpha))
            if relprec < M:
                raise ValueError("result not determined to high enough precision")
        try:
            return self.parent().base_ring()(alpha)
        except ValueError:
            return alpha

    cpdef ModuleElement _rmul_(self, RingElement _left):
        """
        Scalar multiplication.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
        """
        return self._lmul_(_left)

    def diagonal_valuation(self, p=None):
        """
        Returns the largest `m` so that this distribution lies in `Fil^m`.

        INPUT:

        - ``p`` -- (default: None) a positive integral prime

        OUTPUT:

        - the largest integer `m` so that `p^m` divides the `0`-th
          moment, `p^{m-1}` divides the first moment, etc.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(8, 7, 15)
            sage: v = D([7^(5-i) for i in range(1,5)])
            sage: v
            (O(7^4), O(7^3), O(7^2), O(7))
            sage: v.diagonal_valuation(7)
            4
        """
        if p is None:
            p = self.parent()._p
        n = self.precision_absolute()
        return min([n] + [a + self.moment(a).valuation(p) for a in range(n)])

    def valuation(self, p=None):
        """
        Returns the minimum valuation of any moment.

        INPUT:

        - ``p`` -- (default: None) a positive integral prime

        OUTPUT:

        - 

        .. WARNING::

            Since only finitely many moments are computed, this
            valuation may be larger than the actual valuation of this
            distribution.  Moreover, since distributions are
            normalized so that the top moment has precision 1, this valuation may be smaller than the actual valuation (for example, if the actual valuation is 2)

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(8, 7, 15)
            sage: v = D([7^(5-i) for i in range(1,5)])
            sage: v
            (O(7^4), O(7^3), O(7^2), O(7))
            sage: v.valuation(7)
            1
        """
        r"""
        Returns the highest power of `p` which divides all moments of the distribution
        """
        if p is None:
            p = self.parent()._p
        n = self.precision_absolute()
        return min([self.moment(a).valuation(p) for a in range(n)])

    def specialize(self, new_base_ring=None):
        """
        Returns the image of this overconvergent distribution under
        the canonical projection from distributions of weight k to
        Sym^k.

        INPUT:

        - ``new_base_ring`` -- (default: None) a ring giving the
          desired base ring of the result.

        OUTPUT:

        - An element of Sym^k(K), where K is the specified base ring.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
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
        r"""
        Lifts a distribution or element of Sym^k to an overconvergent distribution.

        INPUT:

        - ``p`` -- (default: None) a positive integral prime.  If None
          then p must be available in the parent.

        - ``M`` -- (default: None) a positive integer giving the
          desired number of moments.

        - ``new_base_ring`` -- (default: None) a ring giving the
          desired base ring of the result.

        OUTPUT:

        - An overconvergent distribution with `M` moments whose image
          under the specialization map is this element.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
        """
        V = self.parent().lift(p, M, new_base_ring)
        k = V._k
        p = V.prime()
        M = V.precision_cap()
        R = V.base_ring()
        moments = [R(self.moment(j) * (-1)**j / binomial(k, j)) for j in range(k+1)]
        zero = R(0)
        moments.extend([zero] * (M - k - 1))
        mu = V(moments)
        #val = mu.valuation()
        #if val < 0:
        #    # This seems unnatural
        #    print "scaling by %s^%s to keep things integral"%(p, -val)
        #    mu *= p**(-val)
        return mu

    def _is_malformed(self):
        n = self.precision_absolute()
        for i in range(n):
            if self.moment(i).precision_absolute() < n - i:
                return True
        return False

    def act_right(self,gamma):
        r"""
        The image of this element under the right action by a
        `2 \times 2` matrix.

        INPUT:

        - ``gamma`` -- the matrix by which to act

        OUTPUT:

        - ``self | gamma``

        .. NOTE::

            You may also just use multiplication ``self * gamma``.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
        """
        return self.parent()._act(self, gamma)

cdef class Dist_vector(Dist):
    r"""
    A distribution is stored as a vector whose `j`-th entry is the `j`-th moment of the distribution.

    The `j`-th entry is stored modulo `p^(N-j)` where `N` is the total number of moments.
    (This is the accuracy that is maintained after acting by `\Gamma_0(p)`.)

    INPUTS:

    - ``moments`` -- the list of moments.  If ``check == False`` it
      must be a vector in the appropriate approximation module.

    - ``parent`` -- a :class:`distributions.Distributions_class` or
      :class:`distributions.Symk_class` instance

    - ``check`` -- (default: True) boolean, whether to validate input

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.distributions import Distributions
    """
    def __init__(self,moments, parent, check=True):
        """
        Initialization.

        TESTS::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
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

    def __reduce__(self):
        return (self.__class__,(self.moments,self.parent(),False))

    cdef Dist_vector _new_c(self):
        r"""
        Creates an empty distribution.

        OUTPUT:

        - A distribution with no moments.  The moments are then filled
          in by the calling function.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
        """
        cdef Dist_vector ans = PY_NEW(Dist_vector)
        ans._parent = self._parent
        return ans

    def _repr_(self):
        r"""
        String representation.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
        """
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

            sage: D = Symk(0); d = D(4/3); d
            4/3
            sage: QQ(d)
            4/3

        We get a TypeError if there is more than 1 moment::

            sage: D = Symk(1); d = D([1,2]); d
            (1, 2)
            sage: QQ(d)
            Traceback (most recent call last):
            ...
            TypeError: k must be 0
        """
        if len(self.moments) == 1:
            return QQ(self.moments[0])
        raise TypeError, "k must be 0"

    def moment(self, n):
        r"""
        Returns the `n`-th moment.

        INPUT:

        - ``n`` -- an integer or slice, to be passed on to moments.

        OUTPUT:

        - the `n`-th moment, or a list of moments in the case that `n`
          is a slice.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        r"""
        Returns the `n`-th moment
        """
        return self.moments[n]

    cpdef ModuleElement _add_(self, ModuleElement _right):
        r"""
        Sum of two distributions.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
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
        r"""
        Difference of two distributions.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
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
        r"""
        Scalar product of a distribution with a ring element that coerces into the base ring.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        cdef Dist_vector ans = self._new_c()
        ans.moments = self.moments * right
        return ans

    def precision_absolute(self):
        r"""
        Returns the precision of this distribution.

        The precision is just the number of moments stored, which is
        also k+1 in the case of Sym^k(R).  For overconvergent
        distributions, the precision is the integer `m` so that the
        sequence of moments is known modulo `Fil^m`.

        OUTPUT:

        - An integer giving the number of moments.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        r"""
        Returns the number of moments of the distribution
        """
        return len(self.moments)

    cdef int _cmp_c_impl(left, Element right) except -2:
        r"""
        Comparison.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        return cmp(left.moments, right.moments)

    def zero(self):
        r"""
        

        OUTPUT:

        - 

        """
        cdef Dist_vector ans = self._new_c()
        ans.moments = 0 * self.moments # could make this faster
        return ans

    cpdef normalize(self):
        r"""
        Normalize by reducing modulo `Fil^N`, where `N` is the number of moments.

        If the parent is Symk, then normalize has no effect.  If the
        parent is a space of distributions, then normalize reduces the
        `i`-th moment modulo `p^{N-i}`.

        OUTPUT:

        - this distribtion, after normalizing.

        .. WARNING::

            This function modifies the distribution in place as well as returning it.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        p = self.parent()._p
        if not self.parent().is_symk(): # non-classical
            if self.valuation() < 0:
                verbose("Negative valuation!")
                verbose("%s"%(self.moments))
            #assert self.valuation() >= 0, "moments not integral in normalization"
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
        Only hold on to `M` moments.

        INPUT:

        - ``M`` -- a positive integer less than the precision of this
          distribution.

        OUTPUT:

        - a new distribution with `M` moments equal to the first `M`
          moments of this distribution.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        assert M<=self.precision_absolute(),"not enough moments"

        cdef Dist_vector ans = self._new_c()
        ans.moments = self.moments[:M]
        return ans

    def solve_diff_eqn(self):
        r"""
        Solves the difference equation.

        See Theorem 4.5 and Lemma 4.4 of [PS].

        OUTPUT:

        - a distribution v so that self = v | Delta, where Delta = [1, 1; 0, 1] - 1.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
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
    r"""
    A class for distributions implemented using a C array of longs.

    INPUT:

    - ``moments`` -- the list of moments.  If ``check == False`` it
      must be a vector in the appropriate approximation module.

    - ``parent`` -- a :class:`distributions.Distributions_class` or
      :class:`distributions.Symk_class` instance

    - ``check`` -- (default: True) boolean, whether to validate input

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
    """
    def __init__(self, moments, parent, check=True):
        """
        Initialization.

        TESTS::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
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
        r"""
        

        OUTPUT:

        - 

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        cdef Dist_long ans = PY_NEW(Dist_long)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        return ans

    def _repr_(self):
        r"""
        

        OUTPUT:

        - 

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        self.normalize()
        if self.prec == 1:
            return repr(self.moments[0])
        else:
            return "(" + ", ".join([repr(self.moments[i]) for i in range(self.prec)]) + ")"

    cdef int quasi_normalize(self) except -1:
        r"""
        

        OUTPUT:

        - 

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        cdef int i
        for i in range(self.prec):
            if self.moments[i] > overflow:
                self.moments[i] = self.moments[i] % self.prime_pow.small_powers[self.prec-i]
            elif self.moments[i] < underflow:
                self.moments[i] = self.moments[i] % self.prime_pow.small_powers[self.prec-i]
                self.moments[i] += self.prime_pow.small_powers[self.prec-i]

    cpdef normalize(self):
        r"""
        

        OUTPUT:

        - 

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        cdef int i
        for i in range(self.prec):
            if self.moments[i] < 0:
                self.moments[i] = self.moments[i] % self.prime_pow.small_powers[self.prec-i]
                self.moments[i] += self.prime_pow.small_powers[self.prec-i]
            elif self.moments[i] >= self.prime_pow.small_powers[self.prec-i]:
                self.moments[i] = self.moments[i] % self.prime_pow.small_powers[self.prec-i]
        return self

    def moment(self, _n):
        r"""
        

        INPUT:

        - ``_n`` -- an integer or slice giving an index into the
          moments.

        OUTPUT:

        - 

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        if isinstance(_n, slice):
            a, b, c = _n.indices(self.prec)
            return [self.moment(i) for i in range(a, b, c)]
        cdef int n = _n
        if n < 0:
            n += self.prec
        if n < 0 or n >= self.prec:
            raise IndexError("list index out of range")
        return self.moments[n]

    cpdef ModuleElement _add_(self, ModuleElement _right):
        r"""
        

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
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
        r"""
        

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
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
        r"""
        

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
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
        r"""
        

        OUTPUT:

        - 

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        return self.prec

    cdef int _cmp_c_impl(left, Element _right) except -2:
        r"""
        

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        cdef int i
        cdef Dist_long right = _right
        for i in range(left.prec):
            if left.moments[i] < right.moments[i]:
                return -1
            if left.moments[i] > right.moments[i]:
                return 1
        return 0

    def zero(self):
        r"""
        

        OUTPUT:

        - 

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        cdef Dist_long ans = self._new_c()
        ans.prec = self.prec
        cdef int i
        for i in range(self.prec):
            ans.moments[i] = 0
        return ans

    def reduce_precision(self, M):
        r"""
        

        INPUT:

        - ``M`` -- a positive integer less than the precision of this
          distribution.

        OUTPUT:

        - 

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        if M > self.prec: raise ValueError("not enough moments")
        if M < 0: raise ValueError("precision must be non-negative")
        cdef Dist_long ans = self._new_c()
        ans.prec = M
        cdef int i
        for i in range(ans.prec):
            ans.moments[i] = self.moments[i]
        return ans

    def solve_diff_eqn(self):
        r"""
        

        OUTPUT:

        - 

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
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
    r"""
    

    INPUT:

    - ``Dk`` -- a space of distributions

    - ``character`` -- data specifying a Dirichlet character to apply
      to the top right corner, and a power of the determinant by which
      to scale.  See the documentation of
      :class:`sage.modular.pollack_stevens.distributions.Distributions_factory`
      for more details.

    - ``tuplegen`` -- a callable object that turns matrices into 4-tuples.

    - ``on_left`` -- whether this action should be on the left.

    OUTPUT:

    - 

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
    """
    def __init__(self, Dk, character, tuplegen, on_left):
        r"""
        Initialization.

        TESTS::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        self._k = Dk._k
        if self._k < 0: raise ValueError("k must not be negative")
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
        self._symk = Dk.is_symk()
        self._actmat = {}
        self._maxprecs = {}
        Action.__init__(self, M2ZSpace, Dk, on_left, operator.mul)

    def clear_cache(self):
        r"""
        

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        self._actmat = {}
        self._maxprecs = {}

    cpdef acting_matrix(self, g, M):
        r"""
        

        INPUT:

        - ``g`` -- an instance of
          :class:`sage.matrices.matrix_integer_2x2.Matrix_integer_2x2`

        - ``M`` -- a positive integer giving the precision at which
          ``g`` should act.

        OUTPUT:

        - An `M \times M` matrix so that the action of `g` on a
          distribution with `M` moments is given by a vector-matrix
          multiplication.

        .. NOTE::

            This function caches its results.  To clear the cache use
            :meth:`clear_cache`.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        g = M2Z(g)
        g.set_immutable()
        if not self._maxprecs.has_key(g):
            A = self._compute_acting_matrix(g, M)
            self._actmat[g] = {M:A}
            self._maxprecs[g] = M
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
        r"""
        

        INPUT:

        - ``a``, ``b``, ``c``, ``d`` -- integers, playing the role of
          the corresponding entries of the `2 \times 2` matrix that is
          acting.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        if a*d == b*c:
            raise ValueError("zero determinant")
        if not self._symk:
            if self._p.divides(a):
                raise ValueError("p divides a")
            if not self._Np.divides(c):
                raise ValueError("Np does not divide c")

    cpdef _compute_acting_matrix(self, g, M):
        r"""
        

        INPUT:

        - ``g`` -- an instance of
          :class:`sage.matrices.matrix_integer_2x2.Matrix_integer_2x2`

        - ``M`` -- a positive integer giving the precision at which
          ``g`` should act.

        OUTPUT:

        - 

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        """
        Forms a large M x M matrix say G such that if v is the vector of
        moments of a distribution mu, then v*G is the vector of moments of
        mu|[a,b;c,d]
        """
        raise NotImplementedError

cdef class WeightKAction_vector(WeightKAction):
    cpdef _compute_acting_matrix(self, g, M):
        r"""
        

        INPUT:

        - ``g`` -- an instance of
          :class:`sage.matrices.matrix_integer_2x2.Matrix_integer_2x2`

        - ``M`` -- a positive integer giving the precision at which
          ``g`` should act.

        OUTPUT:

        - 

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        #tim = verbose("Starting")
        a, b, c, d = self._tuplegen(g)
        self._check_mat(a, b, c, d)
        k = self._k
        if self._symk:
            base_ring = QQ
        else:
            base_ring = Zmod(self._p**M)
        R = PowerSeriesRing(base_ring, 'y', default_prec = M)
        y = R.gen()
        #tim = verbose("Checked, made R",tim)
        # special case for small precision, large weight
        scale = (b+d*y)/(a+c*y)
        t = (a+c*y)**k # will already have precision M
        if self._character is not None:
            t *= self._character(a, b, c, d)
        cdef Matrix B = matrix(base_ring,M,M)
        cdef long row, col
        #tim = verbose("Made matrix",tim)
        for col in range(M):
            for row in range(M):
                B.set_unsafe(row, col, t[row])
            t *= scale
        #verbose("Finished loop",tim)
        # the changering here is annoying, but otherwise we have to change ring each time we multiply
        return B.change_ring(self.codomain().base_ring())

    cpdef _call_(self, _v, g):
        r"""
        

        INPUT:

        - ``_v`` -- a :class:`Dist_vector` instance, the distribution
          on which to act.

        - ``g`` -- a
          :class:`sage.matrix.matrix_integer_2x2.Matrix_integer_2x2`
          instance, the `2 \times 2` matrix that is acting.

        OUTPUT:

        - 

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        # if g is a matrix it needs to be immutable
        # hashing on arithmetic_subgroup_elements is by str
        cdef Dist_vector v = <Dist_vector?>_v
        cdef Dist_vector ans = v._new_c()
        try:
            g.set_immutable()
        except AttributeError:
            pass
        ans.moments = v.moments * self.acting_matrix(g, len(v.moments))
        return ans

cdef inline long mymod(long a, unsigned long pM):
    """
    Returns the remainder ``a % pM``.

    INPUT:

    - ``a`` -- a long

    - ``pM`` -- an unsigned long

    OUPUT:

    - ``a % pM`` as a positive integer.
    """
    a = a % pM
    if a < 0:
        a += pM
    return a

cdef class SimpleMat(SageObject):
    r"""
    A simple class emulating a square matrix that holds its values as
    a C array of longs.

    INPUT:

    - ``M`` -- a positive integer, the dimension of the matrix

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
    """
    def __cinit__(self, unsigned long M):
        r"""
        Memory initialization.

        TESTS::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        self._inited = False
        self.M = M
        self._mat = <long*>sage_malloc(M*M*sizeof(long))
        if self._mat == NULL:
            raise MemoryError
        self._inited = True

    def __getitem__(self, i):
        r"""
        

        INPUT:

        - ``i`` -- a tuple containing two slices, each from `0` to `M'` for some `M' < M`

        OUTPUT:

        - A new SimpleMat of size `M'` with the top left `M' \times
          M'` block of values copied over.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
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
        r"""
        Deallocation.

        TESTS::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        sage_free(self._mat)

cdef class WeightKAction_long(WeightKAction):
    cpdef _compute_acting_matrix(self, g, _M):
        r"""
        

        INPUT:

        - ``g`` -- an instance of
          :class:`sage.matrices.matrix_integer_2x2.Matrix_integer_2x2`

        - ``_M`` -- a positive integer giving the precision at which
          ``g`` should act.

        OUTPUT:

        - A :class:`SimpleMat` that gives the action of ``g`` at
          precision ``_M`` in the sense that the moments of the result
          are obtained from the moments of the input by a
          vector-matrix multiplication.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
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
        r"""
        Application of the action.

        INPUT:

        - ``_v`` -- a :class:`Dist_long` instance, the distribution on
          which to act.

        - ``g`` -- a
          :class:`sage.matrix.matrix_integer_2x2.Matrix_integer_2x2`
          instance, the `2 \times 2` matrix that is acting.

        OUTPUT:

        - The image of ``_v`` under the action of ``g``.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
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
    r"""
    

    INPUT:

    - 

    OUTPUT:

    - 

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
    """
    def __init__(self, Dk, on_left):
        """
        Initialization.

        TESTS::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
        """
        Action.__init__(self, ZZ, Dk, on_left, operator.mul)

    cpdef _call_(self, a, b):
        """
        Application of the action.

        INPUT:

        - ``a``, ``b`` -- a :class:`Dist` or scalar, in either order.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: pass # XXX FIX THIS
        """
        if PY_TYPE_CHECK(a, Dist):
            return (<Dist>a)._lmul_(b)
        else:
            return (<Dist?>b)._lmul_(a)
