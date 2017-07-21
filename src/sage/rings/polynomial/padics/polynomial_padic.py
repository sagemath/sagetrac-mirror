"""
Base class for generic `p`-adic polynomials

This provides common functionality for all `p`-adic polynomials, such
as printing and factoring.

AUTHORS:

- Jeroen Demeyer (2013-11-22): initial version, split off from other
  files, made Polynomial_padic the common base class for all p-adic
  polynomials.
- Sebastian Pauli and Brian Sinclair (2017-07-21): omtree, factor for extensions,
  functionality for Eisenstein polynomials: ramification polynomial/polygon, 
  residual polynomials, monge reduction

"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed.math@gmail.com>
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
#       Copyright (C) 2013 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#       Copyright (C) 2017 Sebastian Pauli <s_pauli@uncg.edu>
#       Copyright (C) 2017 Brian Sinclair <brian.sinclair@alumni.uncg.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from six.moves import range
from builtins import zip

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_domain
from sage.rings.padics.precision_error import PrecisionError
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.padics.omtree.omtree import OMTree
from sage.rings.infinity import infinity
from sage.misc.cachefunc import cached_method
from sage.functions.other import binomial
from sage.structure.factorization import Factorization
from sage.rings.integer_ring import ZZ

class Polynomial_padic(Polynomial_generic_domain):
    r"""
    A polynomial over a `p`-adic ring.

    INPUT:

    - ``parent`` -- a polynomial ring over a `p`-adic ring

    - ``is_gen`` -- whether this is the generator of the polynomial ring
      (default: ``False``)

    .. NOTE::

        In contrast to :class:`polynomial_padic_generic.Polynomial_padic_generic`
        (which inherits from this class), this class is meant as a base class
        for implementations which provide their own handling of the polynomial
        data.

    EXAMPLES::

        sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
        sage: R.<x> = Zp(3)[] # indirect doctest
        sage: isinstance(x, Polynomial_padic)
        True
    """

    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        r"""
        Initialization.

        TESTS::

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R.<x> = Zp(3)[]
            sage: type(Polynomial_padic(R))
            <class 'sage.rings.polynomial.padics.polynomial_padic.Polynomial_padic'>

        """
        Polynomial_generic_domain.__init__(self, parent, is_gen=is_gen)

    def _repr(self, name=None):
        r"""
        EXAMPLES::

            sage: R.<w> = PolynomialRing(Zp(5, prec=5, type = 'capped-abs', print_mode = 'val-unit'))
            sage: f = 24 + R(4/3)*w + w^4
            sage: f._repr()
            '(1 + O(5^5))*w^4 + (O(5^5))*w^3 + (O(5^5))*w^2 + (1043 + O(5^5))*w + (24 + O(5^5))'
            sage: f._repr(name='z')
            '(1 + O(5^5))*z^4 + (O(5^5))*z^3 + (O(5^5))*z^2 + (1043 + O(5^5))*z + (24 + O(5^5))'

        TESTS::

            sage: k = Qp(5,10)
            sage: R.<x> = k[]
            sage: f = R([k(0,-3), 0, k(0,-1)]); f
            (O(5^-1))*x^2 + (O(5^-3))
            sage: f + f
            (O(5^-1))*x^2 + (O(5^-3))

        AUTHOR:

        - David Roe (2007-03-03), based on Polynomial_generic_dense._repr()
        """
        s = ""
        coeffs = self.list()
        if name is None:
            name = self.parent().variable_name()
        for n in reversed(range(len(coeffs))):
            x = coeffs[n]
            if x.valuation() != infinity:
                if s:
                    s += " + "
                x = "(%s)"%repr(x)
                if n > 1:
                    var = "*%s^%s"%(name,n)
                elif n==1:
                    var = "*%s"%name
                else:
                    var = ""
                s += (x + var)
        return s or "0"

    def quo_rem(self, right):
        """
        Returns the quotient and remainder of division by right

        EXAMPLES:
            sage: Kx.<x> = PolynomialRing(Zp(7))
            sage: (x^3+7*x+1).quo_rem(x-7)
            ((1 + O(7^20))*x^2 + (7 + O(7^21))*x + (7 + 7^2 + O(7^21)),
             (1 + 7^2 + 7^3 + O(7^20)))
        """
        return self._quo_rem_naive(right)

    def _quo_rem_naive(self, right):
        """
        Naive quotient with remainder operating on padic polynomials as their coefficient lists

        EXAMPLES:
            sage: Kx.<x> = PolynomialRing(ZpFM(7,20))
            sage: (x^3+7*x+1)._quo_rem_naive(x-7)
            ((1 + O(7^20))*x^2 + (7 + O(7^20))*x + (7 + 7^2 + O(7^20)),
             (1 + 7^2 + 7^3 + O(7^20)))
        """
        if right == 0:
            raise ZeroDivisionError, "cannot divide by a polynomial indistinguishable from 0"
        P = self.parent()
        F = list(self); G = list(right);
        fdeg = self.degree()
        gdeg = right.degree()
        Q = [0 for i in range(fdeg-gdeg+1)]
        j=1
        while fdeg >= gdeg:
            a = F[-j]
            if a!=0:
                for i in range(fdeg-gdeg,fdeg+1):
                    F[i] -= a * G[i-(fdeg-gdeg)]
                Q[fdeg-gdeg] += a
            j+=1; fdeg-=1;
        return (P(Q), P(F))

    @cached_method
    def content(self):
        """
        Compute the content of this polynomial.

        OUTPUT:

        If this is the zero polynomial, return the constant coefficient.
        Otherwise, since the content is only defined up to a unit, return the
        content as `\pi^k` with maximal precision where `k` is the minimal
        valuation of any of the coefficients.

        EXAMPLES::

            sage: K = Zp(13,7)
            sage: R.<t> = K[]
            sage: f = 13^7*t^3 + K(169,4)*t - 13^4
            sage: f.content()
            13^2 + O(13^9)
            sage: R(0).content()
            0
            sage: f = R(K(0,3)); f
            (O(13^3))
            sage: f.content()
            O(13^3)

            sage: P.<x> = ZZ[]
            sage: f = x + 2
            sage: f.content()
            1
            sage: fp = f.change_ring(pAdicRing(2, 10))
            sage: fp
            (1 + O(2^10))*x + (2 + O(2^11))
            sage: fp.content()
            1 + O(2^10)
            sage: (2*fp).content()
            2 + O(2^11)

        Over a field it would be sufficient to return only zero or one, as the
        content is only defined up to multiplication with a unit. However, we
        return `\pi^k` where `k` is the minimal valuation of any coefficient::

            sage: K = Qp(13,7)
            sage: R.<t> = K[]
            sage: f = 13^7*t^3 + K(169,4)*t - 13^-4
            sage: f.content()
            13^-4 + O(13^3)
            sage: f = R.zero()
            sage: f.content()
            0
            sage: f = R(K(0,3))
            sage: f.content()
            O(13^3)
            sage: f = 13*t^3 + K(0,1)*t
            sage: f.content()
            13 + O(13^8)

        """
        if self.is_zero():
            return self[0]
        else:
            try:
               val = self.valuation(val_of_var=0)
            except TypeError:
               val = min([c.valuation() for c in self])
            return self.base_ring()(self.base_ring().uniformizer_pow(val))


    def normalized(self,monic_if_possible=False):
        r"""
        Return a polynomial and divisor such that polynomial = self // divisor and
        the polynomial is such that 0 is the lowest valuation amongst the coefficients.

        INPUT::

        - ``monic_if_possible`` Use the leading coefficient as the scaling divisor if possible.

        OUTPUT::

            Normalized polynomial and the scaling divisor

        EXAMPLES::
            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R = ZpFM(3,30,print_mode='terse'); Rx.<x> = R[]
            sage: f = 9*x^9+3*x^2+27
            sage: f.normalized()
            ((3 + O(3^30))*x^9 + ... + (9 + O(3^30)), 3 + O(3^30))

        Here we see the unit factored out of the leading term if possible::

            sage: g = 6*x^9+3*x^2+27
            sage: g.normalized()
            ((2 + O(3^30))*x^9 + ... + (1 + O(3^30))*x^2 + (0 + O(3^30))*x + (9 + O(3^30)), 3 + O(3^30))
            sage: g.normalized(monic_if_possible=True)
            ((1 + O(3^30))*x^9 + ... + (102945566047325 + O(3^30))*x^2 + (0 + O(3^30))*x + (102945566047329 + O(3^30)), 6 + O(3^30))

        """
        normalize_by = self.content()
        self_normal = self.parent()([c>>normalize_by.valuation() for c in self])
        if monic_if_possible and self_normal.leading_coefficient().valuation() == 0:
            # if valuation of leading coefficient of normalized polynomial will be 0
            # normalize by leading coefficient
            normalize_by = normalize_by*self_normal.leading_coefficient()
            self_normal = self.parent()(self_normal/self_normal.leading_coefficient())
        return self_normal,normalize_by

    @cached_method
    def is_eisenstein(self):
        r"""
        Return True if the polynomial self is Eisenstein

        EXAMPLES::
            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R.<x> = ZpFM(3,7)[]
            sage: f = x^2 - x + 15
            sage: f.is_eisenstein()
            False
            sage: f = x^2 + 9*x + 15
            sage: f.is_eisenstein()
            True
        """
        coeffs = self.coefficients(sparse=False)
        n      = self.degree()
        return (coeffs[0].valuation() == 1
                and coeffs[n].valuation() == 0
                and all(coeffs[x].valuation() >= 1 for x in range(1,n)))

    def deformed_eisenstein(self, m, theta, trunc):
        r"""
        Return the "deformed" Eisenstein polynomial, ie. the minimal polynomial of the uniformizer `q` such that `q + theta*q^(m+1) = pi`.
        Note that we also have that `q = pi - theta*pi^(m+1) + O(pi^(m+2))`

        EXAMPLES::

        We deform an Eisenstein polynomial and check the relationship between `q` and `pi`.

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R = ZpFM(3,20,print_mode='terse'); Rx.<x> = R[]
            sage: f = x^9+9*x^2+3
            sage: f.is_eisenstein()
            True
            sage: m = 2; theta = 1; trunc = 180
            sage: g = f.deformed_eisenstein(m,theta,trunc)
            sage: g
            (1 + O(3^20))*x^9 + (3205887336 + O(3^20))*x^8 + (634107690 + O(3^20))*x^7 + (1030936815 + O(3^20))*x^6 + (1721931291 + O(3^20))*x^5 + (2532448143 + O(3^20))*x^4 + (1806953859 + O(3^20))*x^3 + (3022253919 + O(3^20))*x^2 + (460349730 + O(3^20))*x + (3286601949 + O(3^20))
            sage: S.<q> = R.ext(g)
            sage: pi = q + theta*q^(m+1)
            sage: Sy.<y> = S[]
            sage: h = f.move(Sy)
            sage: h(pi).is_zero()
            True

        AUTHORS:

        - Maurizio Monge (2014-11-14): initial version
        - Sebastian Pauli and Brian Sinclair (2017-07-20): bang to sage 8.0, documentation

        """
        if not self.is_eisenstein():
            raise ValueError("this function can only deform Eisenstein polynomials")

        f = self
        n = f.degree()
        x = f.variables()[0]

        # truncation level, because x^trunc = 0 mod piK^prec, as v(x) = 1/n
        g = f(x + theta*x**(m+1)).truncate(trunc) # IMPROVE-ME: use f_0 instead of x^n

        prev = (0, n)
        while(g.degree() > n):
            g_coeffs = g.coefficients(sparse=False)

            #where f has terms we want to kick out
            extra_range = range(n+1, g.degree()+1)

            cur_prec = min([(g_coeffs[i].valuation())*n + i for i in extra_range])
            if cur_prec >= trunc:
                break

            # lexicographically minimal pair (val(f_i), i), for i in the extra range
            min_val, min_idx = min([(g_coeffs[i].valuation(), i) for i in extra_range])

            # the bad terms should be going away...
            assert(prev < (min_val, min_idx))
            prev = (min_val, min_idx)

            # subtract, from g, (badmononial/x^n)*g
            g = g - g_coeffs[min_idx] * g.shift(min_idx-n).truncate(trunc)

        return g.truncate(n+1).monic()

    def move(self,newRx):
        r"""
        Move the polynomial `self` over `R` to a polynomial in `newRx`, the polynomial ring of `newR` which was obtained `R` by change

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: zp = ZpFM(3,7)
            sage: zpx.<x> = zp[]
            sage: f = x^2 - x + 15
            sage: f
            (1 + O(3^7))*x^2 + (2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + O(3^7))*x + (2*3 + 3^2 + O(3^7))
            sage: zp100 = zp.change(prec=100, type="fixed-mod")
            sage: zp100y.<y> = zp100[]
            sage: g = f.move(zp100y)
            sage: g
            (1 + O(3^100))*y^2 + (2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + O(3^100))*y + (2*3 + 3^2 + O(3^100))
            sage: g.move(zpx)
            (1 + O(3^7))*x^2 + (2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + O(3^7))*x + (2*3 + 3^2 + O(3^7))
        """
        if self.is_zero():
            return newRx(0)
        if newRx.base() == newRx.base().base():
            return newRx([newRx.base()(a) for a in self])
        else:
            return newRx([newRx.base()(a.polynomial().list()) if not a.is_zero() else newRx.base()(0) for a in self])

    @cached_method
    def omtree(self, check_squarefree=True, check_field=True):
        r"""
        The OM (Okutsu-Montes/Ore-Mac Lane) tree of the polynomial `self`, which must be monic and square free with integral coefficients.

        AUTHORS:

        - Brian Sinclair and Sebastian Pauli (2017-07-18): initial version

        EXAMPLES::

        We construct the OM tree of a polynomial and extract some of the information it
        yields about the extensions generated by the irreducible factors of the polynomial.

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
            sage: k = Zp(2, prec = 500)
            sage: kx.<x> = k[]
            sage: f = (x^32+16)*((x^2+x+1)^16+80)
            sage: om=f.omtree(check_squarefree=False)
            sage: om.depths()
            [3, 3]
            sage: om.ramification_indices()
            [32, 16]
            sage: om.inertia_degrees()
            [1, 2]
        """
        if not self.is_monic():
            raise Error("OM tree is defined for monic polynomials only")

        if check_field and self.parent().base().is_field():
            raise Error("OM tree is defined for polynomials over valuation rings only")

        if check_squarefree:
            absprec = min([x.precision_absolute() for x in self])
            if self.discriminant().valuation() >= absprec:
                raise PrecisionError("OM tree is not well-defined since the discriminant is zero up to the requested p-adic precision")

        R = self.parent().base()
        RFM = R.change(type="fixed-mod", field=False)
        #Phi = self.change_ring(RFM)
        RFMy = PolynomialRing(RFM,names='y')
        Phi = self.move(RFMy)

        return OMTree(Phi)


    def ramification_polynomial(self):
        r"""
        Returns the ramification polynomial of `self`.

        The ramification polynomial is `self(a*x+a)/a^n` where `a` is a root of `self` and `n` is the degree of self.

        EXAMPLES::
            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R = ZpFM(3,20,print_mode='terse'); Rx.<x> = R[]
            sage: f = x^9+9*x^2+3
            sage: f.is_eisenstein()
            True
            sage: rho = f.ramification_polynomial()
            sage: rho
            (1 + O(a^180))*y^9 + (9 + O(a^180))*y^8 + (36 + O(a^180))*y^7 + (84 + O(a^180))*y^6 + (126 + O(a^180))*y^5 + (126 + O(a^180))*y^4 + (84 + O(a^180))*y^3 + (36 + 729*a + 3486784398*a^2 + 2324522943*a^4 + 3486784374*a^6 + 81*a^8 + O(a^180))*y^2 + (9 + 1458*a + 3486784395*a^2 + 2324522952*a^4 + 3486784347*a^6 + 162*a^8 + O(a^180))*y + (0 + O(a^180))
            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: rp = NewtonPolygon([(i,rho[i].valuation()) for i in range(1,rho.degree()+1)])
            sage: rp == f.ramification_polygon()
            True

        AUTHORS:

        - Sebastian Pauli and Brian Sinclair (2017-07-19): initial version

        """
        if not self.is_eisenstein():
            raise ValueError("the polynomial self must be Eisenstein")

        Rx = self.parent()
        R = Rx.base()
        S = R.ext(self, names='a'); a = S.uniformizer()
        Sy = PolynomialRing(S, names='y')
        n = self.degree()
        r0 = [0]+[sum([(binomial(k,i)*S(self[k]))>>(n-k) for k in range(i,n+1)]) for i in range(1,n)]+[1]
        return Sy(r0)

    @cached_method
    def ramification_polygon(self):
        r"""
        Returns the ramification polygon of `self`.

        The ramification polygon is the Newton polygon of self(a*x+a)/a^n where a is a root of `self` and `n` is the degree of `self`.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R = ZpFM(3,20,print_mode='terse'); Rx.<x> = R[]
            sage: f = x^27+3*x^24+3*x^18+3*x^9+9*x^3+9*x^3+6
            sage: rp = f.ramification_polygon()
            sage: rp.vertices()
            [(1, 51), (3, 24), (9, 9), (27, 0)]
            sage: rp.slopes(repetition=False)
            [-27/2, -5/2, -1/2]

            sage: R = ZpFM(3,20,print_mode='terse'); Rx.<x> = R[]
            sage: f = x^108+3*x^24+3*x^18+3*x^9+9*x^3+9*x^3+6
            sage: rp = f.ramification_polygon()
            sage: rp.vertices()
            [(1, 132), (3, 24), (9, 9), (27, 0), (108, 0)]
            sage: rp.slopes(repetition=False)
            [-54, -5/2, -1/2, 0]

        AUTHORS:

        - Brian Sinclair and Sebastian Pauli (2017-07-19): initial version

        """
        if not self.is_eisenstein():
            raise ValueError("the polynomial self must be Eisenstein")

        from sage.geometry.newton_polygon import NewtonPolygon

        # First we find the ordinates of points above p^k
        verts = []
        vv = [cc.valuation() for cc in self]
        k = self.base_ring()
        p = k.prime()
        n = self.degree()
        su = n.valuation(p)
        abscissa = 1
        for i in range(su):
            abscissa = p**i
            ordinate = min([n * (k(binomial(kk,abscissa)).valuation() + vv[kk] - 1) + kk for kk in range(abscissa,n+1)])
            verts.append((abscissa,ordinate))

        # Now we add the tame segment
        for i in range(p**su,n):
            if binomial(n,i).valuation(p) == 0:
                verts.append((i,0))

        # Finally the point for the monic leading term
        verts.append((n,0))

        return NewtonPolygon(verts)

    @cached_method
    def ramification_polygon_with_colinear_points(self):
        r"""
        Returns the ramification polygon of `self` as a list of points including all points on segments of the lower convex hull.

        The ramification polygon is the Newton polygon of `self(a*x+a)/a^n` where `a` is a root of self and `n` is the degree of `self`.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R = ZpFM(5,20,print_mode='terse'); Rx.<x> = R[]
            sage: f = x^25 + 20*x^6 + 20*x^5 + 5
            sage: f.ramification_polygon_with_colinear_points()
            [(1, 6), (5, 5), (25, 0)]

        The colinear point (5, 5) is missing on the ramification polygon returned by::

            sage: f.ramification_polygon()
            Finite Newton polygon with 2 vertices: (1, 6), (25, 0)

        When the generated extensions has a tamely ramified subextension::

            sage: f = x^100 + 20*x^6 + 20*x^5 + 5
            sage: f.ramification_polygon_with_colinear_points()
            [(1, 6), (5, 5), (25, 0), (50, 0), (75, 0), (100, 0)]
            sage: f.ramification_polygon()
            Finite Newton polygon with 3 vertices: (1, 6), (25, 0), (100, 0)

        AUTHORS:

        - Brian Sinclair (2017-07-19): initial version

        """
        if not self.is_eisenstein():
            raise ValueError("the polynomial self must be Eisenstein")

        # First we find the ordinates of points above p^k
        verts = []
        vv = [cc.valuation() for cc in self]
        k = self.base_ring()
        p = k.prime()
        n = self.degree()
        su = n.valuation(p)
        abscissa = 1
        for i in range(su):
            abscissa = p**i
            ordinate = min([n * (k(binomial(kk,abscissa)).valuation() + vv[kk] - 1) + kk for kk in range(abscissa,n+1)])
            verts.append((abscissa,ordinate))

        # Now we add the tame segment
        for i in range(p**su,n):
            if binomial(n,i).valuation(p) == 0:
                verts.append((i,ZZ(0)))

        # Finally the point for the monic leading term
        verts.append((n,ZZ(0)))

        # Next we need to take the lower convex hull of these points
        def cross(o, a, b):
            """
            2D cross product of the vectors oa and ob.
            """
            return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

        lower = [verts[0]]
        segments = []
        for i in range(1,len(verts)):
            # We check cross < 0 since we want to retain points on the boundary.
            while len(lower) >= 2 and cross(lower[-2], lower[-1], verts[i]) < 0:
                lower.pop()
            lower.append(verts[i])
        if len(lower) <= 1:
            raise ValueError("Not enough vertices")
        return lower

    def hasse_herbrand(self,m):
        r"""
        Returns `n` times the (generalized) Hasse-Herbrand function of `self` evaluated at `m`.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R = ZpFM(2,200,print_mode='terse'); Rx.<x> = R[]
            sage: f = x^16+2
            sage: f.ramification_polygon()
            Finite Newton polygon with 5 vertices: (1, 64), (2, 48), (4, 32), (8, 16), (16, 0)

        We evaluate the Hasse-Herbrand function at various integers.

            sage: [f.hasse_herbrand(m) for m in range(18)]
            [0, 16, 32, 40, 48, 52, 56, 60, 64, 66, 68, 70, 72, 74, 76, 78, 80, 81]

        Now a different degree:

            sage: f = x^80+2
            sage: f.ramification_polygon()
            Finite Newton polygon with 6 vertices: (1, 320), (2, 240), (4, 160), (8, 80), (16, 0), (80, 0)
            sage: [f.hasse_herbrand(m) for m in range(16)]
            [0, 16, 32, 48, 64, 80, 96, 112, 128, 144, 160, 168, 176, 184, 192, 200]

        AUTHORS:

        - Brian Sinclair (2017-07-20): initial version
        """
        return min([pt[1]+m*pt[0] for pt in self.ramification_polygon().vertices()])

    @cached_method
    def residual_polynomials(self):
        """
        Returns a list of the residual polynomials of the ramification polynomial of an Eisenstein polynomial self.

        EXAMPLES::
            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R = ZpFM(3,30); Rx.<x> = R[]
            sage: f = x^9+6*x^3+3
            sage: f.is_eisenstein()
            True
            sage: f.ramification_polygon()
            Finite Newton polygon with 3 vertices: (1, 12), (3, 3), (9, 0)
            sage: f.residual_polynomials()
            [z + 2, z^3 + 1]

        A ramfication polygon with a horizontal segment::

            sage: f = x^90+6*x^3+3
            sage: f.ramification_polygon()
            Finite Newton polygon with 4 vertices: (1, 93), (3, 3), (9, 0), (90, 0)
            sage: f.ramification_polygon().slopes(repetition=False)
            [-45, -1/2, 0]
            sage: f.residual_polynomials()
            [z^2 + 2, z^3 + 1, z^81 + z^72 + 1]

        A ramification polygon with more segments::

            sage: R = ZpFM(2,300,print_mode='terse'); Rx.<x> = R[]
            sage: f = x^16+2
            sage: f.ramification_polynomial()
            (1 + O(a^4800))*y^16 + (16 + O(a^4800))*y^15 + (120 + O(a^4800))*y^14 + (560 + O(a^4800))*y^13 + (1820 + O(a^4800))*y^12 + (4368 + O(a^4800))*y^11 + (8008 + O(a^4800))*y^10 + (11440 + O(a^4800))*y^9 + (12870 + O(a^4800))*y^8 + (11440 + O(a^4800))*y^7 + (8008 + O(a^4800))*y^6 + (4368 + O(a^4800))*y^5 + (1820 + O(a^4800))*y^4 + (560 + O(a^4800))*y^3 + (120 + O(a^4800))*y^2 + (16 + O(a^4800))*y + (0 + O(a^4800))
            sage: f.ramification_polygon()
            Finite Newton polygon with 5 vertices: (1, 64), (2, 48), (4, 32), (8, 16), (16, 0)
            sage: f.residual_polynomials()
            [z + 1, z^2 + 1, z^4 + 1, z^8 + 1]

        AUTHORS:

        - Sebastian Pauli and Brian Sinclair (2017-07-20): initial version

        """
        if not self.is_eisenstein():
            raise ValueError("the polynomial self must be Eisenstein")

        n = self.degree()
        phi0 = self.constant_coefficient()
        phi00 = phi0>>1
        Rx = self.parent()
        R = Rx.base()
        p = R.prime()
        F = R.residue_class_field()
        Fz = PolynomialRing(F,names='z')
        z = Fz.gen()

        rp = self.ramification_polygon_with_colinear_points()
        respols = []
        j = 0
        while j<len(rp)-1:
            k = j
            p_s_k = rp[k][0]
            slope = (rp[j+1][1]-rp[j][1])/(rp[j+1][0]-rp[j][0])
            e = slope.denominator()
            thispol = Fz(0)
            while True:
                p_s_i = rp[j][0]
                a_i, b_i = rp[j][1].quo_rem(n)
                if b_i == 0:
                    a_i -= 1
                    b_i  = n
                thispol += ((self[b_i]*binomial(b_i,p_s_i)*(-phi00)**(-a_i-1))>>(a_i+1)).residue()*z**((p_s_i-p_s_k)//e)
                if j>=len(rp)-1 or (rp[j+1][1]-rp[j][1])/(rp[j+1][0]-rp[j][0]) != slope:
                    break
                j+=1
            respols.append(thispol)

        return respols

    def residual_polynomial_of_component(self,m):
        r"""
        Return the residual polynomial ``S_m`` of the (`-m`)-component of the ramifation polygon of polynomials self, which must be Eisenstein.

        Let N be the ramification polygon then {(k,w) in N | (`-m`)k + w = min{(`-m`)l+u|(l,u) in N} } is the (`-m`)-component of N.

        INPUT::

            A natural number `m`

        OUTPUT::

            The residual polynomial of the (`-m`)-component of the ramificaton polygon of self.

        EXAMPLES::

        In our first example, we have a polynomial over the degree 2 unramified extension of ``\QQ_2``
        which has a ramification polygon with two segments of integral slope::

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R.<g> = ZqFM(4,30); Rx.<x> = R[]
            sage: f = x^8 + 2*g*x^6 + 4*g*x + 2
            sage: f.ramification_polygon()
            Finite Newton polygon with 3 vertices: (1, 9), (2, 6), (8, 0)
            sage: f.ramification_polygon().slopes(repetition=False)
            [-3, -1]
            sage: f.residual_polynomials()
            [g0*z + g0, z^6 + g0]
            sage: [f.residual_polynomial_of_component(m) for m in range(1,10)]
            [z^8 + g0*z^2, z^2, g0*z^2 + g0*z, z, z, z, z, z, z]

        Here we have a nonic polynomial over ``\QQ_3`` whose ramification polygon has no segemnts
        of integral slope::

            sage: R = ZpFM(3,30); Rx.<x> = R[]
            sage: f = x^9 + 6*x^3 + 3
            sage: f.ramification_polygon()
            Finite Newton polygon with 3 vertices: (1, 12), (3, 3), (9, 0)
            sage: f.ramification_polygon().slopes(repetition=False)
            [-9/2, -1/2]
            sage: f.residual_polynomials()
            [z + 2, z^3 + 1]
            sage: [f.residual_polynomial_of_component(m) for m in range(0,10)]
            [z^9, z^3, z^3, z^3, z^3, z, z, z, z, z]

        Here the ramification polygon has a horizontal segment::

            sage: f = x^90 + 6*x^3 + 3
            sage: f.ramification_polygon()
            Finite Newton polygon with 4 vertices: (1, 93), (3, 3), (9, 0), (90, 0)
            sage: f.ramification_polygon().slopes(repetition=False)
            [-45, -1/2, 0]
            sage: f.residual_polynomials()
            [z^2 + 2, z^3 + 1, z^81 + z^72 + 1]
            sage: [f.residual_polynomial_of_component(m) for m in range(0,10)]
            [z^90 + z^81 + z^9, z^3, z^3, z^3, z^3, z^3, z^3, z^3, z^3, z^3]

        AUTHORS:

        - Sebastian Pauli and Brian Sinclair (2017-07-20): initial version

        REFERENCES:

        [PS17] S. Pauli and B. Sinclair, "Enumerating Extensions of (pi)-adic Fields with Given Invariants", International Journal of Number Theory
          (2017)
        """

        if not self.is_eisenstein():
            raise ValueError("residual polynomials are only defined for Eisenstein polynomials")

        Rx = self.parent()
        R = Rx.base()
        F = R.residue_class_field()
        Fz = PolynomialRing(F,names='z')
        z = Fz.gen()

        rp = self.ramification_polygon()

        if -m in rp.slopes(repetition=False):
            i = rp.slopes(repetition=False).index(-m)
            return self.residual_polynomials()[i]*z**rp.vertices()[i][0]
        else:
            L = [v[1]+v[0]*m for v in rp.vertices()]
            mini = min(L)
            mindex = L.index(mini)
            return z**rp.vertices()[mindex][0]

    def has_residual_polynomial_class(self,A):
        r"""
        Checks whether a list `A` of polynomials is in the same residual polynomial class as the residual polynomials of `self`.

        INPUT::

            A list `A` of polynomials over the residual class field of the coefficient ring of the polynomial `self`.

        OUTPUT::

            True if the polynomials in `A` are in thee same residual polynomial class as the residual polynomials of self.

        EXAMPLES::

        Clearly the residual polynomials of `f` are in the residual polynomial class of f.

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic; 
            sage: R = ZpFM(3,30); Rx.<x> = R[]
            sage: F = R.residue_class_field(); Fz.<z>=F[] 
            sage: f = x^9+6*x^3+3; 
            sage: A = f.residual_polynomials()
            sage: A
            [z + 2, z^3 + 1]
            sage: f.has_residual_polynomial_class(A)
            True

        In the following the polynomial g generates an extension isomorphic to the extension generated by f.
        So the residual polynomials of g are in the residual polynomial class of g.

            sage: g = f(2*x)/512
            sage: g.residual_polynomials()
            [2*z + 2, z^3 + 2]
            sage: g.has_residual_polynomial_class(A)
            True

        AUTHORS:

        - Sebastian Pauli and Brian Sinclair (2017-07-21): initial version

        REFERENCES:

        [PS17] S. Pauli and B. Sinclair, "Enumerating Extensions of (\pi)-adic Fields with Given Invariants",
        International Journal of Number Theory (2017)

        """

        if not self.is_eisenstein():
            raise ValueError("residual polynomials are only definesd for Eisenstein polynomials")

        B = self.residual_polynomials()

        if len(A) != len(B):
            return False

        for i in range(0,len(B)):
            if B[i].degree() != A[i].degree():
                return False

        slopes = self.ramification_polygon().slopes()

        Rx = self.parent()
        R = Rx.base()
        F = R.residue_class_field()
        Fz = PolynomialRing(F,names='z')
        z = Fz.gen()

        h = [-s.numerator() for s in slopes]

        for delta in F:
            if delta !=0:
                g=[delta**(-h[-1]*B[-1].degree())]
                for i in range(len(B)-2,-1,-1):
                    g = [g[0]*delta**(-h[i]*B[i].degree())] + g
                Bdelta =  [g[i]*B[i](delta**h[i]*z) for i in range(0,len(B))]
                if A == Bdelta:
                    return True
        return False

    def monge_reduce(self):
        r"""
        Return the Monge-reduced polynomial that generates an extensions isomorphic to the extensions generated by the Eisenstein polynomial `self`.

        When fixing set of representatives for the classes of elements of the residue class field occuring in the Monge reduction, 
        the Monge-reduced polynomials are unique.  We make the following choices.

        If the coefficient ring of `self` is unramified, we choose the representatives of the classes of elements of the residue class field 
        such that their absolute value is minimal.  Otherwise we order the representatives lexicographically.
        Note that this representation depends on the generating polynomial of the unramified part of the extension,

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.polynomial_padic import Polynomial_padic
            sage: R = ZpFM(3,30); Rx.<x> = R[]
            sage: f = x^9+249*x^3+486*x+30
            sage: g = f.monge_reduce()
            sage: g
            (1 + O(3^30))*x^9 + ... + (2*3^2 + O(3^30))*x^6 + ... + (2*3 + 2*3^2 + O(3^30))*x^3 + ... + (3 + O(3^30))

        We now create the extension S of R generated by g and examine the factorization of f over S.
        As the polynomial f has a linear factor over S (and deg(f)=deg(g)).

            sage: S.<a> = R.ext(g)
            sage: Sy.<y> = S[]
            sage: fS = Sy(f)
            sage: fS.omtree().degrees_of_factors()
            [1, 2, 6]

        The Monge-reduction of a polynomial generating a tamely ramified extension::

            sage: f = x^20+249*x^3+486*x+30
            sage: g = f.monge_reduce()
            sage: g
            (1 + O(3^30))*x^20 + ... + (3 + O(3^30))

        The Monge-reduction of a polynomial generating a tamely ramified extensionof large degree::

            sage: f = x^90+249*x^81+486*x^18+30
            sage: g = f.monge_reduce()
            sage: g
            (1 + O(3^30))*x^90 + ... + (2*3 + O(3^30))*x^81 + ... + (3 + 3^3 + O(3^30))

        We use Monge reduction to verify that two polynomials generate isomorphic extensions

            sage: R = ZpFM(5,20); Rx.<x> = R[]
            sage: f = x^25+15625*x^4+5
            sage: g = x^25+5
            sage: f.monge_reduce() == g.monge_reduce()
            True

        Monge-reduction over an unramified extensions::

            sage: R.<g> = ZqFM(4,30); Rx.<x> = R[]
            sage: f = x^8 + 66*g*x^6 + 132*g*x + 258
            sage: g = f.monge_reduce()
            sage: g
            (1 + O(2^30))*x^8 + (O(2^30))*x^7 + (g*2 + O(2^30))*x^6 + ... + (g*2^2 + O(2^30))*x + (2 + O(2^30))

        AUTHORS:

        - Sebastian Pauli and Brian Sinclair (2017-07-20): initial version

        REFERENCES:

        [Mon24] M. Monge, "A family of Eisenstein polynomials
          generating totally ramified extensions, identification of extensions and
          construction of class fields." International Journal of Number Theory
          (2014): 1-29.
        """

        if not self.is_eisenstein():
            raise ValueError("only Eisenstein polynomials can be Monge-reduced")

        f = self
        n = f.degree()
        RT = f.parent()
        T = RT.gen()
        R = RT.base()
        p = R.prime()
        F = R.residue_class_field()
        Fz = PolynomialRing(F,names='z')
        z = Fz.gen()
        r = F.degree()

        def canonical_representative_mult(alpha, modset):
        # terrible, stop reading here
            if F.is_prime_field():
                bset = [ZZ(alpha*s) for s in modset]
            else:
                bset = [[ZZ(b) for b in (alpha*s).polynomial()] for s in modset]
            bset.sort()
            return F(bset[0])

        def canonical_representative_add(alpha, modset):
        # terrible, stop reading here
            if F.is_prime_field():
                bset = [ZZ(alpha+s) for s in modset]
            else:
                bset = [[ZZ(b) for b in (alpha+s).polynomial()] for s in modset]
            bset.sort()
            return F(bset[0])

        def solve_naive(funct,gamma):
        # solve poly(t)=gamma
        # terrible
            for a in F:
                if funct(a)==gamma:
                    return a
            raise Error("No solution found")

        phi0 = f.constant_coefficient()
        eta = (phi0 >> 1).residue()
        alpha = eta

        # reduction step 0 -- taking care of the constant coefficient

        if n == p**n.valuation(p):
            beta = 1
        else:
            beta = canonical_representative_mult(eta,set([a**n for a in F if not a.is_zero()]))

        theta = solve_naive(z**n, beta/alpha)
        Theta = R(theta)
        # f = theta**n*f(theta^-1*T)
        f = RT([f[i]*Theta**(n-i) for i in range(0,n+1)])

        # other reduction steps

        def f_ij(f, lev):
            i = lev % n
            j = (n - i + lev) // n
            fij = F(f[i].padded_list(j+1)[-1])
            return fij, i, j

        for m in range(1,(-min(self.ramification_polygon().slopes())).ceil()):
            alpha, i, j = f_ij(f,m)
            Sm = f.residual_polynomial_of_component(m)
            beta = canonical_representative_add(alpha,[eta**j*Sm(a) for a in F])
            theta = solve_naive(eta**j*Sm,alpha-beta)
            f = f.deformed_eisenstein(m, R(theta), n*R.precision_cap())

        # Find the last break in the Hasse-Herbrand function of self
        hhslope = n
        m = 1
        hhm = 0
        while hhslope > 1:
            lrb = hhm
            hhm = self.hasse_herbrand(m)
            hhslope = hhm - lrb
            m += 1

        # If f generates a tamely ramified extension, then lrb == 0, which
        # would clear all coefficients.  Set to 1, it preserves f_{0,1}.
        lrb = max(lrb,1)

        # All coefficients lexicographically beyond
        #    (lrb mod n, (n-(lrb mod n) + lrb)/n)
        # can be set to zero.
        i = lrb % n
        j = (n - i + lrb) / n
        f = RT([cc % R.uniformizer_pow(j+1) for cc in f.list()[:i]] + [cc % R.uniformizer_pow(j) for cc in f.list()[i:n]] + [R.one()])
        return f

    def factor(self,check_squarefree=True,algorithm="auto"):
        r"""
        Return the factorization of this polynomial.

        Different implementations of algorithms can be selected
        - algorithm="auto" chooses an appropriate algorithm
        - algorithm="pari" use the implementation of the Round 4 algorithm in Pari.  This only works for polynomials over ``\ZZ_p`` and ``\QQ_p``
        - algorithm="om" use the implementation of an OM algorithm in sage

        OUTPUT:

        - A Sage :class:`Factorization`.

        EXAMPLES::

        See the irreducibility of ``x^32+16`` in ``\ZZ_2[x]``::

            sage: len((ZpFM(2,285)['x'](x^32+16)).factor())
            1

        Some times we can find a factoprization although the discriminant (to the give precision) is 0.

        We use an OM algorithm for factoring.

            sage: f = ZpFM(2,50,'terse')['x']( (x^32+16)*(x^32+16+2^16*x^2)+2^34 )
            sage: len(f.factor(algorithm="om",check_squarefree=False)) # long time (5.0s)
            2

        Monic and non-monic examples over ``\QQ_p``::

            sage: R.<t> = PolynomialRing(Qp(3,3,print_mode='terse',print_pos=False))
            sage: pol = t^8 - 1
            sage: for p,e in pol.factor():
            ....:     print("{} {}".format(e, p))
            1 (1 + O(3^3))*t + (1 + O(3^3))
            1 (1 + O(3^3))*t + (-1 + O(3^3))
            1 (1 + O(3^3))*t^2 + (5 + O(3^3))*t + (-1 + O(3^3))
            1 (1 + O(3^3))*t^2 + (-5 + O(3^3))*t + (-1 + O(3^3))
            1 (1 + O(3^3))*t^2 + (0 + O(3^3))*t + (1 + O(3^3))
            sage: R.<t> = PolynomialRing(Qp(5,6,print_mode='terse',print_pos=False))
            sage: pol = 100 * (5*t - 1) * (t - 5)
            sage: pol
            (500 + O(5^9))*t^2 + (-2600 + O(5^8))*t + (500 + O(5^9))
            sage: pol.factor()
            (500 + O(5^9)) * ((1 + O(5^5))*t + (-1/5 + O(5^5))) * ((1 + O(5^6))*t + (-5 + O(5^6)))
            sage: pol.factor().value()
            (500 + O(5^8))*t^2 + (-2600 + O(5^8))*t + (500 + O(5^8))

        The same factorization over ``\ZZ_p``. In this case, the "unit"
        part is a ``p``-adic unit and the power of ``p`` is considered to be
        a factor::

            sage: R.<t> = PolynomialRing(Zp(5,6,print_mode='terse',print_pos=False))
            sage: pol = 100 * (5*t - 1) * (t - 5)
            sage: pol
            (500 + O(5^9))*t^2 + (-2600 + O(5^8))*t + (500 + O(5^9))
            sage: pol.factor()
            (4 + O(5^6)) * ((5 + O(5^7)))^2 * ((1 + O(5^6))*t + (-5 + O(5^6))) * ((5 + O(5^6))*t + (-1 + O(5^6)))
            sage: pol.factor().value()
            (500 + O(5^8))*t^2 + (-2600 + O(5^8))*t + (500 + O(5^8))

        In the following example, the discriminant is zero, so the ``p``-adic
        factorization is not well defined::

            sage: factor(t^2)
            Traceback (most recent call last):
            ...
            PrecisionError: p-adic factorization not well-defined since the discriminant is zero up to the requested p-adic precision

        More examples over ``\ZZ_p``::

            sage: R.<w> = PolynomialRing(Zp(5, prec=6, type = 'capped-abs', print_mode = 'val-unit'))
            sage: f = w^5-1
            sage: f.factor()
            ((1 + O(5^6))*w + (3124 + O(5^6))) * ((1 + O(5^6))*w^4 + (12501 + O(5^6))*w^3 + (9376 + O(5^6))*w^2 + (6251 + O(5^6))*w + (3126 + O(5^6)))

        Over an unramified extension::

            sage: R.<c> = ZqFM(125, prec = 500)
            sage: Rz.<z>=R[]
            sage: g = (z^6+2)^25+5
            sage: Fg = g.factor()
            sage: [f[0].degree() for f in Fg]
            [50, 50, 50]

        Non-monic over an unramified extension::

            sage: R.<c> = Zq(125, prec = 500)
            sage: Rz.<z>=R[]
            sage: g = 25*z^15+5
            sage: Fg = g.factor()
            sage: [f[0].degree() for f in Fg]
            [15]

        Over a totally ramified extension::

            sage: R=Zp(5,500)
            sage: S.<x>=R[]
            sage: W.<w>=R.ext(x^5+25*x+5)
            sage: Wz.<z>=W[]
            sage: a = (z^5+25*z+5)*(z^5+w)+625;
            sage: Fa = a.factor()
            sage: [f[0].degree() for f in Fa]
            [1, 4, 5]

        See :trac:`4038`::

            sage: E = EllipticCurve('37a1')
            sage: K =Qp(7,10)
            sage: EK = E.base_extend(K)
            sage: E = EllipticCurve('37a1')
            sage: K = Qp(7,10)
            sage: EK = E.base_extend(K)
            sage: g = EK.division_polynomial_0(3)
            sage: g.factor()
            (3 + O(7^10)) * ((1 + O(7^10))*x + (1 + 2*7 + 4*7^2 + 2*7^3 + 5*7^4 + 7^5 + 5*7^6 + 3*7^7 + 5*7^8 + 3*7^9 + O(7^10))) * ((1 + O(7^10))*x^3 + (6 + 4*7 + 2*7^2 + 4*7^3 + 7^4 + 5*7^5 + 7^6 + 3*7^7 + 7^8 + 3*7^9 + O(7^10))*x^2 + (6 + 3*7 + 5*7^2 + 2*7^4 + 7^5 + 7^6 + 2*7^8 + 3*7^9 + O(7^10))*x + (2 + 5*7 + 4*7^2 + 2*7^3 + 6*7^4 + 3*7^5 + 7^6 + 4*7^7 + O(7^10)))

        TESTS:

        Check that :trac:`13293` is fixed::

            sage: R.<T> = Qp(3)[]
            sage: f = 1926*T^2 + 312*T + 387
            sage: f.factor()
            (3^2 + 2*3^3 + 2*3^4 + 3^5 + 2*3^6 + O(3^22)) * ((1 + O(3^19))*T + (2*3^-1 + 3 + 3^2 + 2*3^5 + 2*3^6 + 2*3^7 + 3^8 + 3^9 + 2*3^11 + 3^15 + 3^17 + O(3^19))) * ((1 + O(3^20))*T + (2*3 + 3^2 + 3^3 + 3^5 + 2*3^6 + 2*3^7 + 3^8 + 3^10 + 3^11 + 2*3^12 + 2*3^14 + 2*3^15 + 2*3^17 + 2*3^18 + O(3^20)))

        AUTHORS:

        - Jeroen Demeyer (2007)
        - Julian Rueth and David Roe (2013)
        - Sebastian Pauli and Brian Sinclair (2017-07-19): factorization over extensions, non monic, algorithm choice
        """
        def _pari_padic_factorization_to_sage(G, R, leading_coeff):
            r"""
            Given a PARI factorization matrix `G` representing a factorization
            of some polynomial in the `p`-adic polynomial ring `R`,
            return the corresponding Sage factorization. All factors in `G`
            are assumed to have content 1 (this is how PARI returns its
            factorizations).

            INPUT:

            - ``G`` -- PARI factorization matrix, returned by ``factorpadic``.

            - ``R`` -- polynomial ring to be used as parent ring of the factors

            - ``leading_coeff`` -- leading coefficient of the polynomial which
              was factored. This can belong to any ring which can be coerced
              into ``R.base_ring()``.

            OUTPUT:

            - A Sage :class:`Factorization`.

            """
            B = R.base_ring()
            p = B.prime()
            leading_coeff = B(leading_coeff)

            pols = [R(f, absprec=f.padicprec(p)) for f in G[0]]
            exps = [int(e) for e in G[1]]

            # Determine unit part (which is discarded by PARI)
            if B.is_field():
                # When the base ring is a field, we normalize
                # the irreducible factors so they have leading
                # coefficient 1.
                for i in range(len(pols)):
                    lc = pols[i].leading_coefficient()
                    lc = lc.lift_to_precision()  # Ensure we don't lose precision
                    pols[i] *= ~lc
            else:
                # When the base ring is not a field, we normalize
                # the irreducible factors so that the leading term
                # is a power of p.
                c, leading_coeff = leading_coeff.val_unit()
                for i in range(len(pols)):
                    v, upart = pols[i].leading_coefficient().val_unit()
                    upart = upart.lift_to_precision()  # Ensure we don't lose precision
                    pols[i] *= ~upart
                    c -= exps[i] * v
                if c:
                    # Add factor p^c
                    pols.append(R(p))
                    exps.append(c)

            return Factorization(zip(pols, exps), leading_coeff)

        if self == 0:
            raise ArithmeticError("factorization of {!r} is not defined".format(self))

        Rx = self.parent()
        R = Rx.base()

        self_normal, normalize_by = self.normalized(monic_if_possible=True)

        absprec = min([x.precision_absolute() for x in self_normal])
        if check_squarefree and self_normal.discriminant().valuation() >= absprec:
            raise PrecisionError("p-adic factorization not well-defined since the discriminant is zero up to the requested p-adic precision")

        # Factor with pari
        if algorithm=="pari" or (algorithm=="auto" and R == R.base()):
            G = self_normal.__pari__().factorpadic(self.base_ring().prime(), absprec)
            return _pari_padic_factorization_to_sage(G, self.parent(), self.leading_coefficient())

        # The polynomial self_normal is not monic, so we perform a monic transform
        # in the method _factor_non_monic and re-call this factor method.
        if not self_normal.is_monic() and (algorithm=="auto" or algorithm=="om"):
            return self_normal._factor_non_monic(normalize_by,check_squarefree=check_squarefree,algorithm=algorithm)

        # Factoring with OM is straightforward if monic
        if self_normal.is_monic() and (algorithm=="auto" or algorithm=="om"):
            om = self_normal.omtree(check_squarefree=check_squarefree,check_field=False)
            omf = om.factorization()
            return Factorization([[(gm[0]).move(Rx),gm[1]] for gm in omf],normalize_by)

        raise NotImplementedError("No factorization method was selected for the given input")

    def _factor_non_monic(self,normalize_by=None,check_squarefree=True,algorithm="auto"):
        r"""
        Factor non-monic polynomials over valuation rings.

        The valuation of the content of the polynomial has to be 0.

        EXAMPLES::

            sage: R.<x> = ZpFM(3,7)[]
            sage: f = 4*x^2 + 14*x + 8
            sage: f._factor_non_monic()
            (1 + 3 + O(3^7)) * ((1 + O(3^7))*x^2 + (2 + 2*3 + 3^2 + 3^3 + 3^4 + 3^5 + 3^6 + O(3^7))*x + (2 + O(3^7)))

            sage: g = 3*x^2 + 8
            sage: g._factor_non_monic()
            (3 + O(3^7))*x^2 + (O(3^7))*x + (2 + 2*3 + O(3^7))

            sage: h = 3*x^2 + 15
            sage: h._factor_non_monic()
            (3 + O(3^7)) * ((1 + O(3^7))*x + (2 + 2*3^2 + 3^4 + 3^6 + O(3^7))) * ((1 + O(3^7))*x + (1 + 2*3 + 2*3^3 + 3^4 + 2*3^5 + 3^6 + O(3^7)))
        """

        def make_monic(f):
            r"""
            Transform a polynomial over a local ring into a monic one for
            factoring.

            INPUT:

            - ``f`` -- a polynomial over a local ring
            """
            Kx = f.parent()
            K = Kx.base()
            pi = K.uniformizer()

            if f.leading_coefficient() == K(1):
                return f,K(1),0,0

            minval = min([a.valuation() for a in f])
            divi = pi ** minval
            g = Kx([a // divi for a in f])
            multval = -minval

            lead = g.leading_coefficient()
            vlead = lead.valuation()
            exp = -1
            expmult = g.degree()*exp-vlead
            if not g.degree() <= 0:
                while expmult < 0 or min([g[i].valuation()-exp*i+expmult for i in range(1,g.degree()+1)]) < 0:
                    exp += 1
                    expmult = g.degree()*exp-vlead
            else:
                exp = 0
                expmult = -vlead

            g = g * pi ** expmult
            g = Kx([g.coefficients(sparse=False)[i] // pi**(exp*i) for i in range(g.degree()+1)])
            multval = multval + expmult

            uni = g.leading_coefficient()
            g = Kx([a // uni for a in g.coefficients(sparse=False)])

            return g,uni,multval,exp

        def undo_monic(f,exp,h):
            r"""
            Transform back to a monic polynomial
            """
            return h(f.parent().base().uniformizer() ** exp * f.parent().gen())

        def dist_den(facts,multval):
            r"""
            Distribute the (power of the uniformizer) denominator over the factors from the factorization.

            INPUT:

            - ``facts`` -- list of factors
            - ``multval`` -- power of pi in the denominator of the constant
                             in front of our factorization
            """
            if len(facts) == 0:
                return []
            Kx = facts[0].parent()
            K = Kx.base_ring()
            Kpi = facts[0].parent().base_ring().uniformizer()
            newfacts = []
            for fact in facts:
                m = min([a.valuation() for a in fact])
                newfacts.append(Kx([a // Kpi**m for a in fact]))
                multval -= m
            if multval != 0:
                raise StandardError, "Power of pi did not distribute properly over factors of non-monic polynomial."
            return newfacts

        # non monic factor
        f = self
        if normalize_by == None:
            f, normalize_by = f.normalized(monic_if_possible=True)

        Kx = f.parent()
        g,uni,multval,exp = make_monic(f)
        if multval > 0:
            L = f.base_ring().change(prec=f.base_ring().precision_cap()+2*multval,type='fixed-mod',field=False)
            Ly = PolynomialRing(L,names='y')
            f = f.move(Ly)
            g,uni,multval,exp = make_monic(f)
        facts = [fact[0] for fact in g.factor(check_squarefree=check_squarefree,algorithm=algorithm)]
        facts = [undo_monic(f,exp,fact) for fact in facts]
        facts = dist_den(facts,multval)
        if multval > 0:
            uni = Ly(uni).move(Kx)
            facts = [fact.move(Kx) for fact in facts]
        if multval < 0:
            return Factorization([(f.base_ring().uniformizer(),-multval)] + [(fact,1) for fact in facts],unit=uni*normalize_by)
        return Factorization([(fact,1) for fact in facts],unit=uni*normalize_by)
