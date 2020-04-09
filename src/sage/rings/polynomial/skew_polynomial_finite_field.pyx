r"""
Univariate Dense Skew Polynomials over Finite Fields

This module provides the :class:`~sage.rings.polynomial.skew_polynomial_finite_field.SkewPolynomial_finite_field_dense`
which constructs a single univariate skew polynomial over a finite field equipped with the Frobenius
Endomorphism.

AUTHOR::

- Xavier Caruso (2012-06-29): initial version

- Arpit Merchant (2016-08-04): improved docstrings, fixed doctests and refactored classes and methods

"""

#############################################################################
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

#from cysignals.signals cimport sig_on, sig_off
def sig_on(): pass
def sig_off(): pass

import copy
import cysignals
from sage.matrix.constructor import zero_matrix
from sage.rings.ring cimport Ring
from sage.matrix.matrix_dense cimport Matrix_dense
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.all import ZZ
from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.rings.integer cimport Integer
from sage.structure.element cimport RingElement
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.skew_polynomial_element cimport SkewPolynomial_generic_dense
from sage.rings.polynomial.skew_polynomial_finite_order cimport SkewPolynomial_finite_order_dense

from sage.combinat.permutation import Permutation, Permutations
from sage.combinat.partition import Partition
from sage.structure.factorization import Factorization
from sage.misc.mrange import xmrange_iter
from sage.arith.misc import factorial
from sage.combinat.q_analogues import q_jordan


cdef class SkewPolynomial_finite_field_dense(SkewPolynomial_finite_order_dense):
    def reduced_norm_factor(self):
        """
        Return the reduced norm of this polynomial
        factorized in the centre.

        EXAMPLES:

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]

            sage: a = (x^2 + 1) * (x+3)
            sage: a.reduced_norm_factor()
            ((x^3) + 3) * ((x^3) + 2)^2
        """
        if self._norm_factor is None:
            self._norm_factor = self.reduced_norm().factor()
        return self._norm_factor

    def is_irreducible(self):
        """
        Return true if this skew polynomial is irreducible.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]

            sage: a = x^2 + t*x + 1
            sage: a.is_irreducible()
            False
            sage: a.factor()
            (x + 4*t^2 + 4*t + 1) * (x + 3*t + 2)

            sage: a = x^2 + t*x + t + 1
            sage: a.is_irreducible()
            True
            sage: a.factor()
            x^2 + t*x + t + 1

        Skew polynomials of degree `1` are of course irreducible::

            sage: a = x + t
            sage: a.is_irreducible()
            True

        A random irreducible skew polynomial is irreducible::

            sage: a = S.random_irreducible(degree=4,monic=True); a   # random
            x^4 + (t + 1)*x^3 + (3*t^2 + 2*t + 3)*x^2 + 3*t*x + 3*t
            sage: a.is_irreducible()
            True

        By convention, constant skew polynomials are not irreducible::

            sage: S(1).is_irreducible()
            False
            sage: S(0).is_irreducible()
            False
        """
        return self.reduced_norm().is_irreducible()


    def type(self,N):
        """
        INPUT:

        -  ``N`` -- an irreducible polynomial in the
           center of the underlying skew polynomial ring

        OUTPUT:

        The `N`-type of this skew polynomial

        .. NOTE::

            The result is cached.

        DEFINITION:

        The `N`-type of a skew polynomial `a` is the Partition
        `(t_0, t_1, t_2, ...)` defined by

        .. MATH::

            t_0 + \cdots + t_i = \frac{\deg gcd(a,N^i)}{\deg N}

        where `\deg N` is the degree of `N` considered as an
        element in the center.

        This notion has an important mathematic interest because
        it corresponds to the Jordan type of the `N`-typical part
        of the associated Galois representation.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: Z = S.center(); x3 = Z.gen()

            sage: a = x^4 + x^3 + (4*t^2 + 4)*x^2 + (t^2 + 2)*x + 4*t^2
            sage: N = x3^2 + x3 + 1
            sage: a.type(N)
            [1]
            sage: N = x3 + 1
            sage: a.type(N)
            [2]

            sage: a = x^3 + (3*t^2 + 1)*x^2 + (3*t^2 + t + 1)*x + t + 1
            sage: N = x3 + 1
            sage: a.type(N)
            [2, 1]

        If `N` does not divide the reduced map of `a`, the type
        is empty::

            sage: N = x3 + 2
            sage: a.type(N)
            []

        If `a = N`, the type is just `[r]` where `r` is the order
        of the twist map ``Frob``::

            sage: N = x3^2 + x3 + 1
            sage: S(N).type(N)
            [3]

        `N` must be irreducible::

            sage: N = (x3 + 1) * (x3 + 2)
            sage: a.type(N)
            Traceback (most recent call last):
            ...
            ValueError: N is not irreducible
        """
        try:
            return self._types[N]
        except (KeyError, TypeError):
            if not N.is_irreducible():
                raise ValueError("N is not irreducible")
            skew_ring = self._parent
            if self._norm_factor is None:
                m = -1
            else:
                i = [ n for n,_ in self._norm_factor ].index(N)
                m = self._norm_factor[i][1]
            NS = skew_ring(N)
            type = [ ]
            degN = N.degree()
            while True:
                d = self.right_gcd(NS)
                deg = d.degree()/degN
                if deg == 0:
                    break
                if m >= 0:
                    if deg == 1:
                        type += m * [1]
                        break
                    m -= deg
                self = self // d
                type.append(deg)
            # type = Partition(type)
            if self._types is None:
                self._types = { N: type }
            else:
                self._types[N] = type
            return type


    # Finding divisors
    # ----------------

    cdef SkewPolynomial_finite_field_dense _rdivisor_c(P, N):
        """
        cython procedure computing an irreducible monic right divisor
        of `P` whose reduced norm is `N`

        .. WARNING::

            `N` needs to be an irreducible factor of the
            reduced norm of `P`. This function does not check
            this (and his behaviour is not defined if the
            require property doesn't hold).
        """
        cdef skew_ring = P._parent
        cdef Py_ssize_t d = N.degree()
        cdef Py_ssize_t e = P.degree()/d
        cdef SkewPolynomial_finite_field_dense D
        if e == 1:
            D = <SkewPolynomial_finite_field_dense>P._new_c(list(P._coeffs), skew_ring)
            D = D.right_monic()
            return D

        center = N.parent()
        E = center.quo(N)
        PE = PolynomialRing(E, name='T')
        cdef Integer exp
        if skew_ring.characteristic() != 2:
            exp = Integer((E.cardinality()-1)/2)
        cdef SkewPolynomial_finite_field_dense NS = <SkewPolynomial_finite_field_dense>skew_ring(N)
        cdef SkewPolynomial_finite_field_dense Q = <SkewPolynomial_finite_field_dense>(NS // P)
        cdef SkewPolynomial_finite_field_dense R, X
        cdef Matrix_dense M = <Matrix_dense?>MatrixSpace(E,e,e)(0)
        cdef Matrix_dense V = <Matrix_dense?>MatrixSpace(E,e,1)(0)
        cdef Matrix_dense W
        cdef Py_ssize_t i, j, t, r = skew_ring._order
        cdef Polynomial dd, xx, yy, zz

        while True:
            R = <SkewPolynomial_finite_field_dense>skew_ring.random_element((e*r-1,e*r-1))
            R = Q*R
            X = <SkewPolynomial_finite_field_dense>Q._new_c(Q._coeffs[:],Q._parent)
            for j from 0 <= j < e:
                for i from 0 <= i < e:
                    coeffs = [skew_ring._retraction(X[t*r+i]) for t in range(d)]
                    value = E(coeffs)
                    M.set_unsafe(i, j, value)
                X = (R*X) % NS
            for i from 0 <= i < e:
                V.set_unsafe(i, 0, E([skew_ring._retraction(X[t*r+i]) for t in range(d)]))
            W = M._solve_right_nonsingular_square(V)
            if M*W != V:
                skew_ring._new_retraction_map()
                continue
            xx = PE(W.list()+[E(-1)])
            if skew_ring.characteristic() == 2:
                yy = PE.gen()
                zz = PE.gen()
                for i from 1 <= i < d:
                    zz = (zz*zz) % xx
                    yy += zz
                dd = xx.gcd(yy)
                if dd.degree() != 1: continue
            else:
                yy = PE.gen().__pow__(exp,xx) - 1
                dd = xx.gcd(yy)
                if dd.degree() != 1:
                    yy += 2
                    dd = xx.gcd(yy)
                    if dd.degree() != 1: continue
            D = P.right_gcd(R + skew_ring(center((dd[0]/dd[1]).list())))
            if D.degree() == 0:
                continue
            D = D.right_monic()
            return D


    def irreducible_divisor(self,right=True,distribution=None):
        """
        INPUT:

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        -  ``distribution`` -- None (default) or ``uniform``

           - None: no particular specification

           - ``uniform``: the returned irreducible divisor is
             uniformly distributed

        .. NOTE::

            ``uniform`` is a little bit slower.

        OUTPUT:

        -  an irreducible monic ``side`` divisor of ``self``

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^6 + 3*t*x^5 + (3*t + 1)*x^3 + (4*t^2 + 3*t + 4)*x^2 + (t^2 + 2)*x + 4*t^2 + 3*t + 3

            sage: dr = a.irreducible_divisor(); dr  # random
            x^3 + (2*t^2 + t + 4)*x^2 + (4*t + 1)*x + 4*t^2 + t + 1
            sage: a.is_divisible_by(dr)
            True

            sage: dl = a.irreducible_divisor(side=Left); dl  # random
            x^3 + (2*t^2 + t + 1)*x^2 + (4*t^2 + 3*t + 3)*x + 4*t^2 + 2*t + 1
            sage: a.is_divisible_by(dl,side=Left)
            True

        Right divisors are cached. Hence, if we ask again for a
        right divisor, we will get the same answer::

            sage: a.irreducible_divisor()  # random
            x^3 + (2*t^2 + t + 4)*x^2 + (4*t + 1)*x + 4*t^2 + t + 1

        However the algorithm is probabilistic. Hence, if we first
        reinitialiaze `a`, we may get a different answer::

            sage: a = x^6 + 3*t*x^5 + (3*t + 1)*x^3 + (4*t^2 + 3*t + 4)*x^2 + (t^2 + 2)*x + 4*t^2 + 3*t + 3
            sage: a.irreducible_divisor()  # random
            x^3 + (t^2 + 3*t + 4)*x^2 + (t + 2)*x + 4*t^2 + t + 1

        We can also generate uniformly distributed irreducible monic
        divisors as follows::

            sage: a.irreducible_divisor(distribution="uniform")  # random
            x^3 + (4*t + 2)*x^2 + (2*t^2 + 2*t + 2)*x + 2*t^2 + 2
            sage: a.irreducible_divisor(distribution="uniform")  # random
            x^3 + (t^2 + 2)*x^2 + (3*t^2 + 1)*x + 4*t^2 + 2*t
            sage: a.irreducible_divisor(distribution="uniform")  # random
            x^3 + x^2 + (4*t^2 + 2*t + 4)*x + t^2 + 3

        By convention, the zero skew polynomial has no irreducible
        divisor:

            sage: S(0).irreducible_divisor()
            Traceback (most recent call last):
            ...
            ValueError: 0 has no irreducible divisor
        """
        if self.is_zero():
            raise ValueError("0 has no irreducible divisor")
        if not (distribution is None or distribution == "uniform"):
            raise ValueError("distribution must be None or 'uniform'")
        if distribution == "uniform":
            skew_ring = self._parent
            center = skew_ring.center()
            cardcenter = center.base_ring().cardinality()
            gencenter = center.gen()
            count = [ ]
            total = 0
            F = self.reduced_norm_factor()
            for n,_ in F:
                if n == gencenter:
                    total += 1
                else:
                    degn = n.degree()
                    P = self.right_gcd(skew_ring(n))
                    m = P.degree()/degn
                    cardL = cardcenter**degn
                    total += (cardL**m - 1)/(cardL - 1)
                count.append(total)
            if total == 0:
                raise ValueError("No irreducible divisor having given reduced norm")
            random = ZZ.random_element(total)
            for i in range(len(F)):
                if random < count[i]:
                    N = F[i][0]
                    break
        else:
            N = self.reduced_norm_factor()[0][0]
        return self.irreducible_divisor_with_norm(N, right=right, distribution=distribution)


    def irreducible_divisor_with_norm(self, N, right=True, distribution=None):
        """
        INPUT:

        -  ``N`` -- an irreducible polynomial in the center
           of the underlying skew polynomial ring

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        -  ``distribution`` -- None (default) or ``uniform``

           - None: no particular specification

           - ``uniform``: the returned irreducible divisor is
             uniformly distributed

        .. NOTE::

            ``uniform`` is a little bit slower.

        OUTPUT:

        -  an irreducible monic ``side`` divisor of ``self``
           whose reduced norm is similar to `N` (i.e. `N` times
           a unit).

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: Z = S.center(); x3 = Z.gen()
            sage: a = x^6 + 3*x^3 + 2

            sage: d1 = a.irreducible_divisor_with_norm(x3+1); d1   # random
            x + t^2 + 3*t
            sage: a.is_divisible_by(d1)
            True
            sage: d1.reduced_norm()
            (x^3) + 1

            sage: d2 = a.irreducible_divisor_with_norm(x3+2); d2   # random
            x + 2*t^2 + 3*t + 2
            sage: a.is_divisible_by(d2)
            True
            sage: d2.reduced_norm()
            (x^3) + 2

            sage: d3 = a.irreducible_divisor_with_norm(x3+3)
            Traceback (most recent call last):
            ...
            ValueError: No irreducible divisor having given reduced norm

        We can also generate uniformly distributed irreducible monic
        divisors as follows::

            sage: a.irreducible_divisor_with_norm(x3+1,distribution="uniform")   # random
            x + 3*t^2 + 3*t + 1
            sage: a.irreducible_divisor_with_norm(x3+1,distribution="uniform")   # random
            x + 1
            sage: a.irreducible_divisor_with_norm(x3+1,distribution="uniform")   # random
            x + 2*t^2 + 4*t
        """
        cdef SkewPolynomial_finite_field_dense cP1
        if self.is_zero():
            raise "No irreducible divisor having given reduced norm"
        skew_ring = self._parent
        center = skew_ring.center()
        try:
            N = center(N)
        except TypeError:
            raise TypeError("N must be a polynomial in the center")
        cardcenter = center.base_ring().cardinality()
        gencenter = center.gen()

        if N == gencenter:
            if self[0] == 0:
                return skew_ring.gen()
            else:
                raise ValueError("No irreducible divisor having given reduced norm")

        D = None
        try:
            D = self._rdivisors[N]
        except (KeyError, TypeError):
            if N.is_irreducible():
                cP1 = <SkewPolynomial_finite_field_dense>self.right_gcd(self._parent(N))
                cN = N
                if cP1.degree() > 0:
                    D = cP1._rdivisor_c(cN)
            if self._rdivisors is None:
                self._rdivisors = { N: D }
            else:
                self._rdivisors[N] = D
            distribution = ""
        if D is None:
            raise ValueError("No irreducible divisor having given reduced norm")

        NS = self._parent(N)
        degN = N.degree()
        if right:
            if distribution == "uniform":
                P1 = self.right_gcd(NS)
                if P1.degree() != degN:
                    Q1 = NS // P1
                    deg = P1.degree()-1
                    while True:
                        R = Q1*skew_ring.random_element((deg,deg))
                        if P1.right_gcd(R) == 1:
                            break
                    D = P1.right_gcd(D*R)
            return D
        else:
            deg = NS.degree() - 1
            P1 = self.left_gcd(NS)
            while True:
                if distribution == "uniform":
                    while True:
                        R = skew_ring.random_element((deg,deg))
                        if NS.right_gcd(R) == 1:
                            break
                    D = NS.right_gcd(D*R)
                Dp = NS // D
                LDp = P1.right_gcd(Dp)
                LD = P1 // LDp
                if LD.degree() == degN:
                    return LD
                distribution = "uniform"


    def irreducible_divisors(self,right=True):
        """
        INPUT:

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        OUTPUT:

        An iterator over all irreducible monic ``side`` divisors
        of this skew polynomial

        EXAMPLES:

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^4 + 2*t*x^3 + 3*t^2*x^2 + (t^2 + t + 1)*x + 4*t + 3
            sage: iter = a.irreducible_divisors(); iter
            <generator object at 0x...>
            sage: iter.next()   # random
            x + 2*t^2 + 4*t + 4
            sage: iter.next()   # random
            x + 3*t^2 + 4*t + 1

        We can use this function to build the list of all monic
        irreducible divisors of `a`::

            sage: rightdiv = [ d for d in a.irreducible_divisors() ]
            sage: leftdiv = [ d for d in a.irreducible_divisors(side=Left) ]

        We do some checks::

            sage: len(rightdiv) == a.count_irreducible_divisors()
            True
            sage: len(rightdiv) == len(Set(rightdiv))  # check no duplicates
            True
            sage: for d in rightdiv:
            ...       if not a.is_divisible_by(d):
            ...           print "Found %s which is not a right divisor" % d
            ...       elif not d.is_irreducible():
            ...           print "Found %s which is not irreducible" % d

            sage: len(leftdiv) == a.count_irreducible_divisors(side=Left)
            True
            sage: len(leftdiv) == len(Set(leftdiv))  # check no duplicates
            True
            sage: for d in leftdiv:
            ...       if not a.is_divisible_by(d,side=Left):
            ...           print "Found %s which is not a left divisor" % d
            ...       elif not d.is_irreducible():
            ...           print "Found %s which is not irreducible" % d

        Note that left divisors and right divisors differ::

            sage: Set(rightdiv) == Set(leftdiv)
            False

        Note that the algorithm is probabilistic. As a consequence, if we
        build again the list of right monic irreducible divisors of `a`, we
        may get a different ordering::

            sage: rightdiv2 = [ d for d in a.irreducible_divisors() ]
            sage: rightdiv == rightdiv2
            False
            sage: Set(rightdiv) == Set(rightdiv2)
            True
        """
        return self._irreducible_divisors(right)


    def _irreducible_divisors(self, right):
        """
        Return an iterator over all irreducible monic
        divisors of this skew polynomial.

        Do not use this function. Use instead
        ``self.irreducible_divisors()``.
        """
        if self.is_zero():
            return
        skew_ring = self._parent
        center = skew_ring.center()
        kfixed = center.base_ring()
        F = self.reduced_norm_factor()
        for N,_ in F:
            if N == center.gen():
                yield skew_ring.gen()
                continue
            degN = N.degree()
            NS = skew_ring(N)
            if right:
                P = self.right_gcd(NS)
            else:
                P = self.left_gcd(NS)
            m = P.degree()/degN
            if m == 1:
                yield P
                continue
            degrandom = P.degree() - 1
            if right:
                Q,_ = NS.right_quo_rem(P)
                P1 = self.irreducible_divisor_with_norm(N,right=right)
                Q1,_ = P.right_quo_rem(P1)
            else:
                Q,_ = NS.left_quo_rem(P)
                P1 = self.irreducible_divisor_with_norm(N,right=right)
                Q1,_ = P.left_quo_rem(P1)
            while True:
                R = skew_ring.random_element((degrandom,degrandom))
                if right:
                    _, g = (R*Q).left_quo_rem(P)
                    if g.left_gcd(P) != 1: continue
                else:
                    _, g = (R*Q).right_quo_rem(P)
                    if g.right_gcd(P) != 1: continue
                L = Q1
                V = L
                for i in range(1,m):
                    if right:
                        L = (g*L).left_gcd(P)
                        V = V.left_gcd(L)
                    else:
                        L = (L*g).right_gcd(P)
                        V = V.right_gcd(L)
                if V == 1: break
            rng = xmrange_iter([kfixed]*degN,center)
            for i in range(m):
                for pol in xmrange_iter([rng]*i):
                    f = skew_ring(1)
                    for j in range(i):
                        coeff = pol.pop()
                        f = (g*f+coeff).rem(P, right=not right)
                    if right:
                        d = (f*Q1).left_gcd(P)
                        d, _ = P.left_quo_rem(d)
                    else:
                        d = (Q1*f).right_gcd(P)
                        d, _ = P.right_quo_rem(d)
                    yield d


    def count_irreducible_divisors(self,right=True):
        """
        INPUT:

        -  ``side`` -- ``Left`` or ``Right`` (default: Right)

        OUTPUT:

        The number of irreducible monic ``side`` divisors of
        this skew polynomial.

        .. NOTE::

            Actually, one can prove that there are always as
            many left irreducible monic divisors as right
            irreducible monic divisors.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]

        We illustrate that a skew polynomial may have a number of irreducible
        divisors greater than its degree.

            sage: a = x^4 + (4*t + 3)*x^3 + t^2*x^2 + (4*t^2 + 3*t)*x + 3*t
            sage: a.count_irreducible_divisors()
            12
            sage: a.count_irreducible_divisors(side=Left)
            12

        We illustrate that an irreducible polynomial in the center have
        in general a lot of irreducible divisors in the skew polynomial
        ring::

            sage: Z = S.center(); x3 = Z.gen()
            sage: N = x3^5 + 4*x3^4 + 4*x3^2 + 4*x3 + 3; N
            (x^3)^5 + 4*(x^3)^4 + 4*(x^3)^2 + 4*(x^3) + 3
            sage: N.is_irreducible()
            True
            sage: S(N).count_irreducible_divisors()
            9768751
        """
        if self.is_zero():
            return 0
        skew_ring = self.parent()
        cardcenter = skew_ring.center().base_ring().cardinality()
        gencenter = skew_ring.center().gen()
        F = self.reduced_norm_factor()
        val = self.valuation()
        self >>= val
        count = 0
        if val > 0:
            count = 1
        for N,_ in F:
            if N == gencenter:
                continue
            degN = N.degree()
            if right:
                P = self.right_gcd(skew_ring(N))
            else:
                P = self.left_gcd(skew_ring(N))
            m = P.degree()/degN
            cardL = cardcenter**degN
            count += (cardL**m - 1)/(cardL - 1)
        return count


    # Finding factorizations
    # ----------------------

    cdef _factor_c(self):
        """
        Compute a factorization of ``self``
        """
        cdef skew_ring = self._parent
        cdef SkewPolynomial_finite_field_dense poly = self.right_monic()
        cdef list a = poly._coeffs
        cdef Py_ssize_t val = 0
        if len(a) < 0:
            return Factorization([], sort=False, unit=skew_ring.zero_element())
        while a[0].is_zero():
            del a[0]
            val += 1

        cdef Py_ssize_t degQ, degrandom, m, mP, i
        cdef N
        cdef list factors = [ (skew_ring.gen(), val) ]
        cdef SkewPolynomial_finite_field_dense P, Q, P1, NS, g, right, Pn
        cdef RingElement unit = <RingElement>self.leading_coefficient()
        cdef Polynomial gencenter = skew_ring.center().gen()
        cdef Py_ssize_t p = skew_ring.characteristic()
        cdef F = self.reduced_norm_factor()

        for N, m in F:
            if N == gencenter:
                continue
            degN = N.degree()
            if poly.degree() == degN:
                factors.append((poly, 1))
                break
            NS = <SkewPolynomial_finite_field_dense>skew_ring(N)
            P1 = None
            while True:
                P = <SkewPolynomial_finite_field_dense>poly.right_gcd(NS)
                P = P.right_monic()
                mP = P.degree() / degN
                if mP == 0: break
                if mP == 1:
                    factors.append((P,1))
                    poly = poly // P
                    for i from 1 <= i < m:
                        if poly.degree() == degN:
                            factors.append((poly,1))
                            break
                        P = poly.right_gcd(NS).right_monic()
                        factors.append((P, 1))
                        poly = poly // P
                    break
                if P1 is None:
                    P1 = P._rdivisor_c(N)
                Q = <SkewPolynomial_finite_field_dense>NS._new_c(NS._coeffs[:], NS._parent)
                Q = P1 * (Q // P)
                factors.append((P1, 1))
                right = <SkewPolynomial_finite_field_dense>P1._new_c(P1._coeffs[:], P1._parent)
                m -= (mP-1)
                degrandom = P.degree()
                while mP > 2:
                    while True:
                        g = <SkewPolynomial_finite_field_dense>skew_ring.random_element((degrandom, degrandom))
                        g = (Q*g).right_gcd(P)
                        Pn = right._left_lcm_cofactor(g)
                        if Pn.degree() == degN: break
                    Pn = Pn.right_monic()
                    factors.append((Pn, 1))
                    right = Pn * right
                    degrandom -= degN
                    mP -= 1
                poly = poly // right
                P1, _ = P.right_quo_rem(right)
        factors.reverse()
        return Factorization(factors, sort=False, unit=unit)


    cdef _factor_uniform_c(self):
        """
        Compute a uniformly distrbuted factorization of ``self``
        """
        skew_ring = self._parent
        cdef Integer cardE, cardcenter = skew_ring.center().base_ring().cardinality()
        cdef gencenter = skew_ring.center().gen()
        cdef SkewPolynomial_finite_field_dense gen = <SkewPolynomial_finite_field_dense>skew_ring.gen()

        cdef list factorsN = [ ]
        cdef dict dict_divisor = { }
        cdef dict dict_type = { }
        cdef dict dict_right = { }
        cdef Py_ssize_t m
        cdef list type

        for N,m in self.reduced_norm_factor():
            factorsN += m * [N]
            if N == gencenter: continue
            type = list(self.type(N))
            dict_type[N] = type
            if type[0] > 1:
                dict_divisor[N] = self.irreducible_divisor_with_norm(N)
                dict_right[N] = skew_ring(1)
        cdef list indices = list(Permutations(len(factorsN)).random_element())

        cdef RingElement unit = self.leading_coefficient()
        cdef SkewPolynomial_finite_field_dense left = self._new_c(self._coeffs[:],skew_ring)
        left = left.right_monic()
        cdef SkewPolynomial_finite_field_dense right = <SkewPolynomial_finite_field_dense>skew_ring(1)
        cdef SkewPolynomial_finite_field_dense L, R
        cdef SkewPolynomial_finite_field_dense NS, P, Q, D, D1, D2, d
        cdef list factors = [ ]
        cdef list maxtype
        cdef Py_ssize_t i, j, degN, deg
        cdef count, maxcount

        for i in indices:
            N = factorsN[i-1]
            if N == gencenter:
                D1 = gen
            else:
                type = dict_type[N]
                NS = skew_ring(N)
                P = left.right_gcd(NS)
                if type[0] == 1:
                    D1 = P
                else:
                    R = right._new_c(right._coeffs[:],skew_ring)
                    R = R // dict_right[N]
                    D = R._left_lcm_cofactor(dict_divisor[N])
                    maxtype = list(type)
                    maxtype[-1] -= 1
                    degN = N.degree()
                    cardE = cardcenter ** degN
                    maxcount = q_jordan(Partition(maxtype),cardE)
                    Q = NS // P
                    deg = P.degree()-1
                    while 1:
                        while 1:
                            R = <SkewPolynomial_finite_field_dense>skew_ring.random_element((deg,deg))
                            R = Q*R
                            if P.right_gcd(R).degree() == 0:
                                break
                        D1 = P.right_gcd(D*R).right_monic()

                        L = left._new_c(list(left._coeffs),skew_ring)
                        L = L // D1
                        degN = N.degree()
                        for j in range(len(type)):
                            if type[j] == 1:
                                newtype = type[:-1]
                                break
                            d = L.right_gcd(NS).right_monic()
                            deg = d.degree() / degN
                            if deg < type[j]:
                                newtype = type[:]
                                newtype[j] = deg
                                break
                            L = L // d
                        count = q_jordan(Partition(newtype),cardE)
                        if ZZ.random_element(maxcount) < count:
                            break
                    dict_type[N] = newtype

                    D2 = D._new_c(list(D._coeffs),skew_ring)
                    D2 = D2.right_monic()
                    while D2 == D1:
                        while True:
                            R = <SkewPolynomial_finite_field_dense>skew_ring.random_element((deg,deg))
                            R = Q*R
                            if P.right_gcd(R).degree() == 0:
                                break
                        D2 = P.right_gcd(D*R).right_monic()
                    dict_divisor[N] = D1.left_lcm(D2)
            factors.append((D1,1))
            left = left // D1
            right = D1 * right
            dict_right[N] = right._new_c(list(right._coeffs),skew_ring)

        factors.reverse()
        return Factorization(factors,sort=False,unit=unit)


    def factor(self,distribution=None):
        """
        Return a factorization of this skew polynomial.

        INPUT:

        -  ``distribution`` -- None (default) or ``uniform``

           - None: no particular specification

           - ``uniform``: the returned factorization is uniformly
             distributed among all possible factorizations

        .. NOTE::

            ``uniform`` is a little bit slower.

        OUTPUT:

        -  ``Factorization`` -- a factorization of self as a
           product of a unit and a product of irreducible monic
           factors

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^3 + (t^2 + 4*t + 2)*x^2 + (3*t + 3)*x + t^2 + 1
            sage: F = a.factor(); F  # random
            (x + t^2 + 4) * (x + t + 3) * (x + t)
            sage: F.value() == a
            True

        The result of the factorization is cached. Hence, if we try
        again to factor `a`, we will get the same answer::

            sage: a.factor()  # random
            (x + t^2 + 4) * (x + t + 3) * (x + t)

        However, the algorithm is probabilistic. Hence if we first
        reinitialiaze `a`, we may get a different answer::

            sage: a = x^3 + (t^2 + 4*t + 2)*x^2 + (3*t + 3)*x + t^2 + 1
            sage: F = a.factor(); F   # random
            (x + t^2 + t + 2) * (x + 2*t^2 + t + 4) * (x + t)
            sage: F.value() == a
            True

        There is no guarantee on the distribution of the factorizations
        we get that way. (For this particular `a` for example, we get the
        uniform distribution on the subset of all factorizations ending
        by the factor `x + t`.)

        If we rather want uniform distribution among all factorizations,
        we need to specify it as follows::

            sage: a.factor(distribution="uniform")   # random
            (x + t^2 + 4) * (x + t) * (x + t + 3)
            sage: a.factor(distribution="uniform")   # random
            (x + 2*t^2) * (x + t^2 + t + 1) * (x + t^2 + t + 2)
            sage: a.factor(distribution="uniform")   # random
            (x + 2*t^2 + 3*t) * (x + 4*t + 2) * (x + 2*t + 2)

        By convention, the zero skew polynomial has no factorization:

            sage: S(0).factor()
            Traceback (most recent call last):
            ...
            ValueError: factorization of 0 not defined
        """
        if not (distribution is None or distribution == "uniform"):
            raise ValueError("distribution must be None or 'uniform'")
        if self.is_zero():
            raise ValueError("factorization of 0 not defined")
        sig_on()
        if distribution is None:
            if self._factorization is None:
                self._factorization = self._factor_c()
            F = self._factorization
        else:
            F = self._factor_uniform_c()
            if self._factorization is None:
                self._factorization = F
        sig_off()
        return F


    def count_factorizations(self):
        """
        Return the number of factorizations (as a product of a
        unit and a product of irreducible monic factors) of this
        skew polynomial.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^4 + (4*t + 3)*x^3 + t^2*x^2 + (4*t^2 + 3*t)*x + 3*t
            sage: a.count_factorizations()
            216

        We illustrate that an irreducible polynomial in the center have
        in general a lot of distinct factorizations in the skew polynomial
        ring::

            sage: Z = S.center(); x3 = Z.gen()
            sage: N = x3^5 + 4*x3^4 + 4*x3^2 + 4*x3 + 3; N
            (x^3)^5 + 4*(x^3)^4 + 4*(x^3)^2 + 4*(x^3) + 3
            sage: N.is_irreducible()
            True
            sage: S(N).count_factorizations()
            30537115626
        """
        if self.is_zero():
            raise ValueError("factorization of 0 not defined")
        cardcenter = self._parent.center().base_ring().cardinality()
        gencenter = self._parent.center().gen()
        F = self.reduced_norm_factor()
        summ = 0
        count = 1
        for N,m in F:
            summ += m
            if m == 1: continue
            if N != gencenter:
                count *= q_jordan(self.type(N),cardcenter**N.degree())
            count /= factorial(m)
        return count * factorial(summ)


    # Not optimized:
    # many calls to reduced_norm, reduced_norm_factor,_rdivisor_c, which are slow
    def factorizations(self):
        """
        Return an iterator over all factorizations (as a product
        of a unit and a product of irreducible monic factors) of
        this skew polynomial.

        .. NOTE::

            The algorithm is probabilistic. As a consequence, if
            we execute two times with the same input we can get
            the list of all factorizations in two differents orders.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: a = x^3 + (t^2 + 1)*x^2 + (2*t + 3)*x + t^2 + t + 2
            sage: iter = a.factorizations(); iter
            <generator object at 0x...>
            sage: iter.next()   # random
            (x + 3*t^2 + 4*t) * (x + 2*t^2) * (x + 4*t^2 + 4*t + 2)
            sage: iter.next()   # random
            (x + 3*t^2 + 4*t) * (x + 3*t^2 + 2*t + 2) * (x + 4*t^2 + t + 2)

        We can use this function to build the list of factorizations
        of `a`::

            sage: factorizations = [ F for F in a.factorizations() ]

        We do some checks::

            sage: len(factorizations) == a.count_factorizations()
            True
            sage: len(factorizations) == len(Set(factorizations))  # check no duplicates
            True
            sage: for F in factorizations:
            ...       if F.value() != a:
            ...           print "Found %s which is not a correct factorization" % d
            ...           continue
            ...       for d,_ in F:
            ...           if not d.is_irreducible():
            ...               print "Found %s which is not a correct factorization" % d
        """
        if self.is_zero():
            raise ValueError("factorization of 0 not defined")
        def factorizations_rec(P):
            if P.is_irreducible():
                yield [ (P,1) ]
            else:
                for div in P._irreducible_divisors(True):
                    poly = self // div
                    # Here, we should update poly._norm, poly._norm_factor, poly._rdivisors
                    for factors in P._factorizations_rec():
                        factors.append((div,1))
                        yield factors
        for factors in factorizations_rec(self):
            yield Factorization(factors,sort=False,unit=unit)

