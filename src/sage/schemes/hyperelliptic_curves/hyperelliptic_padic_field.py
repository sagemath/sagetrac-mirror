"""
Hyperelliptic curves over a padic field.
"""

#*****************************************************************************
#  Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import hyperelliptic_generic

from sage.rings.all import (PowerSeriesRing, PolynomialRing, ZZ, QQ, O,
                            pAdicField, GF, RR, RationalField, Infinity)
from sage.misc.functional import log, cyclotomic_polynomial
from sage.modules.free_module import VectorSpace
from sage.matrix.constructor import matrix
from sage.modules.all import vector
from sage.matrix.constructor import identity_matrix
from sage.rings.integer import Integer


class HyperellipticCurve_padic_field(hyperelliptic_generic.HyperellipticCurve_generic):

# The functions below were prototyped at the 2007 Arizona Winter School by
# Robert Bradshaw and Ralf Gerkmann, working with Miljan Brakovevic and
# Kiran Kedlaya
# All of the below is with respect to the Monsky Washnitzer cohomology.

    def local_analytic_interpolation(self, P, Q, prec = None):
        """
        For points $P$, $Q$ in the same residue disc,
        this constructs an interpolation from $P$ to $Q$
        (in homogeneous coordinates) in a power series in
        the local parameter $t$, with precision equal to
        the $p$-adic precision of the underlying ring.

        INPUT:

        - P and Q points on self in the same residue disc

        OUTPUT:

        Returns a point $X(t) = ( x(t) : y(t) : z(t) )$ such that

            (1) $X(0) = P$ and $X(1) = Q$ if $P, Q$ are not in the infinite disc
            (2) $X(P[0]^g}/P[1]) = P$ and $X(Q[0]^g/Q[1]) = Q$ if $P, Q$ are in the infinite disc

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)

        A non-Weierstrass disc::

            sage: P = HK(0,3)
            sage: Q = HK(5, 3 + 3*5^2 + 2*5^3 + 3*5^4 + 2*5^5 + 2*5^6 + 3*5^7 + O(5^8))
            sage: x, y, z = HK.local_analytic_interpolation(P, Q)
            sage: x(0) == P[0], x(1) == Q[0], y(0) == P[1], y.polynomial()(1) == Q[1]
            (True, True, True, True)

        A finite Weierstrass disc::

            sage: P = HK.lift_x(1 + 2*5^2)
            sage: Q = HK.lift_x(1 + 3*5^2)
            sage: x, y, z = HK.local_analytic_interpolation(P, Q)
            sage: x(0) == P[0], x.polynomial()(1) == Q[0], y(0) == P[1], y(1) == Q[1]
            (True, True, True, True)

        The infinite disc::

            sage: P = HK.lift_x(5^-2)
            sage: Q = HK.lift_x(4*5^-2)
            sage: x, y,z = HK.local_analytic_interpolation(P, Q)
            sage: x = x/z
            sage: y = y/z
            sage: x(P[0]/P[1]) == P[0]
            True
            sage: x(Q[0]/Q[1]) == Q[0]
            True
            sage: y(P[0]/P[1]) == P[1]
            True
            sage: y(Q[0]/Q[1]) == Q[1]
            True

        An error if points are not in the same disc::

            sage: x, y,z = HK.local_analytic_interpolation(P,HK(1, 0))
            Traceback (most recent call last):
            ...
            ValueError: (5^-2 + O(5^6) : 5^-3 + 4*5^2 + 5^3 + 3*5^4 + O(5^5) : 1 + O(5^8)) and (1 + O(5^8) : 0 : 1 + O(5^8)) are not in the same residue disc

        AUTHORS:

        - Robert Bradshaw (2007-03)
        - Jennifer Balakrishnan (2010-02)
        """
        if prec is None:
            prec = self.base_ring().precision_cap()
        if not self.is_same_disc(P, Q):
            raise ValueError("%s and %s are not in the same residue disc" % (P, Q))
        disc = self.residue_disc(P)
        t = PowerSeriesRing(self.base_ring(), 't', prec).gen(0)
        if disc == self.change_ring(self.base_ring().residue_field())(0, 1, 0):
            x, y = self.local_coordinates_at_infinity(2*prec)
            g = self.genus()
            return (x*t**(2*g+1), y*t**(2*g+1), t**(2*g+1))
        if disc[1] != 0:
            x = P[0]+t*(Q[0]-P[0])
            pts = self.lift_x(x, all=True)
            if pts[0][1][0] == P[1]:
                return pts[0]
            else:
                return pts[1]
        else:
            S = self.find_char_zero_weier_point(P)
            x, y = self.local_coord(S)
            a = P[1]
            b = Q[1] - P[1]
            y = a + b*t
            x = x.polynomial()(y).add_bigoh(x.prec())
            return (x, y, 1)

    def weierstrass_points(self):
        """
        Return the Weierstrass points of self defined over self.base_ring()

        That is, the point at infinity and those points in the support
        of the divisor of $y$

        EXAMPLES::

            sage: K = pAdicField(11, 5)
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
            sage: C.weierstrass_points()
            [(0 : 1 + O(11^5) : 0), (7 + 10*11 + 4*11^3 + O(11^5) : 0 : 1 + O(11^5))]
        """
        f, h = self.hyperelliptic_polynomials()
        if h != 0:
            raise NotImplementedError()
        return [self((0, 1, 0))] + [self((x, 0, 1)) for x in f.roots()]

    def is_in_weierstrass_disc(self, P):
        """
        Checks if $P$ is in a Weierstrass disc

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: P = HK(0,3)
            sage: HK.is_in_weierstrass_disc(P)
            False
            sage: Q = HK(0, 1, 0)
            sage: HK.is_in_weierstrass_disc(Q)
            True
            sage: S = HK(1, 0)
            sage: HK.is_in_weierstrass_disc(S)
            True
            sage: T = HK.lift_x(1+3*5^2); T
            (1 + 3*5^2 + O(5^8) : 2*5 + 4*5^3 + 3*5^4 + 5^5 + 3*5^6 + O(5^7) : 1 + O(5^8))
            sage: HK.is_in_weierstrass_disc(T)
            True

        AUTHOR:

        - Jennifer Balakrishnan (2010-02)
        """
        if (P[1].valuation() == 0 and P != self(0, 1, 0)):
            return False
        else:
            return True

    def is_weierstrass(self, P):
        """
        Checks if $P$ is a Weierstrass point (i.e., fixed by the hyperelliptic involution)

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: P = HK(0,3)
            sage: HK.is_weierstrass(P)
            False
            sage: Q = HK(0, 1, 0)
            sage: HK.is_weierstrass(Q)
            True
            sage: S = HK(1, 0)
            sage: HK.is_weierstrass(S)
            True
            sage: T = HK.lift_x(1+3*5^2); T
            (1 + 3*5^2 + O(5^8) : 2*5 + 4*5^3 + 3*5^4 + 5^5 + 3*5^6 + O(5^7) : 1 + O(5^8))
            sage: HK.is_weierstrass(T)
            False

        AUTHOR:

        - Jennifer Balakrishnan (2010-02)
        """
        if (P[1] == 0 or P[2] == 0):
            return True
        else:
            return False

    def find_char_zero_weier_point(self, Q):
        """
        Given $Q$ a point on self in a Weierstrass disc, finds the
        center of the Weierstrass disc (if defined over self.base_ring())

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: P = HK.lift_x(1 + 2*5^2)
            sage: Q = HK.lift_x(5^-2)
            sage: S = HK(1, 0)
            sage: T = HK(0, 1, 0)
            sage: HK.find_char_zero_weier_point(P)
            (1 + O(5^8) : 0 : 1 + O(5^8))
            sage: HK.find_char_zero_weier_point(Q)
            (0 : 1 + O(5^8) : 0)
            sage: HK.find_char_zero_weier_point(S)
            (1 + O(5^8) : 0 : 1 + O(5^8))
            sage: HK.find_char_zero_weier_point(T)
            (0 : 1 + O(5^8) : 0)

        AUTHOR:

        - Jennifer Balakrishnan
        """
        if self.is_in_weierstrass_disc(Q) is False:
            raise ValueError("%s is not in a Weierstrass disc" % Q)
        points = self.weierstrass_points()
        for P in points:
            if self.is_same_disc(P, Q):
                return P

    def residue_disc(self, P):
        """
        Gives the residue disc of $P$

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: P = HK.lift_x(1 + 2*5^2)
            sage: HK.residue_disc(P)
            (1 : 0 : 1)
            sage: Q = HK(0,3)
            sage: HK.residue_disc(Q)
            (0 : 3 : 1)
            sage: S = HK.lift_x(5^-2)
            sage: HK.residue_disc(S)
            (0 : 1 : 0)
            sage: T = HK(0, 1, 0)
            sage: HK.residue_disc(T)
            (0 : 1 : 0)

        AUTHOR:

        - Jennifer Balakrishnan
        """
        xPv = P[0].valuation()
        yPv = P[1].valuation()
        F = self.base_ring().residue_field()
        HF = self.change_ring(F)
        if P == self(0, 1, 0):
            return HF(0, 1, 0)
        elif yPv > 0:
            if xPv > 0:
                return HF(0, 0, 1)
            if xPv == 0:
                return HF(P[0].list()[0], 0, 1)
        elif yPv ==0:
            if xPv > 0:
                return HF(0, P[1].list()[0], 1)
            if xPv == 0:
                return HF(P[0].list()[0], P[1].list()[0], 1)
        else:
            return HF(0, 1, 0)

    def is_same_disc(self, P, Q):
        """
        Checks if $P, Q$ are in same residue disc

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: P = HK.lift_x(1 + 2*5^2)
            sage: Q = HK.lift_x(5^-2)
            sage: S = HK(1, 0)
            sage: HK.is_same_disc(P, Q)
            False
            sage: HK.is_same_disc(P, S)
            True
            sage: HK.is_same_disc(Q, S)
            False
        """
        if self.residue_disc(P) == self.residue_disc(Q):
            return True
        else:
            return False

    def tiny_integrals(self, F, P, Q):
        r"""
        Evaluate the integrals of $f_i dx/2y$ from $P$ to $Q$ for each $f_i$ in $F$
        by formally integrating a power series in a local parameter $t$

        $P$ and $Q$ MUST be in the same residue disc for this result to make sense.

        INPUT:

        - F a list of functions $f_i$
        - P a point on self
        - Q a point on self (in the same residue disc as P)

        OUTPUT:

        The integrals $\int_P^Q f_i dx/2y$

        EXAMPLES::

            sage: K = pAdicField(17, 5)
            sage: E = EllipticCurve(K, [-31/3, -2501/108]) # 11a
            sage: P = E(K(14/3), K(11/2))
            sage: TP = E.teichmuller(P);
            sage: x, y = E.monsky_washnitzer_gens()
            sage: E.tiny_integrals([1,x], P, TP) == E.tiny_integrals_on_basis(P,TP)
            True

        ::

            sage: K = pAdicField(11, 5)
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
            sage: P = C.lift_x(11^(-2))
            sage: Q = C.lift_x(3*11^(-2))
            sage: C.tiny_integrals([1], P, Q)
            (3*11^3 + 7*11^4 + 4*11^5 + 7*11^6 + 5*11^7 + O(11^8))

        Note that this fails if the points are not in the same residue disc::

            sage: S = C(0, 1/4)
            sage: C.tiny_integrals([1,x,x^2,x^3], P, S)
            Traceback (most recent call last):
            ...
            ValueError: (11^-2 + O(11^3) : 11^-5 + 8*11^-2 + O(11^0) : 1 + O(11^5)) and (0 : 3 + 8*11 + 2*11^2 + 8*11^3 + 2*11^4 + O(11^5) : 1 + O(11^5)) are not in the same residue disc

        """
        x, y, z = self.local_analytic_interpolation(P, Q)  
        # homogeneous coordinates
        x = x/z
        y = y/z
        dt = x.derivative() / (2*y)
        integrals = []
        g = self.genus()
        for f in F:
            try:
                f_dt = f(x, y)*dt
            except TypeError:   # if f is a constant, not callable
                f_dt = f*dt
            if x.valuation() != -2:
                I = sum([f_dt[n]/(n+1) for n in xrange(f_dt.degree()+1)]) # \int_0^1 f dt
            else:
                If_dt = f_dt.integral().laurent_polynomial()
                I = If_dt(Q[0]**g/Q[1]) - If_dt(P[0]**g/P[1])
            integrals.append(I)
        return vector(integrals)

    def tiny_integrals_on_basis(self, P, Q):
        r"""
        Evaluate the integrals $\{\int_P^Q x^i dx/2y \}_{i=0}^{2g-1}$
        by formally integrating a power series in a local parameter $t$.
        $P$ and $Q$ MUST be in the same residue disc for this result to make sense.

        INPUT:

        - P a point on self
        - Q a point on self (in the same residue disc as P)

        OUTPUT:

        The integrals $\{\int_P^Q x^i dx/2y \}_{i=0}^{2g-1}$

        EXAMPLES::

            sage: K = pAdicField(17, 5)
            sage: E = EllipticCurve(K, [-31/3, -2501/108]) # 11a
            sage: P = E(K(14/3), K(11/2))
            sage: TP = E.teichmuller(P);
            sage: E.tiny_integrals_on_basis(P, TP)
            (17 + 14*17^2 + 17^3 + 8*17^4 + O(17^5), 16*17 + 5*17^2 + 8*17^3 + 14*17^4 + O(17^5))

        ::

            sage: K = pAdicField(11, 5)
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
            sage: P = C.lift_x(11^(-2))
            sage: Q = C.lift_x(3*11^(-2))
            sage: C.tiny_integrals_on_basis(P, Q)
            (3*11^3 + 7*11^4 + 4*11^5 + 7*11^6 + 5*11^7 + O(11^8), 3*11 + 10*11^2 + 8*11^3 + 9*11^4 + 7*11^5 + O(11^6), 4*11^-1 + 2 + 6*11 + 6*11^2 + 7*11^3 + O(11^4), 11^-3 + 6*11^-2 + 2*11^-1 + 2 + O(11^2))


        Note that this fails if the points are not in the same residue disc::

            sage: S = C(0, 1/4)
            sage: C.tiny_integrals_on_basis(P,S)
            Traceback (most recent call last):
            ...
            ValueError: (11^-2 + O(11^3) : 11^-5 + 8*11^-2 + O(11^0) : 1 + O(11^5)) and (0 : 3 + 8*11 + 2*11^2 + 8*11^3 + 2*11^4 + O(11^5) : 1 + O(11^5)) are not in the same residue disc

        """
        if P == Q:
            V = VectorSpace(self.base_ring(), 2*self.genus())
            return V(0)
        R = PolynomialRing(self.base_ring(), ['x', 'y'])
        x, y = R.gens()
        return self.tiny_integrals([x**i for i in range(2*self.genus())], P, Q)

    def teichmuller(self, P):
        r"""
        Find a Teichm\:uller point in the same residue class of $P$.

        Because this lift of frobenius acts as $x \mapsto x^p$,
        take the Teichmuller lift of $x$ and then find a matching $y$
        from that.

        EXAMPLES::

            sage: K = pAdicField(7, 5)
            sage: E = EllipticCurve(K, [-31/3, -2501/108]) # 11a
            sage: P = E(K(14/3), K(11/2))
            sage: E.frobenius(P) == P
            False
            sage: TP = E.teichmuller(P); TP
            (0 : 2 + 3*7 + 3*7^2 + 3*7^4 + O(7^5) : 1 + O(7^5))
            sage: E.frobenius(TP) == TP
            True
            sage: (TP[0] - P[0]).valuation() > 0, (TP[1] - P[1]).valuation() > 0
            (True, True)
        """
        K = P[0].parent()
        x = K.teichmuller(P[0])
        pts = self.lift_x(x, all=True)
        p = K.prime()
        if (pts[0][1] - P[1]).valuation() > 0:
            return pts[0]
        else:
            return pts[1]

    def coleman_integrals_on_basis(self, P, Q, algorithm=None):
        r"""
        Computes the Coleman integrals $\{\int_P^Q x^i dx/2y \}_{i=0}^{2g-1}$

        INPUT:

        - P point on self
        - Q point on self
        - algorithm (optional) = None (uses Frobenius) or teichmuller (uses Teichmuller points)

        OUTPUT:

        the Coleman integrals $\{\int_P^Q x^i dx/2y \}_{i=0}^{2g-1}$

        EXAMPLES::

            sage: K = pAdicField(11, 5)
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
            sage: P = C.lift_x(2)
            sage: Q = C.lift_x(3)
            sage: C.coleman_integrals_on_basis(P, Q)
            (10*11 + 6*11^3 + 2*11^4 + O(11^5), 11 + 9*11^2 + 7*11^3 + 9*11^4 + O(11^5), 3 + 10*11 + 5*11^2 + 9*11^3 + 4*11^4 + O(11^5), 3 + 11 + 5*11^2 + 4*11^4 + O(11^5))
            sage: C.coleman_integrals_on_basis(P, Q, algorithm='teichmuller')
            (10*11 + 6*11^3 + 2*11^4 + O(11^5), 11 + 9*11^2 + 7*11^3 + 9*11^4 + O(11^5), 3 + 10*11 + 5*11^2 + 9*11^3 + 4*11^4 + O(11^5), 3 + 11 + 5*11^2 + 4*11^4 + O(11^5))

        ::

            sage: K = pAdicField(11,5)
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
            sage: P = C.lift_x(11^(-2))
            sage: Q = C.lift_x(3*11^(-2))
            sage: C.coleman_integrals_on_basis(P, Q)
            (3*11^3 + 7*11^4 + 4*11^5 + 7*11^6 + 5*11^7 + O(11^8), 3*11 + 10*11^2 + 8*11^3 + 9*11^4 + 7*11^5 + O(11^6), 4*11^-1 + 2 + 6*11 + 6*11^2 + 7*11^3 + O(11^4), 11^-3 + 6*11^-2 + 2*11^-1 + 2 + O(11^2))

        ::

            sage: R = C(0, 1/4)
            sage: a = C.coleman_integrals_on_basis(P, R)  # long time (7s on sage.math, 2011)
            sage: b = C.coleman_integrals_on_basis(R, Q)  # long time (9s on sage.math, 2011)
            sage: c = C.coleman_integrals_on_basis(P, Q)  # long time
            sage: a+b == c  # long time
            True

        ::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: S = HK(1, 0)
            sage: P = HK(0,3)
            sage: T = HK(0, 1, 0)
            sage: Q = HK.lift_x(5^-2)
            sage: R = HK.lift_x(4*5^-2)
            sage: HK.coleman_integrals_on_basis(S, P)
            (2*5^2 + 5^4 + 5^5 + 3*5^6 + 3*5^7 + 2*5^8 + O(5^9), 5 + 2*5^2 + 4*5^3 + 2*5^4 + 3*5^6 + 4*5^7 + 2*5^8 + O(5^9))
            sage: HK.coleman_integrals_on_basis(T, P)
            (2*5^2 + 5^4 + 5^5 + 3*5^6 + 3*5^7 + 2*5^8 + O(5^9), 5 + 2*5^2 + 4*5^3 + 2*5^4 + 3*5^6 + 4*5^7 + 2*5^8 + O(5^9))
            sage: HK.coleman_integrals_on_basis(P, S) == -HK.coleman_integrals_on_basis(S, P)
            True
            sage: HK.coleman_integrals_on_basis(S, Q)
            (4*5 + 4*5^2 + 4*5^3 + O(5^4), 5^-1 + O(5^3))
            sage: HK.coleman_integrals_on_basis(Q, R)
            (4*5 + 2*5^2 + 2*5^3 + 2*5^4 + 5^5 + 5^6 + 5^7 + 3*5^8 + O(5^9), 2*5^-1 + 4 + 4*5 + 4*5^2 + 4*5^3 + 2*5^4 + 3*5^5 + 2*5^6 + O(5^7))
            sage: HK.coleman_integrals_on_basis(S, R) == HK.coleman_integrals_on_basis(S,Q) + HK.coleman_integrals_on_basis(Q, R) 
            True
            sage: HK.coleman_integrals_on_basis(T, T)
            (0, 0)
            sage: HK.coleman_integrals_on_basis(S, T)
            (0, 0)

        AUTHORS:

        - Robert Bradshaw (2007-03): non-Weierstrass points
        - Jennifer Balakrishnan and Robert Bradshaw (2010-02): Weierstrass points
        """
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        from sage.misc.profiler import Profiler
        prof = Profiler()
        prof("setup")
        K = self.base_ring()
        p = K.prime()
        prec = K.precision_cap()
        g = self.genus()
        dim = 2*g
        V = VectorSpace(K, dim)
        #if P or Q is Weierstrass, use the Frobenius algorithm
        if self.is_weierstrass(P):
            if self.is_weierstrass(Q):
                return V(0)
            else:
                PP = None
                QQ = Q
                TP = None
                TQ = self.frobenius(Q)
        elif self.is_weierstrass(Q):
            PP = P
            QQ = None
            TQ = None
            TP = self.frobenius(P)
        elif self.is_same_disc(P, Q):
            return self.tiny_integrals_on_basis(P, Q)
        elif algorithm == 'teichmuller':
            prof("teichmuller")
            PP = TP = self.teichmuller(P)
            QQ = TQ = self.teichmuller(Q)
            evalP, evalQ = TP, TQ
        else:
            prof("frobPQ")
            TP = self.frobenius(P)
            TQ = self.frobenius(Q)
            PP, QQ = P, Q
        prof("tiny integrals")
        if TP is None:
            P_to_TP = V(0)
        else:
            if not(TP is None):
                TPv = (TP[0]**g/TP[1]).valuation()
                xTPv = TP[0].valuation()
            else:
                xTPv = TPv = +Infinity
            if not(TQ is None):
                TQv = (TQ[0]**g/TQ[1]).valuation()
                xTQv = TQ[0].valuation()
            else:
                xTQv = TQv = +Infinity
            offset = (2*g-1)*max(TPv, TQv)
            if offset == +Infinity:
                offset = (2*g-1)*min(TPv,TQv)
            if (offset > prec and (xTPv <0 or xTQv <0) and (self.residue_disc(P) == self.change_ring(GF(p))(0, 1, 0) or self.residue_disc(Q) == self.change_ring(GF(p))(0, 1, 0))):
                newprec = offset + prec
                K = pAdicField(p, newprec)
                A = PolynomialRing(RationalField(),'x')
                f = A(self.hyperelliptic_polynomials()[0])
                from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
                self = HyperellipticCurve(f).change_ring(K)
                xP = P[0]
                xPv = xP.valuation()
                xPnew = K(sum(xP.list()[i]*p**(xPv + i) for i in range(len(xP.list()))))
                PP = P = self.lift_x(xPnew)
                TP = self.frobenius(P)
                xQ = Q[0]
                xQv = xQ.valuation()
                xQnew = K(sum(xQ.list()[i]*p**(xQv + i) for i in range(len(xQ.list()))))
                QQ = Q = self.lift_x(xQnew)
                TQ = self.frobenius(Q)
                V = VectorSpace(K, dim)
            P_to_TP = V(self.tiny_integrals_on_basis(P, TP))
        if TQ is None:
            TQ_to_Q = V(0)
        else:
            TQ_to_Q = V(self.tiny_integrals_on_basis(TQ, Q))
        prof("mw calc")
        try:
            M_frob, forms = self._frob_calc
        except AttributeError:
            M_frob, forms = self._frob_calc = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(self)
        prof("eval f")
        R = forms[0].base_ring()
        try:
            prof("eval f %s"%R)
            if PP is None:
                L = [-f(R(QQ[0]), R(QQ[1])) for f in forms]  ##changed
            elif QQ is None:
                L = [f(R(PP[0]), R(PP[1])) for f in forms]
            else:
                L = [f(R(PP[0]), R(PP[1])) - f(R(QQ[0]), R(QQ[1])) for f in forms]
        except ValueError:
            prof("changing rings")
            forms = [f.change_ring(self.base_ring()) for f in forms]
            prof("eval f %s"%self.base_ring())
            if PP is None:
                L = [-f(QQ[0], QQ[1]) for f in forms]  ##changed
            elif QQ is None:
                L = [f(PP[0], PP[1]) for f in forms]
            else:
                L = [f(PP[0], PP[1]) - f(QQ[0], QQ[1]) for f in forms]
        b = V(L)
        if PP is None:
            b -= TQ_to_Q
        elif QQ is None:
            b -= P_to_TP
        elif algorithm != 'teichmuller':
            b -= P_to_TP + TQ_to_Q
        prof("lin alg")
        M_sys = matrix(K, M_frob).transpose() - 1
        TP_to_TQ = M_sys**(-1) * b
        prof("done")
#        print prof
        if algorithm == 'teichmuller':
            return P_to_TP + TP_to_TQ + TQ_to_Q
        else:
            return TP_to_TQ

    coleman_integrals_on_basis_hyperelliptic = coleman_integrals_on_basis


#    def invariant_differential(self):
#        """
#        Returns the invariant differential $dx/2y$ on self
#
#        EXAMPLES::
#
#            sage: R.<x> = QQ['x']
#            sage: H = HyperellipticCurve(x^3+1)
#            sage: K = Qp(5,8)
#            sage: HK = H.change_ring(K)
#            sage: w = HK.invariant_differential(); w
#            (((1+O(5^8)))*1) dx/2y
#
#        ::
#
#            sage: K = pAdicField(11, 6)
#            sage: x = polygen(K)
#            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
#            sage: C.invariant_differential()
#            (((1+O(11^6)))*1) dx/2y
#
#        """
#        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
#        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self)
#        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(S)
#        return MW.invariant_differential()

    def coleman_integral(self, w, P, Q, algorithm = 'None'):
        r"""
        Returns the Coleman integral `\int_P^Q w`

        INPUT:

        - w differential (if one of P, Q is Weierstrass, w must be odd)
        - P point on self
        - Q point on self
        - algorithm (optional) = None (uses Frobenius) or teichmuller (uses Teichmuller points)

        OUTPUT:

        the Coleman integral $\int_P^Q w$

        EXAMPLES:

        Example of Leprevost from Kiran Kedlaya
        The first two should be zero as $(P-Q) = 30(P-Q)$ in the Jacobian
        and $dx/2y$ and $x dx/2y$ are holomorphic. ::

            sage: K = pAdicField(11, 6)
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
            sage: P = C(-1, 1); P1 = C(-1, -1)
            sage: Q = C(0, 1/4); Q1 = C(0, -1/4)
            sage: x, y = C.monsky_washnitzer_gens()
            sage: w = C.invariant_differential()
            sage: w.coleman_integral(P, Q)
            O(11^6)
            sage: C.coleman_integral(x*w, P, Q)
            O(11^6)
            sage: C.coleman_integral(x^2*w, P, Q)
            7*11 + 6*11^2 + 3*11^3 + 11^4 + 5*11^5 + O(11^6)

        ::

            sage: p = 71; m = 4
            sage: K = pAdicField(p, m)
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 + 33/16*x^4 + 3/4*x^3 + 3/8*x^2 - 1/4*x + 1/16)
            sage: P = C(-1, 1); P1 = C(-1, -1)
            sage: Q = C(0, 1/4); Q1 = C(0, -1/4)
            sage: x, y = C.monsky_washnitzer_gens()
            sage: w = C.invariant_differential()
            sage: w.integrate(P, Q), (x*w).integrate(P, Q)
            (O(71^4), O(71^4))
            sage: R, R1 = C.lift_x(4, all=True)
            sage: w.integrate(P, R)
            21*71 + 67*71^2 + 27*71^3 + O(71^4)
            sage: w.integrate(P, R) + w.integrate(P1, R1)
            O(71^4)

        A simple example, integrating dx::

            sage: R.<x> = QQ['x']
            sage: E= HyperellipticCurve(x^3-4*x+4)
            sage: K = Qp(5, 10)
            sage: EK = E.change_ring(K)
            sage: P = EK(2, 2)
            sage: Q = EK.teichmuller(P)
            sage: x, y = EK.monsky_washnitzer_gens()
            sage: EK.coleman_integral(x.diff(), P, Q)
            5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: Q[0] - P[0]
            5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)

        Yet another example::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x*(x-1)*(x+9))
            sage: K = Qp(7, 10)
            sage: HK = H.change_ring(K)
            sage: import sage.schemes.elliptic_curves.monsky_washnitzer as mw
            sage: M_frob, forms = mw.matrix_of_frobenius_hyperelliptic(HK)
            sage: w = HK.invariant_differential()
            sage: x, y = HK.monsky_washnitzer_gens()
            sage: f = forms[0]
            sage: S = HK(9,36)
            sage: Q = HK.teichmuller(S)
            sage: P = HK(-1,4)
            sage: b = x*w*w._coeff.parent()(f)
            sage: HK.coleman_integral(b, P, Q)
            7 + 7^2 + 4*7^3 + 5*7^4 + 3*7^5 + 7^6 + 5*7^7 + 3*7^8 + 4*7^9 + 4*7^10 + O(7^11)

        ::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3+1)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: w = HK.invariant_differential()
            sage: P = HK(0, 1)
            sage: Q = HK.lift_x(5)
            sage: x, y = HK.monsky_washnitzer_gens()
            sage: (2*y*w).coleman_integral(P, Q)
            5 + O(5^9)
            sage: xloc,yloc,zloc = HK.local_analytic_interpolation(P, Q)
            sage: I2 = (xloc.derivative()/(2*yloc)).integral()
            sage: I2.polynomial()(1) - I2(0)
            3*5 + 2*5^2 + 2*5^3 + 5^4 + 4*5^6 + 5^7 + O(5^9)
            sage: HK.coleman_integral(w, P, Q)
            3*5 + 2*5^2 + 2*5^3 + 5^4 + 4*5^6 + 5^7 + O(5^9)

        Integrals involving Weierstrass points::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: S = HK(1, 0)
            sage: P = HK(0,3)
            sage: negP = HK(0, -3)
            sage: T = HK(0, 1, 0)
            sage: w = HK.invariant_differential()
            sage: x, y = HK.monsky_washnitzer_gens()
            sage: HK.coleman_integral(w*x^3,S,T)
            0
            sage: HK.coleman_integral(w*x^3,T,S)
            0
            sage: HK.coleman_integral(w, S, P)
            2*5^2 + 5^4 + 5^5 + 3*5^6 + 3*5^7 + 2*5^8 + O(5^9)
            sage: HK.coleman_integral(w, T, P)
            2*5^2 + 5^4 + 5^5 + 3*5^6 + 3*5^7 + 2*5^8 + O(5^9)
            sage: HK.coleman_integral(w*x^3, T, P)
            5^2 + 2*5^3 + 3*5^6 + 3*5^7 + O(5^8)
            sage: HK.coleman_integral(w*x^3, S, P)
            5^2 + 2*5^3 + 3*5^6 + 3*5^7 + O(5^8)
            sage: HK.coleman_integral(w, P, negP, algorithm='teichmuller')
            5^2 + 4*5^3 + 2*5^4 + 2*5^5 + 3*5^6 + 2*5^7 + 4*5^8 + O(5^9)
            sage: HK.coleman_integral(w, P, negP)
            5^2 + 4*5^3 + 2*5^4 + 2*5^5 + 3*5^6 + 2*5^7 + 4*5^8 + O(5^9)

        AUTHORS:

        - Robert Bradshaw (2007-03)
        - Kiran Kedlaya (2008-05)
        - Jennifer Balakrishnan (2010-02)
        """
        # TODO: implement Jacobians and show the relationship directly
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        K = self.base_ring()
        prec = K.precision_cap()
        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self, K)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(S)
        w = MW(w)
        f, vec = w.reduce_fast()
        basis_values = self.coleman_integrals_on_basis(P, Q, algorithm)
        dim = len(basis_values)
        x, y = self.local_coordinates_at_infinity(2*prec)
        if self.is_weierstrass(P):
            if self.is_weierstrass(Q):
                return 0
            elif f == 0:
                return sum([vec[i] * basis_values[i] for i in range(dim)])
            elif w._coeff(x, -y)*x.derivative()/(-2*y)+w._coeff(x, y)*x.derivative()/(2*y) == 0:
                return self.coleman_integral(w,self(Q[0], -Q[1]), self(Q[0], Q[1]), algorithm)/2
            else:
                raise ValueError("The differential is not odd: use coleman_integral_from_weierstrass_via_boundary")

        elif self.is_weierstrass(Q):
            if f == 0:
                return sum([vec[i] * basis_values[i] for i in range(dim)])
            elif w._coeff(x, -y)*x.derivative()/(-2*y)+w._coeff(x, y)*x.derivative()/(2*y) == 0:
                return -self.coleman_integral(w,self(P[0], -P[1]), self(P[0], P[1]), algorithm)/2
            else:
                raise ValueError("The differential is not odd: use coleman_integral_from_weierstrass_via_boundary")
        else:
            return f(Q[0], Q[1]) - f(P[0], P[1]) + sum([vec[i] * basis_values[i] for i in range(dim)]) # this is just a dot product...

    def frobenius(self, P=None):
        """
        Returns the $p$-th power lift of Frobenius of $P$

        EXAMPLES::

            sage: K = Qp(11, 5)
            sage: R.<x> = K[]
            sage: E = HyperellipticCurve(x^5 - 21*x - 20)
            sage: P = E.lift_x(2)
            sage: E.frobenius(P)
            (2 + 10*11 + 5*11^2 + 11^3 + O(11^5) : 5 + 9*11 + 2*11^2 + 2*11^3 + O(11^5) : 1 + O(11^5))
            sage: Q = E.teichmuller(P); Q
            (2 + 10*11 + 4*11^2 + 9*11^3 + 11^4 + O(11^5) : 5 + 9*11 + 6*11^2 + 11^3 + 6*11^4 + O(11^5) : 1 + O(11^5))
            sage: E.frobenius(Q)
            (2 + 10*11 + 4*11^2 + 9*11^3 + 11^4 + O(11^5) : 5 + 9*11 + 6*11^2 + 11^3 + 6*11^4 + O(11^5) : 1 + O(11^5))

        ::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: Q = H(0, 0)
            sage: u,v = H.local_coord(Q, prec=100)
            sage: K = Qp(11,5)
            sage: L.<a> = K.extension(x^20-11)
            sage: HL = H.change_ring(L)
            sage: S = HL(u(a),v(a))
            sage: HL.frobenius(S)
            (8*a^22 + 10*a^42 + 4*a^44 + 2*a^46 + 9*a^48 + 8*a^50 + a^52 + 7*a^54 +
            7*a^56 + 5*a^58 + 9*a^62 + 5*a^64 + a^66 + 6*a^68 + a^70 + 6*a^74 +
            2*a^76 + 2*a^78 + 4*a^82 + 5*a^84 + 2*a^86 + 7*a^88 + a^90 + 6*a^92 +
            a^96 + 5*a^98 + 2*a^102 + 2*a^106 + 6*a^108 + 8*a^110 + 3*a^112 +
            a^114 + 8*a^116 + 10*a^118 + 3*a^120 + O(a^122) :
            a^11 + 7*a^33 + 7*a^35 + 4*a^37 + 6*a^39 + 9*a^41 + 8*a^43 + 8*a^45 +
            a^47 + 7*a^51 + 4*a^53 + 5*a^55 + a^57 + 7*a^59 + 5*a^61 + 9*a^63 +
            4*a^65 + 10*a^69 + 3*a^71 + 2*a^73 + 9*a^75 + 10*a^77 + 6*a^79 +
            10*a^81 + 7*a^85 + a^87 + 4*a^89 + 8*a^91 + a^93 + 8*a^95 + 2*a^97 +
            7*a^99 + a^101 + 3*a^103 + 6*a^105 + 7*a^107 + 4*a^109 + O(a^111) :
            1 + O(a^100))

        AUTHORS:

        - Robert Bradshaw and Jennifer Balakrishnan (2010-02)
        """
        try:
            _frob = self._frob
        except AttributeError:
            K = self.base_ring()
            p = K.prime()
            x = K['x'].gen(0)

            f, f2 = self.hyperelliptic_polynomials()
            if f2 != 0:
                raise NotImplementedError("Curve must be in weierstrass normal form.")
            h = (f(x**p) - f**p)

            def _frob(P):
                if P == self(0, 1, 0):
                    return P
                x0 = P[0]
                y0 = P[1]
                try:
                    uN = (1 + h(x0)/y0**(2*p)).sqrt()
                    yres = y0**p * uN
                    xres = x0**p
                    if (yres-y0).valuation() == 0:
                        yres=-yres
                    return self.point([xres, yres, K(1)])
                except (TypeError, NotImplementedError):
                    uN2 = 1 + h(x0)/y0**(2*p)
                    #yfrob2 = f(x)
                    c = uN2.list()[0]
                    v = uN2.valuation()
                    a = uN2.parent().gen()
                    uN = self.newton_sqrt(uN2, c.sqrt()*a**(v//2),
                                          K.precision_cap())
                    yres = y0**p * uN
                    xres = x0**p
                    if (yres - y0).valuation() == 0:
                        yres = -yres
                    try:
                        return self(xres,yres)
                    except ValueError:
                        return self._curve_over_ram_extn(xres, yres)

            self._frob = _frob

        if P is None:
            return _frob
        else:
            return _frob(P)

    def newton_sqrt(self, f, x0, prec):
        r"""
        Takes the square root of the power series $f$ by Newton's method

        NOTE:

        this function should eventually be moved to $p$-adic power series ring

        INPUT:

        - f power series wtih coefficients in $\Q_p$ or an extension
        - x0 seeds the Newton iteration
        - prec precision

        OUTPUT:

        the square root of $f$

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: Q = H(0, 0)
            sage: u,v = H.local_coord(Q, prec=100)
            sage: K = Qp(11,5)
            sage: HK = H.change_ring(K)
            sage: L.<a> = K.extension(x^20-11)
            sage: HL = H.change_ring(L)
            sage: S = HL(u(a),v(a))
            sage: f = H.hyperelliptic_polynomials()[0]
            sage: y = HK.newton_sqrt( f(u(a)^11), a^11,5)
            sage: y^2 - f(u(a)^11)
            O(a^122)

        AUTHOR:

        - Jennifer Balakrishnan
        """
        z = x0
        try:
            x = f.parent().variable_name()
            if x!='a' :  #this is to distinguish between extensions of Qp that are finite vs. not
                S = f.base_ring()[[x]]
                x = S.gen()
        except ValueError:
            pass
        z = x0
        loop_prec = (log(RR(prec))/log(RR(2))).ceil()
        for i in range(loop_prec):
            z = (z+f/z)/2
        try:
            return z + O(x**prec)
        except (NameError, ArithmeticError, TypeError):
            return z

    def curve_over_ram_extn(self, deg):
        r"""
        Returns self over $\Q_p(p^(1/deg))$

        INPUT:

        - deg: the degree of the ramified extension

        OUTPUT:

        self over $\Q_p(p^(1/deg))$

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: K = Qp(11,5)
            sage: HK = H.change_ring(K)
            sage: HL = HK.curve_over_ram_extn(2)
            sage: HL
            Hyperelliptic Curve over Eisenstein Extension of 11-adic Field with capped relative precision 5 in a defined by (1 + O(11^5))*x^2 + (O(11^6))*x + (10*11 + 10*11^2 + 10*11^3 + 10*11^4 + 10*11^5 + O(11^6)) defined by (1 + O(a^10))*y^2 = (1 + O(a^10))*x^5 + (10 + 8*a^2 + 10*a^4 + 10*a^6 + 10*a^8 + O(a^10))*x^3 + (7 + a^2 + O(a^10))*x^2 + (7 + 3*a^2 + O(a^10))*x

        AUTHOR:

        - Jennifer Balakrishnan
        """
        from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
        K = self.base_ring()
        p = K.prime()
        A = PolynomialRing(QQ,'x')
        x = A.gen()
        J = K.extension(x**deg-p, names='a')
        pol = self.hyperelliptic_polynomials()[0]
        H = HyperellipticCurve(A(pol))
        HJ = H.change_ring(J)
        self._curve_over_ram_extn = HJ
        self._curve_over_ram_extn._curve_over_Qp = self
        return HJ

    def get_boundary_point(self, curve_over_extn, P):
        """
        Given self over an extension field, find a point in the disc of $P$ near the boundary

        INPUT:

        - curve_over_extn: self over a totally ramified extension
        - P: Weierstrass point

        OUTPUT:

        a point in the disc of $P$ near the boundary

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(3,6)
            sage: HK = H.change_ring(K)
            sage: P = HK(1, 0)
            sage: J.<a> = K.extension(x^30-3)
            sage: HJ  = H.change_ring(J)
            sage: S = HK.get_boundary_point(HJ, P)
            sage: S
            (1 + 2*a^2 + 2*a^6 + 2*a^18 + a^32 + a^34 + a^36 + 2*a^38 + 2*a^40 + a^42 + 2*a^44 + a^48 + 2*a^50 + 2*a^52 + a^54 + a^56 + 2*a^60 + 2*a^62 + a^70 + 2*a^72 + a^76 + 2*a^78 + a^82 + a^88 + a^96 + 2*a^98 + 2*a^102 + a^104 + 2*a^106 + a^108 + 2*a^110 + a^112 + 2*a^116 + a^126 + 2*a^130 + 2*a^132 + a^144 + 2*a^148 + 2*a^150 + a^152 + 2*a^154 + a^162 + a^164 + a^166 + a^168 + a^170 + a^176 + a^178 + O(a^180) : a + O(a^180) : 1 + O(a^180))

        AUTHOR:

        - Jennifer Balakrishnan
        """
        J = curve_over_extn.base_ring()
        a = J.gen()
        prec2 = J.precision_cap()
        x, y = self.local_coord(P, prec2)
        return curve_over_extn(x(a), y(a))

    def P_to_S(self, P, S):
        r"""
        Given a finite Weierstrass point $P$ and a point $S$
        in the same disc, computes the Coleman integrals $\{\int_P^S x^i dx/2y \}_{i=0}^{2g-1}$

        INPUT:

        - P: finite Weierstrass point
        - S: point in disc of P

        OUTPUT:

        Coleman integrals $\{\int_P^S x^i dx/2y \}_{i=0}^{2g-1}$

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,4)
            sage: HK = H.change_ring(K)
            sage: P = HK(1, 0)
            sage: HJ = HK.curve_over_ram_extn(10)
            sage: S = HK.get_boundary_point(HJ, P)
            sage: HK.P_to_S(P, S)
            (2*a + 4*a^3 + 2*a^11 + 4*a^13 + 2*a^17 + 2*a^19 + a^21 + 4*a^23 + a^25 + 2*a^27 + 2*a^29 + 3*a^31 + 4*a^33 + O(a^35), a^-5 + 2*a + 2*a^3 + a^7 + 3*a^11 + a^13 + 3*a^15 + 3*a^17 + 2*a^19 + 4*a^21 + 4*a^23 + 4*a^25 + 2*a^27 + a^29 + a^31 + O(a^33))

        AUTHOR:

        - Jennifer Balakrishnan
        """
        prec = self.base_ring().precision_cap()
        deg = (S[0]).parent().defining_polynomial().degree()
        prec2= prec*deg
        x, y = self.local_coord(P, prec2)
        g = self.genus()
        integrals = [((x**k*x.derivative()/(2*y)).integral()) for k in range(2*g)]
        val = [I(S[1]) for I in integrals]
        return vector(val)

    def coleman_integral_P_to_S(self, w, P, S):
        r"""
        Given a finite Weierstrass point $P$ and a point $S$
        in the same disc, computes the Coleman integral $\int_P^S w$

        INPUT:

        - w: differential
        - P: Weierstrass point
        - S: point in the same disc of P (S is defined over an extension of $\Q_p$; coordinates
          of S are given in terms of uniformizer $a$)

        OUTPUT:

        Coleman integral $\int_P^S w$ in terms of $a$

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,4)
            sage: HK = H.change_ring(K)
            sage: P = HK(1, 0)
            sage: J.<a> = K.extension(x^10-5)
            sage: HJ  = H.change_ring(J)
            sage: S = HK.get_boundary_point(HJ, P)
            sage: x, y = HK.monsky_washnitzer_gens()
            sage: S[0]-P[0] == HK.coleman_integral_P_to_S(x.diff(), P, S)
            True
            sage: HK.coleman_integral_P_to_S(HK.invariant_differential(), P, S) == HK.P_to_S(P,S)[0]
            True

        AUTHOR:

        - Jennifer Balakrishnan
        """
        prec = self.base_ring().precision_cap()
        deg = S[0].parent().defining_polynomial().degree()
        prec2= prec*deg
        x, y = self.local_coord(P, prec2)
        g = self.genus()
        int_sing = (w.coeff()(x, y)*x.derivative()/(2*y)).integral()
        int_sing_a = int_sing(S[1])
        return int_sing_a

    def S_to_Q(self, S, Q):
        r"""
        Given $S$ a point on self over an extension field, computes the
        Coleman integrals $\{\int_S^Q x^i dx/2y \}_{i=0}^{2g-1}$

        **one should be able to feed $S, Q$ into coleman_integral,
        but currently that segfaults**

        INPUT:

        - S: a point with coordinates in an extension of $\Q_p$ (with unif. a)
        - Q: a non-Weierstrass point defined over $\Q_p$

        OUTPUT:

        the Coleman integrals $\{\int_S^Q x^i dx/2y \}_{i=0}^{2g-1}$ in terms of $a$

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,6)
            sage: HK = H.change_ring(K)
            sage: J.<a> = K.extension(x^20-5)
            sage: HJ  = H.change_ring(J)
            sage: w = HK.invariant_differential()
            sage: x, y = HK.monsky_washnitzer_gens()
            sage: P = HK(1, 0)
            sage: Q = HK(0,3)
            sage: S = HK.get_boundary_point(HJ, P)
            sage: P_to_S = HK.P_to_S(P, S)
            sage: S_to_Q = HJ.S_to_Q(S, Q)
            sage: P_to_S + S_to_Q
            (2*a^40 + a^80 + a^100 + O(a^105), a^20 + 2*a^40 + 4*a^60 + 2*a^80 + O(a^103))
            sage: HK.coleman_integrals_on_basis(P, Q)
            (2*5^2 + 5^4 + 5^5 + 3*5^6 + O(5^7), 5 + 2*5^2 + 4*5^3 + 2*5^4 + 5^6 + O(5^7))

        AUTHOR:

        - Jennifer Balakrishnan
        """
        FS = self.frobenius(S)
        FS = (FS[0], FS[1])
        FQ = self.frobenius(Q)
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        try:
            M_frob, forms = self._frob_calc
        except AttributeError:
            M_frob, forms = self._frob_calc = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(self)
        try:
            HJ = self._curve_over_ram_extn
            K = HJ.base_ring()
        except AttributeError:
            HJ = S.scheme()
            K = self.base_ring()
        g = self.genus()
        prec2 = K.precision_cap()
        p = K.prime()
        dim = 2*g
        V = VectorSpace(K, dim)
        if S == FS:
            S_to_FS = V(dim*[0])
        else:
            P = self(ZZ(FS[0][0]),ZZ(FS[1][0]))
            x, y = self.local_coord(P, prec2)
            integrals = [(x**i*x.derivative()/(2*y)).integral() for i in range(dim)]
            S_to_FS = vector([I.polynomial()(FS[1]) - I.polynomial()(S[1]) for I in integrals])
        if HJ(Q[0], Q[1]) == HJ(FQ):
            FQ_to_Q = V(dim*[0])
        else:
            FQ_to_Q = V(self.tiny_integrals_on_basis(FQ, Q))
        try:
            L = [f(K(S[0]), K(S[1])) - f(K(Q[0]), K(Q[1])) for f in forms]
        except ValueError:
            forms = [f.change_ring(K) for f in forms]
            L = [f(S[0], S[1]) - f(Q[0], Q[1]) for f in forms]
        b = V(L)
        M_sys = matrix(K, M_frob).transpose() - 1
        B = (~M_sys)
        v = [B.list()[i].valuation() for i in range(len(B.list()))]
        vv= min(v)
        B = (p**(-vv)*B).change_ring(K)
        B = p**(vv)*B
        return B*(b-S_to_FS-FQ_to_Q)

    def coleman_integral_S_to_Q(self, w, S, Q):
        r"""
        Computes the Coleman integral `\int_S^Q w`

        **one should be able to feed $S,Q$ into coleman_integral,
        but currently that segfaults**

        INPUT:

        - w -- a differential
        - S -- a point with coordinates in an extension of `\Q_p`
        - Q -- a non-Weierstrass point defined over `\Q_p`

        OUTPUT:

        the Coleman integral `\int_S^Q w` 

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,6)
            sage: HK = H.change_ring(K)
            sage: J.<a> = K.extension(x^20-5)
            sage: HJ  = H.change_ring(J)
            sage: x, y = HK.monsky_washnitzer_gens()
            sage: P = HK(1, 0)
            sage: Q = HK(0,3)
            sage: S = HK.get_boundary_point(HJ, P)
            sage: P_to_S = HK.coleman_integral_P_to_S(y.diff(), P, S)
            sage: S_to_Q = HJ.coleman_integral_S_to_Q(y.diff(), S, Q)
            sage: P_to_S  + S_to_Q
            3 + O(a^119)
            sage: HK.coleman_integral(y.diff(), P, Q)
            3 + O(5^6)

        AUTHOR:

        - Jennifer Balakrishnan
        """
        import sage.schemes.elliptic_curves.monsky_washnitzer as monsky_washnitzer
        K = self.base_ring()
        R = monsky_washnitzer.SpecialHyperellipticQuotientRing(self, K)
        MW = monsky_washnitzer.MonskyWashnitzerDifferentialRing(R)
        w = MW(w)
        f, vec = w.reduce_fast()
        g = self.genus()
        const = f(Q[0], Q[1]) - f(S[0], S[1])
        if vec == vector(2*g*[0]):
            return const
        else:
            basis_values = self.S_to_Q(S, Q)
            dim = len(basis_values)
            dot = sum([vec[i] * basis_values[i] for i in range(dim)])
            return const + dot

    def coleman_integral_from_weierstrass_via_boundary(self, w, P, Q, d):
        r"""
        Computes the Coleman integral $\int_P^Q w$ via a boundary point
        in the disc of $P$, defined over a degree $d$ extension

        INPUT:

        - w -- a differential
        - P -- a Weierstrass point
        - Q -- a non-Weierstrass point
        - d -- degree of extension where coordinates of boundary point lie

        OUTPUT:

        the Coleman integral $\int_P^Q w$, written in terms of the uniformizer
        $a$ of the degree $d$ extension

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: K = Qp(5,6)
            sage: HK = H.change_ring(K)
            sage: P = HK(1, 0)
            sage: Q = HK(0,3)
            sage: x, y = HK.monsky_washnitzer_gens()
            sage: HK.coleman_integral_from_weierstrass_via_boundary(y.diff(), P, Q, 20)
            3 + O(a^119)
            sage: HK.coleman_integral(y.diff(), P, Q)
            3 + O(5^6)
            sage: w = HK.invariant_differential()
            sage: HK.coleman_integral_from_weierstrass_via_boundary(w, P, Q, 20)
            2*a^40 + a^80 + a^100 + O(a^105)
            sage: HK.coleman_integral(w, P, Q)
            2*5^2 + 5^4 + 5^5 + 3*5^6 + O(5^7)

        AUTHOR:

        - Jennifer Balakrishnan
        """
        HJ = self.curve_over_ram_extn(d)
        S = self.get_boundary_point(HJ, P)
        P_to_S = self.coleman_integral_P_to_S(w, P, S)
        S_to_Q = HJ.coleman_integral_S_to_Q(w, S, Q)
        return P_to_S + S_to_Q

#-------------local heights------------------


    def recip_froby(self, x, y, prec=10):
        """
        Given local expansions x(t) and y(t), computes the reciprocal of the Frobenius of y

        EXAMPLES::
        """
        f = self.hyperelliptic_polynomials()[0]
        p = self.base_ring().prime()
        s = 1 + sum(self.Bm(i) * ((f(x**p)-f(x)**p)/(f(x)**p))**i
                    for i in range(1, prec))
        return y**(-p) * s

    def froby(self, x, y, prec=10):
        """
        Given local expansions x(t) and y(t), computes the Frobenius of y

        EXAMPLES::
        """
        f = self.hyperelliptic_polynomials()[0]
        p = self.base_ring().prime()
        s = 1 + sum(self.Bp(i) * ((f(x**p)-f(x)**p)/(f(x)**p))**i
                    for i in range(1, prec))
        return y**p * s

    def sum_of_local_symbols(self, divisor, prec=20):
        """
        For $w$ a differential with given residue divisor and $w_0,...,
        w_{2g-1}$ the basis of de Rham cohomology, computes 
        $\{\sum_P <w,w_i>_P\}_{i=0,...,2g}$, where the sum is taken over all
        points $P$ in the divisor as well as all weierstrass points.

        NOTE: Right now, assuming that divisor = (P)-(Q)
        
        EXAMPLES:
            sage: K = pAdicField(11, 10)
            sage: R.<x> = PolynomialRing(K)
            sage: C = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: P = C(1,6)
            sage: Q = C(-2, 12)
        
        
        AUTHOR:

        - Jennifer Balakrishnan (2007-12)
        """
        x, y = self.monsky_washnitzer_gens()
        w = self.invariant_differential()
        g = self.genus()
        P = divisor[0][1]
        Q = divisor[1][1]
        local = vector([divisor[1][0]*self.coleman_integral(w*x**j, P, Q)
                        for j in range(2*g)])
        return local

    def differential_log(self, divisor, prec=40):
        r"""
        Given the hyperelliptic curve `C`, computes the log of a
        differential with given residue divisor lies in `H^1_dR(C)`

        This is Psi(w)

        if g = 1, W is spanned by x^(g+1) dx/y, ... x^(2g-1) dx/y else
        W is unit root subspace, given by Frob^n (x^(g+1) dx/y), ...,
        Frob^n (x^(2g-1) dx/y)

        EXAMPLES::

            sage: K = pAdicField(11,5)
            sage: R.<x> = PolynomialRing(K)
            sage: C = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: P = C(1,6)
            sage: Q = C(-2, 12)
            sage: C.differential_log([(1, P), (-1, Q)])
            (11^-1 + 4 + 10*11 + 3*11^2 + O(11^3), 5*11^-1 + 7 + 11 + 9*11^2 + O(11^3), 3 + 10*11 + 8*11^2 + 7*11^3 + O(11^4), 6*11 + 10*11^2 + 4*11^3 + O(11^4))

        AUTHOR:

        - Jennifer Balakrishnan (2007-12)
        """
        A = self.cup_product_matrix(prec)
        # A = self._cpm
        v = self.sum_of_local_symbols(divisor, prec)
        g = self.genus()
        if g == 1:
            self._subspace = identity_matrix(2)
            return A**(-1)*v
        else:
            from sage.schemes.elliptic_curves.monsky_washnitzer import matrix_of_frobenius_hyperelliptic
            try:
                M_frob, forms = self._frob_calc
            except AttributeError:
                M_frob, forms = self._frob_calc = matrix_of_frobenius_hyperelliptic(self)
            default_prec = self.base_ring().precision_cap()
            X = (M_frob**(default_prec)).matrix_from_columns([g, 2*g-1]).transpose().list()
            I = identity_matrix(2*g).matrix_from_columns([0, g-1]).transpose().list()
            V = matrix(2*g, 2*g, I + X).transpose()
            self._subspace = V
            return V**(-1)*(A**(-1)*v)

    def differential_log_projection(self, divisor, prec=20):
        """
        The component of differential_log that lies in W
        Have to specify W in differential_log

        EXAMPLES:
            sage: K = pAdicField(11,5)
            sage: R.<x> = PolynomialRing(K)
            sage: C = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: P = C(1,6)
            sage: Q = C(-2, 12)


        AUTHOR:

        - Jennifer Balakrishnan (2008-01)
        """
        g = self.genus()
        v = self.differential_log(divisor, prec)
        return vector([0]*g + [v[i] for i in range(g, 2*g)])

    def differential_log_holomorphic(self, divisor, prec=20):
        """
        The holomorphic component of differential_log.

        aka ``eta''

        EXAMPLES:
            sage: K = pAdicField(11,5)
            sage: R.<x> = PolynomialRing(K)
            sage: C = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: P = C(1,6)
            sage: Q = C(-2, 12)


        AUTHOR:

        - Jennifer Balakrishnan (2008-01)
        """
        #print "differential_log_holomorphic"
        g = self.genus()
        v = self.differential_log(divisor, prec)
        w = vector([v[i] for i in range(g)]+[0]*g)
        return w

    def square_root_extension(self, num):
        """
        takes a square root in an p-adic extension
        (might not return the right one)
        """
        #print "square_root_extension"
        p = self.base_ring().prime()
        #print num
        if num.valuation() == 0:
            c = num.list()[0]
            i = 1
            while i <p//2:
                if (i**2)%p == c:
                    break
                i=i+1
        ###Newton's method for square root
        j = 0
        while j <100:       #change prec later
            root = (i+num/i)/2
            if i == root:
                break
            i = root
            #print i
            j = j+1
        if i**2 == num:
            #print "found a root!"
            return i

    def find_pth_root_point(self, P, all=False):
        r"""
        Given P=(a, b), finds a point P'=(a', b') over Qp(a^(1/p) such that
        a'^p = a

        if all is ``True``, find all pth roots

        """
        #print "find_pth_root_point"
        K = self.base_ring()
        p = K.prime()
        xP = P[0]

        g = self.hyperelliptic_polynomials()[0]
        ###working over the appropriate field
        R = QQ['x']
        x = R.gen()
        if xP**p == xP:
            f = cyclotomic_polynomial(p, var='y')
            J = K.extension(f(x+1), names='a')
        else:
            #print "not cyclotomic"
            J = K.extension((x+QQ(xP))**p-QQ(xP), names='a')
        a = J.gen()
        HJ = self.change_ring(J)

        ###find the pth roots of x(P)

        if xP**p == xP:
            xPfracpow = (1+a) * xP
        else:
            xPfracpow = a + xP
        if g(xPfracpow) == 0:
            return HJ(xPfracpow, 0)
        yPfracpow = HJ.square_root_extension(g(xPfracpow))
        Pfracpow = HJ(xPfracpow,yPfracpow)
        P = HJ(P[0], P[1])
        if P[0] == Pfracpow[0] :
            xnew = (a+xP)*Pfracpow[0]
            ynew = HJ.square_root_extension(g(xnew))
            Pfracpow = HJ(xnew,ynew)
        if ((Pfracpow[1]).list()[0] == (P[1]).list()[0]):
            point =  Pfracpow
        else:
            point = HJ(Pfracpow[0], -Pfracpow[1])
        if all is False:
            #print point[0]**p == P[0]
            return point
        else:
            if xP**p != xP:
                print "sorry, we can only print all of the roots when the extension is cyclotomic."
                return point
            else:
                pts = [point]
                xs = [(a+1)**i*(pts[0][0]) for i in range(1, p-1)]
                ys = [HJ.square_root_extension(g(x)) for x in xs]
                ynew = []
                for y in ys:
                    if y.list()[0] != (P[1]).list()[0] :
                        ynew = ynew + [-y]
                    else:
                        ynew = ynew + [y]
                pts += [HJ(xs[i],ynew[i]) for i in range(len(xs))]
                return pts

    def local_analytic_interpolation_cyclotomic(self, P, Q, prec=30):
        """
        Given P and x(Q), with P,Q
        in the same residue disc and P defined over Qp,
        this computes the local analytic interpolation
        between P,Q

        USE: for non-weierstrass points
        """
        R = self.base_ring()[['t']]
        t = R.gen()
        x, y = self.local_coord(P, prec)      # figure out precision here
        return x((Q-P[0])*t), y((Q-P[0])*t)

    def sum_of_local_symbols_extension(self, divisor, P, Q, extension=False,
                                       prec=80):
        """
        Computes the vector of local symbols (<w,w_i>_A)_{i=0...2g-1}, where w is a differential form with residue divisor "divisor"

        if extension is True, creates the appropriate field extension and curve over that extension

        P is fixed point for all computations to link constants of integration
        Q is the residue disc where all the stuff should be happening

        .. TODO: merge with sum_of_local_symbols, make more modular
        """
        p = self.base_ring().prime()
        #if (Q[0]**p==Q[0]):
        #    cyc = True
        #else:
        #    cyc = False

        wstrass = bool(Q[1] == 0)

        if extension:
            A = self.find_pth_root_point(Q)
        else:
            A = Q

        g = self.hyperelliptic_polynomials()[0]
        gen = self.genus()

        ###working over the appropriate field

        cyc = False
        if extension:
            R = QQ['x']
            x = R.gen()
            if cyc:
                f = cyclotomic_polynomial(p)
                J = self.base_ring().extension(f(x+1), names='a')
            else:
                if Q[0]**p == Q[0]:
                    f = (x+Q[0])**p - Q[0]
                    d = f.degree()
                    ff = sum(f.list()[i]*x**(i-1) for i in range(1, d+1))
                    J = self.base_ring().extension(ff, names='a')
                else:
                    J = self.base_ring().extension((x+Q[0])**p-Q[0], names='a')
            a = J.gen()
            H = self.change_ring(J)
        else:
            H = self

        x, y = H.local_coord(A, prec=30)     #worry about prec later
        ###formal antiderivative of w_i


        I2 = vector(J,[0]*2*gen)
        ###if working over an extension, need tiny integral + coleman integral over Qp
        if extension is True:
            xx, yy = H.local_analytic_interpolation_cyclotomic(Q, A[0], prec)  # this will change when it's weiestrass
            print xx(0) == Q[0]
            print yy(0) == Q[1]
            print xx(1) == A[0]
            print yy(1) == A[1]
            Q_to_A = [(xx.derivative()*xx**i/(2*yy)).integral() for i in range(2*gen)] #changed 1210   ##changed 03/04
            I = vector([f(1)-f(0) for f in Q_to_A])  #changed 1210

            #print "I = %s"%I

            P = H(P[0], P[1])
        ###plus a Coleman integral to offset the constant if P, A aren't in the same residue disc #this will change when it's weierstrass
            xm, ym = self.monsky_washnitzer_gens()
            omega = self.invariant_differential()

            if ((P[0]).list()[0] != (A[0]).list()[0] or (P[1]).list()[0] != (A[1]).list()[0]):
                I2 = vector([self.coleman_integral(omega*xm**i, divisor[0][1],
                                                   divisor[1][1]) for i in range(2*gen)])
        else:
            I = vector([self.coleman_integral(omega*xm**i, divisor[0][1], A) for i in range(2*gen)]) #weierstrass case
        print "I = %s"%I
        print "I2 = %s"%I2
        w  = self.frob_diff_nw(divisor, x, y, prec=10)

        v = [(w*(I[n]+I2[n])) for n in range(2*gen)]

        v = [f.residue() for f in v]
        #v = [(I[n])+(I2[n])+Fw_i[n] for n in range(2*gen)]
        #print "v = %s (before trace)"%v
        try:
            return vector([f.trace() for f in v])
        except AttributeError:
            return vector(self.base_ring(), v)
        #return (w*(I+I3+Fw_i)).residue_at_0().trace()

    def is_good(self, P, R):
        """
        checks if P is good wrt R

        EXAMPLES::
        """
        # S = div2[1][1]
        if P[0] == R[0] and P[1] == -R[1]:
            return True
        else:
            return not self.is_neg_disc(P, R)

    def diff(self, divisor, x, y, tiny=False, alpha=False):
        """
        ###needs to be fixed to account for neg discs
        writing differential with residue divisor "divisor" in terms of x, y
        (an interpolation usually)

        alpha=True: truncates diffP to get rid of meaningless terms

        EXAMPLES::
        """
        P = divisor[0][1]
        Q = divisor[1][1]
        a, b, nn = P
        c, d, nn = Q
        g = self.genus()
        f = self.hyperelliptic_polynomials()[0]
        if (a == c and b == -d):
            return b*x.derivative()/(y*(x-a))

        elif (tiny is False or self.is_good(P, (x(Integer(0)), y(Integer(0))))):
            forP = (y+b)/(x-a)
        else:
            hP = ((f(x)-f(a))/(x-a)).truncate(2*g+1)
            if hP.list()[-1] == 0:
                print "bug with division"
                print "this only works for genus 1 curves"
                if g > 1:
                    return "sorry, this only works for genus 1 right now"
                else:
                    hP = x**2+x*a+a**2 + f.list()[2]*(x+a)+f.list()[1]
            forP = hP/(y-b)
            #print "forP = %s"%forP
            #print "old forP = %s"%((y+b)/(x-a))
        if tiny==False or self.is_good(Q, (x(Integer(0)),y(Integer(0)))):
            forQ = (y+d)/(x-c)
        else:
            print "is bad and recomputing diff"
            hQ = ((f(x)-f(c))/(x-c)).truncate(2*g+1)
            if hQ.list()[-1] == 0:
                print "bug with division"
                print "this only works for genus 1 curves"
                if g > 1:
                    return "sorry, this only works for genus 1 right now"
                else:
                    hQ = x**2+x*c+c**2 + f.list()[2]*(x+c)+f.list()[1]
            forQ =hQ/(y-d)
            #print "forQ = %s"%forQ
            #print "old forQ = %s"%((y+d)/(x-c))
        #print "forP = %s"%forP
        #print "forQ = %s"%forQ

        #forP = (y+b)/(x-a)
        #forQ = (y+d)/(x-c)
        if alpha:
            #print forP.list()
            #print "this is the length: %s"%len(forP.list())
            i = 0
            while i < len(forP.list())-3:
                if (forP.list()[i] == 0 and forP.list()[i+1] == 0 and forP.list()[i+2] == 0):
                    forP = forP.truncate(i)
                    break
                else:
                    i += 1
            #print "forP = %s"%forP
            i = 0
            while i < len(forQ.list())-3:
                if (forQ.list()[i]==0 and forQ.list()[i+1]==0 and forQ.list()[i+2]==0):
                    forQ = forQ.truncate(i)
                    break
                else:
                    i += 1
            #print "forQ = %s"%forQ
        #return derivative(x)/(2*y)*((y+b)/(x-a)-(y+d)/(x-c))
        return x.derivative()/(2*y)*(forP-forQ)

    def frob_diff_nw(self, divisor, x, y , prec=7):
        """
        ###needs to be fixed to account for negative discs
        Action of Frobenius on differential for nonweierstrass point (x, y are local coordinates of that point)

        EXAMPLES::
        """
        p = self.base_ring().prime()
        P = divisor[0][1]
        Q = divisor[1][1]
        a, b, nn = P
        c, d, nn = Q
        if (a == c and b == -d):
            return p*b*x**(p-1)*x.derivative()*self.recip_froby(x, y, prec)/(x**p-a)
        phiy = self.froby(x, y, prec)
        forP = (phiy + b)/(x**p-a)
        forQ = (phiy + d)/(x**p-c)
        return (p/2)*x**(p-1)*x.derivative()*self.recip_froby(x, y, prec)*(forP - forQ)

    def frob_diff_wstrass(self, divisor, x, y, prec=10):
        r"""
        Returns the action of Frobenius on the differential associated to
        the divisor `(P) - (Q)` with respect to local coordinates of a
        Weierstrass point

        .. NOTE::

            This agrees with :meth:`frob_diff_w_alt` for Weierstrass points
            (but this is faster).

        This will not work for the prime `p = 2`.

        This is `phi^*w`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x*(x-1)*(x+9))
            sage: K = Qp(7, 10)
            sage: HK = H.change_ring(K)
            sage: P = HK(-1,4)
            sage: Pprime = HK(25/16, 195/64)

        AUTHOR:

        - Jennifer Balakrishnan (2008-02)
        """
        p = self.base_ring().prime()
        a, b, z = divisor[0][1]
        c, d, zz = divisor[1][1]
        return p*x**(p-1)/(2*(x**p-a)*(x**p-c))*(a-c+((b-d)*x**p+a*d-b*c)*self.recip_froby(x, y, prec))*x.derivative()

    def finite_weierstrass_points(self):
        """
        Gives a list of finite Weierstrass points of self

        .. TODO:: 

            This should be moved to hyperelliptic_generic

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x*(x-1)*(x+9))
        """
        f = self.hyperelliptic_polynomials()[0]
        return [j*[(i, 0)] for i, j in f.roots()]

    def alpha(self, divisor, prec=5, all=False):
        """
        Returns alpha = phi^*w - p*w in the various local coordinates

        here w has residue divisor (P)-(Q) and

        if all = True: P, Q, pth roots, then weierstrass
        else: just finite weierstrass

        #alpha[0] is at P
        #alpha[1] is at Q
        #alpha[2] is at zeta_p x(P)^(1/p)
        #alpha[3] is at zeta_p x(Q)^(1/p)
        have to re-index the following:

        alpha[4] is at finite weierstrass_1
        ...
        alpha[2g+4] is at finite weierstrass_{2g+1}
        alpha[2g+5] is at infnity

        So we have the following consistency checks (since res in each
        disc is supposed to be 0):
        (?? this also works when supp(divisor) is just in 1 disc)

        alpha[0].residue() + alpha[2].residue_at_0().trace() = 0
        alpha[1].residue() + alpha[3].residue_at_0().trace() = 0
        alpha[4].residue() = ... = alpha[2g+5].residue() = 0

        TESTS::

            sage: R.<x> = QQ['x']
            sage: K = Qp(7,6)
            sage: H = HyperellipticCurve(x*(x-1)*(x+9))
            sage: HK = H.change_ring(K)
            sage: P = HK(-1, 4)
            sage: Q = HK(9, -36)
            sage: a = HK.alpha([(1, P), (-1, Q)])
            sage: (a[2]).residue().trace()+(a[0]).residue()
            O(7^5)
            sage: (a[3]).residue().trace()+(a[1]).residue()
            O(7^4)
            sage: [(a[i]).residue() for i in range(4,8)]
            [0, 0, 0, 0]
        """
        g = self.genus()
        p = self.base_ring().prime()
        prec = self.base_ring().precision_cap()
        # this overrides what the user inputs, maybe change?
        if all:
            D1 = [divisor[0][1], divisor[1][1]]
        try:
            D2 = self._pth_roots
        except AttributeError:
            D2 = []    # what to do here is unclear to me ???
        # a vector of phi^*w
        frob = []
        diff = []
        if all:
            for i in range(2):
                x, y = self.local_coord(D1[i],20)
                frob += [self.frob_diff_nw(divisor, x, y, prec)]
                diff += [self.diff(divisor, x, y)]
            for i in range(2):
                x, y = self.local_coord(D2[i],20)
                frob += [self.frob_diff_nw(divisor, x, y, prec)]
                diff += [self.diff(divisor, x, y)]
        r = self._fwstrass
        for R in r:
            if R[0] != 0:
                x, y = self.local_coord(R, 2*p*prec - p - 3 + 2*g-1)
                # this seems to be the min for prec=5
                frob += [self.frob_diff_wstrass(divisor, x, y, prec)]
                diff += [self.diff(divisor, x, y)]
        return [frob[i] - p*diff[i] for i in range(len(frob))]

    def psi_frob_diff(self, divisor, prec=10):
        """
        This computes Psi(Frob(diff)), where diff is the form
        corresponding with residue divisor "divisor"

        #1210 - changed prec from 5 to 10

        #this works again as is, 12/10!

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x*(x-1)*(x+9))
            sage: K = Qp(7, 10)
            sage: HK = H.change_ring(K)
            sage: P = HK(-1, 4)
            sage: Q = HK(9, -36)
            sage: Pprime = HK(-1, -4)
            sage: Qprime = HK(9, 36)
            sage: HK.init_height([(1, P), (-1, Q)],[(1, Pprime),(-1, Qprime)], 10)
            sage: HK.psi_frob_diff([(1, P), (-1, Q)])
            (2*7 + 2*7^2 + 2*7^3 + 5*7^4 + 5*7^5 + 4*7^6 + 2*7^7 + 7^8 + O(7^9), 4*7 + 2*7^3 + 5*7^7 + 3*7^8 + O(7^9))
            sage: HK.frob_psi_diff([(1, P), (-1, Q)])
            (2*7 + 2*7^2 + 2*7^3 + 5*7^4 + 5*7^5 + 4*7^6 + 2*7^7 + 7^8 + 3*7^9 + O(7^10), 4*7 + 2*7^3 + 5*7^7 + 3*7^8 + 5*7^9 + O(7^10))

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x*(x-1)*(x+9))
            sage: K = Qp(7, 10)
            sage: HK = H.change_ring(K)
            sage: P = HK(-1, 4)
            sage: Q = HK(9, -36)
            sage: Pprime = HK(-1, -4)
            sage: Qprime = HK(9, 36)
            sage: HK.init_height([(1,Q),(-1,Qprime)],[(1, P),(-1, Pprime)], 10)
            sage: HK.psi_frob_diff([(1,Q),(-1,Qprime)])
            (5*7 + 2*7^2 + 4*7^5 + 7^6 + 7^7 + O(7^8), 7 + 5*7^2 + 2*7^3 + 4*7^4 + 3*7^5 + 6*7^6 + 6*7^7 + O(7^8))
            sage: HK.frob_psi_diff([(1,Q),(-1,Qprime)])
            (5*7 + 2*7^2 + 4*7^5 + 7^6 + 7^7 + 3*7^9 + O(7^10), 7 + 5*7^2 + 2*7^3 + 4*7^4 + 3*7^5 + 6*7^6 + 6*7^7 + 2*7^8 + 4*7^9 + O(7^10))

        Does not work for::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x*(x-1)*(x+9))
            sage: K = Qp(7, 10)
            sage: HK = H.change_ring(K)
            sage: P = HK(-1, 4)
            sage: Q = HK(9, -36)
            sage: Pprime = HK(-1, -4)
            sage: Qprime = HK(9, 36)
            sage: HK.init_height([(1, P),(-1, Pprime)],[(1,Q),(-1,Qprime)], 10)
            sage: HK.psi_frob_diff([(1, P),(-1, Pprime)])
            sage: HK.frob_psi_diff([(1, P),(-1, Pprime)])

        """
        # cup = self._cpm
        cup = self.cup_product_matrix()  # prec
        g = self.genus()
        p = self.base_ring().prime()
        P = divisor[0][1]
        Q = divisor[1][1]

        local_at_P = self.sum_of_local_symbols_extension([(1, P), (-1, Q)],
                                                         P, P, True)
        #print "local_at_P = %s"%local_at_P
        local_at_Q = self.sum_of_local_symbols_extension([(1, P), (-1, Q)],
                                                         P, Q, True)
        local = local_at_P + local_at_Q

        #print "local = %s"%local
        r = self._fwstrass
        local_wstrass=vector(self.base_ring(), [0]*2*g)
        for R in r:
            if R[0] != 0:
                x, y = self.local_coord(R, 20*p)  # this seems to be the min for prec=5
                frob = self.frob_diff_wstrass(divisor, x, y, prec)
                local_wstrass = local_wstrass + vector(self.base_ring(), [(frob* ((x**j*x.derivative()/(2*y) ).integral()) ).residue() for j in range(2*g)])
        #print "local_wstrass = %s"%local_wstrass
        return cup**(-1)*(local + local_wstrass)

    def frob_psi_diff(self, divisor, prec=20):
        """
        Computes Frobenius of Psi(diff), where diff has residue divisor
        "divisor"

        works!
        (see above for consistency check)

        EXAMPLES::
        """
        from sage.schemes.elliptic_curves.monsky_washnitzer import matrix_of_frobenius_hyperelliptic
        M_frob, forms = self._frob_calc = matrix_of_frobenius_hyperelliptic(self)
        if divisor == self._div1:
            #print "using cached value of log(div1)"
            psiw = self._diff_log_div1
        elif divisor == self._div2:
            #print "using cashed value of log(div2)"
            psiw = self._diff_log_div2
        else:
            #print "recomputing log because something's wrong"
            psiw = self.differential_log(divisor, prec)
        V = self._subspace
        return M_frob*V*psiw

    def psi_alpha(self, divisor, prec=20):
        """
        Computes Psi(alpha)= Psi(phi^*w-p*w) as phi^*(Psi(w))-p*Psi(w)

        EXAMPLES::
        """
        frob_psiw = self.frob_psi_diff(divisor, prec)
        p = self.base_ring().prime()
        #psiw = self.differential_log(divisor, prec)
        if divisor == self._div1:
            psiw = self._diff_log_div1
        elif divisor == self._div2:
            psiw = self._diff_log_div2
        else:
            psiw = self.differential_log(divisor, prec)
        V = self._subspace
        return frob_psiw - p*V*psiw

    def res_alpha_int_beta(self, divisor1, divisor2, prec=5):
        """
        returns sum(res(alpha*(int(beta) + c))) where c is the right constant of integration
        and the sum is over non-Weierstrass poles of alpha

        EXAMPLES::
        """
    #    print "res_alpha_int_beta"
        p = self.base_ring().prime()
        g = self.genus()
        pth_roots = self._pth_roots
        P = divisor1[0][1]
        ####integrate beta from P to P_i, take the trace
        x, y = self.local_analytic_interpolation_cyclotomic(P, pth_roots[0][0], p*prec+2*g+1)    # will have to change prec
        betaP = self.diff(divisor2, x, y, True, True)
        int_betaP = (betaP).integral()
        print "int_betaP = %s" % int_betaP
        I1 = int_betaP(1)-int_betaP(0)
        print "I1 = %s" % I1
        if I1 != 0:
            I1 = I1.trace()
        else:
            I1 = 0
        print "I1 via trace = %s" % I1
        Q = divisor1[1][1]
        if self.is_same_disc(P,Q):
            ##if in same disc, then trace(\int(beta, P,Q_i))
            xx, yy = self.local_analytic_interpolation_cyclotomic(P, pth_roots[1][0],2*p*prec)
            betaQ = self.diff(divisor2, xx, yy, True, True)
            int_betaQ = (betaQ).integral()
            I2 = int_betaQ(1)-int_betaQ(0)

            ##also have to integrate beta from P to Q
            xx2, yy2, z = self.local_analytic_interpolation(P, Q, 2*p*prec)  ##prev: lai2
            fix = self.diff(divisor2, xx2, yy2, True, True)
            fix = (fix).integral()
            fix = fix(1)-fix(0)
            xA, yA = self.local_coord(Q, prec=30)
            alphaQ = self.frob_diff_nw(divisor1, xA, yA, prec=10) - p*self.diff(divisor1, xA, yA)
            const = alphaQ.list()[0]
            #print "const = %s"%const
            print "same disc"
            #print "I2 = %s"%(-I2.trace() + const*fix)
            return I1 - I2.trace() +  const*fix

        else:

            xx, yy = self.local_analytic_interpolation_cyclotomic(Q, pth_roots[1][0],2*p*prec)       ###will have to change prec
            betaQ = self.diff(divisor2, xx, yy,True, True)
            int_betaQ = (betaQ).integral()
            #print "int_betaQ = %s"%int_betaQ
            I2 = int_betaQ(1)-int_betaQ(0)
            print "I2 = %s"%I2

            if I2 != 0:
                I2 = I2.trace()
            else:
                I2 = 0
            print "I2 trace = %s"%I2.trace()
            return I1-I2

    def res(self, divisor1, divisor2, prec=12):
        """
        #changed default of prec = 5 to prec = 12
        returns sum(res(alpha*integral(beta)))
        (but need to sum over all of the right residue discs)
        alpha[0] is at P
        alpha[1] is at Q
        alpha[2] is at zeta_p x(P)^(1/p)
        alpha[3] is at zeta_p x(Q)^(1/p)
        alpha[4] is at finite weierstrass_1
        ...
        alpha[2g+4] is at finite weierstrass_{2g+1}
        alpha[2g+5] is at infinity

        EXAMPLES::
        """
        from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
        p = self.base_ring().prime()
        #print "prec = %s"%prec
        A = self.alpha(divisor1, prec)
        ###need to cache values here
        #div3 = [self.find_pth_root_point(divisor1[i][1]) for i in range(2)]
        div3 = self._pth_roots
        #A = A.list()
        B = []
        r = self._fwstrass
        prec = self.base_ring().precision_cap()
        g = self.genus()
        for R in r:
            if R[0] != 0:
                x, y = self.local_coord(R, 2*p*prec - p -3 + 2*g-1)
                B += [self.diff(divisor2, x, y)]
        print "for A"
        print [f.valuation() for f in A]
        print [f.degree() for f in A]
        print "for B"
        print [f.valuation() for f in B]
        dd = [f.degree() for f in B]
        print dd
        Anew = []
        for i in range(len(A)):
            v = B[i].valuation()
            try:
                Anew += [A[i].truncate_neg(v-max(dd))]
            except AttributeError:
                Anew = [A[i]]
        res = [Anew[i]*((B[i]).integral()) for i in range(len(Anew))]
        res = [theta.residue() for theta in res]
        t = []
        for f in res:
            try:
                t = t + [f.trace()]
            except AttributeError:
                t = t + [f]
        return sum(t)

    def psiA_cup_psiB(self, divisor1, divisor2, prec=5):
        """
        Returns the cup product of psiA and psiB

        JSB (2008-05)

        EXAMPLES::
        """
        psiA = self.psi_alpha(divisor1)
        V = self._subspace
        psiA = psiA
        print "psiA wrt standard basis is %s" % psiA

        psiB = self._diff_log_div2
        psiB = V*psiB
        print "psiB wrt standard basis is = %s" % psiB
        return self.cup(psiA, psiB)

    def eta_integral(self, divisor1, divisor2, prec=5):
        """
        Integral of eta

        EXAMPLES::
        """
        #coeffs = self.differential_log_holomorphic(divisor1, prec)
        coeffs = self._diff_log_hol_div1
        int = self.coleman_integrals_on_basis(divisor2[1][1], divisor2[0][1])
        return coeffs*int

    def is_neg_disc(self, P, Q):
        """
        checks if P is in disc(-Q)

        EXAMPLES::
        """
        K = self.base_ring()
        if (P[0]).parent() == (Q[0]).parent():
            if (P[0].list()[0] == Q[0].list()[0] and K(P[1].list()[0])==K((-Q[1]).list()[0])):
                return True
            else:
                return False
        else:
            if (P[0].list()[0] == Q[0].list()[0].list()[0] and K(P[1].list()[0])==K((-Q[1]).list()[0].list()[0])):
                return True
            else:
                return False

    def omega_integral(self, divisor1, divisor2, prec=5):
        """
        EXAMPLES::
        """
        p = self.base_ring().prime()
        P = divisor1[0][1]
        Q = divisor1[1][1]
        R = divisor2[0][1]
        S = divisor2[1][1]
        f = self.hyperelliptic_polynomials()[0]
        if self.is_same_disc(R, S):
            x, y, z = self.local_analytic_interpolation(S, R, 5*prec)
            int_diff = self.diff(divisor1, x, y, True, False).integral()
            I = int_diff(Integer(1)) - int_diff(Integer(0))

        FR = self.frobenius(R)
        if R != FR:
            x, y, z = self.local_analytic_interpolation(R, FR, 5*prec)
            R_to_FR = self.diff(divisor1, x, y, tiny=True).integral()
            R_to_FR = R_to_FR(Integer(1)) - R_to_FR(Integer(0))
        else:
            R_to_FR = 0

        FS = self.frobenius(S)
        if S != FS:
            x, y, z = self.local_analytic_interpolation(FS, S, 5*prec)
            FS_to_S = self.diff(divisor1, x, y, tiny=True).integral()
            FS_to_S = FS_to_S(Integer(1)) - FS_to_S(Integer(0))
        else:
            FS_to_S = 0

        res = self.res(divisor1, divisor2, prec)
        ab = self.res_alpha_int_beta(divisor1, divisor2, prec)
        cups = self.psiA_cup_psiB(divisor1, divisor2, prec)

        return (cups + res + ab - FS_to_S - R_to_FR)/(1-p)

    def height(self, divisor1, divisor2, prec=5):
        """
        The p-part of the Coleman-Gross height pairing of divisor1 and
        divisor2

        If self has ordinary reduction at self.base_ring().prime(),
        the height pairing is symmetric.

        GENUS 1 EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x*(x-1)*(x+9))
            sage: K = Qp(7, 10)
            sage: HK = H.change_ring(K)
            sage: P = HK(9,36)
            sage: Q = HK.teichmuller(P)
            sage: Pprime = HK(-4, 10)
            sage: Qprime = HK.teichmuller(Pprime)
            sage: HK.height([(1, P),(-1,Q)],[(1, Pprime),(-1, Qprime)], 10)
            2*7^2 + 5*7^3 + 7^4 + 7^5 + 2*7^6 + 3*7^7 + 7^8 + 3*7^9 + O(7^10)
            sage: HK.height([(1, Pprime),(-1,Qprime)],[(1, P),(-1, Q)], 10)
            2*7^2 + 5*7^3 + 7^4 + 7^5 + 2*7^6 + 3*7^7 + 7^8 + 3*7^9 + O(7^10)

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x*(x-1)*(x+9))
            sage: K = Qp(7, 10)
            sage: HK = H.change_ring(K)
            sage: P = HK(-1, 4)
            sage: Q = HK(-1, -4)
            sage: R = HK(-4, -10)
            sage: S = HK(-4, 10)
            sage: Pprime = HK(25/16, 195/64)
            sage: Qprime = HK(25/16, -195/64)

        Test that h_7(P-Q, R-S) + h_7(P-Q, S-Pprime) = h_7(P-Q,R-Pprime)::

            sage: HK.height([(1, P),(-1,Q)],[(1,R),(-1,S)],9)
            6*7 + 5*7^2 + 2*7^3 + 4*7^4 + 7^5 + 3*7^6 + 7^7 + 4*7^9 + O(7^10)
            sage: HK.height([(1, P),(-1,Q)],[(1,S),(-1, Pprime)],9)
            4*7 + 2*7^2 + 3*7^3 + 6*7^4 + 5*7^5 + 4*7^6 + 6*7^7 + 2*7^8 + 5*7^9 + O(7^10)
            sage: HK.height([(1, P),(-1,Q)],[(1,R),(-1, Pprime)],9)
            3*7 + 7^2 + 6*7^3 + 3*7^4 + 7^6 + 7^7 + 3*7^8 + 2*7^9 + O(7^10)

            sage: 6*7 + 5*7^2 + 2*7^3 + 4*7^4 + 7^5 + 3*7^6 + O(7^7)+4*7 + 2*7^2 + 3*7^3 + 6*7^4 + 5*7^5 + 4*7^6 + 6*7^7 + 2*7^8 + 5*7^9 + O(7^10)
            3*7 + 7^2 + 6*7^3 + 3*7^4 + 7^6 + O(7^7)

        Test that h_7(Pprime-P, R-S) = h_7(Q-Qprime, R-S), where (Pprime)-(P) ~ (Q)-(Qprime)::

            sage: HK.height([(1,Pprime),(-1,P)],[(1,R),(-1,S)],9)
            3*7 + 7^3 + 7^4 + 7^5 + 2*7^6 + 2*7^7 + 7^8 + O(7^10)

            sage: HK.height([(1,Q),(-1,Qprime)],[(1,R),(-1,S)],9)
            3*7 + 7^3 + 7^4 + 7^5 + 2*7^6 + 2*7^7 + 7^8 + O(7^10)

        GENUS 2 EXAMPLES
        (with respect to W the unit root subspace)::

            sage: R.<x> = PolynomialRing(pAdicField(11, 10))
            sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: P = H(-4,24)
            sage: Pprime = H.teichmuller(P)
            sage: Q = H(5,30)
            sage: Qprime = H.teichmuller(Q)
            sage: H.height([(1,Q),(-1,Qprime)],[(1,P),(-1,Pprime)], 10)
            6*11^2 + 9*11^3 + 4*11^4 + 2*11^5 + 6*11^6 + 4*11^7 + 6*11^8 + 11^9 + O(11^10)
            sage: H.height([(1,P),(-1,Pprime)],[(1,Q),(-1,Qprime)], 10)
            6*11^2 + 9*11^3 + 4*11^4 + 2*11^5 + 6*11^6 + 4*11^7 + 6*11^8 + 11^9 + O(11^10)

            sage: R.<x> = PolynomialRing(pAdicField(11, 10))
            sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: P = H(-4,24)
            sage: Pprime = H(-4, -24)
            sage: Q = H(5,30)
            sage: Qprime = H(5, -30)
            sage: H.height([(1,P),(-1,Pprime)],[(1,Q),(-1,Qprime)], 10)
            6*11^-1 + 10 + 7*11 + 6*11^2 + 3*11^3 + 7*11^4 + 7*11^5 + 11^6 + O(11^8)
            sage: H.height([(1,Q),(-1,Qprime)],[(1,P),(-1,Pprime)], 10)
            6*11^-1 + 10 + 7*11 + 6*11^2 + 3*11^3 + 7*11^4 + 7*11^5 + 11^6 + O(11^8)

            sage: R.<x> = Qp(11, 10)['x']
            sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: P = H(-4,24)
            sage: Q = H(5,30)
            sage: R = H(1,6)
            sage: S = H(-2, 12)
            sage: H.height([(1,P),(-1,Q)],[(1,R),(-1,S)], 10)
            7*11^-1 + 3 + 8*11 + 6*11^2 + 5*11^3 + 7*11^4 + 3*11^5 + 9*11^6 + 6*11^7 + O(11^8)
            sage: H.height([(1,R),(-1,S)],[(1,P),(-1,Q)], 10)
            7*11^-1 + 3 + 8*11 + 6*11^2 + 5*11^3 + 7*11^4 + 3*11^5 + 9*11^6 + 6*11^7 + O(11^8)

        """
        self.init_height(divisor1, divisor2, prec)
        omega = self.omega_integral(divisor1, divisor2, prec)
        eta = self.eta_integral(divisor1, divisor2, prec)
        return omega - eta

    def is_ordinary(self):
        """
        Determine if self.base_ring().prime() is a prime of ordinary reduction

        EXAMPLES::

            sage: R.<x> = PolynomialRing(pAdicField(11, 10))
            sage: H = HyperellipticCurve(x^5-23*x^3+18*x^2+40*x)
            sage: H.is_ordinary()
            True
        """
        try:
            M_frob, forms = self._frob_calc
        except AttributeError:
            from sage.schemes.elliptic_curves.monsky_washnitzer import matrix_of_frobenius_hyperelliptic
            M_frob, forms = self._frob_calc = matrix_of_frobenius_hyperelliptic(self)
        return bool(M_frob.charpoly().list()[self.genus()].valuation() == 0)

    def init_height(self, divisor1, divisor2, prec=5):
        """
        initializes and caches certain quantities for height computation so as to avoid repeats

        EXAMPLES::
        """
        from sage.schemes.elliptic_curves.monsky_washnitzer import matrix_of_frobenius_hyperelliptic
        g = self.genus()
        K = self.base_ring()
        p = K.prime()
        try:
            M_frob, forms = self._frob_calc
        except AttributeError:
            M_frob, forms = self._frob_calc = matrix_of_frobenius_hyperelliptic(self)
        A = matrix(K, M_frob).transpose() - 1
        m = max(a.valuation() for a in A.list())
        self._prec = max(m, prec)
        print "Current working precision is %s. If you would like a final answer with %s guaranteed digits of precision, increase working precision to %s." % (prec, prec, self._prec+prec)
        self._fwstrass = self.finite_weierstrass_points()
        self._div1 = divisor1
        self._div2 = divisor2
        self._cpm = self.cup_product_matrix()  # prec
        self._pth_roots = [self.find_pth_root_point(divisor1[i][1]) for i in range(2)]

        self._diff_log_div1 = self.differential_log(divisor1)  # prec
        self._diff_log_div2 = self.differential_log(divisor2)
        self._diff_log_hol_div1 =  vector([self._diff_log_div1[i] for i in range(g)]+[0]*g)
