# -*- coding: utf-8 -*-
r"""
Points on elliptic curves

The base class ``EllipticCurvePoint_field``, derived from
``AdditiveGroupElement``, provides support for points on elliptic
curves defined over general fields.  The derived classes
``EllipticCurvePoint_number_field`` and
``EllipticCurvePoint_finite_field`` provide further support for point
on curves defined over number fields (including the rational field
`\QQ`) and over finite fields.

The class ``EllipticCurvePoint``, which is based on
``SchemeMorphism_point_projective_ring``, currently has little extra
functionality.

EXAMPLES:

An example over `\QQ`::

    sage: E = EllipticCurve('389a1')
    sage: P = E(-1,1); P
    (-1 : 1 : 1)
    sage: Q = E(0,-1); Q
    (0 : -1 : 1)
    sage: P+Q
    (4 : 8 : 1)
    sage: P-Q
    (1 : 0 : 1)
    sage: 3*P-5*Q
    (328/361 : -2800/6859 : 1)

An example over a number field::

    sage: K.<i> = QuadraticField(-1)
    sage: E = EllipticCurve(K,[1,0,0,0,-1])
    sage: P = E(0,i); P
    (0 : i : 1)
    sage: P.order()
    +Infinity
    sage: 101*P-100*P == P
    True

An example over a finite field::

    sage: K.<a> = GF((101,3))
    sage: E = EllipticCurve(K,[1,0,0,0,-1])
    sage: P = E(40*a^2 + 69*a + 84 , 58*a^2 + 73*a + 45)
    sage: P.order()
    1032210
    sage: E.cardinality()
    1032210

Arithmetic with a point over an extension of a finite field::

    sage: k.<a> = GF((5,2))
    sage: E = EllipticCurve(k,[1,0]); E
    Elliptic Curve defined by y^2 = x^3 + x over Finite Field in a of size 5^2
    sage: P = E([a,2*a+4])
    sage: 5*P
    (2*a + 3 : 2*a : 1)
    sage: P*5
    (2*a + 3 : 2*a : 1)
    sage: P + P + P + P + P
    (2*a + 3 : 2*a : 1)

::

    sage: F = Zmod(3)
    sage: E = EllipticCurve(F,[1,0]);
    sage: P = E([2,1])
    sage: import sys
    sage: n = sys.maxsize
    sage: P*(n+1)-P*n == P
    True

Arithmetic over `\ZZ/N\ZZ` with composite `N` is supported.  When an
operation tries to invert a non-invertible element, a
ZeroDivisionError is raised and a factorization of the modulus appears
in the error message::

    sage: N = 1715761513
    sage: E = EllipticCurve(Integers(N),[3,-13])
    sage: P = E(2,1)
    sage: LCM([2..60])*P
    Traceback (most recent call last):
    ...
    ZeroDivisionError: Inverse of 1520944668 does not exist (characteristic = 1715761513 = 26927*63719)

AUTHORS:

- William Stein (2005) -- Initial version

- Robert Bradshaw et al....

- John Cremona (Feb 2008) -- Point counting and group structure for
  non-prime fields, Frobenius endomorphism and order, elliptic logs

- John Cremona (Aug 2008) -- Introduced ``EllipticCurvePoint_number_field`` class

- Tobias Nagel, Michael Mardaus, John Cremona (Dec 2008) -- `p`-adic elliptic logarithm over `\QQ`

- David Hansen (Jan 2009) -- Added ``weil_pairing`` function to ``EllipticCurvePoint_finite_field`` class

- Mariah Lenox (March 2011) -- Added ``tate_pairing`` and ``ate_pairing``
  functions to ``EllipticCurvePoint_finite_field`` class

"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import math

import sage.plot.all as plot

from sage.rings.padics.factory import Qp
from sage.rings.padics.precision_error import PrecisionError

import sage.rings.all as rings
import sage.rings.abc
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
import sage.groups.generic as generic
from sage.libs.pari.all import pari, PariError
from cypari2.pari_instance import prec_words_to_bits
from sage.structure.sequence import Sequence
from sage.structure.richcmp import richcmp

from sage.structure.coerce_actions import IntegerMulAction

from sage.schemes.curves.projective_curve import Hasse_bounds
from sage.schemes.projective.projective_point import (SchemeMorphism_point_projective_ring,
                                                      SchemeMorphism_point_abelian_variety_field)
from sage.schemes.generic.morphism import is_SchemeMorphism

from .constructor import EllipticCurve

oo = rings.infinity       # infinity


class EllipticCurvePoint(SchemeMorphism_point_projective_ring):
    """
    A point on an elliptic curve.
    """
    pass


class EllipticCurvePoint_field(SchemeMorphism_point_abelian_variety_field):
    """
    A point on an elliptic curve over a field.  The point has coordinates
    in the base field.

    EXAMPLES::

        sage: E = EllipticCurve('37a')
        sage: E([0,0])
        (0 : 0 : 1)
        sage: E(0,0)               # brackets are optional
        (0 : 0 : 1)
        sage: E([GF(5)(0), 0])     # entries are coerced
        (0 : 0 : 1)

        sage: E(0.000, 0)
        (0 : 0 : 1)

        sage: E(1,0,0)
        Traceback (most recent call last):
        ...
        TypeError: Coordinates [1, 0, 0] do not define a point on
        Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    ::

        sage: E = EllipticCurve([0,0,1,-1,0])
        sage: S = E(QQ); S
        Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

        sage: K.<i>=NumberField(x^2+1)
        sage: E=EllipticCurve(K,[0,1,0,-160,308])
        sage: P=E(26,-120)
        sage: Q=E(2+12*i,-36+48*i)
        sage: P.order() == Q.order() == 4  # long time (3s)
        True
        sage: 2*P==2*Q
        False

    ::

        sage: K.<t>=FractionField(PolynomialRing(QQ,'t'))
        sage: E=EllipticCurve([0,0,0,0,t^2])
        sage: P=E(0,t)
        sage: P,2*P,3*P
        ((0 : t : 1), (0 : -t : 1), (0 : 1 : 0))

    TESTS::

        sage: loads(S.dumps()) == S
        True
        sage: E = EllipticCurve('37a')
        sage: P = E(0,0); P
        (0 : 0 : 1)
        sage: loads(P.dumps()) == P
        True
        sage: T = 100*P
        sage: loads(T.dumps()) == T
        True

    Test pickling an elliptic curve that has known points on it::

        sage: e = EllipticCurve([0, 0, 1, -1, 0]); g = e.gens(); loads(dumps(e)) == e
        True

    Test that the refactoring from :trac:`14711` did preserve the behaviour
    of domain and codomain::

        sage: E=EllipticCurve(QQ,[1,1])
        sage: P=E(0,1)
        sage: P.domain()
        Spectrum of Rational Field
        sage: K.<a>=NumberField(x^2-3,'a')
        sage: P=E.base_extend(K)(1,a)
        sage: P.domain()
        Spectrum of Number Field in a with defining polynomial x^2 - 3
        sage: P.codomain()
        Elliptic Curve defined by y^2 = x^3 + x + 1 over Number Field in a with defining polynomial x^2 - 3
        sage: P.codomain() == P.curve()
        True
    """
    def __init__(self, curve, v, check=True):
        """
        Constructor for a point on an elliptic curve.

        INPUT:

        - curve -- an elliptic curve
        - v -- data determining a point (another point, the integer
                 0, or a tuple of coordinates)

        EXAMPLES::

            sage: E = EllipticCurve('43a')
            sage: P = E([2, -4, 2]); P
            (1 : -2 : 1)
            sage: P == E([1,-2])
            True
            sage: P = E(0); P
            (0 : 1 : 0)
            sage: P=E(2, -4, 2); P
            (1 : -2 : 1)
        """
        point_homset = curve.point_homset()
        if is_SchemeMorphism(v) or isinstance(v, EllipticCurvePoint_field):
            v = list(v)
        elif v == 0:
            # some of the code assumes that E(0) has integral entries
            # regardless of the base ring...
            # R = self.base_ring()
            # v = (R.zero(),R.one(),R.zero())
            v = (0, 1, 0)
        if check:
            # mostly from SchemeMorphism_point_projective_field
            d = point_homset.codomain().ambient_space().ngens()
            if not isinstance(v, (list, tuple)):
                raise TypeError("Argument v (= %s) must be a scheme point, list, or tuple." % str(v))
            if len(v) != d and len(v) != d-1:
                raise TypeError("v (=%s) must have %s components" % (v, d))
            v = Sequence(v, point_homset.value_ring())
            if len(v) == d-1:     # very common special case
                v.append(v.universe()(1))

            n = len(v)
            all_zero = True
            for i in range(n):
                c = v[n-1-i]
                if c:
                    all_zero = False
                    if c == 1:
                        break
                    for j in range(n-i):
                        v[j] /= c
                    break
            if all_zero:
                raise ValueError("%s does not define a valid point "
                                 "since all entries are 0" % repr(v))

            x, y, z = v
            if z == 0:
                test = x
            else:
                a1, a2, a3, a4, a6 = curve.ainvs()
                test = y**2 + (a1*x+a3)*y - (((x+a2)*x+a4)*x+a6)
            if not test == 0:
                raise TypeError("Coordinates %s do not define a point on %s" % (list(v), curve))

        SchemeMorphism_point_abelian_variety_field.__init__(self, point_homset, v, check=False)
        # AdditiveGroupElement.__init__(self, point_homset)

    def _repr_(self):
        """
        Return a string representation of this point.

        EXAMPLES::

            sage: E = EllipticCurve('39a')
            sage: P = E([-2, 1, 1])
            sage: P._repr_()
            '(-2 : 1 : 1)'
        """
        return self.codomain().ambient_space()._repr_generic_point(self._coords)

    def _latex_(self):
        """
        Return a LaTeX representation of this point.

        EXAMPLES::

            sage: E = EllipticCurve('40a')
            sage: P = E([3, 0])
            sage: P._latex_()
            '\\left(3 : 0 : 1\\right)'
        """
        return self.codomain().ambient_space()._latex_generic_point(self._coords)

    def __getitem__(self, n):
        """
        Return the n'th coordinate of this point.

        EXAMPLES::

            sage: E = EllipticCurve('42a')
            sage: P = E([-17, -51, 17])
            sage: [P[i] for i in [2,1,0]]
            [1, -3, -1]
        """
        return self._coords[n]

    def __iter__(self):
        """
        Return the coordinates of this point as a list.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: list(E([0,0]))
            [0, 0, 1]
        """
        return iter(self._coords)

    def __tuple__(self):
        """
        Return the coordinates of this point as a tuple.

        EXAMPLES::

            sage: E = EllipticCurve('44a')
            sage: P = E([1, -2, 1])
            sage: P.__tuple__()
            (1, -2, 1)
        """
        return tuple(self._coords)  # Warning: _coords is a list!

    def _richcmp_(self, other, op):
        """
        Comparison function for points to allow sorting and equality testing.

        EXAMPLES::

            sage: E = EllipticCurve('45a')
            sage: P = E([2, -1, 1])
            sage: P == E(0)
            False
            sage: P+P == E(0)
            True
        """
        if not isinstance(other, EllipticCurvePoint_field):
            try:
                other = self.codomain().ambient_space()(other)
            except TypeError:
                return NotImplemented
        return richcmp(self._coords, other._coords, op)

    def __pari__(self):
        r"""
        Converts this point to PARI format.

        EXAMPLES::

            sage: E = EllipticCurve([0,0,0,3,0])
            sage: O = E(0)
            sage: P = E.point([1,2])
            sage: O.__pari__()
            [0]
            sage: P.__pari__()
            [1, 2]

        The following implicitly calls O.__pari__() and P.__pari__()::

            sage: pari(E).elladd(O,P)
            [1, 2]

        TESTS:

        Try the same over a finite field::

            sage: E = EllipticCurve(GF(11), [0,0,0,3,0])
            sage: O = E(0)
            sage: P = E.point([1,2])
            sage: O.__pari__()
            [0]
            sage: P.__pari__()
            [Mod(1, 11), Mod(2, 11)]

        We no longer need to explicitly call ``pari(O)`` and ``pari(P)``
        after :trac:`11868`::

            sage: pari(E).elladd(O, P)
            [Mod(1, 11), Mod(2, 11)]
        """
        if self[2]:
            return pari([self[0]/self[2], self[1]/self[2]])
        else:
            return pari([0])

    def scheme(self):
        """
        Return the scheme of this point, i.e., the curve it is on.
        This is synonymous with :meth:`curve` which is perhaps more
        intuitive.

        EXAMPLES::

            sage: E=EllipticCurve(QQ,[1,1])
            sage: P=E(0,1)
            sage: P.scheme()
            Elliptic Curve defined by y^2 = x^3 + x + 1 over Rational Field
            sage: P.scheme() == P.curve()
            True
            sage: K.<a>=NumberField(x^2-3,'a')
            sage: P=E.base_extend(K)(1,a)
            sage: P.scheme()
            Elliptic Curve defined by y^2 = x^3 + x + 1 over Number Field in a with defining polynomial x^2 - 3
        """
        # The following text is just not true: it applies to the class
        # EllipticCurvePoint, which appears to be never used, but does
        # not apply to EllipticCurvePoint_field which is simply derived
        # from AdditiveGroupElement.
        #
        # "Technically, points on curves in Sage are scheme maps from
        #  the domain Spec(F) where F is the base field of the curve to
        #  the codomain which is the curve.  See also domain() and
        #  codomain()."

        return self.codomain()

    def order(self):
        r"""
        Return the order of this point on the elliptic curve.

        If the point is zero, returns 1, otherwise raise a
        NotImplementedError.

        For curves over number fields and finite fields, see below.

        .. NOTE::

            :meth:`additive_order` is a synonym for :meth:`order`

        EXAMPLES::

            sage: K.<t>=FractionField(PolynomialRing(QQ,'t'))
            sage: E=EllipticCurve([0,0,0,-t^2,0])
            sage: P=E(t,0)
            sage: P.order()
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of order of a point not implemented over general fields.
            sage: E(0).additive_order()
            1
            sage: E(0).order() == 1
            True
        """
        if self.is_zero():
            return Integer(1)
        raise NotImplementedError("Computation of order of a point "
                                  "not implemented over general fields.")

    additive_order = order

    def curve(self):
        """
        Return the curve that this point is on.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: P = E([-1,1])
            sage: P.curve()
            Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field
        """
        return self.scheme()

    def __bool__(self):
        """
        Return ``True`` if this is not the zero point on the curve.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: P = E(0); P
            (0 : 1 : 0)
            sage: P.is_zero()
            True
            sage: P = E.gens()[0]
            sage: P.is_zero()
            False
        """
        return bool(self[2])

    

    def has_finite_order(self):
        """
        Return ``True`` if this point has finite additive order as an
        element of the group of points on this curve.

        For fields other than number fields and finite fields, this is
        NotImplemented unless self.is_zero().

        EXAMPLES::

            sage: K.<t>=FractionField(PolynomialRing(QQ,'t'))
            sage: E=EllipticCurve([0,0,0,-t^2,0])
            sage: P = E(0)
            sage: P.has_finite_order()
            True
            sage: P=E(t,0)
            sage: P.has_finite_order()
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of order of a point not implemented over general fields.
            sage: (2*P).is_zero()
            True
        """
        if self.is_zero():
            return True
        return self.order() != oo

    is_finite_order = has_finite_order  # for backward compatibility

    def has_infinite_order(self):
        """
        Return True if this point has infinite additive order as an element
        of the group of points on this curve.

        For fields other than number fields and finite fields, this is
        NotImplemented unless self.is_zero().

        EXAMPLES::

            sage: K.<t>=FractionField(PolynomialRing(QQ,'t'))
            sage: E=EllipticCurve([0,0,0,-t^2,0])
            sage: P = E(0)
            sage: P.has_infinite_order()
            False
            sage: P=E(t,0)
            sage: P.has_infinite_order()
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of order of a point not implemented over general fields.
            sage: (2*P).is_zero()
            True
        """
        if self.is_zero():
            return False
        return self.order() == oo

    def plot(self, **args):
        """
        Plot this point on an elliptic curve.

        INPUT:

        - ``**args`` -- all arguments get passed directly onto the point
          plotting function.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: P = E([-1,1])
            sage: P.plot(pointsize=30, rgbcolor=(1,0,0))
            Graphics object consisting of 1 graphics primitive
        """
        if self.is_zero():
            return plot.text("$\\infty$", (-3, 3), **args)
        else:
            return plot.point((self[0], self[1]), **args)

    def _add_(self, right):
        """
        Add self to right.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: P = E([-1,1]); Q = E([0,0])
            sage: P + Q
            (1 : 0 : 1)
            sage: P._add_(Q) == P + Q
            True

        Example to show that bug :trac:`4820` is fixed::

            sage: [type(c) for c in 2*EllipticCurve('37a1').gen(0)]
            [<... 'sage.rings.rational.Rational'>,
            <... 'sage.rings.rational.Rational'>,
            <... 'sage.rings.rational.Rational'>]

        Checks that :trac:`15964` is fixed::

            sage: N = 1715761513
            sage: E = EllipticCurve(Integers(N),[3,-13])
            sage: P = E(2,1)
            sage: LCM([2..60])*P
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse of 1520944668 does not exist
            (characteristic = 1715761513 = 26927*63719)

            sage: N = 35
            sage: E = EllipticCurve(Integers(N),[5,1])
            sage: P = E(0,1)
            sage: 4*P
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse of 28 does not exist
            (characteristic = 35 = 7*5)
        """
        # Use Prop 7.1.7 of Cohen "A Course in Computational Algebraic
        # Number Theory"
        if self.is_zero():
            return right
        if right.is_zero():
            return self
        E = self.curve()
        a1, a2, a3, a4, a6 = E.ainvs()
        x1, y1 = self[0], self[1]
        x2, y2 = right[0], right[1]
        if x1 == x2 and y1 == -y2 - a1*x2 - a3:
            return E(0)  # point at infinity

        if x1 == x2 and y1 == y2:
            try:
                m = (3*x1*x1 + 2*a2*x1 + a4 - a1*y1) / (2*y1 + a1*x1 + a3)
            except ZeroDivisionError:
                R = E.base_ring()
                if R.is_finite():
                    N = R.characteristic()
                    N1 = N.gcd(Integer(2*y1 + a1*x1 + a3))
                    N2 = N//N1
                    raise ZeroDivisionError("Inverse of %s does not exist (characteristic = %s = %s*%s)" % (2*y1 + a1*x1 + a3, N, N1, N2))
                else:
                    raise ZeroDivisionError("Inverse of %s does not exist" % (2*y1 + a1*x1 + a3))
        else:
            try:
                m = (y1-y2)/(x1-x2)
            except ZeroDivisionError:
                R = E.base_ring()
                if R.is_finite():
                    N = R.characteristic()
                    N1 = N.gcd(Integer(x1-x2))
                    N2 = N//N1
                    raise ZeroDivisionError("Inverse of %s does not exist (characteristic = %s = %s*%s)" % (x1-x2, N, N1, N2))
                else:
                    raise ZeroDivisionError("Inverse of %s does not exist" % (x1-x2))

        x3 = -x1 - x2 - a2 + m*(m+a1)
        y3 = -y1 - a3 - a1*x3 + m*(x1-x3)
        # See trac #4820 for why we need to coerce 1 into the base ring here:
        return E.point([x3, y3, E.base_ring()(1)], check=False)

    def _sub_(self, right):
        """
        Subtract right from  self.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: P = E([-1,1]); Q = E([0,0])
            sage: P - Q
            (4 : 8 : 1)
            sage: P - Q == P._sub_(Q)
            True
            sage: (P - Q) + Q
            (-1 : 1 : 1)
            sage: P
            (-1 : 1 : 1)
        """
        return self + (-right)

    def __neg__(self):
        """
        Return the additive inverse of this point.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: P = E([-1,1])
            sage: Q = -P; Q
            (-1 : -2 : 1)
            sage: Q + P
            (0 : 1 : 0)

        Example to show that bug :trac:`4820` is fixed::

            sage: [type(c) for c in -EllipticCurve('37a1').gen(0)]
            [<... 'sage.rings.rational.Rational'>,
            <... 'sage.rings.rational.Rational'>,
            <... 'sage.rings.rational.Rational'>]
        """
        if self.is_zero():
            return self
        E, x, y = self.curve(), self[0], self[1]
        # See trac #4820 for why we need to coerce 1 into the base ring here:
        return E.point([x, -y - E.a1()*x - E.a3(), E.base_ring()(1)], check=False)

    def xy(self):
        """
        Return the `x` and `y` coordinates of this point, as a 2-tuple.
        If this is the point at infinity a ZeroDivisionError is raised.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: P = E([-1,1])
            sage: P.xy()
            (-1, 1)
            sage: Q = E(0); Q
            (0 : 1 : 0)
            sage: Q.xy()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
        """
        if self[2] == 1:
            return self[0], self[1]
        else:
            return self[0]/self[2], self[1]/self[2]

    def is_divisible_by(self, m):
        """
        Return True if there exists a point `Q` defined over the same
        field as self such that `mQ` == self.

        INPUT:

        - ``m`` -- a positive integer.

        OUTPUT:

        (bool) -- True if there is a solution, else False.

        .. WARNING::

            This function usually triggers the computation of the
            `m`-th division polynomial of the associated elliptic
            curve, which will be expensive if `m` is large, though it
            will be cached for subsequent calls with the same `m`.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: Q = 5*E(0,0); Q
            (-2739/1444 : -77033/54872 : 1)
            sage: Q.is_divisible_by(4)
            False
            sage: Q.is_divisible_by(5)
            True

        A finite field example::

            sage: E = EllipticCurve(GF(101),[23,34])
            sage: E.cardinality().factor()
            2 * 53
            sage: Set([T.order() for T in E.points()])
            {1, 106, 2, 53}
            sage: len([T for T in E.points() if T.is_divisible_by(2)])
            53
            sage: len([T for T in E.points() if T.is_divisible_by(3)])
            106

        TESTS:

        This shows that the bug reported at :trac:`10076` is fixed::

            sage: K = QuadraticField(8,'a')
            sage: E = EllipticCurve([K(0),0,0,-1,0])
            sage: P = E([-1,0])
            sage: P.is_divisible_by(2)
            False
            sage: P.division_points(2)
            []

        Note that it is not sufficient to test that
        ``self.division_points(m,poly_only=True)`` has roots::

            sage: P.division_points(2, poly_only=True).roots()
            [(1/2*a - 1, 1), (-1/2*a - 1, 1)]

            sage: tor = E.torsion_points(); len(tor)
            8
            sage: [T.order() for T in tor]
            [2, 4, 4, 2, 1, 2, 4, 4]
            sage: all(T.is_divisible_by(3) for T in tor)
            True
            sage: sorted(T for T in tor if T.is_divisible_by(2))
            [(0 : 1 : 0), (1 : 0 : 1)]
            sage: sorted(Set([2*T for T in tor]))
            [(0 : 1 : 0), (1 : 0 : 1)]
        """
        # Coerce the input m to an integer
        m = Integer(m)

        # Check for trivial cases of m = 1, -1 and 0.
        if m == 1 or m == -1:
            return True
        if m == 0:
            return self == 0  # then m*self=self for all m!
        m = m.abs()

        # Now the following line would of course be correct, but we
        # work harder to be more efficient:
        # return len(self.division_points(m)) > 0

        P = self

        # If P has finite order n and gcd(m,n)=1 then the result is
        # True.  However, over general fields computing P.order() is
        # not implemented.

        try:
            n = P.order()
            if not n == oo:
                if m.gcd(n) == 1:
                    return True
        except NotImplementedError:
            pass

        P_is_2_torsion = (P == -P)
        g = P.division_points(m, poly_only=True)

        if not P_is_2_torsion:
            # In this case deg(g)=m^2, and each root in K lifts to two
            # points Q,-Q both in E(K), of which exactly one is a
            # solution.  So we just check the existence of roots:
            return len(g.roots()) > 0

        # Now 2*P==0

        if m % 2 == 1:
            return True  # P itself is a solution when m is odd

        # Now m is even and 2*P=0.  Roots of g in K may or may not
        # lift to solutions in E(K), so we fall back to the default.
        # Note that division polynomials are cached so this is not
        # inefficient:

        return len(self.division_points(m)) > 0

    def division_points(self, m, poly_only=False):
        r"""
        Return a list of all points `Q` such that `mQ=P` where `P` = self.

        Only points on the elliptic curve containing self and defined
        over the base field are included.

        INPUT:

        - ``m`` -- a positive integer

        - ``poly_only`` -- bool (default: False); if True return
          polynomial whose roots give all possible `x`-coordinates of
          `m`-th roots of self.

        OUTPUT:

        (list) -- a (possibly empty) list of solutions `Q` to `mQ=P`,
        where `P` = self.

        EXAMPLES:

        We find the five 5-torsion points on an elliptic curve::

            sage: E = EllipticCurve('11a'); E
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: P = E(0); P
            (0 : 1 : 0)
            sage: P.division_points(5)
            [(0 : 1 : 0), (5 : -6 : 1), (5 : 5 : 1), (16 : -61 : 1), (16 : 60 : 1)]

        Note above that 0 is included since [5]*0 = 0.

        We create a curve of rank 1 with no torsion and do a consistency check::

            sage: E = EllipticCurve('11a').quadratic_twist(-7)
            sage: Q = E([44,-270])
            sage: (4*Q).division_points(4)
            [(44 : -270 : 1)]

        We create a curve over a non-prime finite field with group of
        order `18`::

            sage: k.<a> = GF((5,2))
            sage: E = EllipticCurve(k, [1,2+a,3,4*a,2])
            sage: P = E([3,3*a+4])
            sage: factor(E.order())
            2 * 3^2
            sage: P.order()
            9

        We find the `1`-division points as a consistency check -- there
        is just one, of course::

            sage: P.division_points(1)
            [(3 : 3*a + 4 : 1)]

        The point `P` has order coprime to 2 but divisible by 3, so::

            sage: P.division_points(2)
            [(2*a + 1 : 3*a + 4 : 1), (3*a + 1 : a : 1)]

        We check that each of the 2-division points works as claimed::

            sage: [2*Q for Q in P.division_points(2)]
            [(3 : 3*a + 4 : 1), (3 : 3*a + 4 : 1)]

        Some other checks::

            sage: P.division_points(3)
            []
            sage: P.division_points(4)
            [(0 : 3*a + 2 : 1), (1 : 0 : 1)]
            sage: P.division_points(5)
            [(1 : 1 : 1)]

        An example over a number field (see :trac:`3383`)::

            sage: E = EllipticCurve('19a1')
            sage: K.<t> = NumberField(x^9-3*x^8-4*x^7+16*x^6-3*x^5-21*x^4+5*x^3+7*x^2-7*x+1)
            sage: EK = E.base_extend(K)
            sage: E(0).division_points(3)
            [(0 : 1 : 0), (5 : -10 : 1), (5 : 9 : 1)]
            sage: EK(0).division_points(3)
            [(0 : 1 : 0), (5 : 9 : 1), (5 : -10 : 1)]
            sage: E(0).division_points(9)
            [(0 : 1 : 0), (5 : -10 : 1), (5 : 9 : 1)]
            sage: EK(0).division_points(9)
            [(0 : 1 : 0), (5 : 9 : 1), (5 : -10 : 1), (-150/121*t^8 + 414/121*t^7 + 1481/242*t^6 - 2382/121*t^5 - 103/242*t^4 + 629/22*t^3 - 367/242*t^2 - 1307/121*t + 625/121 : 35/484*t^8 - 133/242*t^7 + 445/242*t^6 - 799/242*t^5 + 373/484*t^4 + 113/22*t^3 - 2355/484*t^2 - 753/242*t + 1165/484 : 1), (-150/121*t^8 + 414/121*t^7 + 1481/242*t^6 - 2382/121*t^5 - 103/242*t^4 + 629/22*t^3 - 367/242*t^2 - 1307/121*t + 625/121 : -35/484*t^8 + 133/242*t^7 - 445/242*t^6 + 799/242*t^5 - 373/484*t^4 - 113/22*t^3 + 2355/484*t^2 + 753/242*t - 1649/484 : 1), (-1383/484*t^8 + 970/121*t^7 + 3159/242*t^6 - 5211/121*t^5 + 37/484*t^4 + 654/11*t^3 - 909/484*t^2 - 4831/242*t + 6791/484 : 927/121*t^8 - 5209/242*t^7 - 8187/242*t^6 + 27975/242*t^5 - 1147/242*t^4 - 1729/11*t^3 + 1566/121*t^2 + 12873/242*t - 10871/242 : 1), (-1383/484*t^8 + 970/121*t^7 + 3159/242*t^6 - 5211/121*t^5 + 37/484*t^4 + 654/11*t^3 - 909/484*t^2 - 4831/242*t + 6791/484 : -927/121*t^8 + 5209/242*t^7 + 8187/242*t^6 - 27975/242*t^5 + 1147/242*t^4 + 1729/11*t^3 - 1566/121*t^2 - 12873/242*t + 10629/242 : 1), (-4793/484*t^8 + 6791/242*t^7 + 10727/242*t^6 - 18301/121*t^5 + 2347/484*t^4 + 2293/11*t^3 - 7311/484*t^2 - 17239/242*t + 26767/484 : 30847/484*t^8 - 21789/121*t^7 - 34605/121*t^6 + 117164/121*t^5 - 10633/484*t^4 - 29437/22*t^3 + 39725/484*t^2 + 55428/121*t - 176909/484 : 1), (-4793/484*t^8 + 6791/242*t^7 + 10727/242*t^6 - 18301/121*t^5 + 2347/484*t^4 + 2293/11*t^3 - 7311/484*t^2 - 17239/242*t + 26767/484 : -30847/484*t^8 + 21789/121*t^7 + 34605/121*t^6 - 117164/121*t^5 + 10633/484*t^4 + 29437/22*t^3 - 39725/484*t^2 - 55428/121*t + 176425/484 : 1)]

        TESTS:

        Check that :trac:`24844` is fixed::

            sage: p = next_prime(1000000)
            sage: E = EllipticCurve(GF(p), 123, 456)
            sage: pts = E(0).division_points(3)
            sage: P = pts[1]; P
            (389063 : 124244 : 1)
            sage: P._order
            3

        When we successfully divide a point known to have infinite
        order, the points returned know that they also have infinite
        order::

            sage: E = EllipticCurve([0,0,1,-1,0])
            sage: P = E(-1,0)
            sage: P.order()
            +Infinity
            sage: pts = P.division_points(3);   len(pts)
            1
            sage: [(Q,Q._order) for Q in pts]
            [((0 : -1 : 1), +Infinity)]

        When we successfully divide a point of known finite order `n`,
        the points returned know that they also have finite order `nk`
        for some divisor `k` of `m`::

          sage: E = EllipticCurve([1, 0, 1, -19, 26])
          sage: [(Q,Q._order) for Q in E(0).division_points(12)]
          [((-5 : 2 : 1), 2),
          ((-2 : -7 : 1), 6),
          ((-2 : 8 : 1), 6),
          ((0 : 1 : 0), 1),
          ((1 : -4 : 1), 6),
          ((1 : 2 : 1), 6),
          ((7/4 : -11/8 : 1), 2),
          ((3 : -2 : 1), 2),
          ((4 : -7 : 1), 3),
          ((4 : 2 : 1), 3),
          ((13 : -52 : 1), 6),
          ((13 : 38 : 1), 6)]
          sage: P = E(4,-7)
          sage: P.order()
          3
          sage: [(Q,Q._order) for Q in P.division_points(4)]
          [((-2 : -7 : 1), 6), ((1 : 2 : 1), 6), ((4 : -7 : 1), 3), ((13 : 38 : 1), 6)]
        """
        # Coerce the input m to an integer
        m = Integer(m)
        # Check for trivial cases of m = 1, -1 and 0.
        if m == 1 or m == -1:
            return [self]
        if m == 0:
            if self == 0:  # then every point Q is a solution, but...
                return [self]
            else:
                return []

        # ans will contain the list of division points.
        ans = []

        # We compute a polynomial g whose roots give all possible x
        # coordinates of the m-division points.  The number of
        # solutions (over the algebraic closure) is m^2, assuming that
        # the characteristic does not divide m.

        E = self.curve()
        P = self
        nP = -P
        P_is_2_torsion = (P == nP)

        # If self is the 0, then self is a solution, and the correct
        # poly is the m'th division polynomial
        if P == 0:
            ans.append(P)
            g = E.division_polynomial(m)
        else:
            # The poly g here is 0 at x(Q) iff x(m*Q) = x(P).
            g = E._multiple_x_numerator(m) - P[0]*E._multiple_x_denominator(m)

            # When 2*P=0, then -Q is a solution iff Q is.  For even m,
            # no 2-torsion point is a solution, so that g is the
            # square of a polynomial g1 of degree m^2/2, and each root
            # of g1 leads to a pair of solutions Q, -Q to m*Q=P.  For
            # odd m, P itself is the only 2-torsion solution, so g has
            # the form (x-x(P))*g1(x)^2 where g1 has degree (m^2-1)/2
            # and each root of g1 leads to a pair Q, -Q.

            if P_is_2_torsion:
                if m % 2 == 0:
                    # This computes g.sqrt() which is not implemented
                    g = g.gcd(g.derivative())*g.leading_coefficient().sqrt()

            # When 2*P!=0, then for each solution Q to m*Q=P, -Q is
            # not a solution (and points of order 2 are not
            # solutions).  Hence the roots of g are distinct and each
            # gives rise to precisely one solution Q.

                else:
                    g0 = g.variables()[0] - P[0]
                    g = g // g0
                    g = g.gcd(g.derivative())*g.leading_coefficient().sqrt()
                    g = g0*g

        if poly_only:
            return g

        for x in g.roots(multiplicities=False):
            if E.is_x_coord(x):
                # Make a point on the curve with this x coordinate.
                Q = E.lift_x(x)
                nQ = -Q
                mQ = m*Q
                # if P==-P then Q works iff -Q works, so we include
                # both unless they are equal:
                if P_is_2_torsion:
                    if mQ == P:
                        ans.append(Q)
                        if nQ != Q:
                            ans.append(nQ)
                else:
                    # P is not 2-torsion so at most one of Q, -Q works
                    # and we must try both:
                    if mQ == P:
                        ans.append(Q)
                    elif mQ == nP:
                        ans.append(nQ)

        if not ans:
            return ans

        # set orders of points found when self's order is known:
        if self.is_zero():
            self._order = Integer(1)
        try:
            n = self._order  # do not compute, just use if already known
            if n == oo:
                for Q in ans:
                    Q._order = oo
            else:
                mfac = m.factor()
                for Q in ans:
                    R = n * Q
                    Q._order = n * generic.order_from_multiple(R, m, factorization=mfac, operation='+')
        except AttributeError:  # do nothing about order if self's order unknown
            pass

        # Finally, sort and return
        ans.sort()
        return ans

    def _divide_out(self, p):
        r"""
        Return `(Q,k)` where `p^kQ` == self and `Q` cannot be divided by `p`.

        .. WARNING::

            It is up to the caller to make sure that this does not loop
            endlessly.  It is used in
            ``EllipticCurve_generic._p_primary_torsion_basis()``, when
            self will always have (finite) order which is a power of `p`,
            so that the order of `Q` increases by a factor of `p` at each
            stage.

            Since it will clearly be in danger of looping when
            self.is_zero(), this case is caught, but otherwise caveat
            user.

        EXAMPLES::

            sage: E = EllipticCurve('37a1')
            sage: P = E([0, 0])
            sage: R = 12*P
            sage: R._divide_out(2)
            ((-1 : -1 : 1), 2)
            sage: R._divide_out(3)
            ((2 : -3 : 1), 1)
            sage: R._divide_out(5)
            ((1357/841 : 28888/24389 : 1), 0)
            sage: R._divide_out(12)
            Traceback (most recent call last):
            ...
            ValueError: p (=12) should be prime.
        """
        p = Integer(p)
        if not p.is_prime():
            raise ValueError("p (=%s) should be prime." % p)

        if self.is_zero():
            raise ValueError("self must not be 0.")

        k = 0
        Q = self
        pts = Q.division_points(p)
        while pts:
            Q = pts[0]
            k += 1
            pts = Q.division_points(p)
        return (Q, k)

    def set_order(self, value):
        r"""
        Set the value of self._order to value.

        Use this when you know a priori the order of this point to avoid a
        potentially expensive order calculation.

        INPUT:

        - ``value`` -- positive integer

        OUTPUT:

        ``None``

        EXAMPLES:

        This example illustrates basic usage.

        ::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 6
            sage: G = E(5, 0)
            sage: G.set_order(2)
            sage: 2*G
            (0 : 1 : 0)

        We now give a more interesting case, the NIST-P521 curve. Its
        order is too big to calculate with Sage, and takes a long time
        using other packages, so it is very useful here.

        ::

            sage: p = 2^521 - 1
            sage: prev_proof_state = proof.arithmetic()
            sage: proof.arithmetic(False) # turn off primality checking
            sage: F = GF(p)
            sage: A = p - 3
            sage: B = 1093849038073734274511112390766805569936207598951683748994586394495953116150735016013708737573759623248592132296706313309438452531591012912142327488478985984
            sage: q = 6864797660130609714981900799081393217269435300143305409394463459185543183397655394245057746333217197532963996371363321113864768612440380340372808892707005449
            sage: E = EllipticCurve([F(A), F(B)])
            sage: G = E.random_point()
            sage: G.set_order(q)
            sage: G.order() * G  # This takes practically no time.
            (0 : 1 : 0)
            sage: proof.arithmetic(prev_proof_state) # restore state

        It is an error to pass a `value` equal to `0`::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 6
            sage: G = E.random_point()
            sage: G.set_order(0)
            Traceback (most recent call last):
            ...
            ValueError: Value 0 illegal for point order
            sage: G.set_order(1000)
            Traceback (most recent call last):
            ...
            ValueError: Value 1000 illegal: outside max Hasse bound

        It is also very likely an error to pass a value which is not the actual
        order of this point. How unlikely is determined by the factorization of
        the actual order, and the actual group structure::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 6
            sage: G = E(5, 0)   # G has order 2
            sage: G.set_order(11)
            Traceback (most recent call last):
            ...
            ValueError: Value 11 illegal: 11 * (5 : 0 : 1) is not the identity

        However, set_order can be fooled, though it's not likely in "real cases
        of interest". For instance, the order can be set to a multiple the
        actual order::

            sage: E = EllipticCurve(GF(7), [0, 1]) # This curve has order 6
            sage: G = E(5, 0)   # G has order 2
            sage: G.set_order(8)
            sage: G.order()
            8

        AUTHORS:

         - Mariah Lenox (2011-02-16)
        """
        value = Integer(value)

        E = self.curve()
        q = E.base_ring().order()
        if value <= 0:
            raise ValueError('Value %s illegal for point order' % value)
        low, hi = Hasse_bounds(q)
        if value > hi:
            raise ValueError('Value %s illegal: outside max Hasse bound' % value)
        if value * self != E(0):
            raise ValueError('Value %s illegal: %s * %s is not the identity' % (value, value, self))
        self._order = value

    # #############################  end  ################################

    def _line_(self, R, Q):
        r"""
        Computes the value at `Q` of a straight line through points
        self and `R`.

        INPUT:

        - ``R, Q`` -- points on self.curve() with ``Q`` nonzero.

        OUTPUT:

        An element of the base field self.curve().base_field().
        A ValueError is raised if ``Q`` is zero.

        EXAMPLES::

            sage: F.<a>=GF((2,5))
            sage: E=EllipticCurve(F,[0,0,1,1,1])
            sage: P = E(a^4 + 1, a^3)
            sage: Q = E(a^4, a^4 + a^3)
            sage: O = E(0)
            sage: P._line_(P,-2*P) == 0
            True
            sage: P._line_(Q,-(P+Q)) == 0
            True
            sage: O._line_(O,Q) == F(1)
            True
            sage: P._line_(O,Q) == a^4 - a^4 + 1
            True
            sage: P._line_(13*P,Q) == a^4
            True
            sage: P._line_(P,Q) == a^4 + a^3 + a^2 + 1
            True

        See :trac:`7116`::

            sage: P._line_ (Q,O)
            Traceback (most recent call last):
            ...
            ValueError: Q must be nonzero.

        .. NOTE::

            This function is used in _miller_ algorithm.

        AUTHORS:

        - David Hansen (2009-01-25)
        """
        if Q.is_zero():
            raise ValueError("Q must be nonzero.")

        if self.is_zero() or R.is_zero():
            if self == R:
                return self.curve().base_field().one()
            if self.is_zero():
                return Q[0] - R[0]
            if R.is_zero():
                return Q[0] - self[0]
        elif self != R:
            if self[0] == R[0]:
                return Q[0] - self[0]
            else:
                l = (R[1] - self[1])/(R[0] - self[0])
                return Q[1] - self[1] - l * (Q[0] - self[0])
        else:
            a1, a2, a3, a4, a6 = self.curve().a_invariants()
            numerator = (3*self[0]**2 + 2*a2*self[0] + a4 - a1*self[1])
            denominator = (2*self[1] + a1*self[0] + a3)
            if denominator == 0:
                return Q[0] - self[0]
            else:
                l = numerator/denominator
                return Q[1] - self[1] - l * (Q[0] - self[0])

    def _miller_(self, Q, n):
        r"""
        Return the value at `Q` of the rational function `f_{n,P}`, where the
        divisor of `f_{n,P}` is `n[P]-[nP]-(n-1)[O]`.

        INPUT:

        - ``Q`` -- a nonzero point on self.curve().

        - ``n`` -- an nonzero integer. If `n<0` then return `Q`
                   evaluated at `1/(v_{nP}*f_{n,P})` (used in the ate pairing).

        OUTPUT:

        An element in the base field self.curve().base_field()
        A ValueError is raised if `Q` is zero.

        EXAMPLES::

            sage: F.<a>=GF((2,5))
            sage: E=EllipticCurve(F,[0,0,1,1,1])
            sage: P = E(a^4 + 1, a^3)
            sage: Fx.<b>=GF((2,(4*5)))
            sage: Ex=EllipticCurve(Fx,[0,0,1,1,1])
            sage: phi=Hom(F,Fx)(F.gen().minpoly().roots(Fx)[0][0])
            sage: Px=Ex(phi(P.xy()[0]),phi(P.xy()[1]))
            sage: Qx = Ex(b^19 + b^18 + b^16 + b^12 + b^10 + b^9 + b^8 + b^5 + b^3 + 1, b^18 + b^13 + b^10 + b^8 + b^5 + b^4 + b^3 + b)
            sage: Px._miller_(Qx,41) == b^17 + b^13 + b^12 + b^9 + b^8 + b^6 + b^4 + 1
            True
            sage: Qx._miller_(Px,41) == b^13 + b^10 + b^8 + b^7 + b^6 + b^5
            True
            sage: P._miller_(E(0),41)
            Traceback (most recent call last):
            ...
            ValueError: Q must be nonzero.

        An example of even order::

            sage: F.<a> = GF((19,4))
            sage: E = EllipticCurve(F,[-1,0])
            sage: P = E(15*a^3 + 17*a^2 + 14*a + 13,16*a^3 + 7*a^2 + a + 18)
            sage: Q = E(10*a^3 + 16*a^2 + 4*a + 2, 6*a^3 + 4*a^2 + 3*a + 2)
            sage: x=P.weil_pairing(Q,360)
            sage: x^360 == F(1)
            True

        You can use the _miller_ function on linearly dependent
        points, but with the risk of a dividing with zero::

            sage: Px._miller_(2*Px,41)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero in finite field

        A small example of embedding degree 6::

            sage: q = 401; F = GF(q); a = 146; b = 400; k = 6
            sage: E = EllipticCurve([F(a), F(b)])
            sage: R.<x> = F[]; K.<a> = GF(q^k, modulus=x^6 + 4*x^4 + 115*x^3 + 81*x^2 + 51*x + 3)
            sage: EK = E.base_extend(K)
            sage: P = E([F(338), F(227)])
            sage: Q_x = 333*a^5 + 391*a^4 + 160*a^3 + 335*a^2 + 71*a + 93
            sage: Q_y = 343*a^5 + 273*a^4 + 26*a^3 + 342*a^2 + 340*a + 210
            sage: Q = EK([Q_x, Q_y])
            sage: P._miller_(Q, 127)
            371*a^5 + 39*a^4 + 355*a^3 + 233*a^2 + 20*a + 275

        A series of small examples and small torsions.  We start with
        `n=1`, which is trivial: the function is the constant
        1::

            sage: E = EllipticCurve([GF(7)(0), 2])
            sage: P = E(5, 1); Q = E(0, 3); I = E(0)
            sage: I._miller_(P, 1)
            1
            sage: I._miller_(Q, 1)
            1

        A two-torsion example. In this case `f_{n,P}(Q) = x_Q - x_P`::

            sage: E = EllipticCurve([GF(7)(-1), 0])
            sage: P = E(0,0); Q = E(1, 0);
            sage: P._miller_(P, 2)
            0
            sage: Q._miller_(Q, 2)
            0
            sage: P._miller_(Q, 2)
            1
            sage: Q._miller_(P, 2)
            6

        A three-torsion example::

            sage: E = EllipticCurve([GF(7)(0), 2])
            sage: P = E(5, 1); Q = E(0, 3);
            sage: P._miller_(Q, 3)
            4

        A 4-torsion example::

            sage: E = EllipticCurve([GF(7)(-1), 0])
            sage: P = E(5, 1); Q = E(4, 2)
            sage: P._miller_(Q, 4)
            3

        A 5-torsion example::

            sage: E = EllipticCurve([GF(7)(-1), 4])
            sage: P = E(4, 1); Q = E(6, 5)
            sage: P._miller_(Q, 5)
            1

        A 6-torsion example::

            sage: E = EllipticCurve([GF(7)(3), 1])
            sage: P = E(5, 1); Q = E(3, 3)
            sage: P._miller_(Q, 6)
            5

        An example which is part of an ate pairing calculation. The trace of
        the curve is negative, so it should exercise the `n<0` case in the
        code::

            sage: p = 2017; A = 1; B = 30; r = 29; t = -70; k = 7;
            sage: F = GF(p); R.<x> = F[]
            sage: E = EllipticCurve([F(A), F(B)]); P = E(369, 716)
            sage: K.<a> = GF(p^k, modulus=x^k+2); EK = E.base_extend(K)
            sage: Qx = 1226*a^6 + 1778*a^5 + 660*a^4 + 1791*a^3 + 1750*a^2 + 867*a + 770
            sage: Qy = 1764*a^6 + 198*a^5 + 1206*a^4 + 406*a^3 + 1200*a^2 + 273*a + 1712
            sage: Q = EK(Qx, Qy)
            sage: Q._miller_(P, t-1)
            1311*a^6 + 1362*a^5 + 1177*a^4 + 807*a^3 + 1331*a^2 + 1530*a + 1931

        ALGORITHM:

            Double-and-add. See also [Mil2004]_.

        AUTHORS:

        - David Hansen (2009-01-25)
        - Mariah Lenox (2011-03-07) -- Added more doctests and support for
          negative n.
        """
        if Q.is_zero():
            raise ValueError("Q must be nonzero.")
        if n.is_zero():
            raise ValueError("n must be nonzero.")
        n_is_negative = False
        if n < 0:
            n = n.abs()
            n_is_negative = True

        one = self.curve().base_field().one()
        t = one
        V = self
        nbin = n.bits()
        i = n.nbits() - 2
        while i > -1:
            S = 2*V
            ell = V._line_(V, Q)
            vee = S._line_(-S, Q)
            t = (t**2)*(ell/vee)
            V = S
            if nbin[i] == 1:
                S = V+self
                ell = V._line_(self, Q)
                vee = S._line_(-S, Q)
                t = t*(ell/vee)
                V = S
            i = i-1
        if n_is_negative:
            vee = V._line_(-V, Q)
            t = 1/(t*vee)
        return t

    def weil_pairing(self, Q, n, algorithm=None):
        r"""
        Compute the Weil pairing of this point with another point `Q`
        on the same curve.

        INPUT:

        - ``Q`` -- another point on the same curve as ``self``.

        - ``n`` -- an integer `n` such that `nP = nQ = (0:1:0)`, where
          `P` is ``self``.

        - ``algorithm`` (default: ``None``) -- choices are ``pari``
          and ``sage``. PARI is usually significantly faster, but it
          only works over finite fields. When ``None`` is given, a
          suitable algorithm is chosen automatically.

        OUTPUT:

        An `n`'th root of unity in the base field of the curve.

        EXAMPLES::

            sage: F.<a>=GF((2,5))
            sage: E=EllipticCurve(F,[0,0,1,1,1])
            sage: P = E(a^4 + 1, a^3)
            sage: Fx.<b>=GF((2,4*5))
            sage: Ex=EllipticCurve(Fx,[0,0,1,1,1])
            sage: phi=Hom(F,Fx)(F.gen().minpoly().roots(Fx)[0][0])
            sage: Px=Ex(phi(P.xy()[0]),phi(P.xy()[1]))
            sage: O = Ex(0)
            sage: Qx = Ex(b^19 + b^18 + b^16 + b^12 + b^10 + b^9 + b^8 + b^5 + b^3 + 1, b^18 + b^13 + b^10 + b^8 + b^5 + b^4 + b^3 + b)
            sage: Px.weil_pairing(Qx,41) == b^19 + b^15 + b^9 + b^8 + b^6 + b^4 + b^3 + b^2 + 1
            True
            sage: Px.weil_pairing(17*Px,41) == Fx(1)
            True
            sage: Px.weil_pairing(O,41) == Fx(1)
            True

        An error is raised if either point is not `n`-torsion::

            sage: Px.weil_pairing(O,40)
            Traceback (most recent call last):
            ...
            ValueError: points must both be n-torsion

        A larger example (see :trac:`4964`)::

            sage: P,Q = EllipticCurve(GF((19,4),'a'),[-1,0]).gens()
            sage: P.order(), Q.order()
            (360, 360)
            sage: z = P.weil_pairing(Q,360)
            sage: z.multiplicative_order()
            360

        An example over a number field::

            sage: P,Q = EllipticCurve('11a1').change_ring(CyclotomicField(5)).torsion_subgroup().gens()
            sage: P,Q = (P.element(), Q.element())
            sage: (P.order(),Q.order())
            (5, 5)
            sage: P.weil_pairing(Q,5)
            zeta5^2
            sage: Q.weil_pairing(P,5)
            zeta5^3

        TESTS:

        Check that the original Sage implementation still works::

            sage: GF(65537^2).inject_variables()
            Defining z2
            sage: E = EllipticCurve(GF(65537^2), [0,1])
            sage: P = E(22, 28891)
            sage: Q = E(-93, 40438*z2 + 31573)
            sage: P.weil_pairing(Q, 7282, algorithm='sage')
            19937*z2 + 65384

        Passing an unknown ``algorithm=`` argument should fail::

            sage: P.weil_pairing(Q, 7282, algorithm='_invalid_')
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm

        ALGORITHM:

        - For ``algorithm='pari'``: :pari:`ellweilpairing`.

        - For ``algorithm='sage'``:
          Implemented using Proposition 8 in [Mil2004]_.  The value 1 is
          returned for linearly dependent input points.  This condition
          is caught via a DivisionByZeroError, since the use of a
          discrete logarithm test for linear dependence is much too slow
          for large `n`.

        AUTHORS:

        - David Hansen (2009-01-25)
        - Lorenz Panny (2022): ``algorithm='pari'``
        """
        P = self
        E = P.curve()

        if Q.curve() is not E:
            raise ValueError("points must both be on the same curve")

        if algorithm is None:
            if E.base_field().is_finite():
                algorithm = 'pari'
            else:
                algorithm = 'sage'

        if algorithm == 'pari':
            if pari.ellmul(E,P,n) != [0] or pari.ellmul(E,Q,n) != [0]:
                raise ValueError("points must both be n-torsion")
            return E.base_field()(pari.ellweilpairing(E, P, Q, n))

        if algorithm != 'sage':
            raise ValueError('unknown algorithm')

        # Test if P, Q are both in E[n]
        if n*P or n*Q:
            raise ValueError("points must both be n-torsion")

        one = E.base_field().one()

        # Case where P = Q
        if P == Q:
            return one

        # Case where P = O or Q = O
        if P.is_zero() or Q.is_zero():
            return one

        # The non-trivial case P != Q

        # Note (2010-12-29): The following code block should not be
        # used.  It attempts to reduce the pairing calculation to order
        # d = gcd(|P|,|Q|) instead of order n, but this approach is
        # counterproductive, since calculating |P| and |Q| is much
        # slower than calculating the pairing [BGN05].
        #
        # [BGN05] D. Boneh, E. Goh, and K. Nissim, "Evaluating 2-DNF
        # Formulas on Ciphertexts", TCC 2005, LNCS 3378, pp. 325-341.
        #
        # Reduction to order d = gcd(|P|,|Q|); value is a d'th root of unity
        # try:
        #     nP = P.order()
        # except AttributeError:
        #     nP = generic.order_from_multiple(P,n,operation='+')
        # try:
        #     nQ = Q.order()
        # except AttributeError:
        #     nQ = generic.order_from_multiple(Q,n,operation='+')
        # d = arith.gcd(nP,nQ)
        # if d==1:
        #     return one
        #
        # P = (nP//d)*P # order d
        # Q = (nQ//d)*Q # order d
        # n = d

        try:
            return ((-1)**n.test_bit(0))*(P._miller_(Q, n)/Q._miller_(P, n))
        except ZeroDivisionError:
            return one

    def tate_pairing(self, Q, n, k, q=None):
        r"""
        Return Tate pairing of `n`-torsion point `P = self` and point `Q`.

        The value returned is `f_{n,P}(Q)^e` where `f_{n,P}` is a function with
        divisor `n[P]-n[O].`. This is also known as the "modified Tate
        pairing". It is a well-defined bilinear map.

        INPUT:

        - ``P=self`` -- Elliptic curve point having order n

        - ``Q`` -- Elliptic curve point on same curve as P (can be any order)

        - ``n`` -- positive integer: order of P

        - ``k`` -- positive integer: embedding degree

        - ``q`` -- positive integer: size of base field (the "big"
          field is `GF(q^k)`. `q` needs to be set only if its value
          cannot be deduced.)

        OUTPUT:

        An `n`'th root of unity in the base field self.curve().base_field()

        EXAMPLES:

        A simple example, pairing a point with itself, and pairing a point with
        another rational point::

            sage: p = 103; A = 1; B = 18; E = EllipticCurve(GF(p), [A, B])
            sage: P = E(33, 91); n = P.order(); n
            19
            sage: k = GF(n)(p).multiplicative_order(); k
            6
            sage: P.tate_pairing(P, n, k)
            1
            sage: Q = E(87, 51)
            sage: P.tate_pairing(Q, n, k)
            1
            sage: set_random_seed(35)
            sage: P.tate_pairing(P,n,k)
            1

        We now let Q be a point on the same curve as above, but defined over
        the pairing extension field, and we also demonstrate the bilinearity of
        the pairing::

            sage: K.<a> = GF((p,k))
            sage: EK = E.base_extend(K); P = EK(P)
            sage: Qx = 69*a^5 + 96*a^4 + 22*a^3 + 86*a^2 + 6*a + 35
            sage: Qy = 34*a^5 + 24*a^4 + 16*a^3 + 41*a^2 + 4*a + 40
            sage: Q = EK(Qx, Qy);

        Multiply by cofactor so Q has order n::

            sage: h = 551269674; Q = h*Q
            sage: P = EK(P); P.tate_pairing(Q, n, k)
            24*a^5 + 34*a^4 + 3*a^3 + 69*a^2 + 86*a + 45
            sage: s = Integer(randrange(1,n))
            sage: ans1 = (s*P).tate_pairing(Q, n, k)
            sage: ans2 = P.tate_pairing(s*Q, n, k)
            sage: ans3 = P.tate_pairing(Q, n, k)^s
            sage: ans1 == ans2 == ans3
            True
            sage: (ans1 != 1) and (ans1^n == 1)
            True

        Here is an example of using the Tate pairing to compute the Weil
        pairing (using the same data as above)::

            sage: e = Integer((p^k-1)/n); e
            62844857712
            sage: P.weil_pairing(Q, n)^e
            94*a^5 + 99*a^4 + 29*a^3 + 45*a^2 + 57*a + 34
            sage: P.tate_pairing(Q, n, k) == P._miller_(Q, n)^e
            True
            sage: Q.tate_pairing(P, n, k) == Q._miller_(P, n)^e
            True
            sage: P.tate_pairing(Q, n, k)/Q.tate_pairing(P, n, k)
            94*a^5 + 99*a^4 + 29*a^3 + 45*a^2 + 57*a + 34

        An example where we have to pass the base field size (and we again have
        agreement with the Weil pairing)::

            sage: F.<a>=GF((2,5))
            sage: E=EllipticCurve(F,[0,0,1,1,1])
            sage: P = E(a^4 + 1, a^3)
            sage: Fx.<b>=GF((2,4*5))
            sage: Ex=EllipticCurve(Fx,[0,0,1,1,1])
            sage: phi=Hom(F,Fx)(F.gen().minpoly().roots(Fx)[0][0])
            sage: Px=Ex(phi(P.xy()[0]),phi(P.xy()[1]))
            sage: Qx = Ex(b^19+b^18+b^16+b^12+b^10+b^9+b^8+b^5+b^3+1, b^18+b^13+b^10+b^8+b^5+b^4+b^3+b)
            sage: Px.tate_pairing(Qx, n=41, k=4)
            Traceback (most recent call last):
            ...
            ValueError: Unexpected field degree: set keyword argument q equal to the size of the base field (big field is GF(q^4)).
            sage: num = Px.tate_pairing(Qx, n=41, k=4, q=32); num
            b^19 + b^14 + b^13 + b^12 + b^6 + b^4 + b^3
            sage: den = Qx.tate_pairing(Px, n=41, k=4, q=32); den
            b^19 + b^17 + b^16 + b^15 + b^14 + b^10 + b^6 + b^2 + 1
            sage: e = Integer((32^4-1)/41); e
            25575
            sage: Px.weil_pairing(Qx, 41)^e == num/den
            True

        .. NOTE::

            This function uses Miller's algorithm, followed by a naive
            exponentiation. It does not do anything fancy. In the case
            that there is an issue with `Q` being on one of the lines
            generated in the `r*P` calculation, `Q` is offset by a random
            point `R` and P.tate_pairing(Q+R,n,k)/P.tate_pairing(R,n,k)
            is returned.

        AUTHORS:

        - Mariah Lenox (2011-03-07)
        """
        P = self
        E = P.curve()

        if not Q.curve() is E:
            raise ValueError("Points must both be on the same curve")

        K = E.base_ring()
        d = K.degree()
        if q is None:
            if d == 1:
                q = K.order()
            elif d == k:
                q = K.base_ring().order()
            else:
                raise ValueError("Unexpected field degree: set keyword argument q equal to the size of the base field (big field is GF(q^%s))." % k)

        if n*P != E(0):
            raise ValueError('This point is not of order n=%s' % n)

        # In small cases, or in the case of pairing an element with
        # itself, Q could be on one of the lines in the Miller
        # algorithm. If this happens we try again, with an offset of a
        # random point.
        try:
            ret = self._miller_(Q, n)
            e = Integer((q**k - 1)/n)
            ret = ret**e
        except (ZeroDivisionError, ValueError):
            R = E.random_point()
            ret = self.tate_pairing(Q + R, n, k)/self.tate_pairing(R, n, k)
        return ret

    def ate_pairing(self, Q, n, k, t, q=None):
        r"""
        Return ate pairing of `n`-torsion points `P=self` and `Q`.

        Also known as the `n`-th modified ate pairing. `P` is `GF(q)`-rational,
        and `Q` must be an element of `Ker(\pi-p)`, where `\pi` is the
        `q`-frobenius map (and hence `Q` is `GF(q^k)`-rational).

        INPUT:

        - ``P=self`` -- a point of order `n`, in `ker(\pi-1)`, where
          `\pi` is the `q`-Frobenius map (e.g., `P` is `q-rational`).

        - ``Q`` -- a point of order `n` in `ker(\pi-q)`

        - ``n`` -- the order of `P` and `Q`.

        - ``k`` -- the embedding degree.

        - ``t`` -- the trace of Frobenius of the curve over `GF(q)`.

        - ``q`` -- (default: None) the size of base field (the "big"
          field is `GF(q^k)`). `q` needs to be set only if its value
          cannot be deduced.

        OUTPUT:

        FiniteFieldElement in `GF(q^k)` -- the ate pairing of `P` and `Q`.

        EXAMPLES:

        An example with embedding degree 6::

            sage: p = 7549; A = 0; B = 1; n = 157; k = 6; t = 14
            sage: F = GF(p); E = EllipticCurve(F, [A, B])
            sage: R.<x> = F[]; K.<a> = GF((p,k), modulus=x^k+2)
            sage: EK = E.base_extend(K)
            sage: P = EK(3050, 5371); Q = EK(6908*a^4, 3231*a^3)
            sage: P.ate_pairing(Q, n, k, t)
            6708*a^5 + 4230*a^4 + 4350*a^3 + 2064*a^2 + 4022*a + 6733
            sage: s = Integer(randrange(1, n))
            sage: (s*P).ate_pairing(Q, n, k, t) == P.ate_pairing(s*Q, n, k, t)
            True
            sage: P.ate_pairing(s*Q, n, k, t) == P.ate_pairing(Q, n, k, t)^s
            True

        Another example with embedding degree 7 and positive trace::

            sage: p = 2213; A = 1; B = 49; n = 1093; k = 7; t = 28
            sage: F = GF(p); E = EllipticCurve(F, [A, B])
            sage: R.<x> = F[]; K.<a> = GF((p,k), modulus=x^k+2)
            sage: EK = E.base_extend(K)
            sage: P = EK(1583, 1734)
            sage: Qx = 1729*a^6+1767*a^5+245*a^4+980*a^3+1592*a^2+1883*a+722
            sage: Qy = 1299*a^6+1877*a^5+1030*a^4+1513*a^3+1457*a^2+309*a+1636
            sage: Q = EK(Qx, Qy)
            sage: P.ate_pairing(Q, n, k, t)
            1665*a^6 + 1538*a^5 + 1979*a^4 + 239*a^3 + 2134*a^2 + 2151*a + 654
            sage: s = Integer(randrange(1, n))
            sage: (s*P).ate_pairing(Q, n, k, t) == P.ate_pairing(s*Q, n, k, t)
            True
            sage: P.ate_pairing(s*Q, n, k, t) == P.ate_pairing(Q, n, k, t)^s
            True

        Another example with embedding degree 7 and negative trace::

            sage: p = 2017; A = 1; B = 30; n = 29; k = 7; t = -70
            sage: F = GF(p); E = EllipticCurve(F, [A, B])
            sage: R.<x> = F[]; K.<a> = GF((p,k), modulus=x^k+2)
            sage: EK = E.base_extend(K)
            sage: P = EK(369, 716)
            sage: Qx = 1226*a^6+1778*a^5+660*a^4+1791*a^3+1750*a^2+867*a+770
            sage: Qy = 1764*a^6+198*a^5+1206*a^4+406*a^3+1200*a^2+273*a+1712
            sage: Q = EK(Qx, Qy)
            sage: P.ate_pairing(Q, n, k, t)
            1794*a^6 + 1161*a^5 + 576*a^4 + 488*a^3 + 1950*a^2 + 1905*a + 1315
            sage: s = Integer(randrange(1, n))
            sage: (s*P).ate_pairing(Q, n, k, t) == P.ate_pairing(s*Q, n, k, t)
            True
            sage: P.ate_pairing(s*Q, n, k, t) == P.ate_pairing(Q, n, k, t)^s
            True

        Using the same data, we show that the ate pairing is a power of the
        Tate pairing (see [HSV2006]_ end of section 3.1)::

            sage: c = (k*p^(k-1)).mod(n); T = t - 1
            sage: N = gcd(T^k - 1, p^k - 1)
            sage: s = Integer(N/n)
            sage: L = Integer((T^k - 1)/N)
            sage: M = (L*s*c.inverse_mod(n)).mod(n)
            sage: P.ate_pairing(Q, n, k, t) == Q.tate_pairing(P, n, k)^M
            True

        An example where we have to pass the base field size (and we again have
        agreement with the Tate pairing). Note that though `Px` is not
        `F`-rational, (it is the homomorphic image of an `F`-rational point) it
        is nonetheless in `ker(\pi-1)`, and so is a legitimate input::

            sage: q = 2^5; F.<a>=GF(q)
            sage: n = 41; k = 4; t = -8
            sage: E=EllipticCurve(F,[0,0,1,1,1])
            sage: P = E(a^4 + 1, a^3)
            sage: Fx.<b>=GF(q^k)
            sage: Ex=EllipticCurve(Fx,[0,0,1,1,1])
            sage: phi=Hom(F,Fx)(F.gen().minpoly().roots(Fx)[0][0])
            sage: Px=Ex(phi(P.xy()[0]),phi(P.xy()[1]))
            sage: Qx = Ex(b^19+b^18+b^16+b^12+b^10+b^9+b^8+b^5+b^3+1, b^18+b^13+b^10+b^8+b^5+b^4+b^3+b)
            sage: Qx = Ex(Qx[0]^q, Qx[1]^q) - Qx  # ensure Qx is in ker(pi - q)
            sage: Px.ate_pairing(Qx, n, k, t)
            Traceback (most recent call last):
            ...
            ValueError: Unexpected field degree: set keyword argument q equal to the size of the base field (big field is GF(q^4)).
            sage: Px.ate_pairing(Qx, n, k, t, q)
            b^19 + b^18 + b^17 + b^16 + b^15 + b^14 + b^13 + b^12 + b^11 + b^9 + b^8 + b^5 + b^4 + b^2 + b + 1
            sage: s = Integer(randrange(1, n))
            sage: (s*Px).ate_pairing(Qx, n, k, t, q) == Px.ate_pairing(s*Qx, n, k, t, q)
            True
            sage: Px.ate_pairing(s*Qx, n, k, t, q) == Px.ate_pairing(Qx, n, k, t, q)^s
            True
            sage: c = (k*q^(k-1)).mod(n); T = t - 1
            sage: N = gcd(T^k - 1, q^k - 1)
            sage: s = Integer(N/n)
            sage: L = Integer((T^k - 1)/N)
            sage: M = (L*s*c.inverse_mod(n)).mod(n)
            sage: Px.ate_pairing(Qx, n, k, t, q) == Qx.tate_pairing(Px, n, k, q)^M
            True

        It is an error if `Q` is not in the kernel of `\pi - p`, where `\pi` is
        the Frobenius automorphism::

            sage: p = 29; A = 1; B = 0; n = 5; k = 2; t = 10
            sage: F = GF(p); R.<x> = F[]
            sage: E = EllipticCurve(F, [A, B]);
            sage: K.<a> = GF((p,k), modulus=x^k+2); EK = E.base_extend(K)
            sage: P = EK(13, 8); Q = EK(13, 21)
            sage: P.ate_pairing(Q, n, k, t)
            Traceback (most recent call last):
            ...
            ValueError: Point (13 : 21 : 1) not in Ker(pi - q)

        It is also an error if `P` is not in the kernel os `\pi - 1`::

            sage: p = 29; A = 1; B = 0; n = 5; k = 2; t = 10
            sage: F = GF(p); R.<x> = F[]
            sage: E = EllipticCurve(F, [A, B]);
            sage: K.<a> = GF((p,k), modulus=x^k+2); EK = E.base_extend(K)
            sage: P = EK(14, 10*a); Q = EK(13, 21)
            sage: P.ate_pairing(Q, n, k, t)
            Traceback (most recent call last):
            ...
            ValueError: This point (14 : 10*a : 1) is not in Ker(pi - 1)

        .. NOTE::

            First defined in the paper of [HSV2006]_, the ate pairing
            can be computationally effective in those cases when the
            trace of the curve over the base field is significantly
            smaller than the expected value. This implementation is
            simply Miller's algorithm followed by a naive
            exponentiation, and makes no claims towards efficiency.

        AUTHORS:

        - Mariah Lenox (2011-03-08)
        """
        P = self
        # check for same curve
        E = P.curve()
        O = E(0)
        if not Q.curve() is E:
            raise ValueError("Points must both be on the same curve")

        # set q to be the order of the base field
        K = E.base_ring()
        if q is None:
            d = K.degree()
            if d == k:
                q = K.base_ring().order()
            else:
                raise ValueError("Unexpected field degree: set keyword argument q equal to the size of the base field (big field is GF(q^%s))." % k)

        # check order of P
        if n*P != O:
            raise ValueError('This point %s is not of order n=%s' % (P, n))

        # check for P in kernel pi - 1:
        piP = E(P[0]**q, P[1]**q)
        if piP - P != O:
            raise ValueError('This point %s is not in Ker(pi - 1)' % P)

        # check for Q in kernel pi - q:
        piQ = E(Q[0]**q, Q[1]**q)
        if piQ - q*Q != O:
            raise ValueError('Point %s not in Ker(pi - q)' % Q)

        T = t-1
        ret = Q._miller_(P, T)
        e = Integer((q**k - 1)/n)
        ret = ret**e
        return ret


class EllipticCurvePoint_number_field(EllipticCurvePoint_field):
    """
    A point on an elliptic curve over a number field.

    Most of the functionality is derived from the parent class
    ``EllipticCurvePoint_field``.  In addition we have support for
    orders, heights, reduction modulo primes, and elliptic logarithms.

    EXAMPLES::

        sage: E = EllipticCurve('37a')
        sage: E([0,0])
        (0 : 0 : 1)
        sage: E(0,0)               # brackets are optional
        (0 : 0 : 1)
        sage: E([GF(5)(0), 0])     # entries are coerced
        (0 : 0 : 1)

        sage: E(0.000, 0)
        (0 : 0 : 1)

        sage: E(1,0,0)
        Traceback (most recent call last):
        ...
        TypeError: Coordinates [1, 0, 0] do not define a point on
        Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    ::

        sage: E = EllipticCurve([0,0,1,-1,0])
        sage: S = E(QQ); S
        Abelian group of points on Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    TESTS::

        sage: loads(S.dumps()) == S
        True
        sage: P = E(0,0); P
        (0 : 0 : 1)
        sage: loads(P.dumps()) == P
        True
        sage: T = 100*P
        sage: loads(T.dumps()) == T
        True

    Test pickling an elliptic curve that has known points on it::

        sage: e = EllipticCurve([0, 0, 1, -1, 0]); g = e.gens(); loads(dumps(e)) == e
        True
    """

    def order(self):
        r"""
        Return the order of this point on the elliptic curve.

        If the point has infinite order, returns +Infinity.  For
        curves defined over `\QQ`, we call PARI; over other
        number fields we implement the function here.

        .. NOTE::

            :meth:`additive_order` is a synonym for :meth:`order`

        EXAMPLES::

            sage: E = EllipticCurve([0,0,1,-1,0])
            sage: P = E([0,0]); P
            (0 : 0 : 1)
            sage: P.order()
            +Infinity

        ::

            sage: E = EllipticCurve([0,1])
            sage: P = E([-1,0])
            sage: P.order()
            2
            sage: P.additive_order()
            2
        """
        try:
            return self._order
        except AttributeError:
            pass

        if self.is_zero():
            self._order = Integer(1)
            return self._order

        E = self.curve()

        # First try PARI
        try:
            n = E.pari_curve().ellorder(self)
            if n:
                n = Integer(n)
            else:
                n = oo
            self._order = n
            return n
        except PariError:
            pass

        # Get the torsion order if known, else a bound on (multiple
        # of) the order.  We do not compute the torsion if it is not
        # already known, since computing the bound is faster (and is
        # also cached).

        try:
            N = E._torsion_order
        except AttributeError:
            N = E._torsion_bound()

        # Now self is a torsion point iff it is killed by N:
        if not (N*self).is_zero():
            self._order = oo
            return self._order

        # Finally we find the exact order using the generic code:
        self._order = generic.order_from_multiple(self, N, operation='+')
        return self._order

    additive_order = order

    def has_finite_order(self):
        """
        Return ``True`` iff this point has finite order on the elliptic curve.

        EXAMPLES::

            sage: E = EllipticCurve([0,0,1,-1,0])
            sage: P = E([0,0]); P
            (0 : 0 : 1)
            sage: P.has_finite_order()
            False

        ::

            sage: E = EllipticCurve([0,1])
            sage: P = E([-1,0])
            sage: P.has_finite_order()
            True
        """
        if self.is_zero():
            return True
        return self.order() != oo

    def has_infinite_order(self):
        r"""
        Return True iff this point has infinite order on the elliptic curve.

        EXAMPLES::

            sage: E = EllipticCurve([0,0,1,-1,0])
            sage: P = E([0,0]); P
            (0 : 0 : 1)
            sage: P.has_infinite_order()
            True

        ::

            sage: E = EllipticCurve([0,1])
            sage: P = E([-1,0])
            sage: P.has_infinite_order()
            False
        """
        if self.is_zero():
            return False
        return self.order() == oo

    def is_on_identity_component(self, embedding=None):
        r"""
        Returns True iff this point is on the identity component of
        its curve with respect to a given (real or complex) embedding.

        INPUT:

        - ``self`` -- a point on a curve over any ordered field (e.g. `\QQ`)

        - ``embedding`` -- an embedding from the base_field of the
          point's curve into `\RR` or `\CC`; if ``None`` (the
          default) it uses the first embedding of the base_field into
          `\RR` if any, else the first embedding into `\CC`.

        OUTPUT:

        (bool) -- ``True`` iff the point is on the identity component of
        the curve.  (If the point is zero then the result is True.)

        EXAMPLES:

        For `K=\QQ` there is no need to specify an embedding::

            sage: E=EllipticCurve('5077a1')
            sage: [E.lift_x(x).is_on_identity_component() for x in srange(-3,5)]
            [False, False, False, False, False, True, True, True]

        An example over a field with two real embeddings::

            sage: L.<a> = QuadraticField(2)
            sage: E=EllipticCurve(L,[0,1,0,a,a])
            sage: P=E(-1,0)
            sage: [P.is_on_identity_component(e) for e in L.embeddings(RR)]
            [False, True]

        We can check this as follows::

            sage: [e(E.discriminant())>0 for e in L.embeddings(RR)]
            [True, False]
            sage: e = L.embeddings(RR)[0]
            sage: E1 = EllipticCurve(RR,[e(ai) for ai in E.ainvs()])
            sage: e1,e2,e3 = E1.two_division_polynomial().roots(RR,multiplicities=False)
            sage: e1 < e2 < e3 and e(P[0]) < e3
            True
        """
        if self.is_zero():       # trivial case
            return True

        e = embedding
        # It is also trivially true if we have a complex embedding
        if e is not None:
            if not isinstance(e.codomain(), sage.rings.abc.RealField):
                return True

        # find a suitable embedding if none was supplied:
        E = self.curve()
        K = E.base_field()
        if e is None:
            try:
                e = K.embeddings(rings.RealField())[0]
            except IndexError:
                e = K.embeddings(rings.ComplexField())[0]

        # If there is only one component, the result is True:
        if not isinstance(e.codomain(), sage.rings.abc.RealField):  # complex embedding
            return True
        if e(E.discriminant()) < 0:  # only one component
            return True

        # Now we have a real embedding and two components and have to work:
        gx = E.two_division_polynomial()
        gxd = gx.derivative()
        gxdd = gxd.derivative()
        return (e(gxd(self[0])) > 0 and e(gxdd(self[0])) > 0)

    def has_good_reduction(self, P=None):
        r"""
        Returns True iff this point has good reduction modulo a prime.

        INPUT:

        - ``P`` -- a prime of the base_field of the point's curve, or
          ``None`` (default)

        OUTPUT:

        (bool) If a prime `P` of the base field is specified, returns
        True iff the point has good reduction at `P`; otherwise,
        return true if the point has god reduction at all primes in
        the support of the discriminant of this model.

        EXAMPLES::

            sage: E = EllipticCurve('990e1')
            sage: P = E.gen(0); P
            (15 : 51 : 1)
            sage: [E.has_good_reduction(p) for p in [2,3,5,7]]
            [False, False, False, True]
            sage: [P.has_good_reduction(p) for p in [2,3,5,7]]
            [True, False, True, True]
            sage: [E.tamagawa_exponent(p) for p in [2,3,5,7]]
            [2, 2, 1, 1]
            sage: [(2*P).has_good_reduction(p) for p in [2,3,5,7]]
            [True, True, True, True]
            sage: P.has_good_reduction()
            False
            sage: (2*P).has_good_reduction()
            True
            sage: (3*P).has_good_reduction()
            False

        ::

            sage: K.<i> = NumberField(x^2+1)
            sage: E = EllipticCurve(K,[0,1,0,-160,308])
            sage: P = E(26,-120)
            sage: E.discriminant().support()
            [Fractional ideal (i + 1),
            Fractional ideal (-i - 2),
            Fractional ideal (2*i + 1),
            Fractional ideal (3)]
            sage: [E.tamagawa_exponent(p) for p in E.discriminant().support()]
            [1, 4, 4, 4]
            sage: P.has_good_reduction()
            False
            sage: (2*P).has_good_reduction()
            False
            sage: (4*P).has_good_reduction()
            True

        TESTS:

        An example showing that :trac:`8498` is fixed::

            sage: E = EllipticCurve('11a1')
            sage: K.<t> = NumberField(x^2+47)
            sage: EK = E.base_extend(K)
            sage: T = EK(5,5)
            sage: P = EK(-2, -1/2*t - 1/2)
            sage: p = K.ideal(11)
            sage: T.has_good_reduction(p)
            False
            sage: P.has_good_reduction(p)
            True
        """
        if self.is_zero():       # trivial case
            return True

        E = self.curve()
        if P is None:
            return all(self.has_good_reduction(Pr)
                       for Pr in E.discriminant().support())
        K = E.base_field()
        from sage.schemes.elliptic_curves.ell_local_data import check_prime
        P = check_prime(K, P)

        # If the curve has good reduction at P, the result is True:
        t = E.local_data(P).bad_reduction_type()
        if t is None:
            return True

        # Make sure the curve is integral and locally minimal at P:
        Emin = E.local_minimal_model(P)
        urst = E.isomorphism_to(Emin)
        Q = urst(self)

        # Scale the homogeneous coordinates of the point to be primitive:
        xyz = list(Q)
        e = min([c.valuation(P) for c in xyz])
        if e != 0:
            if K is rings.QQ:
                pi = P
            else:
                pi = K.uniformizer(P)
            pie = pi**e
            xyz = [c/pie for c in xyz]

        # Evaluate the partial derivatives at the point to see if they
        # are zero mod P

        # See #8498: sometimes evaluating F's derivatives at xyz
        # returns a constant polynomial instead of a constant

        F = Emin.defining_polynomial()
        for v in F.variables():
            c = (F.derivative(v))(xyz)
            try:
                val = c.valuation(P)
            except AttributeError:
                val = c.constant_coefficient().valuation(P)
            if val == 0:
                return True
        return False

    def reduction(self, p):
        """
        This finds the reduction of a point `P` on the elliptic curve
        modulo the prime `p`.

        INPUT:

        - ``self`` -- A point on an elliptic curve.

        - ``p`` -- a prime number

        OUTPUT:

        The point reduced to be a point on the elliptic curve modulo `p`.

        EXAMPLES::

            sage: E = EllipticCurve([1,2,3,4,0])
            sage: P = E(0,0)
            sage: P.reduction(5)
            (0 : 0 : 1)
            sage: Q = E(98,931)
            sage: Q.reduction(5)
            (3 : 1 : 1)
            sage: Q.reduction(5).curve() == E.reduction(5)
            True

        ::

            sage: F.<a> = NumberField(x^2+5)
            sage: E = EllipticCurve(F,[1,2,3,4,0])
            sage: Q = E(98,931)
            sage: Q.reduction(a)
            (3 : 1 : 1)
            sage: Q.reduction(11)
            (10 : 7 : 1)

        ::

            sage: F.<a> = NumberField(x^3+x^2+1)
            sage: E = EllipticCurve(F,[a,2])
            sage: P = E(a,1)
            sage: P.reduction(F.ideal(5))
            (abar : 1 : 1)
            sage: P.reduction(F.ideal(a^2-4*a-2))
            (abar : 1 : 1)
        """
        P = self
        E = P.curve()
        return E.reduction(p)(P)

    def height(self, precision=None, normalised=True, algorithm='pari'):
        r"""Return the Néron-Tate canonical height of the point.

        INPUT:

        - ``self`` -- a point on an elliptic curve over a number field
          `K`.

        - ``precision`` -- positive integer, or None (default). The
          precision in bits of the result. If None, the default real
          precision is used.

        - ``normalised`` -- boolean. If True (default), the height is
          normalised to be invariant under extension of `K`. If False,
          return this normalised height multiplied by the degree of
          `K`.

        - ``algorithm`` -- string: either 'pari' (default) or 'sage'.
          If 'pari' and the base field is `\QQ`, use the PARI library
          function; otherwise use the Sage implementation.

        OUTPUT:

        The rational number 0, or a non-negative real number.

        There are two normalisations used in the literature, one of
        which is double the other. We use the larger of the two, which
        is the one appropriate for the BSD conjecture. This is
        consistent with [Cre1997]_ and double that of [Sil2009]_.

        See :wikipedia:`Néron-Tate height`.

        .. NOTE::

            The correct height to use for the regulator in the BSD
            formula is the non-normalised height.

        EXAMPLES::

            sage: E = EllipticCurve('11a'); E
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: P = E([5,5]); P
            (5 : 5 : 1)
            sage: P.height()
            0
            sage: Q = 5*P
            sage: Q.height()
            0

        ::

            sage: E = EllipticCurve('37a'); E
            Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: P = E([0,0])
            sage: P.height()
            0.0511114082399688
            sage: P.order()
            +Infinity
            sage: E.regulator()
            0.0511114082399688...

            sage: def naive_height(P):
            ....:     return log(RR(max(abs(P[0].numerator()), abs(P[0].denominator()))))
            sage: for n in [1..10]:
            ....:     print(naive_height(2^n*P)/4^n)
            0.000000000000000
            0.0433216987849966
            0.0502949347635656
            0.0511006335618645
            0.0511007834799612
            0.0511013666152466
            0.0511034199907743
            0.0511106492906471
            0.0511114081541082
            0.0511114081541180

        ::

            sage: E = EllipticCurve('4602a1'); E
            Elliptic Curve defined by y^2 + x*y  = x^3 + x^2 - 37746035*x - 89296920339 over Rational Field
            sage: x = 77985922458974949246858229195945103471590
            sage: y = 19575260230015313702261379022151675961965157108920263594545223
            sage: d = 2254020761884782243
            sage: E([ x / d^2,  y / d^3 ]).height()
            86.7406561381275

        ::

            sage: E = EllipticCurve([17, -60, -120, 0, 0]); E
            Elliptic Curve defined by y^2 + 17*x*y - 120*y = x^3 - 60*x^2 over Rational Field
            sage: E([30, -90]).height()
            0

            sage: E = EllipticCurve('389a1'); E
            Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field
            sage: [P,Q] = [E(-1,1),E(0,-1)]
            sage: P.height(precision=100)
            0.68666708330558658572355210295
            sage: (3*Q).height(precision=100)/Q.height(precision=100)
            9.0000000000000000000000000000
            sage: _.parent()
            Real Field with 100 bits of precision

        Canonical heights over number fields are implemented as well::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-2)
            sage: E = EllipticCurve([a, 4]); E
            Elliptic Curve defined by y^2 = x^3 + a*x + 4 over Number Field in a with defining polynomial x^3 - 2
            sage: P = E((0,2))
            sage: P.height()
            0.810463096585925
            sage: P.height(precision=100)
            0.81046309658592536863991810577
            sage: P.height(precision=200)
            0.81046309658592536863991810576865158896130286417155832378086
            sage: (2*P).height() / P.height()
            4.00000000000000
            sage: (100*P).height() / P.height()
            10000.0000000000

        Setting normalised=False multiplies the height by the degree of `K`::

            sage: E = EllipticCurve('37a')
            sage: P = E([0,0])
            sage: P.height()
            0.0511114082399688
            sage: P.height(normalised=False)
            0.0511114082399688
            sage: K.<z> = CyclotomicField(5)
            sage: EK = E.change_ring(K)
            sage: PK = EK([0,0])
            sage: PK.height()
            0.0511114082399688
            sage: PK.height(normalised=False)
            0.204445632959875

        Some consistency checks::

            sage: E = EllipticCurve('5077a1')
            sage: P = E([-2,3,1])
            sage: P.height()
            1.36857250535393

            sage: EK = E.change_ring(QuadraticField(-3,'a'))
            sage: PK = EK([-2,3,1])
            sage: PK.height()
            1.36857250535393

            sage: K.<i> = NumberField(x^2+1)
            sage: E = EllipticCurve(K, [0,0,4,6*i,0])
            sage: Q = E.lift_x(-9/4); Q
            (-9/4 : -27/8*i : 1)
            sage: Q.height()
            2.69518560017909
            sage: (15*Q).height() / Q.height()
            225.000000000000

            sage: E = EllipticCurve('37a')
            sage: P = E([0,-1])
            sage: P.height()
            0.0511114082399688
            sage: K.<a> = QuadraticField(-7)
            sage: ED = E.quadratic_twist(-7)
            sage: Q = E.isomorphism_to(ED.change_ring(K))(P); Q
            (0 : -7/2*a - 1/2 : 1)
            sage: Q.height()
            0.0511114082399688
            sage: Q.height(precision=100)
            0.051111408239968840235886099757

        An example to show that the bug at :trac:`5252` is fixed::

            sage: E = EllipticCurve([1, -1, 1, -2063758701246626370773726978, 32838647793306133075103747085833809114881])
            sage: P = E([-30987785091199, 258909576181697016447])
            sage: P.height()
            25.8603170675462
            sage: P.height(precision=100)
            25.860317067546190743868840741
            sage: P.height(precision=250)
            25.860317067546190743868840740735110323098872903844416215577171041783572513
            sage: P.height(precision=500)
            25.8603170675461907438688407407351103230988729038444162155771710417835725129551130570889813281792157278507639909972112856019190236125362914195452321720

            sage: P.height(precision=100) == P.non_archimedean_local_height(prec=100)+P.archimedean_local_height(prec=100)
            True

        An example to show that the bug at :trac:`8319` is fixed (correct height when the curve is not minimal)::

            sage: E = EllipticCurve([-5580472329446114952805505804593498080000,-157339733785368110382973689903536054787700497223306368000000])
            sage: xP = 204885147732879546487576840131729064308289385547094673627174585676211859152978311600/23625501907057948132262217188983681204856907657753178415430361
            sage: P = E.lift_x(xP)
            sage: P.height()
            157.432598516754
            sage: Q = 2*P
            sage: Q.height() # long time (4s)
            629.730394067016
            sage: Q.height()-4*P.height() # long time
            0.000000000000000

        An example to show that the bug at :trac:`12509` is fixed (precision issues)::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2-x-1)
            sage: v = [0, a + 1, 1, 28665*a - 46382, 2797026*a - 4525688]
            sage: E = EllipticCurve(v)
            sage: P = E([72*a - 509/5,  -682/25*a - 434/25])
            sage: P.height()
            1.38877711688727
            sage: (2*P).height()/P.height()
            4.00000000000000
            sage: (2*P).height(precision=100)/P.height(precision=100)
            4.0000000000000000000000000000
            sage: (2*P).height(precision=1000)/P.height(precision=1000)
            4.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

        This shows that the bug reported at :trac:`13951` has been fixed::

            sage: E = EllipticCurve([0,17])
            sage: P1 = E(2,5)
            sage: P1.height()
            1.06248137652528
            sage: F = E.change_ring(QuadraticField(-3,'a'))
            sage: P2 = F([2,5])
            sage: P2.height()
            1.06248137652528
        """
        if self.has_finite_order():
            return rings.QQ(0)

        E = self.curve()
        K = E.base_ring()

        if precision is None:
            precision = rings.RealField().precision()

        known_prec = -1
        try:
            height = self.__height
            known_prec = height.prec()
            if known_prec > precision:
                height = rings.RealField(precision)(height)
        except AttributeError:
            pass

        if known_prec < precision:
            if algorithm == 'pari' and K is rings.QQ:
                Emin = E.minimal_model()
                iso = E.isomorphism_to(Emin)
                P = iso(self)
                h = Emin.pari_curve().ellheight(P, precision=precision)
                height = rings.RealField(precision)(h)
            else:
                height = (self.non_archimedean_local_height(prec=precision)
                          + self.archimedean_local_height(prec=precision))

        # The cached height is the one that is independent of the base field.
        self.__height = height
        if not normalised:
            height *= K.degree()
        return height

    def archimedean_local_height(self, v=None, prec=None, weighted=False):
        """
        Compute the local height of self at the archimedean place `v`.

        INPUT:

        - ``self`` -- a point on an elliptic curve over a number field
          `K`.

        - ``v`` -- a real or complex embedding of K, or None (default).
          If `v` is a real or complex embedding, return the local
          height of self at `v`. If `v` is None, return the total
          archimedean contribution to the global height.

        - ``prec`` -- integer, or None (default). The precision of the
          computation. If None, the precision is deduced from v.

        - ``weighted`` -- boolean. If False (default), the height is
          normalised to be invariant under extension of `K`. If True,
          return this normalised height multiplied by the local degree
          if `v` is a single place, or by the degree of `K` if `v` is
          None.

        OUTPUT:

        A real number. The normalisation is twice that in Silverman's
        paper [Sil1988]_. Note that this local height depends on the
        model of the curve.

        ALGORITHM:

        See [Sil1988]_, Section 4.

        EXAMPLES:

        Examples 1, 2, and 3 from [Sil1988]_::

            sage: K.<a> = QuadraticField(-2)
            sage: E = EllipticCurve(K, [0,-1,1,0,0]); E
            Elliptic Curve defined by y^2 + y = x^3 + (-1)*x^2 over Number Field in a with defining polynomial x^2 + 2 with a = 1.414213562373095?*I
            sage: P = E.lift_x(2+a); P
            (a + 2 : 2*a + 1 : 1)
            sage: P.archimedean_local_height(K.places(prec=170)[0]) / 2
            0.45754773287523276736211210741423654346576029814695

            sage: K.<i> = NumberField(x^2+1)
            sage: E = EllipticCurve(K, [0,0,4,6*i,0]); E
            Elliptic Curve defined by y^2 + 4*y = x^3 + 6*i*x over Number Field in i with defining polynomial x^2 + 1
            sage: P = E((0,0))
            sage: P.archimedean_local_height(K.places()[0]) / 2
            0.510184995162373

            sage: Q = E.lift_x(-9/4); Q
            (-9/4 : -27/8*i : 1)
            sage: Q.archimedean_local_height(K.places()[0]) / 2
            0.654445619529600

        An example over the rational numbers::

            sage: E = EllipticCurve([0, 0, 0, -36, 0])
            sage: P = E([-3, 9])
            sage: P.archimedean_local_height()
            1.98723816350773

        Local heights of torsion points can be non-zero (unlike the
        global height)::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve([0, 0, 0, K(1), 0])
            sage: P = E(i, 0)
            sage: P.archimedean_local_height()
            0.346573590279973

        TESTS:

        See :trac:`12509`::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2-x-1)
            sage: v = [0, a + 1, 1, 28665*a - 46382, 2797026*a - 4525688]
            sage: E = EllipticCurve(v)
            sage: P = E([72*a - 509/5,  -682/25*a - 434/25])
            sage: P.archimedean_local_height()
            -0.220660795546828

        See :trac:`19276`::

            sage: K.<a> = NumberField(x^2-x-104)
            sage: E = EllipticCurve([1, a - 1, 1, -816765673272*a - 7931030674178, 1478955604013312315*a + 14361086227143654561])
            sage: P = E(5393511/49*a + 52372721/49 , -33896210324/343*a - 329141996591/343 )
            sage: P.height()
            0.974232017827741

        See :trac:`29966`::

            sage: K.<a> = NumberField(x^3 - x^2 - 6*x + 2)
            sage: E = EllipticCurve([1, -a^2 + 2*a + 4, 0, -6056450500590472699700624*a^2 - 11239394326797569935861742*a + 4241549693833829432516231, 1904879037869682826729875958079326124520*a^2 + 3535022146945771697732350459284777382011*a - 1334055169621036218710397707677347972626])
            sage: P = E([1033399668533*a^2 + 1917754693229*a - 723726883800 , 12536493059202326563*a^2 + 23264879148900575548*a - 8779756111574815918 , 1])
            sage: P.height()
            0.297318833424763
            sage: (2*P).height() / P.height()
            4.00000000000000
            sage: P.height(200)
            0.29731883342476341806143743594519935578696537745294661858984
            sage: (2*P).height(200) / P.height(200)
            4.0000000000000000000000000000000000000000000000000000000000
        """
        from sage.rings.number_field.number_field import refine_embedding
        from sage.all import RealField, ComplexField, Infinity

        E = self.curve()
        K = E.base_ring()

        if v is None:

            if prec is None:
                prec = 53
            if K is rings.QQ:
                v = K.embeddings(rings.RR)[0]
                h = self.archimedean_local_height(v, prec+10)
            else:
                r1, r2 = K.signature()
                pl = K.places()
                h = (sum(self.archimedean_local_height(pl[i], prec+10, weighted=False)
                         for i in range(r1))
                     + 2 * sum(self.archimedean_local_height(pl[i], prec+10, weighted=False)
                               for i in range(r1, r1 + r2)))
                if not weighted:
                    h /= K.degree()
            return RealField(prec)(h)

        prec_v = v.codomain().prec()
        if prec is None:
            prec = prec_v
        if K is rings.QQ:
            v = K.embeddings(RealField())[0]
        v_inf = refine_embedding(v, Infinity)
        v_is_real = v_inf(K.gen()).imag().is_zero()

        # Find a suitable working precision.  See trac#29666 for an
        # example where 100 extra bits is not enough when the
        # discriminant is ~1e-92, but this code uses a working
        # precision of 333 and gets it right.

        D = v_inf(E.discriminant())

        if D.abs().real_number(RealField()).round():
            extra_prec = 100
        else:  # then |D| is small
            extra_prec = 10 + (1/D).abs().real_number(RealField()).round().nbits()

        working_prec = prec + extra_prec
        RC = RealField(working_prec) if v_is_real else ComplexField(working_prec)
        # print("Using working precision {}, |D| = {}".format(working_prec, RC(D).abs()))

        # NB We risk losing much precision if we compute the embedding
        # of K into RR or CC to some precision and then apply that to
        # elements of K.  Instead we map elements of K into AA or Qbar
        # (with infinite precision) and then trim back to RR or CC.

        x = RC(v_inf(self[0]))
        b2, b4, b6, b8 = [RC(v_inf(b)) for b in E.b_invariants()]

        # The following comes from Silverman Theorem 4.2.  Silverman
        # uses decimal precision d, so his term (5/3)d =
        # (5/3)*(log(2)/log(10))*prec = 0.5017*prec, which we round
        # up.  The rest of the expression was wrongly transcribed in
        # Sage versions <5.6 (see #12509).

        H = max(RC(4).abs(), b2.abs(), 2*b4.abs(), 2*b6.abs(), b8.abs())
        absdisc = RC(v_inf(E.discriminant())).abs()
        adl3 = 0 if absdisc >= 1 else absdisc.log() / 3
        nterms = int(math.ceil(0.51*working_prec + 0.5 + 0.75 * (7 + 4*H.log()/3 - adl3).log()))

        b2p = b2 - 12
        b4p = b4 - b2 + 6
        b6p = b6 - 2*b4 + b2 - 4
        b8p = b8 - 3*b6 + 3*b4 - b2 + 3

        fz = lambda T: 1 - T**2 * (b4 + T*(2*b6 + T*b8))
        fzp = lambda T: 1 - T**2 * (b4p + T*(2*b6p + T*b8p))
        fw = lambda T: T*(4 + T*(b2 + T*(2*b4 + T*b6)))
        fwp = lambda T: T*(4 + T*(b2p + T*(2*b4p + T*b6p)))

        if abs(x) >= .5:
            t = 1/x
            beta = True
        else:
            t = 1/(x+1)
            beta = False
        lam = -t.abs().log()
        mu = 0
        four_to_n = rings.QQ(1)

        for n in range(nterms):
            if beta:
                w = fw(t)
                z = fz(t)
                if abs(w) <= 2 * z.abs():
                    mu += four_to_n * z.abs().log()
                    t = w/z
                else:
                    mu += four_to_n * (z+w).abs().log()
                    t = w/(z+w)
                    beta = not beta
            else:
                w = fwp(t)
                z = fzp(t)
                if abs(w) <= 2 * z.abs():
                    mu += four_to_n * z.abs().log()
                    t = w/z
                else:
                    mu += four_to_n * (z-w).abs().log()
                    t = w/(z-w)
                    beta = not beta
            four_to_n >>= 2

        h = RealField(prec)(lam + mu/4)
        if weighted and not v_is_real:
            h *= 2
        return h

    def non_archimedean_local_height(self, v=None, prec=None,
                                     weighted=False, is_minimal=None):
        """
        Compute the local height of self at the non-archimedean place `v`.

        INPUT:

        - ``self`` -- a point on an elliptic curve over a number field
          `K`.

        - ``v`` -- a non-archimedean place of `K`, or None (default).
          If `v` is a non-archimedean place, return the local height
          of self at `v`. If `v` is None, return the total
          non-archimedean contribution to the global height.

        - ``prec`` -- integer, or None (default). The precision of the
          computation. If None, the height is returned symbolically.

        - ``weighted`` -- boolean. If False (default), the height is
          normalised to be invariant under extension of `K`. If True,
          return this normalised height multiplied by the local degree
          if `v` is a single place, or by the degree of `K` if `v` is
          None.

        OUTPUT:

        A real number. The normalisation is twice that in Silverman's
        paper [Sil1988]_. Note that this local height depends on the
        model of the curve.

        ALGORITHM:

        See [Sil1988]_, Section 5.

        EXAMPLES:

        Examples 2 and 3 from [Sil1988]_::

            sage: K.<i> = NumberField(x^2+1)
            sage: E = EllipticCurve(K, [0,0,4,6*i,0]); E
            Elliptic Curve defined by y^2 + 4*y = x^3 + 6*i*x over Number Field in i with defining polynomial x^2 + 1
            sage: P = E((0,0))
            sage: P.non_archimedean_local_height(K.ideal(i+1))
            -1/2*log(2)
            sage: P.non_archimedean_local_height(K.ideal(3))
            0
            sage: P.non_archimedean_local_height(K.ideal(1-2*i))
            0

            sage: Q = E.lift_x(-9/4); Q
            (-9/4 : -27/8*i : 1)
            sage: Q.non_archimedean_local_height(K.ideal(1+i))
            2*log(2)
            sage: Q.non_archimedean_local_height(K.ideal(3))
            0
            sage: Q.non_archimedean_local_height(K.ideal(1-2*i))
            0
            sage: Q.non_archimedean_local_height()
            2*log(2)

        An example over the rational numbers::

            sage: E = EllipticCurve([0, 0, 0, -36, 0])
            sage: P = E([-3, 9])
            sage: P.non_archimedean_local_height()
            -log(3)

        Local heights of torsion points can be non-zero (unlike the
        global height)::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve([0, 0, 0, K(1), 0])
            sage: P = E(i, 0)
            sage: P.non_archimedean_local_height()
            -1/2*log(2)

        TESTS::

            sage: Q.non_archimedean_local_height(prec=100)
            1.3862943611198906188344642429
            sage: (3*Q).non_archimedean_local_height()
            1/2*log(75923153929839865104)

            sage: F.<a> = NumberField(x^4 + 2*x^3 + 19*x^2 + 18*x + 288)
            sage: F.ring_of_integers().basis()
            [1, 5/6*a^3 + 1/6*a, 1/6*a^3 + 1/6*a^2, a^3]
            sage: F.class_number()
            12
            sage: E = EllipticCurve('37a').change_ring(F)
            sage: P = E((-a^2/6 - a/6 - 1, a)); P
            (-1/6*a^2 - 1/6*a - 1 : a : 1)
            sage: P[0].is_integral()
            True
            sage: P.non_archimedean_local_height()
            0

        This shows that the bug reported at :trac:`13951` has been fixed::

            sage: E = EllipticCurve([0,17])
            sage: P = E(2,5)
            sage: P.non_archimedean_local_height(2)
            -2/3*log(2)
        """
        if prec:
            log = lambda x: rings.RealField(prec)(x).log()
        else:
            from sage.functions.log import log

        if v is None:
            D = self.curve().discriminant()
            K = self.curve().base_ring()
            if K is rings.QQ:
                factorD = D.factor()
                if self[0] == 0:
                    c = 1
                else:
                    c = self[0].denominator()
                # The last sum is for bad primes that divide c where
                # the model is not minimal.
                h = (log(c)
                     + sum(self.non_archimedean_local_height(p, prec, weighted=True, is_minimal=(e < 12))
                           for p, e in factorD if not p.divides(c))
                     + sum(self.non_archimedean_local_height(p, prec, weighted=True)
                           - c.valuation(p) * log(p)
                           for p, e in factorD if e >= 12 and p.divides(c)))
            else:
                factorD = K.factor(D)
                if self[0] == 0:
                    c = K.ideal(1)
                else:
                    c = K.ideal(self[0]).denominator()
                # The last sum is for bad primes that divide c where
                # the model is not minimal.
                h = (log(c.norm())
                     + sum(self.non_archimedean_local_height(v, prec, weighted=True, is_minimal=(e < 12))
                           for v, e in factorD if not v.divides(c))
                     + sum(self.non_archimedean_local_height(v, prec, weighted=True)
                           - c.valuation(v) * log(v.norm())
                           for v, e in factorD if e >= 12 and v.divides(c)))
                if not weighted:
                    h /= K.degree()
            return h

        if is_minimal:
            E = self.curve()
            P = self
            offset = ZZ.zero()
        else:
            E = self.curve().local_minimal_model(v)
            P = self.curve().isomorphism_to(E)(self)
            # Silverman's normalization is not invariant under change of model,
            # but it all cancels out in the global height.
            offset = (self.curve().discriminant()/E.discriminant()).valuation(v)

        a1, a2, a3, a4, a6 = E.a_invariants()
        b2, b4, b6, b8 = E.b_invariants()
        c4 = E.c4()
        x, y = P.xy()
        D = E.discriminant()
        N = D.valuation(v)
        A = (3*x**2 + 2*a2*x + a4 - a1*y).valuation(v)
        B = (2*y+a1*x+a3).valuation(v)
        C = (3*x**4 + b2*x**3 + 3*b4*x**2 + 3*b6*x + b8).valuation(v)
        if A <= 0 or B <= 0:
            r = max(0, -x.valuation(v))
        elif c4.valuation(v) == 0:
            n = min(B, N/2)
            r = -n*(N-n)/N
        elif C >= 3*B:
            r = -2*B/3
        else:
            r = -C/4
        r -= offset/6
        if not r:
            return rings.QQ.zero()
        else:
            if E.base_ring() is rings.QQ:
                Nv = Integer(v)
            else:
                Nv = v.norm()
                if not weighted:
                    r = r / (v.ramification_index() * v.residue_class_degree())
            return r * log(Nv)

    def elliptic_logarithm(self, embedding=None, precision=100,
                           algorithm='pari'):
        r"""
        Return the elliptic logarithm of this elliptic curve point.

        An embedding of the base field into `\RR` or `\CC` (with
        arbitrary precision) may be given; otherwise the first real
        embedding is used (with the specified precision) if any, else
        the first complex embedding.

        INPUT:

        - ``embedding``: an embedding of the base field into `\RR` or `\CC`

        - ``precision``: a positive integer (default 100) setting the
          number of bits of precision for the computation

        - ``algorithm``: either 'pari' (default for real embeddings)
          to use PARI's :pari:`ellpointtoz`, or 'sage' for a native
          implementation.  Ignored for complex embeddings.

        ALGORITHM:

        See [Coh1993]_ for the case of real embeddings,
        and Cremona, J.E. and Thongjunthug, T. 2010 for the complex
        case.

        AUTHORS:

        - Michael Mardaus (2008-07),
        - Tobias Nagel (2008-07) -- original version from [Coh1993]_.
        - John Cremona (2008-07) -- revision following eclib code.
        - John Cremona (2010-03) -- implementation for complex embeddings.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: E.discriminant() > 0
            True
            sage: P = E([-1,1])
            sage: P.is_on_identity_component ()
            False
            sage: P.elliptic_logarithm (precision=96)
            0.4793482501902193161295330101 + 0.985868850775824102211203849...*I
            sage: Q=E([3,5])
            sage: Q.is_on_identity_component()
            True
            sage: Q.elliptic_logarithm (precision=96)
            1.931128271542559442488585220

        An example with negative discriminant, and a torsion point::

            sage: E = EllipticCurve('11a1')
            sage: E.discriminant() < 0
            True
            sage: P = E([16,-61])
            sage: P.elliptic_logarithm(precision=70)
            0.25384186085591068434
            sage: E.period_lattice().real_period(prec=70) / P.elliptic_logarithm(precision=70)
            5.0000000000000000000

        A larger example.  The default algorithm uses PARI and makes
        sure the result has the requested precision::

            sage: E = EllipticCurve([1, 0, 1, -85357462, 303528987048]) #18074g1
            sage: P = E([4458713781401/835903744, -64466909836503771/24167649046528, 1])
            sage: P.elliptic_logarithm()  # 100 bits
            0.27656204014107061464076203097

        The native algorithm 'sage' used to have trouble with
        precision in this example, but no longer::

            sage: P.elliptic_logarithm(algorithm='sage')  # 100 bits
            0.27656204014107061464076203097

        This shows that the bug reported at :trac:`4901` has been fixed::

            sage: E = EllipticCurve("4390c2")
            sage: P = E(683762969925/44944,-565388972095220019/9528128)
            sage: P.elliptic_logarithm()
            0.00025638725886520225353198932529
            sage: P.elliptic_logarithm(precision=64)
            0.000256387258865202254
            sage: P.elliptic_logarithm(precision=65)
            0.0002563872588652022535
            sage: P.elliptic_logarithm(precision=128)
            0.00025638725886520225353198932528666427412
            sage: P.elliptic_logarithm(precision=129)
            0.00025638725886520225353198932528666427412
            sage: P.elliptic_logarithm(precision=256)
            0.0002563872588652022535319893252866642741168388008346370015005142128009610936373
            sage: P.elliptic_logarithm(precision=257)
            0.00025638725886520225353198932528666427411683880083463700150051421280096109363730

        Examples over number fields::

            sage: K.<a> = NumberField(x^3-2)
            sage: embs = K.embeddings(CC)
            sage: E = EllipticCurve([0,1,0,a,a])
            sage: Ls = [E.period_lattice(e) for e in embs]
            sage: [L.real_flag for L in Ls]
            [0, 0, -1]
            sage: P = E(-1,0)  # order 2
            sage: [L.elliptic_logarithm(P) for L in Ls]
            [-1.73964256006716 - 1.07861534489191*I, -0.363756518406398 - 1.50699412135253*I, 1.90726488608927]

            sage: E = EllipticCurve([-a^2 - a - 1, a^2 + a])
            sage: Ls = [E.period_lattice(e) for e in embs]
            sage: pts = [E(2*a^2 - a - 1 , -2*a^2 - 2*a + 6 ), E(-2/3*a^2 - 1/3 , -4/3*a - 2/3 ), E(5/4*a^2 - 1/2*a , -a^2 - 1/4*a + 9/4 ), E(2*a^2 + 3*a + 4 , -7*a^2 - 10*a - 12 )]
            sage: [[L.elliptic_logarithm(P) for P in pts] for L in Ls]
            [[0.250819591818930 - 0.411963479992219*I, -0.290994550611374 - 1.37239400324105*I, -0.693473752205595 - 2.45028458830342*I, -0.151659609775291 - 1.48985406505459*I], [1.33444787667954 - 1.50889756650544*I, 0.792633734249234 - 0.548467043256610*I, 0.390154532655013 + 0.529423541805758*I, 0.931968675085317 - 0.431006981443071*I], [1.14758249500109 + 0.853389664016075*I, 2.59823462472518 + 0.853389664016075*I, 1.75372176444709, 0.303069634723001]]

        ::

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve([0,0,0,9*i-10,21-i])
            sage: emb = K.embeddings(CC)[1]
            sage: L = E.period_lattice(emb)
            sage: P = E(2-i,4+2*i)
            sage: L.elliptic_logarithm(P,prec=100)
            0.70448375537782208460499649302 - 0.79246725643650979858266018068*I
        """
        from sage.rings.number_field.number_field import refine_embedding
        from sage.rings.all import RealField, ComplexField, QQ

        # Check the trivial case:

        C = ComplexField(precision)
        if self.is_zero():
            return C.zero()

        # find a suitable embedding if none was supplied:

        E = self.curve()
        K = E.base_field()
        rational = (K is QQ)
        emb = embedding

        if emb is None:
            emb = K.embeddings(RealField(precision))
            if emb:
                emb = emb[0]
            else:
                emb = K.embeddings(ComplexField(precision))[0]
        else:
            # Get the precision of the supplied embedding
            prec = emb.codomain().precision()
            # if the precision parameter is greater, refine the embedding:
            if precision > prec:
                emb = refine_embedding(emb, precision)

        L = E.period_lattice(emb)

        if algorithm == 'sage' or not isinstance(emb.codomain, sage.rings.abc.RealField):
            return L.elliptic_logarithm(self, precision)

        if algorithm != 'pari':
            raise ValueError("algorithm must be either 'pari' or 'sage'")

        # From now on emb() is a real embedding of K into
        # RealField(precision).  We interface with the PARI library.

        x, y = self.xy()
        if rational:        # work with exact coordinates
            E_work = E
            pt_pari = pari([x, y])
        else:               # use the embedding to get real coordinates
            ai = [emb(a) for a in E.a_invariants()]
            E_work = EllipticCurve(ai)  # defined over RR
            pt_pari = pari([emb(x), emb(y)])
        working_prec = precision
        E_pari = E_work.pari_curve()
        log_pari = E_pari.ellpointtoz(pt_pari, precision=working_prec)

        while prec_words_to_bits(log_pari.precision()) < precision:
            # result is not precise enough, re-compute with double
            # precision. if the base field is not QQ, this
            # requires modifying the precision of the embedding,
            # the curve, and the point
            working_prec = 2*working_prec
            if not rational:
                emb = refine_embedding(emb, working_prec)
                ai = [emb(a) for a in E.a_invariants()]
                E_work = EllipticCurve(ai)  # defined over RR
                pt_pari = pari([emb(x), emb(y)])
            E_pari = E_work.pari_curve()
            log_pari = E_pari.ellpointtoz(pt_pari, precision=working_prec)

        # normalization step
        r, i = C(log_pari)
        wR, wI = L.basis(prec=precision)
        k = (r/wR).floor()
        if k:
            r -= k*wR
        if self.is_on_identity_component(emb):
            return C(r)
        # Now there are two components and P is on the non-identity one
        return C(r)+C(wI/2)

    def padic_elliptic_logarithm(self, p, absprec=20):
        r"""
        Computes the `p`-adic elliptic logarithm of this point.

        INPUT:

        ``p`` - integer: a prime ``absprec`` - integer (default: 20):
        the initial `p`-adic absolute precision of the computation

        OUTPUT:

        The `p`-adic elliptic logarithm of self, with precision ``absprec``.

        AUTHORS:

        - Tobias Nagel
        - Michael Mardaus
        - John Cremona

        ALGORITHM:

        For points in the formal group (i.e. not integral at `p`) we
        take the ``log()`` function from the formal groups module and
        evaluate it at `-x/y`.  Otherwise we first multiply the point
        to get into the formal group, and divide the result
        afterwards.

        .. TODO::

            See comments at :trac:`4805`.  Currently the absolute
            precision of the result may be less than the given value
            of absprec, and error-handling is imperfect.

        EXAMPLES::

            sage: E = EllipticCurve([0,1,1,-2,0])
            sage: E(0).padic_elliptic_logarithm(3)
            0
            sage: P = E(0,0)
            sage: P.padic_elliptic_logarithm(3)
            2 + 2*3 + 3^3 + 2*3^7 + 3^8 + 3^9 + 3^11 + 3^15 + 2*3^17 + 3^18 + O(3^19)
            sage: P.padic_elliptic_logarithm(3).lift()
            660257522
            sage: P = E(-11/9,28/27)
            sage: [(2*P).padic_elliptic_logarithm(p)/P.padic_elliptic_logarithm(p) for p in prime_range(20)]  # long time (3s)
            [2 + O(2^19), 2 + O(3^20), 2 + O(5^19), 2 + O(7^19), 2 + O(11^19), 2 + O(13^19), 2 + O(17^19), 2 + O(19^19)]
            sage: [(3*P).padic_elliptic_logarithm(p)/P.padic_elliptic_logarithm(p) for p in prime_range(12)]  # long time (2s)
            [1 + 2 + O(2^19), 3 + 3^20 + O(3^21), 3 + O(5^19), 3 + O(7^19), 3 + O(11^19)]
            sage: [(5*P).padic_elliptic_logarithm(p)/P.padic_elliptic_logarithm(p) for p in prime_range(12)]  # long time (2s)
            [1 + 2^2 + O(2^19), 2 + 3 + O(3^20), 5 + O(5^19), 5 + O(7^19), 5 + O(11^19)]

        An example which arose during reviewing :trac:`4741`::

            sage: E = EllipticCurve('794a1')
            sage: P = E(-1,2)
            sage: P.padic_elliptic_logarithm(2) # default precision=20
            2^4 + 2^5 + 2^6 + 2^8 + 2^9 + 2^13 + 2^14 + 2^15 + O(2^16)
            sage: P.padic_elliptic_logarithm(2, absprec=30)
            2^4 + 2^5 + 2^6 + 2^8 + 2^9 + 2^13 + 2^14 + 2^15 + 2^22 + 2^23 + 2^24 + O(2^26)
            sage: P.padic_elliptic_logarithm(2, absprec=40)
            2^4 + 2^5 + 2^6 + 2^8 + 2^9 + 2^13 + 2^14 + 2^15 + 2^22 + 2^23 + 2^24 + 2^28 + 2^29 + 2^31 + 2^34 + O(2^35)
        """
        if not p.is_prime():
            raise ValueError('p must be prime')
        debug = False  # True
        if debug:
            print("P=", self, "; p=", p, " with precision ", absprec)
        E = self.curve()
        Q_p = Qp(p, absprec)
        if self.has_finite_order():
            return Q_p(0)
        while True:
            try:
                Ep = E.change_ring(Q_p)
                P = Ep(self)
                x, y = P.xy()
                break
            except (PrecisionError, ArithmeticError, ZeroDivisionError):
                absprec *= 2
                Q_p = Qp(p, absprec)
        if debug:
            print("x,y=", (x, y))
        f = 1   # f will be such that f*P is in the formal group E^1(Q_p)
        if x.valuation() >= 0:   # P is not in E^1
            if not self.has_good_reduction(p):   # P is not in E^0
                n = E.tamagawa_exponent(p)   # n*P has good reduction at p
                if debug:
                    print("Tamagawa exponent = =", n)
                f = n
                P = n*P   # lies in E^0
                if debug:
                    print("P=", P)
                try:
                    x, y = P.xy()
                except ZeroDivisionError:
                    raise ValueError("Insufficient precision in "
                                     "p-adic_elliptic_logarithm()")
                if debug:
                    print("x,y=", (x, y))
            if x.valuation() >= 0:   # P is still not in E^1
                t = E.local_data(p).bad_reduction_type()
                if t is None:
                    m = E.reduction(p).abelian_group().exponent()
                else:
                    m = p - t
                if debug:
                    print("mod p exponent = =", m)
                    # now m*(n*P) reduces to the identity mod p, so is
                    # in E^1(Q_p)
                f *= m
                P = m*P   # lies in E^1
                try:
                    x, y = P.xy()
                except ZeroDivisionError:
                    raise ValueError("Insufficient precision in "
                                     "p-adic_elliptic_logarithm()")
                if debug:
                    print("f=", f)
                    print("x,y=", (x, y))
        vx = x.valuation()
        vy = y.valuation()
        v = vx-vy
        if not (v > 0 and vx == -2*v and vy == -3*v):
            raise ValueError("Insufficient precision in "
                             "p-adic_elliptic_logarithm()")
        try:
            t = -x/y
        except (ZeroDivisionError, PrecisionError):
            raise ValueError("Insufficient precision in "
                             "p-adic_elliptic_logarithm()")
        if debug:
            print("t=", t, ", with valuation ", v)
        phi = Ep.formal().log(prec=1+absprec//v)
        return phi(t)/f


class EllipticCurvePoint_finite_field(EllipticCurvePoint_field):
    r"""
    Class for elliptic curve points over finite fields.
    """
    def _magma_init_(self, magma):
        """
        Return a string representation of self that ``MAGMA`` can
        use for input.

        EXAMPLES::

            sage: E = EllipticCurve(GF(17), [1,-1])
            sage: P = E([13, 4])
            sage: P._magma_init_(magma)                               # optional - magma
            'EllipticCurve([_sage_ref...|GF(17)!0,GF(17)!0,GF(17)!0,GF(17)!1,GF(17)!16])![13,4]'
        """
        E = self.curve()._magma_init_(magma)
        x, y = self.xy()
        return "%s![%s,%s]" % (E, x, y)

    def _acted_upon_(self, other, side):
        r"""
        We implement ``_acted_upon_`` to keep track of cached
        point orders when scalar multiplications are applied.

        EXAMPLES::

            sage: P = EllipticCurve(GF(65537), [2,2]).lift_x(6)
            sage: P.order().factor()
            2^2 * 3 * 37^2
            sage: getattr(74*P, '_order', None)
            222
            sage: getattr(P*4070, '_order', None)
            222
            sage: getattr(506*P*37, '_order', None)
            222
        """
        k = ZZ(other)
        E = self.curve()

        try:
            pariQ = pari.ellmul(E, self, k)
        except PariError:
            pariQ = None

        if pariQ is not None:
            if pariQ == [0]:
                vQ = 0
            else:
                assert len(pariQ) == 2
                vQ = Sequence(tuple(pariQ) + (1,), E.base_field())
            Q = EllipticCurvePoint_finite_field(E, vQ, check=False)

        else:
            Q = IntegerMulAction(ZZ, self.parent())._act_(k, self)

        n = getattr(self, '_order', None)
        if n is not None:
            Q._order = n // n.gcd(k)  # Lagrange's theorem

        return Q

    def discrete_log(self, Q, ord=None):
        r"""
        Return the discrete logarithm of `Q` to base `P` = ``self``,
        that is, an integer `x` such that `xP = Q`.

        A :class:`ValueError` is raised if there is no solution.

        ALGORITHM:

        To compute the actual logarithm, :pari:`elllog` is called.

        However, ``elllog()`` does not guarantee termination if `Q`
        is not a multiple of `P`, so we first need to check subgroup
        membership. This is done as follows:

        - Let `n` denote the order of `P`. First check that `nQ` equals
          the point at infinity (and hence the order of `Q` divides `n`).

        - If the curve order `\#E` has been cached, check whether
          `\gcd(n^2, \#E) = n`. If this holds, the curve has cyclic
          `n`-torsion, hence all points whose order divides `n` must be
          multiples of `P` and we are done.

        - Otherwise (if this test is inconclusive), check that the Weil
          pairing of `P` and `Q` is trivial.

        INPUT:

        - ``Q`` (point) -- another point on the same curve as ``self``.

        OUTPUT:

        (integer) -- The discrete logarithm of `Q` with respect to `P`,
        which is an integer `x` with `0\le x<\mathrm{ord}(P)` such that
        `xP=Q`, if one exists.

        AUTHORS:

        - John Cremona. Adapted to use generic functions 2008-04-05.
        - Lorenz Panny (2022): switch to PARI.

        EXAMPLES::

            sage: F = GF((3,6),'a')
            sage: a = F.gen()
            sage: E = EllipticCurve([0,1,1,a,a])
            sage: E.cardinality()
            762
            sage: P = E.gens()[0]
            sage: Q = 400*P
            sage: P.discrete_log(Q)
            400

        TESTS:

        Some random testing::

            sage: sz = randint(8,32)
            sage: e = randint(1,3)
            sage: p = random_prime(ceil(2**(sz/e)))
            sage: E = EllipticCurve(j=GF((p,e),'a').random_element())
            sage: P = E.random_point()
            sage: Q = randrange(2**999) * P
            sage: x = P.discrete_log(Q)
            sage: x*P == Q
            True

        Doctest deprecation::

            sage: P.discrete_log(Q, ord=P.order())
            doctest:warning
            ...
            DeprecationWarning: The "ord" argument to .discrete_log() is obsolete. ...
        """
        if ord is not None:
            from sage.misc.superseded import deprecation
            deprecation(33121, 'The "ord" argument to .discrete_log() is obsolete. Use the .set_order() method instead.')
            self.set_order(ord)
        if Q not in self.parent():
            raise ValueError('not a point on the same curve')
        n = self.order()
        if n*Q:
            raise ValueError('ECDLog problem has no solution (order of Q does not divide order of P)')
        E = self.curve()
        if hasattr(E, '_order') and E._order.gcd(n**2) == n:
            pass    # cyclic rational n-torsion -> okay
        elif self.weil_pairing(Q, n) != 1:
            raise ValueError('ECDLog problem has no solution (non-trivial Weil pairing)')
        return ZZ(pari.elllog(self.curve(), Q, self, n))

    def has_finite_order(self):
        r"""
        Return ``True`` if this point has finite additive order as an element
        of the group of points on this curve.

        Since the base field is finite, the answer will always be ``True``.

        EXAMPLES::

            sage: E = EllipticCurve(GF(7), [1,3])
            sage: P = E.points()[3]
            sage: P.has_finite_order()
            True
        """
        return True

    def order(self):
        r"""
        Return the order of this point on the elliptic curve.

        ALGORITHM: Use PARI function :pari:`ellorder`.

        .. NOTE::

            :meth:`additive_order` is a synonym for :meth:`order`

        EXAMPLES::

            sage: k.<a> = GF((5,5))
            sage: E = EllipticCurve(k,[2,4]); E
            Elliptic Curve defined by y^2 = x^3 + 2*x + 4 over Finite Field in a of size 5^5
            sage: P = E(3*a^4 + 3*a , 2*a + 1 )
            sage: P.order()
            3227
            sage: Q = E(0,2)
            sage: Q.order()
            7
            sage: Q.additive_order()
            7

        ::

            sage: p=next_prime(2^150)
            sage: E=EllipticCurve(GF(p),[1,1])
            sage: P=E(831623307675610677632782670796608848711856078, 42295786042873366706573292533588638217232964)
            sage: P.order()
            1427247692705959881058262545272474300628281448
            sage: P.order() == E.cardinality()
            True

        The next example has `j(E)=0`::

            sage: p = 33554501
            sage: F.<u> = GF((p,2))
            sage: E = EllipticCurve(F,[0,1])
            sage: E.j_invariant()
            0
            sage: P = E.random_point()
            sage: P.order() # random
            16777251

        Similarly when `j(E)=1728`::

            sage: p = 33554473
            sage: F.<u> = GF((p,2))
            sage: E = EllipticCurve(F,[1,0])
            sage: E.j_invariant()
            1728
            sage: P = E.random_point()
            sage: P.order() # random
            46912611635760

        TESTS:

        Check that the order actually gets cached (:trac:`32786`)::

            sage: E = EllipticCurve(GF(31337), [42,1])
            sage: P = E.lift_x(1)
            sage: hasattr(P, '_order')
            False
            sage: P.order()
            15649
            sage: P._order
            15649

        The curve order should also get cached as a side effect
        of computing a point order::

            sage: E._order
            31298
        """
        try:
            return self._order
        except AttributeError:
            pass

        E = self.curve()

        if getattr(E, '_order', None) is None:
            # The curve order will be computed and cached by PARI during
            # ellorder() anyway. We might as well cache it here too.
            E._order = Integer(E.pari_curve().ellcard())

        self._order = Integer(E.pari_curve().ellorder(self, E._order))
        return self._order

    additive_order = order
