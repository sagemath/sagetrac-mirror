"""
General L-series

AUTHOR:

- William Stein

.. TODO::

    - Symmetric powers (and modular degree -- see trac 9758)
    - Triple product L-functions: Gross-Kudla, Zhang, etc -- see the
      code in triple_prod/triple.py
    - Support L-calc L-function
    - Make it so we use exactly one GP session for *all* of the
      Dokchitser L-functions
    - Tensor products
    - Genus 2 curves, via smalljac and genus2reduction
    - Fast L-series of elliptic curves over number fields (not just
      sqrt(5)), via smalljac
    - Inverse of number_of_coefficients function.
"""

################################################################################
#
# (c) Copyright 2011 William Stein
#
#  This file is part of SAGE
#
#  SAGE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  SAGE is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
################################################################################

import copy
import math
import types

from sage.all import prime_range, cached_method, sqrt, SR, vector
from sage.arith.misc import divisors, primes
from sage.lfunctions.dokchitser import Dokchitser
from sage.misc.all import prod
from sage.misc.mrange import cartesian_product_iterator
from sage.modular.abvar.abvar import is_ModularAbelianVariety
from sage.modular.all import Newform
from sage.modular.dirichlet import is_DirichletCharacter
from sage.modular.modsym.space import is_ModularSymbolsSpace
from sage.rings.all import ZZ, Integer, QQ, O, ComplexField, CDF, infinity as oo
from sage.rings.number_field.number_field_base import is_NumberField
from sage.rings.rational_field import is_RationalField
from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve
from sage.structure.factorization import Factorization
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.categories.monoids import Monoids
import sage.modular.modform.element

from sage.lfunctions.eulerprod_fast import extend_multiplicatively_generic

I = sqrt(-1)


def prec(s):
    """
    Return precision of s, if it has a precision attribute.  Otherwise
    return 53.  This is useful internally in this module.

    EXAMPLES::

        sage: from sage.lfunctions.eulerprod import prec
        sage: prec(ComplexField(100)(1))
        100
        sage: prec(RealField(125)(1))
        125
        sage: prec(1/3)
        53
    """
    if hasattr(s, 'prec'):
        return s.prec()
    return 53


def norm(a):
    """
    Return the norm of a, for a in either a number field or QQ.

    This is a function used internally in this module, mainly because
    elements of QQ and ZZ have no norm method.

    EXAMPLES::

        sage: from sage.lfunctions.eulerprod import norm
        sage: K.<a> = NumberField(x^2-x-1)
        sage: norm(a+5)
        29
        sage: (a+5).norm()
        29
        sage: norm(17)
        17
    """
    try:
        return a.norm()
    except AttributeError:
        return a


def tiny(prec):
    """
    Return a number that we consider tiny to the given precision ``prec``
    in bits.   This is used in various places as "zero" to the given
    precision for various checks, e.g., of correctness of the
    functional equation.
    """
    return max(1e-8, 1.0 / 2**(prec - 1))


def prime_below(P):
    """
    Return the prime in ZZ below the prime P (or element of QQ).

    EXAMPLES::

        sage: from sage.lfunctions.eulerprod import prime_below
        sage: K.<a> = NumberField(x^2-x-1)
        sage: prime_below(K.prime_above(11))
        11
        sage: prime_below(K.prime_above(5))
        5
        sage: prime_below(K.prime_above(3))
        3
        sage: prime_below(7)
        7
    """
    try:
        return P.smallest_integer()
    except AttributeError:
        return ZZ(P)


class LSeriesDerivative(object):
    """
    The formal derivative of an L-series.

    EXAMPLES::

        sage: from sage.lfunctions.eulerprod import LSeries
        sage: L = LSeries('delta')
        sage: L.derivative()
        First derivative of L-function associated to Ramanujan's Delta
        (a weight 12 cusp form)
        sage: L.derivative()(11/2)
        0.125386233743526

    We directly create an instance of the class (users should not need
    to do this)::

        sage: from sage.lfunctions.eulerprod import LSeriesDerivative
        sage: Ld = LSeriesDerivative(L, 2); Ld
        Second derivative of L-function associated to Ramanujan's Delta
        (a weight 12 cusp form)
        sage: type(Ld)
        <class 'sage.lfunctions.eulerprod.LSeriesDerivative'>
    """
    def __init__(self, lseries, k):
        """
        INPUT:

        - lseries -- any LSeries object (derives from LseriesAbstract)
        - k -- positive integer
        """
        k = ZZ(k)
        if k <= 0:
            raise ValueError("k must be a positive integer")
        self._lseries = lseries
        self._k = k

    def __eq__(self, right):
        return (self._lseries, self._k) == (right._lseries, right._k)

    def __ne__(self, right):
        return not(self == right)

    def __call__(self, s):
        """
        Return the value of this derivative at s, which must coerce to a
        complex number.  The computation is to the same precision as s,
        or to 53 bits of precision if s is exact.

        As usual, if s has large imaginary part, then there could be
        substantial precision loss (a warning is printed in that case).

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: L = LSeries(EllipticCurve('389a'))
            sage: L1 = L.derivative(1); L1(1)
            -1.94715429754927e-20
            sage: L1(2)
            0.436337613850735
            sage: L1(I)
            -19.8890471908356 + 31.2633280771869*I
            sage: L2 = L.derivative(2); L2(1)
            1.51863300057685
            sage: L2(I)
            134.536162459604 - 62.6542402272310*I
        """
        return self._lseries._function(prec(s)).derivative(s, self._k)

    def __repr__(self):
        """
        String representation of this derivative of an L-series.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: L = LSeries('delta')
            sage: L.derivative(1).__repr__()
            "First derivative of L-function associated to Ramanujan's Delta (a weight 12 cusp form)"
            sage: L.derivative(2).__repr__()
            "Second derivative of L-function associated to Ramanujan's Delta (a weight 12 cusp form)"
            sage: L.derivative(3).__repr__()
            "Third derivative of L-function associated to Ramanujan's Delta (a weight 12 cusp form)"
            sage: L.derivative(4).__repr__()
            "4-th derivative of L-function associated to Ramanujan's Delta (a weight 12 cusp form)"
            sage: L.derivative(2011).__repr__()
            "2011-th derivative of L-function associated to Ramanujan's Delta (a weight 12 cusp form)"
        """
        k = self._k
        if k == 1:
            kth = 'First'
        elif k == 2:
            kth = 'Second'
        elif k == 3:
            kth = 'Third'
        else:
            kth = '%s-th' % k
        return "%s derivative of %s" % (kth, self._lseries)

    def derivative(self, k=1):
        """
        Return the k-th derivative of this derivative object.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: f = Newforms(43,2,names='a')[1]; f
            q + a1*q^2 - a1*q^3 + (-a1 + 2)*q^5 + O(q^6)
            sage: L = LSeries(f); L1 = L.derivative()
            sage: L1(1)
            0.331674007376949
            sage: L(1)
            0.620539857407845
            sage: L = LSeries(f); L1 = L.derivative(); L1
            First derivative of L-series of a degree 2 newform of level 43
            and weight 2
            sage: L1.derivative()
            Second derivative of L-series of a degree 2 newform of level 43
            and weight 2
            sage: L1.derivative(3)
            4-th derivative of L-series of a degree 2 newform of level 43
            and weight 2
        """
        if k == 0:
            return self
        return LSeriesDerivative(self._lseries, self._k + k)


class LSeriesParentClass(Parent):
    def __init__(self):
        Parent.__init__(self, category=Monoids().Commutative())

    def __contains__(self, x):
        return isinstance(x, (LSeriesAbstract, LSeriesProduct))

    def __call__(self, x):
        """
        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: L = LSeries('zeta')
            sage: P = L.parent(); P
            All L-series objects and their products
            sage: P(L) is L
            True

            sage: P(L^3)
            (Riemann Zeta function viewed as an L-series)^3
            sage: (L^3).parent()
            All L-series objects and their products

        You can make the L-series attached to an object by coercing
        that object into the parent::

            sage: P(EllipticCurve('11a'))
            L-series of Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field

        We also allow coercing in 1, because this is very useful for
        formal factorizations and products::

            sage: P(1)
            1

        Other numbers do not coerce in, of course::

            sage: P(2)
            Traceback (most recent call last):
            ...
            TypeError
        """
        if isinstance(x, LSeriesAbstract):
            return x
        elif isinstance(x, LSeriesProduct):
            return x
        elif x == 1:
            return x
        else:
            try:
                return LSeries(x)
            except NotImplementedError:
                raise TypeError

    def __repr__(self):
        """
        Return string representation of this parent object.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeriesParent
            sage: LSeriesParent.__repr__()
            'All L-series objects and their products'
        """
        return "All L-series objects and their products"

    def __eq__(self, right):
        """
        Return equality if right is an instance of
        LSeriesParentClass.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeriesParent, LSeriesParentClass
            sage: LSeriesParent == LSeriesParentClass()
            True
        """
        if isinstance(right, LSeriesParentClass):
            return True
        return False


LSeriesParent = LSeriesParentClass()


class LSeriesAbstract(Element):
    r"""
    L-series defined by an Euler product.

    The parameters that define the 'shape' of the L-series are:

             conductor, hodge_numbers, weight, epsilon, poles, base_field

    Let gamma(s) = prod( Gamma((s+h)/2) for h in hodge_numbers ).
    Denote this L-series by L(s), and let `L^*(s) = A^s \gamma(s) L(s)`, where
    `A = \sqrt(N)/\pi^{d/2}`, where d = len(hodge_numbers) and N = conductor.
    Then the functional equation is

                  Lstar(s) = epsilon * Lstar(weight - s).

    To actually use this class we create a derived class that in
    addition implements a method _local_factor(P), that takes as input
    a prime ideal P of K=base_field, and returns a polynomial, which
    is typically the reversed characteristic polynomial of Frobenius
    at P of Gal(Kbar/K) acting on the maximal unramified quotient of
    some Galois representation.  This class automatically computes the
    Dirichlet series coefficients `a_n` from the local factors of the
    `L`-function.

    The derived class may optionally -- and purely as an optimization
    -- define a method self._precompute_local_factors(bound,
    prec=oo), which is typically called before

        [_local_factor(P) for P with norm(P) < bound]

    is called in the course of various computations.
    """
    def __init__(self,
                 conductor,
                 hodge_numbers,
                 weight,
                 epsilon,
                 poles,
                 residues,
                 base_field,
                 is_selfdual=True,
                 prec=53):
        """
        INPUT:
            - ``conductor`` -- list or number (in a subset of the
              positive real numbers); if the conductor is a list, then
              each conductor is tried in order (to the precision prec
              below) until we find one that works.
            - ``hodge_numbers`` -- list of numbers (in a subring of
              the complex numbers)
            - ``weight`` -- number (in a subset of the positive real numbers)
            - ``epsilon`` -- number (in a subring of the complex numbers)
            - ``poles`` -- list of numbers (in subring of complex
              numbers); poles of the *completed* L-function
            - ``residues`` -- list of residues at each pole given in
              poles or string "automatic"
            - ``base_field`` -- QQ or a number field; local L-factors
              correspond to nonzero prime ideals of this field.
            - ``is_selfdual`` -- bool (default: True)
            - ``prec`` -- integer (default: 53); precision to use when
              trying to figure out parameters using the functional
              equation

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeriesAbstract
            sage: L = LSeriesAbstract(conductor=1, hodge_numbers=[0], weight=1,
            ....:     epsilon=1, poles=[1], residues=[-1], base_field=QQ)
            sage: type(L)
            <class 'sage.lfunctions.eulerprod.LSeriesAbstract'>
            sage: L._conductor
            1
            sage: L._hodge_numbers
            [0]
            sage: L._weight
            1
            sage: L._epsilon
            1
            sage: L._poles
            [1]
            sage: L._residues
            [-1]
            sage: L._base_field
            Rational Field
            sage: L
            Euler Product L-series with conductor 1, Hodge numbers [0],
            weight 1, epsilon 1, poles [1], residues [-1] over Rational Field
        """
        self._anlist = {oo: []}

        (self._conductor, self._hodge_numbers, self._weight, self._epsilon,
         self._poles, self._residues, self._base_field, self._is_selfdual) = (
            conductor, hodge_numbers, weight, epsilon, poles, residues,
            base_field, is_selfdual)

        # the following parameters allow for specifying a list of
        # possibilities:
        #
        #    conductor -- list of possible conductors
        #    hodge_numbers -- give a list of lists
        #    weight -- list of possible weights
        #    poles (residues must be the default 'automatic')
        #    epsilon -- list of possible epsilon's.

        # 1. Figure out for which parameters we have multiple options
        # 2. Run through them until checking each until one is found
        #    that works.
        v = []
        if isinstance(conductor, list):
            v.append('_conductor')
        if hodge_numbers and isinstance(hodge_numbers[0], list):
            v.append('_hodge_numbers')
        if isinstance(weight, list):
            v.append('_weight')
        if poles and isinstance(poles[0], list):
            v.append('_poles')
        if isinstance(epsilon, list):
            v.append('_epsilon')

        if v:
            found_params = False
            for X in cartesian_product_iterator([getattr(self, attr)
                                                 for attr in v]):
                kwds = dict((v[i], X[i]) for i in range(len(v)))
                if self._is_valid_parameters(prec=prec, save=True, **kwds):
                    found_params = True
                    break
            if not found_params:
                raise RuntimeError("no choice of values for %s works" % (', '.join(v)))
        Element.__init__(self, LSeriesParent)

    def _is_valid_parameters(self, prec=53, save=True, **kwds):
        valid = False
        try:
            old = [(k, getattr(self, k)) for k in kwds]
            for k, v in kwds.items():
                setattr(self, k, v)
            self._function.clear_cache()
            self._function(prec=prec)
            try:
                self._function(prec=prec)
                valid = True
            except RuntimeError:
                pass
        finally:
            if not save:
                for k, v in old:
                    setattr(self, k, v)
            return valid

    def __eq__(self, right):
        if self is right:
            return True
        return all(getattr(self, a)() == getattr(right, a)()
                   for a in ['degree', 'weight', 'conductor',
                             'epsilon', 'base_field'])

    def __ne__(self, right):
        return not(self == right)

    def __hash__(self):
        return hash(tuple(getattr(self, a)()
                          for a in ['degree', 'weight', 'conductor',
                                    'epsilon', 'base_field']))

    def parent(self):
        """
        Return parent of this L-series, which is the collection of all
        L-series and their products.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: L = LSeries('delta');  P = L.parent(); P
            All L-series objects and their products
            sage: L in P
            True
        """
        return LSeriesParent

    def __pow__(self, n):
        """
        Return the n-th power of this L-series, where n can be any integer.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: L = LSeries('delta');
            sage: L3 = L^3; L3
            (L-function associated to Ramanujan's Delta
            (a weight 12 cusp form))^3
            sage: L3(1)
            0.0000524870430366548
            sage: L(1)^3
            0.0000524870430366548
            sage: M = L^(-3); M(1)
            19052.3211471761
            sage: L(1)^(-3)
            19052.3211471761

        Higher precision::

            sage: M = L^(-3); M(RealField(100)(1))
            19052.321147176093380952680193
            sage: L(RealField(100)(1))^(-3)
            19052.321147176093380952680193

        Special case -- 0th power -- is not allowed::

            sage: L^0
            Traceback (most recent call last):
            ...
            ValueError: product must be nonempty
        """
        return LSeriesProduct([(self, ZZ(n))])

    def _mul_(self, right):
        """
        Multiply two L-series, or an L-series times a formal product
        of L-series.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: d = LSeries('delta'); z = LSeries('zeta')
            sage: d * z
            (Riemann Zeta function viewed as an L-series) *
            (L-function associated to Ramanujan's Delta (a weight 12 cusp form))
            sage: d * (d * z)
            (Riemann Zeta function viewed as an L-series) *
            (L-function associated to Ramanujan's Delta (a weight 12 cusp form))^2
        """
        if isinstance(right, LSeriesAbstract):
            return LSeriesProduct([(self, 1), (right, 1)])
        elif isinstance(right, LSeriesProduct):
            return right * self
        raise TypeError

    def _div_(self, right):
        """
        Divide two L-series or formal L-series products.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: d = LSeries('delta'); z = LSeries('zeta')
            sage: d / z
            (Riemann Zeta function viewed as an L-series)^-1 *
            (L-function associated to Ramanujan's Delta (a weight 12 cusp form))
            sage: d / (z^3)
            (Riemann Zeta function viewed as an L-series)^-3 *
            (L-function associated to Ramanujan's Delta (a weight 12 cusp form))
        """
        if isinstance(right, LSeriesAbstract):
            return LSeriesProduct([(self, 1), (right, -1)])
        elif isinstance(right, LSeriesProduct):
            F = Factorization([(self, 1)]) / right._factorization
            return LSeriesProduct(F)
        raise TypeError

    def conductor(self):
        """
        Return the conductor of this L-series.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: LSeries('zeta').conductor()
            1
            sage: LSeries('delta').conductor()
            1
            sage: LSeries(EllipticCurve('11a')).conductor()
            11
            sage: LSeries(Newforms(33)[0]).conductor()
            33
            sage: LSeries(DirichletGroup(37).0).conductor()
            37
            sage: LSeries(kronecker_character(7)).conductor()
            28
            sage: kronecker_character(7).conductor()
            28
            sage: L = LSeries(EllipticCurve('11a3').base_extend(QQ[sqrt(2)]), prec=5)
            sage: L.conductor().factor()
            2^6 * 11^2
        """
        return self._conductor

    def hodge_numbers(self):
        """
        Return the Hodge numbers of this L-series.

        These define the local Gamma factors.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: LSeries('zeta').hodge_numbers()
            [0]
            sage: LSeries('delta').hodge_numbers()
            [0, 1]
            sage: LSeries(EllipticCurve('11a')).hodge_numbers()
            [0, 1]
            sage: LSeries(Newforms(43,names='a')[1]).hodge_numbers()
            [0, 1]
            sage: LSeries(DirichletGroup(37).0).hodge_numbers()
            [1]
            sage: LSeries(EllipticCurve(QQ[sqrt(-1)],[1,2]), prec=5).hodge_numbers() # long time
            [0, 0, 1, 1]
        """
        return list(self._hodge_numbers)

    def weight(self):
        """
        Return the weight of this L-series.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: LSeries('zeta').weight()
            1
            sage: LSeries('delta').weight()
            12
            sage: LSeries(EllipticCurve('389a')).weight()
            2
            sage: LSeries(Newforms(43,names='a')[1]).weight()
            2
            sage: LSeries(Newforms(6,4)[0]).weight()
            4
            sage: LSeries(DirichletGroup(37).0).weight()
            1
            sage: L = LSeries(EllipticCurve('11a3').base_extend(QQ[sqrt(2)]),prec=5); L.weight()
            2
        """
        return self._weight

    def poles(self):
        """
        Poles of the *completed* L-function with the extra Gamma
        factors included.

        WARNING: These are not just the poles of ``self``.

        OUTPUT:

        - list of numbers

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: LSeries('zeta').poles()
            [0, 1]
            sage: LSeries('delta').poles()
            []
        """
        return list(self._poles)

    def residues(self, prec=oo):
        """
        Residues of the *completed* L-function at each pole.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: LSeries('zeta').residues()
            [9, 8]
            sage: LSeries('delta').residues()
            []
            sage: v = LSeries('zeta').residues()
            sage: v.append(10)
            sage: LSeries('zeta').residues()
            [9, 8]

        The residues of the Dedekind Zeta function of a field are
        dynamically computed::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: L = LSeries(K); L
            Dedekind Zeta function of Number Field in a
            with defining polynomial x^2 + 1

        If you just call residues you get back that they are
        automatically computed::

            sage: L.residues()
            'automatic'

        But if you call with a specific precision, they are computed
        using that precision::

            sage: L.residues(prec=53)
            [-0.886226925452758]
            sage: L.residues(prec=200)
            [-0.88622692545275801364908374167057259139877472806119356410690]
        """
        if self._residues == 'automatic':
            if prec is oo:
                return self._residues
            else:
                C = ComplexField(prec)
                self._function(prec=prec).__check_init()
                return [C(a)
                        for a in self._function(prec=prec)._gp_call_inst('Lresidues')]
        else:
            return list(self._residues)

    def base_field(self):
        """
        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: LSeries('zeta').base_field()
            Rational Field
            sage: L = LSeries(EllipticCurve('11a3').base_extend(QQ[sqrt(2)]), prec=5); L.base_field()
            Number Field in sqrt2 with defining polynomial x^2 - 2
        """
        return self._base_field

    def epsilon(self, prec=oo):
        """
        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: LSeries('zeta').epsilon()
            1
            sage: LSeries('delta').epsilon()
            1
            sage: LSeries(EllipticCurve('389a')).epsilon()
            1
            sage: LSeries(EllipticCurve('37a')).epsilon()
            -1
            sage: LSeries(Newforms(389,names='a')[1]).epsilon()
            -1
            sage: LSeries(Newforms(6,4)[0]).epsilon()
            1
            sage: LSeries(DirichletGroup(7).0).epsilon()
            1/7*I*sqrt(7)*((((((e^(2/21*I*pi) + 1)*e^(2/21*I*pi) + 1)*e^(1/21*I*pi) - 1)*e^(1/21*I*pi) - 1)*e^(2/21*I*pi) - 1)*e^(1/21*I*pi) - 1)*e^(1/21*I*pi)
            sage: LSeries(DirichletGroup(7).0).epsilon(prec=53)
            0.386513572759156 + 0.922283718859307*I

        In this example, the epsilon factor is computed when the curve
        is created. The prec parameter determines the floating point precision
        used in computing the epsilon factor::

            sage: L = LSeries(EllipticCurve('11a3').base_extend(QQ[sqrt(2)]), prec=5); L.epsilon()
            -1

        Here is extra confirmation that the rank is really odd over
        the quadratic field::

            sage: EllipticCurve('11a').quadratic_twist(2).rank()
            1

        We can compute with the L-series too::

            sage: L(RealField(5)(2))
            0.53

        For L-functions of newforms with nontrivial character, the
        epsilon factor is harder to find (we do not have a good
        algorithm implemented to find it) and might not even be 1 or
        -1, so it is set to 'solve'.  In this case, the functional
        equation is used to determine the solution.::

            sage: f = Newforms(kronecker_character_upside_down(7),3)[0]; f
            q - 3*q^2 + 5*q^4 + O(q^6)
            sage: L = LSeries(f)
            sage: L.epsilon()
            'solve'
            sage: L(3/2)
            0.332981771482934
            sage: L.epsilon()
            1

        Here is an example with nontrivial character::

            sage: f = Newforms(DirichletGroup(7).0, 5, names='a')[0]; f
            q + a0*q^2 + ((zeta6 - 2)*a0 - zeta6 - 1)*q^3 + (-4*zeta6*a0 + 2*zeta6 - 2)*q^4 + ((4*zeta6 - 2)*a0 + 9*zeta6 - 18)*q^5 + O(q^6)
            sage: L = LSeries(f)

        First trying to evaluate with the default (53 bits) of
        precision fails, since for some reason (that I do not
        understand) the program is unable to find a valid epsilon
        factor::

            sage: L(0)
            Traceback (most recent call last):
            ...
            RuntimeError: unable to determine epsilon from functional equation working to precision 53, since we get epsilon=0.806362085925390 - 0.00491051026156280*I, which is not sufficiently close to 1

        However, when we evaluate to 100 bits of precision it works::

            sage: L(RealField(100)(0))
            0.00000000000000000000000000000

        The epsilon factor is *not* known to infinite precision::

            sage: L.epsilon()
            'solve'

        But it is now known to 100 bits of precision, and here it is::

            sage: L.epsilon(100)
            0.42563106101692403875896879406 - 0.90489678963824790765479396740*I

        When we try to compute to higher precision, again Sage solves
        for the epsilon factor numerically::

            sage: L(RealField(150)(1))
            0.26128389551787271923496480408992971337929665 - 0.29870133769674001421149135036267324347896657*I

        And now it is known to 150 bits of precision.  Notice that
        this is consistent with the value found above, and has
        absolute value (very close to) 1.

            sage: L.epsilon(150)
            0.42563106101692403875896879406038776338921622 - 0.90489678963824790765479396740501409301704122*I
            sage: abs(L.epsilon(150))
            1.0000000000000000000000000000000000000000000
        """
        if self._epsilon == 'solve':
            return 'solve'
        if hasattr(self._epsilon, 'prec') and self._epsilon.prec() < prec:
            return 'solve'
        if prec is not oo:
            C = ComplexField(prec)
            if isinstance(self._epsilon, list):
                return [C(x) for x in self._epsilon]
            return C(self._epsilon)
        return self._epsilon

    def is_selfdual(self):
        """
        Return whether this L-series is self dual.

        EXAMPLES:

        Many L-series are self dual::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: LSeries('zeta').is_selfdual()
            True
            sage: LSeries('delta').is_selfdual()
            True
            sage: LSeries(Newforms(6,4)[0]).is_selfdual()
            True

        Nonquadratic characters have non-self dual L-series::

            sage: LSeries(DirichletGroup(7).0).is_selfdual()
            False
            sage: LSeries(kronecker_character(7)).is_selfdual()
            True

        Newforms with non-quadratic characters also have non-self dual L-seris::

            sage: L = LSeries(Newforms(DirichletGroup(7).0, 5, names='a')[0])
            sage: L.is_selfdual()
            False
        """
        return self._is_selfdual

    @cached_method
    def degree(self):
        """
        Return the degree of this L-function.

        This is by definition the number of Gamma factors (e.g., the
        number of Hodge numbers) divided by the degree of the base
        field.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: LSeries('zeta').degree()
            1
            sage: LSeries(DirichletGroup(5).0).degree()
            1
            sage: LSeries(EllipticCurve('11a')).degree()
            2

        The L-series attached to this modular symbols space of
        dimension 2 is a product of 2 degree 2 L-series, hence has
        degree 4::

            sage: M = ModularSymbols(43,2,sign=1).cuspidal_subspace()[1]; M.dimension()
            2
            sage: L = LSeries(M); L
            L-series attached to Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field
            sage: L.factor()
            (L-series of a degree 2 newform of level 43 and weight 2) *
            (L-series of a degree 2 newform of level 43 and weight 2)
            sage: L.degree()
            4

            sage: x = var('x'); K.<a> = NumberField(x^2-x-1)
            sage: LSeries(EllipticCurve([0,-a,a,0,0])).degree()
            2
        """
        n = len(self.hodge_numbers())
        d = self.base_field().degree()
        assert n % d == 0, "degree of base field must divide the number of Hodge numbers"
        return n // d

    def twist(self, chi, conductor=None, epsilon=None, prec=53):
        r"""
        Return the quadratic twist of this L-series by the character chi, which
        must be a character of self.base_field().

        Thus chi should take as input prime ideals (or primes) of the
        ring of integers of the base field, and output something that
        can be coerced to the complex numbers.

        INPUT:

        - `\chi` -- 1-dimensional character

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: E = EllipticCurve('11a')
            sage: L = LSeries(E); L
            L-series of Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: L3 = L.twist(DirichletGroup(3).0); L3
            Twist of L-series of Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field by Dirichlet character modulo 3 of conductor 3 mapping 2 |--> -1
            sage: L3._chi
            Dirichlet character modulo 3 of conductor 3 mapping 2 |--> -1
            sage: L3._L
            L-series of Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: L3(1)
            1.68449633297548
            sage: F = E.quadratic_twist(-3)
            sage: L3.conductor()
            99
            sage: F.conductor()
            99
            sage: F.lseries()(1)
            1.68449633297548
            sage: L3.anlist(20)
            [0, 1, 2, 0, 2, -1, 0, -2, 0, 0, -2, -1, 0, 4, -4, 0, -4, 2, 0, 0, -2]
            sage: F.anlist(20)
            [0, 1, 2, 0, 2, -1, 0, -2, 0, 0, -2, -1, 0, 4, -4, 0, -4, 2, 0, 0, -2]
            sage: L3.anlist(1000) == F.anlist(1000)
            True
            sage: L3.local_factor(11)
            T + 1

        A higher degree twist::

            sage: L = LSeries(EllipticCurve('11a'))
            sage: L5 = L.twist(DirichletGroup(5).0); L5
            Twist of L-series of Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field by Dirichlet character modulo 5 of conductor 5 mapping 2 |--> zeta4
            sage: L5(1)
            1.28009593569230 - 0.681843202124309*I
            sage: L5.epsilon()
            'solve'
            sage: L5.epsilon(53)
            0.989989082587826 + 0.141143956147310*I
            sage: L5.conductor()
            275
            sage: L5.taylor_series(center=1, degree=3)
            1.28009593569230 - 0.681843202124309*I + (-0.536450338806282 + 0.166075270978779*I)*z + (0.123743053129226 + 0.320802890011298*I)*z^2 + O(z^3)


        WARNING!! Twisting is not implemented in full generality when
        the conductors are not coprime.  One case where we run into
        trouble is when twisting lowers the level of a newform.  Below
        we take the form of level 11 and weight 2, twist it by the
        character chi of conductor 3 to get a form of level 99.  Then
        we take the L-series of the level 99 form, and twist that by
        chi, which should be the L-series attached to the form of
        level 11.  Unfortunately, our code for working out the local
        L-factors doesn't succeed in this case, hence the local factor
        is wrong, so the functional equation is not satisfied.

            sage: f = Newform('11a'); f
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6)
            sage: L = LSeries(f)
            sage: chi = DirichletGroup(3).0
            sage: Lc = L.twist(chi); Lc
            Twist of L-series of a degree 1 newform of level 11 and weight 2 by Dirichlet character modulo 3 of conductor 3 mapping 2 |--> -1
            sage: Lc.anlist(20)
            [0, 1, 2, 0, 2, -1, 0, -2, 0, 0, -2, -1, 0, 4, -4, 0, -4, 2, 0, 0, -2]
            sage: g = Newform('99d'); g
            q + 2*q^2 + 2*q^4 - q^5 + O(q^6)
            sage: list(g.qexp(20))
            [0, 1, 2, 0, 2, -1, 0, -2, 0, 0, -2, -1, 0, 4, -4, 0, -4, 2]
            sage: Lt = Lc.twist(chi, conductor=11)
            Traceback (most recent call last):
            ...
            RuntimeError: no choice of values for _epsilon works
            sage: Lt = Lc.twist(chi,conductor=11, epsilon=1)
            sage: Lt(1)
            Traceback (most recent call last):
            ...
            RuntimeError: invalid L-series parameters: functional equation not satisfied

        This is because the local factor is wrong::

            sage: Lt.local_factor(3)
            1
            sage: L.local_factor(3)
            3*T^2 + T + 1
        """
        return LSeriesTwist(self, chi=chi, conductor=conductor,
                            epsilon=epsilon, prec=prec)

    @cached_method
    def local_factor(self, P, prec=oo):
        """
        Return the local factor of the L-function at the prime P of
        self._base_field.  The result is cached.

        INPUT:

        - a prime P of the ring of integers of the base_field
        - prec -- oo or positive integer (bits of precision)

        OUTPUT:

        - a polynomial, e.g., something like "1-a*T+p*T^2".

        EXAMPLES:

        You must overload this in the derived class::

            sage: from sage.lfunctions.eulerprod import LSeriesAbstract
            sage: L = LSeriesAbstract(conductor=1, hodge_numbers=[0], weight=1, epsilon=1, poles=[1], residues=[-1], base_field=QQ)
            sage: L.local_factor(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: must be implemented in the derived class
        """
        return self._local_factor(P, prec=prec)

    def _local_factor(self, P, prec=oo):
        """
        Compute local factor at prime P.  This must be overwritten in the derived class.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeriesAbstract
            sage: L = LSeriesAbstract(conductor=1, hodge_numbers=[0], weight=1, epsilon=1, poles=[1], residues=[-1], base_field=QQ)
            sage: L.local_factor(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: must be implemented in the derived class
        """
        raise NotImplementedError("must be implemented in the derived class")

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeriesAbstract
            sage: L = LSeriesAbstract(conductor=1, hodge_numbers=[0], weight=1, epsilon=1, poles=[1], residues=[-1], base_field=QQ)
            sage: L.__repr__()
            'Euler Product L-series with conductor 1, Hodge numbers [0], weight 1, epsilon 1, poles [1], residues [-1] over Rational Field'
        """
        return "Euler Product L-series with conductor %s, Hodge numbers %s, weight %s, epsilon %s, poles %s, residues %s over %s" % (
            self._conductor, self._hodge_numbers, self._weight,
            self.epsilon(), self._poles, self._residues, self._base_field)

    def _precompute_local_factors(self, bound, prec):
        """
        Derived classes may use this as a 'hint' that _local_factors
        will soon get called for primes of norm less than the bound.

        In the base class, this is a no-op, and it is not necessary to
        overload this class.

        INPUT:

        - ``bound`` -- integer
        - ``prec`` -- integer

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeriesAbstract
            sage: L = LSeriesAbstract(conductor=1, hodge_numbers=[0], weight=1, epsilon=1, poles=[1], residues=[-1], base_field=QQ)
            sage: L._precompute_local_factors(100, 53)
        """
        pass

    def _primes_above(self, p):
        """
        Return the primes of the ring of integers of the base field above the integer p.

        INPUT:

        - p -- prime integer (no type checking necessarily done)

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeriesAbstract
            sage: L = LSeriesAbstract(conductor=1, hodge_numbers=[0], weight=1, epsilon=1, poles=[1], residues=[-1], base_field=QQ)
            sage: L._primes_above(3)
            [3]
            sage: L = LSeriesAbstract(conductor=1, hodge_numbers=[0], weight=1, epsilon=1, poles=[1], residues=[-1], base_field=QQ[sqrt(-1)])
            sage: L._primes_above(5)
            [Fractional ideal (-I - 2), Fractional ideal (I - 2)]
            sage: L._primes_above(3)
            [Fractional ideal (3)]
        """
        K = self._base_field
        if is_RationalField(K):
            return [p]
        else:
            return K.primes_above(p)

    def zeros(self, n):
        """
        Return the imaginary parts of the first n nontrivial zeros of
        this L-function on the critical line in the upper half plane,
        as 32-bit real numbers.

        INPUT:

        - n -- nonnegative integer

        EXAMPLES::

        """
        return self._lcalc().zeros(n)

    def _lcalc(self):
        """
        Return Rubinstein Lcalc object attached to this L-series.

        This is useful both for evaluating the L-series, especially
        when the imaginary part is large, and for computing the zeros
        in the critical strip.

        EXAMPLES::
        """
        # this will require using Rubinstein's L-calc
        raise NotImplementedError

    def anlist(self, bound, prec=oo):
        """
        Return list `v` of Dirichlet series coefficients `a_n` for `n`
        up to and including ``bound``, where `v[n] = a_n`.

        If at least one of the coefficients is ambiguous, e.g., the
        :meth:`local_factor` method returns a list of possibilities,
        then this function instead returns a generator over all
        possible `a_n` lists.

        In particular, we include v[0]=0 as a convenient place holder
        to avoid having `v[n-1] = a_n`.

        INPUT:

        - ``bound`` -- nonnegative integer

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries; L = LSeries('zeta')
            sage: L.anlist(30)
            [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: from sage.lfunctions.eulerprod import LSeries; L = LSeries(EllipticCurve('11a'))
            sage: L.anlist(30)
            [0, 1, -2, -1, 2, 1, 2, -2, 0, -2, -2, 1, -2, 4, 4, -1, -4, -2, 4, 0, 2, 2, -2, -1, 0, -4, -8, 5, -4, 0, 2]
            sage: K.<a> = NumberField(x^2-x-1); L = LSeries(EllipticCurve([0,-a,a,0,0]))
            sage: L.anlist(30)
            [0, 1, 0, 0, -2, -1, 0, 0, 0, -4, 0, 3, 0, 0, 0, 0, 0, 0, 0, 5, 2, 0, 0, 0, 0, -4, 0, 0, 0, 11, 0]
        """
        # First check if we know anlist to infinite bit precision up
        # to given bound:
        if len(self._anlist[oo]) > bound:
            if prec is oo:
                # request it to infinite precision
                return self._anlist[oo][:bound + 1]
            else:
                # finite precision request
                C = ComplexField(prec)
                return [C(a) for a in self._anlist[oo]]

        if prec is not oo:
            # check if numerically computed already to at least this precision
            t = [(precision, cfs) for precision, cfs in self._anlist.items()
                 if precision >= prec and len(cfs) > bound]
            if t:
                C = ComplexField(prec)
                return [C(a) for a in t[0][1][:bound + 1]]

        self._precompute_local_factors(bound + 1, prec=prec)
        compute_anlist_multiple = False
        LF = []
        for p in prime_range(bound + 1):
            lf = []
            for P in self._primes_above(p):
                if norm(P) <= bound:
                    F = self._local_factor(P, prec)
                    if isinstance(F, list):
                        compute_anlist_multiple = True
                    lf.append(F)
            LF.append((p, lf))

        coefficients = self._compute_anlist(LF, bound, prec)
        if not compute_anlist_multiple:
            coefficients = list(coefficients)[0]
            # save answer in cache
            self._anlist[prec] = coefficients
        return coefficients

    def _compute_anlist(self, LF, bound, prec):
        """
        Iterator over possible anlists, given LF, bound, and prec.

        INPUT:

        - ``LF`` -- list of pairs (p, [local factors (or lists of them) at primes over p])
        - ``bound`` -- positive integer
        - ``prec`` -- positive integer (bits of precision)
        """
        coefficients = [0, 1] + [0] * (bound - 1)

        for i, (p, v) in enumerate(LF):
            if v:
                # Check for the possibility of several different
                # choices of Euler factor at a given prime.  If this
                # happens, we switch gears.
                some_list = False
                for j, lf in enumerate(v):
                    if isinstance(lf, list):
                        some_list = True
                        for f in list(lf):
                            LF0 = copy.deepcopy(LF)
                            LF0[i][1][j] = f
                            for z in self._compute_anlist(LF0, bound, prec):
                                yield z
                if some_list:
                    return
                # Not several factors -- just compute the a_{p^r} up
                # to the required bound:
                f = prod(v)
                accuracy_p = int(math.floor(math.log(bound) / math.log(p))) + 1
                T = f.parent().gen()
                series_p = (f + O(T**accuracy_p))**(-1)
                for j in range(1, accuracy_p):
                    coefficients[p**j] = series_p[j]

        # fill in non-prime power coefficients
        extend_multiplicatively_generic(coefficients)
        yield list(coefficients)

    def _symbolic_(self, R, bound=10, prec=oo):
        """
        Convert self into the symbolic ring as a truncated Dirichlet
        series, including terms up to `n^s` where n=bound.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: L = LSeries(EllipticCurve('37a'))
            sage: SR(L)
            4/10^s + 6/9^s - 1/7^s + 6/6^s - 2/5^s + 2/4^s - 3/3^s - 2/2^s + 1
            sage: L._symbolic_(SR, 20)
            -4/20^s - 12/18^s - 4/16^s + 6/15^s + 2/14^s - 2/13^s - 6/12^s - 5/11^s + 4/10^s + 6/9^s - 1/7^s + 6/6^s - 2/5^s + 2/4^s - 3/3^s - 2/2^s + 1
        """
        s = R.var('s')
        a = self.anlist(bound, prec)
        return sum(a[n] / n**s for n in range(1, bound + 1))

    def __call__(self, s):
        """
        Return value of this L-function at s.

        If s is a real or complex number to prec bits of precision,
        then the result is also computed to prec bits of precision.

        If s has infinite or unknown precision, then the L-value is
        computed to 53 bits of precision.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: L = LSeries(EllipticCurve('37a'))
            sage: L(1)
            0.000000000000000
            sage: L(2)
            0.381575408260711
            sage: z = L(RealField(100)(2)); z
            0.38157540826071121129371040958
            sage: z.prec()
            100
            sage: L(RealField(150)(2))
            0.38157540826071121129371040958008663667709753

        WARNING: There will be precision loss (with a warning) if the imaginary
        part of the input is large::

            sage: L = LSeries('zeta')
            sage: L(1/2 + 40*I)
            verbose -1 (...: dokchitser.py, __call__) Warning: Loss of 14 decimal digits due to cancellation
            0.793046013671137 - 1.04127377821427*I
            sage: L(ComplexField(200)(1/2 + 40*I))
            verbose -1 (...: dokchitser.py, __call__) Warning: Loss of 14 decimal digits due to cancellation
            0.79304495256192867196489258889793696080572220439302833315881 - 1.0412746146510650200518905953910554313275550685861559488384*I

        An example with a pole::

            sage: L = LSeries('zeta')
            sage: L(2)  # good
            1.64493406684823
            sage: L(1)  # a pole!
            Traceback (most recent call last):
            ...
            ZeroDivisionError: pole at 1
        """
        if s in self._poles:
            raise ZeroDivisionError("pole at %s" % s)
        return self._function(prec(s))(s)

    def derivative(self, k=1):
        """
        Return the k-th derivative of self.

        INPUT:

        - k -- (default: 1) nonnegative integer

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: L = LSeries('zeta')
            sage: L.derivative()
            First derivative of Riemann Zeta function viewed as an L-series

        We numerically approximate the derivative at two points and compare
        with evaluating the derivative object::

            sage: eps=1e-10;  (L(2+eps) - L(2))/eps
            -0.937547817159157
            sage: Lp = L.derivative(); Lp(2)
            -0.937548254315844
            sage: eps=1e-10;  (L(2+I+eps) - L(2+I))/eps
            0.0624900131640516 + 0.489033813444451*I
            sage: Lp(2+I)
            0.0624900021906470 + 0.489033591679899*I

        Higher derivatives::

            sage: L.derivative(2)
            Second derivative of Riemann Zeta function viewed as an L-series
            sage: L.derivative(2)(2)
            1.98928023429890
            sage: L.derivative(3)
            Third derivative of Riemann Zeta function viewed as an L-series
            sage: L.derivative(4)
            4-th derivative of Riemann Zeta function viewed as an L-series
            sage: L.derivative(5)
            5-th derivative of Riemann Zeta function viewed as an L-series

        Derivative of derivative::

            sage: L.derivative().derivative()
            Second derivative of Riemann Zeta function viewed as an L-series

        Using the derivative function in Sage works::

            sage: derivative(L)
            First derivative of Riemann Zeta function viewed as an L-series
            sage: derivative(L,2)
            Second derivative of Riemann Zeta function viewed as an L-series
        """
        if k == 0:
            return self
        return LSeriesDerivative(self, k)

    def taylor_series(self, center=None, degree=6, variable='z', prec=53):
        """
        Return the Taylor series expansion of self about the given
        center point to the given degree in the specified variable
        numerically computed to the precision prec in bits.  If the
        center is not specified it defaults to weight / 2.

        INPUT:

        - ``center`` -- None or number that coerces to the complex numbers
        - ``degree`` -- integer
        - ``variable`` -- string or symbolic variable
        - ``prec`` -- positive integer (floating point bits of precision)

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries; L = LSeries('zeta')
            sage: L.taylor_series()
            -1.46035450880959 - 3.92264613920915*z - 8.00417850696433*z^2 - 16.0005515408865*z^3 - 31.9998883216853*z^4 - 64.0000050055172*z^5 + O(z^6)
            sage: RealField(53)(zeta(1/2))
            -1.46035450880959
            sage: L.taylor_series(center=2, degree=4, variable='t', prec=30)
            1.6449341 - 0.93754825*t + 0.99464012*t^2 - 1.0000243*t^3 + O(t^4)
            sage: RealField(30)(zeta(2))
            1.6449341

        """
        if center is None:
            center = ComplexField(prec)(self._weight) / 2
        return self._function(prec).taylor_series(center, degree, variable)

    def analytic_rank(self, tiny=1e-8, prec=53):
        center = ComplexField(prec)(self._weight) / 2
        degree = 4
        while True:
            f = self.taylor_series(center, degree, prec=prec)
            i = 0
            while i < degree:
                if abs(f[i]) > tiny:
                    return i
                i += 1
            degree += 2

    @cached_method
    def _function(self, prec=53, T=1.2):
        """
        Return Dokchitser object that allows for computation of this
        L-series computed to enough terms so that the functional
        equation checks out with the given value of T and precision.

        This is used behind the scenes for evaluation and computation
        of Taylor series.
        """
        eps = self.epsilon(prec)
        return self._dokchitser(prec, eps, T=T)

    def _dokchitser_unitialized(self, prec, epsilon):
        # Create the Dokchitser object
        if epsilon == 'solve':
            eps = 'X'
        else:
            eps = epsilon
        return Dokchitser(conductor=self.conductor(),
                          gammaV=self.hodge_numbers(),
                          weight=self.weight(),
                          eps=eps,
                          poles=self.poles(),
                          residues=self.residues(),
                          prec=prec)

    def number_of_coefficients(self, prec=53, T=1.2):
        """
        Return the number of Dirichlet series coefficients that will
        be needed in order to evaluate this L-series (near the real
        line) to prec bits of precision and have the functional
        equation test pass with the given value of T.

        INPUT:

        - prec -- integer

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: L = LSeries(DirichletGroup(5).0)
            sage: L.number_of_coefficients(20)
            8
            sage: L.number_of_coefficients()
            11
            sage: L.number_of_coefficients(1000)
            43
        """
        return self._dokchitser_unitialized(prec, self.epsilon(prec)).num_coeffs(T)

    def _dokchitser(self, prec, epsilon, T=1.2):
        L = self._dokchitser_unitialized(prec, epsilon)

        # Find out how many coefficients of the Dirichlet series are needed
        # to compute to the requested precision.
        n = L.num_coeffs(T=T)
        # if n >= 500:   # TODO: for debugging only -- remove later
        #    print "num coeffs =", n

        # Compute the Dirichlet series coefficients
        X = self.anlist(n, prec)
        if isinstance(X, types.GeneratorType):
            # Several possible coefficients -- we try them until finding on that works.
            coeff_lists = X
        else:
            # Only one to try
            coeff_lists = [X]

        tiny0 = tiny(prec)
        for coeffs in coeff_lists:
            # Define a string that when evaluated in PARI defines a function
            # a(k), which returns the Dirichlet coefficient a_k.
            s = 'v=[%s]; a(k)=v[k];' % ','.join(str(z) if isinstance(z, (int, Integer)) else z._pari_init_() for z in coeffs[1:])

            # Tell the L-series / PARI about the coefficients.

            if self.is_selfdual():
                L.init_coeffs('a(k)', pari_precode=s)
            else:
                # Have to directly call gp_eval, since case of functional equation having two different
                # (conjugate) L-functions isn't supported in Dokchitser class (yet).
                L._Dokchitser__init = True
                L._gp_eval(s)
                L._gp_call_inst('initLdata', '"a(k)",1.3,"conj(a(k))"')

            if epsilon == 'solve':
                epsvec = L._gp_call_inst("checkfeq")
                L._gp_eval('sgneq = Vec(%s)' % epsvec)
                epsilon = ComplexField(prec)(L._gp_set_inst('sgn', '-sgneq[2]/sgneq[1]'))
                if abs(abs(epsilon) - 1) > tiny0:
                    raise RuntimeError("unable to determine epsilon from functional equation working to precision %s, since we get epsilon=%s, which is not sufficiently close to 1" % (prec, epsilon))
                # 1, -1 are two common special cases, where it is clear what the
                # infinite precision version is.
                if epsilon == 1:
                    self._epsilon = 1
                elif epsilon == -1:
                    self._epsilon = -1
                else:
                    self._epsilon = epsilon

            fe = L.check_functional_equation()
            if abs(fe) <= tiny0:
                # one worked!
                self._anlist[prec] = coeffs
                return L
            else:
                pass

        # They all failed.
        raise RuntimeError("invalid L-series parameters: functional equation not satisfied")

    def check_functional_equation(self, T, prec=53):
        return self._function(prec=prec).check_functional_equation(T)


class LSeriesProductEvaluator(object):
    def __init__(self, factorization, prec):
        self._factorization = factorization
        self._prec = prec

    def __call__(self, s):
        try:
            v = self._functions
        except AttributeError:
            self._functions = [(L._function(self._prec), e)
                               for L, e in self._factorization]
            v = self._functions
        return prod(f(s)**e for f, e in v)


class LSeriesProduct(Element):
    """
    A formal product of L-series.
    """
    def __init__(self, F):
        """
        INPUT:

        - `F` -- list of pairs (L,e) where L is an L-function and e is a nonzero integer.
        """
        if not isinstance(F, Factorization):
            F = Factorization(F, sort=False, simplify=False)
        F.sort(key=lambda L: (L[0]._conductor, L[0]._weight))
        if len(F) == 0:
            raise ValueError("product must be nonempty")
        self._factorization = F
        Element.__init__(self, LSeriesParent)

    def __eq__(self, right):
        # TODO: make work even if right not a product
        return self._factorization == right._factorization

    def __ne__(self, right):
        return not(self == right)

    def is_selfdual(self):
        """
        Return whether every factor of self is self dual.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: L1 = LSeries('zeta'); L1.is_selfdual()
            True
            sage: L2 = LSeries(DirichletGroup(9).0); L2.is_selfdual()
            False
            sage: (L1*L1).is_selfdual()
            True
            sage: (L1*L2).is_selfdual()
            False
            sage: (L2*L2).is_selfdual()
            False
        """
        return all(L.is_selfdual() for L, _ in self._factorization)

    def conductor(self):
        """
        Return the conductor of this product.

        We define the conductor to be the product with multiplicities
        of the conductors of the factors of ``self``.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: J = J0(33); L = LSeries(J)
            sage: L.conductor()
            3993
            sage: L.conductor().factor()
            3 * 11^3
            sage: J.decomposition()
            [
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33),
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33),
            Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)
            ]

        Of course, the conductor as we have defined it need not be an integer::

            sage: L = LSeries(J[0])/LSeries(J[2])^2; L
            (L-series of a degree 1 newform of level 11 and weight 2) * (L-series of a degree 1 newform of level 33 and weight 2)^-2
            sage: L.conductor()
            1/99
        """
        return prod(L.conductor()**e for L, e in self._factorization)

    def __pow__(self, n):
        """
        Return the n-th power of this formal L-series product, where n can be any integer.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: L = LSeries('zeta') * LSeries(DirichletGroup(9).0)
            sage: L(2)
            1.80817853715812 + 0.369298119218816*I
            sage: L = LSeries('zeta') * LSeries(DirichletGroup(9).0)
            sage: L3 = L^3; L3
            (Riemann Zeta function viewed as an L-series)^3 * (L-series attached to Dirichlet character modulo 9 of conductor 9 mapping 2 |--> zeta6)^3
            sage: L3(2)
            5.17205298762567 + 3.57190597873829*I
            sage: CC((zeta(2)*LSeries(DirichletGroup(9).0)(2))^3)
            5.17205298762567 + 3.57190597873829*I
            sage: L^(-1999)
            (Riemann Zeta function viewed as an L-series)^-1999 * (L-series attached to Dirichlet character modulo 9 of conductor 9 mapping 2 |--> zeta6)^-1999
            sage: (L^(-1999)) (2)
            8.90248311986228e-533 - 6.20123089437732e-533*I
        """
        return LSeriesProduct(self._factorization**ZZ(n))

    def _mul_(self, right):
        if isinstance(right, LSeriesAbstract):
            return LSeriesProduct(self._factorization * Factorization([(right, 1)]))
        elif isinstance(right, LSeriesProduct):
            return LSeriesProduct(self._factorization * right._factorization)
        else:
            raise TypeError

    def _div_(self, right):
        if isinstance(right, LSeriesAbstract):
            return LSeriesProduct(self._factorization * Factorization([(right, -1)]))
        elif isinstance(right, LSeriesProduct):
            return LSeriesProduct(self._factorization / right._factorization)
        raise TypeError

    def parent(self):
        return LSeriesParent

    def factor(self):
        return self._factorization

    def hodge_numbers(self):
        """
        Return the Hodge numbers of this product of L-series.
        """

    def degree(self):
        return sum(e * L.degree() for L, e in self._factorization)

    def local_factor(self, P, prec=oo):
        return prod(L.local_factor(P, prec)**e for L, e in self._factorization)

    def __getitem__(self, *args, **kwds):
        return self._factorization.__getitem__(*args, **kwds)

    def __repr__(self):
        return self.factor().__repr__()

    def __call__(self, s):
        return self._function(prec(s))(s)

    def derivative(self, k=1):
        raise NotImplementedError

    def taylor_series(self, center=None, degree=6, variable='z', prec=53):
        """
        Return the Taylor series of self.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: L1 = LSeries('zeta'); L2 = LSeries('delta')
            sage: f = L1 * L2; f
            (Riemann Zeta function viewed as an L-series) * (L-function associated to Ramanujan's Delta (a weight 12 cusp form))
            sage: f.taylor_series(center=2, degree=4, variable='w', prec=30)
            0.24077647 + 0.10066485*w + 0.061553731*w^2 - 0.041923238*w^3 + O(w^4)
            sage: L1.taylor_series(2, 4, 'w', 30) * L2.taylor_series(2, 4, 'w', 30)
            0.24077647 + 0.10066485*w + 0.061553731*w^2 - 0.041923238*w^3 + O(w^4)
            sage: f = L1 / L2; f
            (Riemann Zeta function viewed as an L-series) * (L-function associated to Ramanujan's Delta (a weight 12 cusp form))^-1
            sage: f.taylor_series(center=2, degree=4, variable='w', prec=30)
            11.237843 - 17.508629*w + 21.688182*w^2 - 24.044641*w^3 + O(w^4)
            sage: L1.taylor_series(2, 4, 'w', 30) / L2.taylor_series(2, 4, 'w', 30)
            11.237843 - 17.508629*w + 21.688182*w^2 - 24.044641*w^3 + O(w^4)

        """
        return prod(L.taylor_series(center, degree, variable, prec)**e
                    for L, e in self._factorization)

    def analytic_rank(self, prec=53):
        """
        Return sum of the order of vanishing counted with
        multiplicities of each factor at their center point.

        WARNING: The analytic rank is computed numerically, so is
        definitely not provably correct.

        EXAMPLES:

        We compute the analytic rank of the non-simple 32-dimensional modular abelian variety `J_0(389)`.

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: M = ModularSymbols(389,sign=1).cuspidal_subspace()
            sage: L = LSeries(M); L
            L-series attached to Modular Symbols subspace of dimension 32 of Modular Symbols space of dimension 33 for Gamma_0(389) of weight 2 with sign 1 over Rational Field

        We first attempt computation of the analytic rank with the default of 53 bits precision::

            sage: L.analytic_rank()
            Traceback (most recent call last):
            ...
            RuntimeError: invalid L-series parameters: functional equation not satisfied

        The above failed because trying to compute one of the degree
        20 newforms resulting in some internal error when double
        checking the functional equation.  So we try with slightly more precision::

            sage: L.analytic_rank(70)
            13

        This works, since the factors have dimensions 1,2,3,6,20, and
        the one of degree 1 has rank 2, the ones of degree 2,3,6 have
        rank 2,3,6, respectively, and the one of degree 20 has rank 0::

            sage: 2*1 + 2 + 3 + 6
            13
        """
        return sum(e * L.analytic_rank(prec=prec)
                   for L, e in self._factorization)

    def weight(self):
        """
        Return the weight of this L-series, which is the sum of the weights
        of the factors counted with multiplicity.
        """
        return sum(e * L.weight() for L, e in self._factorization)

    @cached_method
    def _function(self, prec=53):
        """
        Return Dokchitser object that allows for computation of this
        L-series.  This is used behind the scenes for evaluation and
        computation of Taylor series.
        """
        return LSeriesProductEvaluator(self._factorization, prec)


class LSeriesZeta(LSeriesAbstract):
    """
    EXAMPLES::

        sage: from sage.lfunctions.eulerprod import LSeries
        sage: L = LSeries('zeta'); L
        Riemann Zeta function viewed as an L-series
        sage: L(2)
        1.64493406684823
        sage: L(2.000000000000000000000000000000000000000000000000000)
        1.64493406684822643647241516664602518921894990120680
        sage: zeta(2.000000000000000000000000000000000000000000000000000)
        1.64493406684822643647241516664602518921894990120680
        sage: L.local_factor(3)
        -T + 1
        sage: L.local_factor(5)
        -T + 1
        sage: L.anlist(30)
        [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    """
    def __init__(self):
        LSeriesAbstract.__init__(self, conductor=1, hodge_numbers=[0],
                                 weight=1, epsilon=1,
                                 poles=[0, 1], residues=[9, 8], base_field=QQ)
        T = ZZ['T'].gen()
        self._lf = 1 - T

    def _local_factor(self, P, prec):
        return self._lf

    def __repr__(self):
        return "Riemann Zeta function viewed as an L-series"


class LSeriesDelta(LSeriesAbstract):
    """
    EXAMPLES::

        sage: from sage.lfunctions.eulerprod import LSeriesDelta; L = LSeriesDelta()
        sage: L.anlist(10)
        [0, 1, -24, 252, -1472, 4830, -6048, -16744, 84480, -113643, -115920]
        sage: list(delta_qexp(11))
        [0, 1, -24, 252, -1472, 4830, -6048, -16744, 84480, -113643, -115920]
        sage: L.anlist(10^4) == list(delta_qexp(10^4+1))
        True
    """
    def __init__(self):
        LSeriesAbstract.__init__(self, conductor=1, hodge_numbers=[0, 1],
                                 weight=12, epsilon=1,
                                 poles=[], residues=[], base_field=QQ)
        self._T = ZZ['T'].gen()
        self._lf = {}

    def _local_factor(self, P, prec):
        """
        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeriesDelta; L = LSeriesDelta()
            sage: L.local_factor(2)
            2048*T^2 + 24*T + 1
            sage: L._local_factor(11, oo)   # really this is called
            285311670611*T^2 - 534612*T + 1

        The representation is reducible modulo 691::

            sage: L.local_factor(2).factor_mod(691)
            (666) * (T + 387) * (T + 690)
            sage: L.local_factor(3).factor_mod(691)
            (251) * (T + 234) * (T + 690)
            sage: L.local_factor(11).factor_mod(691)
            (468) * (T + 471) * (T + 690)

        because of the 691 here::

            sage: bernoulli(12)
            -691/2730
        """
        try:
            return self._lf[P]
        except KeyError:
            pass
        self._precompute_local_factors(P + 1)
        return self._lf[P]

    def _precompute_local_factors(self, bound, prec=oo):
        """
        Precompute local factors up to the given bound.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeriesDelta; L = LSeriesDelta()
            sage: L._lf
            {}
            sage: L._precompute_local_factors(10)
            sage: L._lf
            {2: 2048*T^2 + 24*T + 1, 3: 177147*T^2 - 252*T + 1, 5: 48828125*T^2 - 4830*T + 1, 7: 1977326743*T^2 + 16744*T + 1}
        """
        from sage.modular.all import delta_qexp
        T = self._T
        T2 = T**2
        f = delta_qexp(bound)
        for p in prime_range(bound):
            if p not in self._lf:
                self._lf[p] = 1 - f[p] * T + (p**11) * T2

    def __repr__(self):
        return "L-function associated to Ramanujan's Delta (a weight 12 cusp form)"


class LSeriesEllipticCurve(LSeriesAbstract):
    def __init__(self, E, prec=53):
        """
        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeriesEllipticCurve
            sage: L = LSeriesEllipticCurve(EllipticCurve('389a'))
            sage: L(2)
            0.360092863578881
        """
        E = E.global_minimal_model()
        self._E = E
        K = E.base_field()
        d = K.degree()
        self._N = E.conductor()
        LSeriesAbstract.__init__(self,
                                 conductor=norm(self._N) * K.discriminant()**2,
                                 hodge_numbers=[0] * d + [1] * d,
                                 weight=2,
                                 epsilon=[1, -1],
                                 poles=[],
                                 residues=[],
                                 base_field=K,
                                 prec=prec)

    def elliptic_curve(self):
        return self._E

    def __repr__(self):
        return "L-series of %s" % self.elliptic_curve()

    def _local_factor(self, P, prec):
        R = ZZ['T']
        T = R.gen()
        q = norm(P)
        p = prime_below(P)
        f = ZZ(q).ord(p)
        if P.divides(self._N):
            a = self._E.local_data(P).bad_reduction_type()
            return 1 - a * (T**f)
        else:
            a = q + 1 - self._E.reduction(P).count_points()
            return 1 - a * (T**f) + q * (T**(2 * f))


class LSeriesEllipticCurveQQ(LSeriesEllipticCurve):
    def __init__(self, E):
        E = E.global_minimal_model()
        self._E = E
        self._N = E.conductor()
        self._lf = {}
        self._T = ZZ['T'].gen()
        LSeriesAbstract.__init__(self, conductor=self._N,
                                 hodge_numbers=[0, 1], weight=2,
                                 epsilon=E.root_number(),
                                 poles=[], residues=[], base_field=QQ)

    def _lf0(self, p):
        a = self._E.ap(p)
        T = self._T
        if self._N % p == 0:
            if self._N % (p * p) == 0:
                return T.parent()(1)
            else:
                return 1 - a * T
        else:
            return 1 - a * T + p * T * T

    def _precompute_local_factors(self, bound, prec=oo):
        for p in primes(bound):
            if p not in self._lf:
                self._lf[p] = self._lf0(p)

    def _local_factor(self, P, prec):
        if P in self._lf:
            return self._lf[P]
        else:
            return self._lf0(P)

    def _primes_above(self, p):
        return [p]


class LSeriesDedekindZeta(LSeriesAbstract):
    """
    EXAMPLES::

        sage: from sage.lfunctions.eulerprod import LSeries
        sage: K.<a> = NumberField(x^3 - 2)
        sage: L = LSeries(K); L
        Dedekind Zeta function of Number Field in a with defining polynomial x^3 - 2
        sage: L(2)
        1.60266326190044
        sage: L.residues()
        'automatic'
        sage: L.residues(prec=53)
        [-4.77632833933856]
        sage: L.residues(prec=100)
        [-4.7763283393385594030639875094]
    """
    def __init__(self, K):
        if not K.is_absolute():
            K = K.absolute_field(names='a')
        self._K = K
        d = K.degree()
        sigma = K.signature()[1]
        LSeriesAbstract.__init__(self,
                                 conductor=abs(K.discriminant()),
                                 hodge_numbers=[0] * (d - sigma) + [1] * sigma,
                                 weight=1,
                                 epsilon=1,
                                 poles=[1],
                                 residues='automatic',
                                 base_field=K,
                                 is_selfdual=True)
        self._T = ZZ['T'].gen()

    def number_field(self):
        return self._K

    def __repr__(self):
        return "Dedekind Zeta function of %s" % self._K

    def _local_factor(self, P, prec):
        T = self._T
        return 1 - T**P.residue_class_degree()


class LSeriesDirichletCharacter(LSeriesAbstract):
    """
    EXAMPLES::

        sage: from sage.lfunctions.eulerprod import LSeries;  L = LSeries(DirichletGroup(5).0)
        sage: L(3)
        0.988191681624057 + 0.0891051883457395*I
    """
    def __init__(self, chi):
        if not chi.is_primitive():
            raise NotImplementedError("chi must be primitive")
        if chi.is_trivial():
            raise NotImplementedError("chi must be nontrivial")
        if chi.base_ring().characteristic() != 0:
            raise ValueError("base ring must have characteristic 0")
        self._chi = chi
        LSeriesAbstract.__init__(self, conductor=chi.conductor(),
                                 hodge_numbers=[1] if chi.is_odd() else [0],
                                 weight=1,
                                 epsilon=None,
                                 poles=[],  # since primitive
                                 residues=[],  # since primitive
                                 base_field=QQ,
                                 is_selfdual=chi.order() <= 2)
        self._T = ZZ['T'].gen()

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries;  L = LSeries(DirichletGroup(3).0)
            sage: L.__repr__()
            'L-series attached to Dirichlet character modulo 3 of conductor 3 mapping 2 |--> -1'
        """
        return "L-series attached to %s" % self._chi

    def character(self):
        return self._chi

    def epsilon(self, prec=oo):
        chi = self._chi
        if prec is oo:
            # answer in symbolic ring
            return (sqrt(-1) * SR(chi.gauss_sum())) / chi.modulus().sqrt()
        else:
            C = ComplexField(prec)
            x = C(chi.modulus()).sqrt() / chi.gauss_sum_numerical(prec=prec)
            if chi.is_odd():
                x *= C.gen()
            return 1 / x

    def _local_factor(self, P, prec):
        a = self._chi(P)
        if prec is not oo:
            a = ComplexField(prec)(a)
        return 1 - a * self._T


class LSeriesModularSymbolsAbstract(LSeriesAbstract):
    def modular_symbols(self):
        return self._M

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeries
            sage: f = Newforms(43,2,names='a')[1]; f
            q + a1*q^2 - a1*q^3 + (-a1 + 2)*q^5 + O(q^6)
            sage: LSeries(f).__repr__()
            'L-series of a degree 2 newform of level 43 and weight 2'
        """
        return "L-series of a degree %s newform of level %s and weight %s" % (self._M.dimension(), self._M.level(), self._M.weight())

    def _precompute_local_factors(self, bound, prec):
        primes = [p for p in prime_range(bound)
                  if p not in self._lf or self._lf[p][0] < prec]
        self._do_precompute(primes, prec)

    def _do_precompute(self, primes, prec):
        E, v = self._M.compact_system_of_eigenvalues(primes)
        if prec == 53:
            C = CDF
        elif prec is oo:
            if v.base_ring() == QQ:
                C = QQ
            else:
                C = CDF
        else:
            C = ComplexField(prec)

        phi = v.base_ring().embeddings(C)[self._conjugate]
        v = vector(C, [phi(a) for a in v])
        aplist = E.change_ring(C) * v
        T = C['T'].gen()
        T2 = T**2
        chi = self._M.character()
        k = self.weight()
        for i in range(len(primes)):
            p = primes[i]
            s = chi(p)
            if s != 0:
                s *= p**(k - 1)
            F = 1 - aplist[i] * T + s * T2
            self._lf[p] = (prec, F)

    def _local_factor(self, P, prec):
        # TODO: ugly -- get rid of all "prec=oo" in whole program --
        # always use oo.
        if P in self._lf and self._lf[P][0] >= prec:
            return self._lf[P][1]
        else:
            self._do_precompute([P], prec)
            return self._lf[P][1]


class LSeriesModularSymbolsNewformGamma0(LSeriesModularSymbolsAbstract):
    def __init__(self, M, conjugate=0, check=True, epsilon=None):
        """
        INPUT:

        - M -- a simple, new, cuspidal modular symbols space with
          sign 1
        - conjugate -- (default: 0), integer between 0 and dim(M)-1
        - check -- (default: True), if True, checks that M is
          simple, new, cuspidal, which can take a very long time,
          depending on how M was defined
        - epsilon -- (default: None), if not None, should be the sign
          in the functional equation, which is -1 or 1.  If this is
          None, then epsilon is computed by computing the sign of
          the main Atkin-Lehner operator on M.  If you have a faster
          way to determine epsilon, use it.

        EXAMPLES::

            sage: from sage.lfunctions.eulerprod import LSeriesModularSymbolsNewformGamma0
            sage: M = ModularSymbols(43,sign=1).cuspidal_subspace()[1]; M
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field
            sage: L0 = LSeriesModularSymbolsNewformGamma0(M,0,check=False,epsilon=1); L0
            L-series of a degree 2 newform of level 43 and weight 2
            sage: L0.taylor_series()
            0.620539857407845 + 0.331674007376949*z - 0.226392184536589*z^2 + 0.0960519649929789*z^3 - 0.00451826124421802*z^4 - 0.0203363026933833*z^5 + O(z^6)
            sage: L1 = LSeriesModularSymbolsNewformGamma0(M,1,check=False,epsilon=1); L1
            L-series of a degree 2 newform of level 43 and weight 2
            sage: L1(1)
            0.921328017272472
            sage: L1.taylor_series()
            0.921328017272472 + 0.492443075089339*z - 0.391019352704047*z^2 + 0.113271812405127*z^3 + 0.0213067052584679*z^4 - 0.0344198080536274*z^5 + O(z^6)
        """
        if M.dimension() == 0:
            raise ValueError("modular symbols space must positive dimension")
        chi = M.character()
        if chi is None or not chi.is_trivial():
            raise ValueError("modular symbols space must have trivial character")
        self._M = M
        N = M.level()

        if check:
            if not M.is_simple():
                raise ValueError("modular symbols space must be simple")
            if not M.is_new():
                raise ValueError("modular symbols space must be new")
            if not M.is_cuspidal():
                raise ValueError("modular symbols space must be cuspidal")

        k = M.weight()
        if epsilon is None:
            w = M.atkin_lehner_operator(N).matrix()
            if w not in [-1, 1]:
                raise ValueError("modular symbols space must have constant Atkin-Lehner operator")
            epsilon = (-1)**(k / 2) * w[0, 0]

        conjugate = ZZ(conjugate)
        if conjugate < 0 or conjugate >= M.dimension():
            raise ValueError("conjugate must a nonnegative integer less than the dimension")
        self._conjugate = conjugate
        self._lf = {}

        LSeriesAbstract.__init__(self,
                                 conductor=N,
                                 hodge_numbers=[0, 1],
                                 weight=k,
                                 epsilon=epsilon,
                                 poles=[],  # since primitive
                                 residues=[],  # since primitive
                                 base_field=QQ)


def _is_valid_modsym_space(M):
    if not is_ModularSymbolsSpace(M):
        raise TypeError("must be a modular symbols space")
    if M.dimension() == 0:
        raise ValueError("modular symbols space must positive dimension")
    if M.character() is None:
        raise ValueError("modular symbols space must have associated character")
    if not M.is_simple():
        raise ValueError("modular symbols space must be simple")
    if not M.is_new():
        raise ValueError("modular symbols space must be new")
    if not M.is_cuspidal():
        raise ValueError("modular symbols space must be cuspidal")


class LSeriesModularSymbolsNewformCharacter(LSeriesModularSymbolsAbstract):
    def __init__(self, M, conjugate=0):
        _is_valid_modsym_space(M)
        chi = M.character()
        self._M = M
        N = M.level()

        # See Remark 5.0.2 in [Diamond-Im] which says: "Let f be a newform of level N
        # which is a common eigenform under all the Hecke operators T_p.  Then
        # w_N(f) = c*fbar, where fbar = sum bar(a_n) q^n and c is a scalar.  The functional
        # equation may be rewritten as   Lambda(s, f) = c * i^k * Lambda(k-s, fbar).
        # That said, this seems hard to compute, so we just solve using the
        # functional equation.
        epsilon = 'solve'

        k = M.weight()
        conjugate = ZZ(conjugate)
        if conjugate < 0 or conjugate >= M.dimension():
            raise ValueError("conjugate must a nonnegative integer "
                             "less than the dimension")
        self._conjugate = conjugate

        LSeriesAbstract.__init__(self,
                                 conductor=N,
                                 hodge_numbers=[0, 1],
                                 weight=k,
                                 epsilon=epsilon,
                                 poles=[],  # since primitive
                                 residues=[],  # since primitive
                                 base_field=QQ,
                                 is_selfdual=chi.order() <= 2)
        self._lf = {}


def _new_modsym_space_with_multiplicity(M):
    """
    Return a simple new modular symbols space N and an integer d such
    that M is isomorphic to `N^d` as a module over the anemic Hecke algebra.

    INPUT:

    - M -- a sign=1 modular simple space for the full Hecke
      algebra (including primes dividing the level) that cannot be
      decomposed further by the Hecke operators.  None of the
      conditions on M are explicitly checked.

    OUTPUT:

    - N -- a simple new modular symbols space
    - d -- a positive integer
    """
    if M.is_new():
        return [(M, 1)]
    raise NotImplementedError


def LSeriesModularSymbolsNewform(M, i=0):
    chi = M.character()
    if chi is None:
        raise NotImplementedError
    elif chi.is_trivial():
        return LSeriesModularSymbolsNewformGamma0(M, i)
    else:
        return LSeriesModularSymbolsNewformCharacter(M, i)


class LSeriesModularSymbolsMotive(LSeriesProduct):
    """
    The product of L-series attached to the modular symbols space M.
    """
    def __init__(self, M):
        self._M = M
        if not is_ModularSymbolsSpace(M):
            raise TypeError("X must be a modular symbols space or have a modular symbols method")
        self._M = M
        D = M.decomposition()
        for A in D:
            _is_valid_modsym_space(A)
        F = []
        for A in D:
            for X in _new_modsym_space_with_multiplicity(A):
                N, d = X
                for i in range(N.dimension()):
                    F.append((LSeriesModularSymbolsNewform(N, i), d))
        LSeriesProduct.__init__(self, F)

    def modular_symbols(self):
        return self._M

    def __repr__(self):
        return "L-series attached to %s" % self._M


class LSeriesModularAbelianVariety(LSeriesProduct):
    """
    The product of L-series attached to the modular abelian variety A.

    EXAMPLES::

        sage: from sage.lfunctions.eulerprod import LSeries
        sage: L = LSeries(J0(54)); L
        L-series attached to Abelian variety J0(54) of dimension 4
        sage: L.factor()
        (L-series of a degree 1 newform of level 27 and weight 2)^2 * (L-series of a degree 1 newform of level 54 and weight 2) * (L-series of a degree 1 newform of level 54 and weight 2)
        sage: L(1)
        0.250717238804658
        sage: L.taylor_series(prec=20)
        0.25072 + 0.59559*z + 0.15099*z^2 - 0.35984*z^3 + 0.056934*z^4 + 0.17184*z^5 + O(z^6)

    Independent check of L(1)::

        sage: prod(EllipticCurve(lbl).lseries()(1) for lbl in ['54a', '54b', '27a', '27a'])
        0.250717238804658

    Different check that totally avoids using Dokchitser::

        sage: prod(EllipticCurve(lbl).lseries().at1(prec=53)[0] for lbl in ['54a', '54b', '27a', '27a'])
        0.250848605530185
    """
    def __init__(self, A):
        self._A = A
        D = A.decomposition()
        F = None
        for A in D:
            # TODO: This is ugly, but I don't know a cleaner way to do it yet.
            # Could be painfully inefficient in general.
            f = Newform(A.newform_label(), names='a')
            M = f.modular_symbols(sign=1)
            d = ZZ(A.dimension() / M.dimension())
            L = LSeriesModularSymbolsMotive(M)**d
            if F is None:
                F = L
            else:
                F *= L
        if F is None:
            raise ValueError("abelian variety must have positive dimension")
        LSeriesProduct.__init__(self, F.factor())

    def abelian_variety(self):
        return self._A

    def __repr__(self):
        return "L-series attached to %s" % self._A


class LSeriesTwist(LSeriesAbstract):
    """
    Twist of an L-series by a character.
    """
    def __init__(self, L, chi, conductor=None, epsilon=None, prec=53):
        """
        INPUT:

        - `L` -- an L-series
        - ``chi`` -- a character of the base field of L
        - ``conductor`` -- None, or a list of conductors to try
        - ``prec`` -- precision to use when trying conductors, if
          conductor is a list
        """
        self._L = L
        self._chi = chi

        if not chi.is_primitive():
            raise ValueError("character must be primitive")

        A = ZZ(L.conductor())
        B = chi.conductor()
        if conductor is None:
            if A.gcd(B) != 1:
                # Make a list of all possible conductors, and let the
                # functional equation figure it out.
                smallest = ZZ(A)
                while smallest.gcd(B) != 1:
                    smallest = smallest // smallest.gcd(B)
                biggest = A * (B**L.degree())
                assert biggest % smallest == 0
                #
                # TODO: improve this using the theorem stated
                # on page 1 of http://wstein.org/papers/padictwist/
                #
                conductor = [smallest * d
                             for d in divisors(biggest // smallest)]
            else:
                conductor = A * (B**L.degree())
        hodge_numbers = L.hodge_numbers()
        weight = L.weight()
        if epsilon is None:
            if L.epsilon() != 'solve':
                if chi.order() <= 2:
                    if A.gcd(B) == 1:
                        epsilon = L.epsilon() * chi(-A)
                    else:
                        epsilon = [L.epsilon(), -L.epsilon()]
                else:
                    epsilon = 'solve'
            else:
                epsilon = 'solve'
        is_selfdual = L.is_selfdual()
        poles = []  # ???  TODO -- no clue here.
        residues = 'automatic'
        base_field = L.base_field()
        LSeriesAbstract.__init__(self, conductor, hodge_numbers,
                                 weight, epsilon,
                                 poles, residues, base_field, is_selfdual, prec)

    def _local_factor(self, P, prec):
        L0 = self._L.local_factor(P, prec)
        chi = self._chi
        T = L0.parent().gen()
        c = chi(P)
        if prec is not oo:
            c = ComplexField(prec)(c)
        return L0(c * T)

    def __repr__(self):
        return "Twist of %s by %s" % (self._L, self._chi)

    def untwisted_lseries(self):
        return self._L

    def twist_character(self):
        return self._chi


def LSeries(X, *args, **kwds):
    """
    Return the L-series of X, where X can be any of the following:

        - elliptic curve over a number field (including QQ)
        - Dirichlet character
        - cuspidal newform
        - new cuspidal modular symbols space -- need not be simple
        - string: 'zeta' (Riemann Zeta function), 'delta'
        - modular elliptic curve attached to Hilbert modular forms space

    For convenience, if L is returned, then L._X is set to X.

    EXAMPLES::

    The Dedekind Zeta function of a number field::

        sage: from sage.lfunctions.eulerprod import LSeries
        sage: K.<a> = NumberField(x^2 + 1)
        sage: L = LSeries(K); L
        Dedekind Zeta function of Number Field in a with defining polynomial x^2 + 1
        sage: L(2)
        1.50670300992299

        sage: K.zeta_coefficients(100) == L.anlist(100)[1:]
        True

        sage: L = LSeries(ModularSymbols(43, weight=2,sign=1).cuspidal_subspace().decomposition()[1])
        sage: L(1)
        0.571720756464112
        sage: L.factor()[0][0](1)
        0.620539857407845
        sage: L.factor()[1][0](1)
        0.921328017272472
        sage: L._X
        Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field
        sage: L = LSeries(ModularSymbols(DirichletGroup(13).0^2, weight=2,sign=1).cuspidal_subspace())
        sage: L(1)
        0.298115272465799 - 0.0402203326076733*I

    The L-series of a modular abelian variety with both new and old parts::

        sage: L = LSeries(J0(33)); L
        L-series attached to Abelian variety J0(33) of dimension 3
        sage: L.factor()
        (L-series of a degree 1 newform of level 11 and weight 2)^2 * (L-series of a degree 1 newform of level 33 and weight 2)
        sage: L.local_factor(2, prec=oo)
        8*T^6 + 12*T^5 + 12*T^4 + 8*T^3 + 6*T^2 + 3*T + 1
        sage: L(1)
        0.0481553138900504

    We check the above computation of L(1) via independent methods (and implementations)::

        sage: prod(EllipticCurve(lbl).lseries().at1(prec=53)[0] for lbl in ['11a', '11a', '33a'])
        0.0481135342926321
        sage: prod(EllipticCurve(lbl).lseries()(1) for lbl in ['11a', '11a', '33a'])
        0.0481553138900504

    A nonsimple new modular symbols space of level 43::

        sage: L = LSeries(ModularSymbols(43,sign=1).cuspidal_subspace())
        sage: L
        L-series attached to Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field
        sage: L(1)
        0.000000000000000
        sage: L.taylor_series()
        0.000000000000000 + 0.196399786632435*z + 0.314922741074845*z^2 - 0.0797083673829092*z^3 - 0.161630566287135*z^4 + 0.123939472976207*z^5 + O(z^6)
        sage: L.factor()
        (L-series of a degree 1 newform of level 43 and weight 2) * (L-series of a degree 2 newform of level 43 and weight 2) * (L-series of a degree 2 newform of level 43 and weight 2)
        sage: L.analytic_rank()
        1
        sage: D = ModularSymbols(43,sign=1).cuspidal_subspace().decomposition()
        sage: L0 = LSeries(D[0]); L1 = LSeries(D[1])
        sage: L0.taylor_series() * L1.taylor_series()
        0.000000000000000 + 0.196399786632435*z + 0.314922741074845*z^2 - 0.0797083673829091*z^3 - 0.161630566287135*z^4 + 0.123939472976207*z^5 + O(z^6)
        sage: L0.factor()
        L-series of a degree 1 newform of level 43 and weight 2
        sage: L1.factor()
        (L-series of a degree 2 newform of level 43 and weight 2) * (L-series of a degree 2 newform of level 43 and weight 2)

    """
    L = _lseries(X, *args, **kwds)
    L._X = X
    return L


def _lseries(X, *args, **kwds):
    """
    Helper function used by LSeries function.
    """
    if is_EllipticCurve(X):
        K = X.base_ring()
        if is_RationalField(K):
            return LSeriesEllipticCurveQQ(X, *args, **kwds)
        else:
            return LSeriesEllipticCurve(X)

    if is_DirichletCharacter(X):
        if X.is_trivial() and X.is_primitive():
            return LSeriesZeta(*args, **kwds)
        else:
            return LSeriesDirichletCharacter(X, *args, **kwds)

    if is_NumberField(X):
        return LSeriesDedekindZeta(X, *args, **kwds)

    if isinstance(X, sage.modular.modform.element.Newform):
        return LSeriesModularSymbolsNewform(X.modular_symbols(sign=1), *args, **kwds)

    if is_ModularSymbolsSpace(X):
        if X.sign() != 1:
            raise NotImplementedError
        return LSeriesModularSymbolsMotive(X, *args, **kwds)

    if isinstance(X, str):
        y = X.lower()
        if y == 'zeta':
            return LSeriesZeta(*args, **kwds)
        elif y == 'delta':
            return LSeriesDelta(*args, **kwds)
        else:
            raise ValueError('unknown L-series "%s"' % y)

    if is_ModularAbelianVariety(X):
        return LSeriesModularAbelianVariety(X, *args, **kwds)

    raise NotImplementedError
