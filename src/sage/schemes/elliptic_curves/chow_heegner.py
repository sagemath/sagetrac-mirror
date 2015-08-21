##############################################################################
#       Copyright (C) 2011 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################
"""
Chow-Heegner Points

This is a package for numerically computing Chow-Heegner points
associated to pairs (E,F) of non-isogenous optimal elliptic curves
over the rational numbers of the same conductor. For a description of
the algorithm, see http://wstein.org/papers/chow_heegner/.

AUTHORS:

- William Stein (November 2011, January 2012)

EXAMPLES::

    sage: P = EllipticCurve('57a').chow_heegner_point(EllipticCurve('57b'))
    sage: P.numerical_approx(deg1=100)
    (1.44444444444... : -1.03703703703... : 1.00000000000000)
    sage: P.point_exact(deg1=100)
    (13/9 : -28/27 : 1)
    sage: P.index(deg1=100)
    8

TESTS::

    sage: P = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b'))
    sage: TestSuite(P).run()
"""

# General Sage imports
from sage.misc.all import verbose
from sage.misc.cachefunc import cached_method
from sage.misc.misc import cputime
from sage.structure.element import parent
from sage.parallel.all import parallel
from sage.rings.all import (QQ, ZZ, CDF, RDF, RR, infinity, xgcd, NumberField, ComplexField)

from sage.rings.complex_field  import is_ComplexField
from sage.rings.rational_field import is_RationalField

from sage.modular.all import Gamma0

# Elliptic curve imports
from constructor import EllipticCurve
from ell_generic import is_EllipticCurve

# Fast Cython functions specifically for this code.
from sage.schemes.elliptic_curves.chow_heegner_fast import (
    cdf_roots_of_rdf_poly, Polynomial_RDF_gsl, ComplexPolynomial)

######################################################################
# Gamma0(N) equivalence of points in the upper half plane
######################################################################

def _slz2_rep_in_fundom_helper(z):
    """
    Return element z2 of the standard fundamental domain for SL2Z
    equivalent to z and a transformation g in SL2Z that maps z to z2,
    modulo rounding errors.  For precision guarantees, see the
    function slz2_rep_in_fundom below.

    INPUT:

    - `z` -- floating point complex number in the upper half plane

    OUTPUT:

    - element of upper half plane
    - 2x2 integer matrix of determinant 1

    EXAMPLES:

    This helper function just runs the raw algorithm, so can yield
    results that are subject to substantial rounding errors::

        sage: from sage.schemes.elliptic_curves.chow_heegner import _slz2_rep_in_fundom_helper
        sage: z = CDF(110.1, 1e-10)
        sage: zz,g = _slz2_rep_in_fundom_helper(z)
        sage: zz
        0.341867712408 + 99999999.6769*I

    The value of zz above is wrong, as we see applying the raw
    algorithm but to higher precision using an interval (so that we
    are aware of the influence of rounding errors)::

        sage: zz,g = _slz2_rep_in_fundom_helper(ComplexIntervalField(100)(z))
        sage: zz
        0.241867713702? + 9.99999996768825706223026?e7*I

    The sl2z_rep_in_fundom function takes care of using intervals
    automatically, and to sufficient precision, but may be much
    slower.
    """
    # We use the standard algorithm: replace z by z-n and z by -1/z
    # until |Re(z)|<=1/2 and |z|>=1.
    from sage.modular.all import SL2Z
    gamma = SL2Z(1)
    S, T = SL2Z.gens()
    change = True
    half = z.real().parent()(1)/2
    while change:
        change = False
        t = z.real()
        if abs(t) > half:
            change = True
            # |t - n| <= 1/2
            # -1/2 <= t-n <= 1/2
            # n - 1/2 <= t < = 1/2+n
            # n <= t + 1/2 <= n + 1, so n = floor(t+1/2)
            n = (t + half).floor()  # avoid rounding error with 0.5
            z -= n
            if hasattr(n, 'center'):
                k = round(n.center())
            else:
                k = n
            gamma *= T**k
        if abs(z) < 1:
            change = True
            z = -1/z
            gamma *= S
    return z, gamma**(-1)

def sl2z_rep_in_fundom(z, eps=None):
    """
    Return element z2 of the standard fundamental domain for SL2Z
    equivalent to z and a transformation g in SL2Z.  The transformation
    g, which is an exact element of SL2Z, will take any sufficiently
    high precision approximation of z to z2 with an error of at most
    eps.

    Since this function uses interval arithmetic to be careful about
    precision, it is slower than _slz2_rep_in_fundom_helper, which
    makes no attempt to deal with precision issues.

    INPUT:

    - `z` -- floating point complex number in the upper half plane
    - `eps` -- a small positive real or None; if None, uses 2**(-prec),
      where prec is the number of bits of precision of z.

    OUTPUT:

    - element of upper half plane
    - 2x2 integer matrix of determinant 1

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.chow_heegner import sl2z_rep_in_fundom
        sage: z = CDF(2.17,.19)
        sage: w, g = sl2z_rep_in_fundom(z)
        sage: w
        0.384615384615 + 2.92307692308*I
        sage: g
        [-3  7]
        [-1  2]
        sage: g.acton(z)
        0.384615384615 + 2.92307692308*I

    An example in which the imaginary part is small, which leads to
    rounding issues::

        sage: z = CDF(110.1, 1e-10)
        sage: zz, g = sl2z_rep_in_fundom(z)
        sage: zz
        0.241867713702 + 99999999.6769*I
        sage: g
        [ -56841 6258194]
        [     10   -1101]
        sage: g.acton(z)         # far from zz due to rounding during action
        -5684.1 + 99999999.6275*I
        sage: g.acton(CIF(z))    # indeed to double precision, everything lost
        0.?e5 + 1.0000000?e8*I
        sage: g.acton(ComplexField(100)(z))   # get correct zz above with higher prec
        0.24186771370190053192231958133 + 9.9999999676882570622302579775e7*I

    Here is a similer, but much more extreme, example::

        sage: z = CC(110.1 + 1e-100*I)
        sage: zz, g = sl2z_rep_in_fundom(z); zz
        -0.499999999999858 + 8.07793566946316e72*I
        sage: g.acton(z)
        -0.499999999999858
        sage: g.acton(ComplexField(10^3)(z))
        -0.49999999999985...000 + 8.07793566946...e72*I

    TESTS:

    The value of eps must be positive and small or a ValueError
    exception is raised::

        sage: from sage.schemes.elliptic_curves.chow_heegner import *
        sage: sl2z_rep_in_fundom(CC(1+I), eps=0)   # must be positive
        Traceback (most recent call last):
        ...
        ValueError: eps must be positive
        sage: sl2z_rep_in_fundom(CC(1+I), eps=-1)  # must be positive
        Traceback (most recent call last):
        ...
        ValueError: eps must be positive
        sage: sl2z_rep_in_fundom(CC(1+I), eps=1)   # must be small
        Traceback (most recent call last):
        ...
        ValueError: prec (=0) must be >= 2 and <= -1.

    TESTS:

    We require that the input is in the upper half plane::

        sage: from sage.schemes.elliptic_curves.chow_heegner import *
        sage: sl2z_rep_in_fundom(CC(1,-1))
        Traceback (most recent call last):
        ...
        ValueError: z must be in the upper half plane
        sage: sl2z_rep_in_fundom(CC(1,-1))
        Traceback (most recent call last):
        ...
        ValueError: z must be in the upper half plane
        sage: sl2z_rep_in_fundom(float(3))
        Traceback (most recent call last):
        ...
        ValueError: z must be in the upper half plane
        sage: sl2z_rep_in_fundom(3)
        Traceback (most recent call last):
        ...
        ValueError: z must be in the upper half plane

    The Python complex type is allowed::

        sage: sl2z_rep_in_fundom(complex(2.17,.19))
        (
                                          [-3  7]
        0.384615384615 + 2.92307692308*I, [-1  2]
        )
    """
    if isinstance(z, complex):
        z = CDF(z)

    if isinstance(z, (float, int, long)) or z.imag() <= 0:
        raise ValueError, "z must be in the upper half plane"

    if eps is None:
        # set default value of eps if it is not given
        eps = RR(2)**(-z.parent().prec())
    else:
        # make sure eps is *positive* and in RR.
        eps = RR(eps)
        if eps <= 0:
            raise ValueError, "eps must be positive"

    # To avoid rounding errors, we use interval arithmetic for the
    # actual calculation (until the precision is sufficient).
    from sage.rings.complex_interval_field import (
        is_ComplexIntervalField, ComplexIntervalField)

    # We start with twice the precision to avoid an extra iteration:
    # the documentation warns that this function may be slow.
    prec = -2*eps.log(2).ceil()
    while True:
        # make the interval version of z to the given precision
        w = ComplexIntervalField(prec)(z)
        # run the algorithm on that interval
        w2, g = _slz2_rep_in_fundom_helper(w)
        if w2.diameter() <= eps:
            # There wasn't too much precision loss, so we go with this result.
            return z.parent()(w2.center()), g
        else:
            # Too much precision loss -- try again with twice as many bits.
            prec *= 2


def canonicalize_sl2z(a, g=None, eps_ratio=2):
    """
    Assume that a = g(z) is in the fundamental domain for SL2Z.
    Adjust a by applying T^(-1) or S so that a is the canonical
    representative in the fundamental domain, so a is not on the right
    edge, and if a is on the unit circle, then it is on the left hand
    side.  Also, modify g so that the relation a = g(z) continues to
    hold.

    Special case: if g is None, just adjust a, ignoring g.

    INPUT:

    - a -- element of fundamental domain (so in ComplexField(prec))
    - g -- None or element of SL2Z
    - eps_ratio -- default 2; used in testing equality of the floating
      point number a with 1 and 1/2. See comment in source code.

    OUTPUT:

    - new a and a matrix g (or None)

    EXAMPLES:

    When a is on the right edge::

        sage: from sage.schemes.elliptic_curves.chow_heegner import canonicalize_sl2z
        sage: g = SL2Z([1,0,0,1])
        sage: canonicalize_sl2z(CDF(1/2, 1), g)
        (
                      [ 1 -1]
        -0.5 + 1.0*I, [ 0  1]
        )

    Don't bother to compute the matrix g::

        sage: canonicalize_sl2z(CDF(1/2, 1))
        (-0.5 + 1.0*I, None)

    A higher precision example::

        sage: canonicalize_sl2z(ComplexField(200)(1/2,1))
        (-0.50000000000000000000000000000000000000000000000000000000000 + 1.0000000000000000000000000000000000000000000000000000000000*I, None)

    When a is on the right side of the unit circle::

        sage: canonicalize_sl2z(CDF(1/sqrt(2),1/sqrt(2)), SL2Z([1,0,0,1]))
        (
                                            [ 0 -1]
        -0.707106781187 + 0.707106781187*I, [ 1  0]
        )

    When a is in the middle (so a does not change)::

        sage: canonicalize_sl2z(CDF(0,1), SL2Z([1,0,0,1]))
        (
               [1 0]
        1.0*I, [0 1]
        )

    Using a bigger eps_ratio makes it so comparison of floating point
    numbers is less sensitive::

        sage: canonicalize_sl2z(CDF(0.500001, 1))
        (0.500001 + 1.0*I, None)
        sage: canonicalize_sl2z(CDF(0.500001, 1), eps_ratio=4)
        (-0.499999 + 1.0*I, None)
    """
    if g is not None:
        S, T = g.parent().gens()
    # There are two cases where points on boundary are identified,
    # as explained in Theorem 1 of Chapter VII of Serre's "A Course
    # in Arithmetic".
    # Here we have to check for equality of a floating point number
    # with 1/2 and with 1.  In each case, exact equality of course
    # almost never happens, so we must use some notion of "eps",
    # and we use eps=2**(-prec//eps_ratio), where prec is the number of bits
    # of precision of the parent, and eps_ratio=2 by default.
    C = a.parent()
    eps = 2**(-C.prec() // eps_ratio)
    if abs(a.real() - C(1)/2)<eps:
        a -= 1
        if g is not None: g = T**(-1)*g
    elif abs(abs(a) - 1)<eps and a.real() > 0:
        # points are sl2z equivalent on boundary of unit circle
        a = -1/a
        if g is not None: g = S*g
    return a, g

def is_sl2z_equivalent(z1, z2, prec):
    """
    Return True if the complex numbers z1 and z2 in the upper half
    plane are equivalent modulo the action of SL_2(Z), at least to the
    given number prec of bits of absolute precision.  Canonical
    representatives equivalent to z1 and z2 are found, and declared
    equal if the absolute value of their difference is less than
    2**(-prec).

    INPUT:

    - `z_1`, `z_2` -- floating point complex numbers in the upper
      half plane
    - ``prec`` -- positive integer (bits of absolute precision)

    OUTPUT:

    - bool

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.chow_heegner import is_sl2z_equivalent

    Trivial special case::

        sage: is_sl2z_equivalent(CDF(5,7), CDF(5,7), 53)
        True

    Apply an element of SL2Z, and note that precision 53 fails due to
    rounding errors, though 40 bits works::

        sage: z1 = CDF(5,7)
        sage: z2 = (-5*z1 - 6)/(11*z1 + 13); z2
        -0.455131242301 + 0.000663318487634*I
        sage: is_sl2z_equivalent(z1, z2, 53)
        False
        sage: is_sl2z_equivalent(z1, z2, 40)
        True

    Here the two elements are not equivalent::

        sage: is_sl2z_equivalent(z1, .5*z2, 20)
        False

    This is the same as above, but to higher precision::

        sage: z1 = ComplexField(200)(5,7)
        sage: z2 = (-5*z1 - 6)/(11*z1 + 13)
        sage: is_sl2z_equivalent(z1, z2, 150)
        True
    """
    w1, _ = sl2z_rep_in_fundom(z1)
    w2, _ = sl2z_rep_in_fundom(z2)
    a1, _ = canonicalize_sl2z(w1)
    a2, _ = canonicalize_sl2z(w2)
    return abs(a1 - a2) < 2**(-prec)

def is_gamma0N_equivalent(z1, z2, N, prec):
    """
    Return True if z1 and z2 are equivalent modulo the action of
    `\\Gamma_0(N)`` to the given number of bits of precision.

    INPUT:

    - `z_1`, `z_2` -- floating point complex numbers in the upper
      half plane
    - `N` -- positive integer
    - ``prec`` -- positive integer (bits of absolute precision)

    OUTPUT:

    - bool

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.chow_heegner import is_gamma0N_equivalent

    Construct a point and its image via the action of a matrix in
    Gamma0(11), and do some basic checks::

        sage: z1 = CDF(5,7)
        sage: z2 = (-5*z1 - 6)/(11*z1 + 13); z2
        -0.455131242301 + 0.000663318487634*I
        sage: is_gamma0N_equivalent(z1, z2, 11, 30)
        True
        sage: is_gamma0N_equivalent(z1, z2, 12, 30)
        False

    Test under image of all generators of Gamma0(12))::

        sage: all(is_gamma0N_equivalent(z1, g.acton(z1), 12, 30) for g in Gamma0(12).gens())
        True

    Test to higher precision using Gamma0(15)::

        sage: z1 = ComplexField(200)(5,7)
        sage: all(is_gamma0N_equivalent(z1, g.acton(z1), 15, 150) for g in Gamma0(15).gens())
        True

    These two points have the property that z1 and z2 are equivalent
    modulo SL2(Z), and 5*z1 and 5*z2 are *also* equivalent modulo
    SL2(Z), but z1 and z2 are not equivalent modulo Gamma0(N).  Note
    that both points are equivalent to I, which has extra
    automorphisms.  This is interesting because generators for the
    modular function field j(z) and j(5*z) agree on z1 and z2, but
    the points are inequivalent (this illustrates that this model
    for the modular curve X0(5) is singular)::

        sage: from sage.schemes.elliptic_curves.chow_heegner import is_gamma0N_equivalent, is_sl2z_equivalent
        sage: z1 = CDF(-2,1)/5
        sage: z2 = CDF(2,1)/5
        sage: is_gamma0N_equivalent(z1,z2,5, 30)
        False
        sage: is_sl2z_equivalent(z1, z2, 30)
        True
        sage: is_sl2z_equivalent(5*z1, 5*z2, 30)
        True
    """
    C = ComplexField(prec)
    eps = RR(2)**(-prec)
    w1, g1 = sl2z_rep_in_fundom(z1, eps)  # g1(z1) = w1 = canonical rep
    w2, g2 = sl2z_rep_in_fundom(z2, eps)  # g2(z2) = w2 = canonical rep
    a1, g1 = canonicalize_sl2z(w1, g1)
    a2, g2 = canonicalize_sl2z(w2, g2)
    if abs(a1 - a2) >= eps:
        # The points are not even sl2z-equivalent, so they can't be
        # Gamma_0(N) equivalent
        return False

    # Now we may assume that g1(z1) = g2(z2), because of the adjustments
    # made above.  We double check.
    # This C2 is purely used for the assert below and nothing else, so this
    # "magic constant" not so evil.
    C2 = ComplexField(2*prec+10)
    assert abs(g1.acton(C2(z1)) - g2.acton(C2(z2))) < eps

    # So now we have z := g1(z1) = g2(z2), both in the standard
    # fundamental domain.
    #
    # The nontrivial elements of PSL2Z with a fixed point z in the
    # standard fundamental domain for the upper half plane are
    # Stab(z), where
    #
    #     * z = i, so Stab(z) generated by S (where S has order 2)
    #     * z = rho = exp(2*pi*i/3) so Stab(z) generated by S*T
    #     * z = -rhobar = exp(pi*i/3) so Stab(z) generated by T*S
    #
    # The elements in PSL2Z that send z1 to z2 are the elements
    # g2^(-1)*A*g1 for A in Stab(z), so we just check if any are in
    # Gamma0(N).

    g2i = g2**(-1)
    g = g2i*g1
    i = C.gen()
    pi = C.pi()
    if g[1,0]%N == 0:
        return True
    S, T = g1.parent().gens()
    if a1 == i:
        return (g2i*S*g1)[1,0]%N == 0
    elif a1 == ((2*pi*i)/3).exp():
        return (g2i*S*T*g1)[1,0]%N == 0 or (g2i*S*T*S*T*g1)[1,0]%N == 0
    return False


class X0NPoint(object):
    """
    An affine point on the modular curve `X_0(N)` represented as a
    floating point complex number.  We use this for representing
    points in the fiber over a point on an elliptic curve.

    TESTS::

        sage: from sage.schemes.elliptic_curves.chow_heegner import X0NPoint
        sage: x = X0NPoint(CDF(1,1), 11, 30)
        sage: TestSuite(x).run()
    """
    def __init__(self, z, N, prec):
        """
        INPUT:

        - `z` -- floating point complex number in the upper half plane
        - `N` -- positive integer
        - `prec` -- positive integer; bits of precision used in comparisons

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import X0NPoint
            sage: x = X0NPoint(CDF(1,1), 11, 30); x
            [1.0000000 + 1.0000000*I]
            sage: type(x)
            <class 'sage.schemes.elliptic_curves.chow_heegner.X0NPoint'>
        """
        self._z = z
        self._N = N
        self._prec = prec
        self._C = ComplexField(prec)

    def z(self):
        """
        Return representative point z in the upper half plane as a
        floating point number.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import X0NPoint
            sage: X0NPoint(CDF(1,1), 11, 30).z()
            1.0 + 1.0*I
        """
        return self._z

    def __cmp__(self, right):
        """
        Compare self and right.  If they have different levels, then
        the one with smaller level is considered smaller.  The points
        are considered equal if the two points are equivalent modulo
        the action of Gamma0(N) to the minimum of the precisions of
        self and right.  If the levels are the same, but the points
        are not equivalent, this function compares the underlying
        floating point representatives as complex numbers.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import X0NPoint
            sage: g = SL2Z([-5,-1,11,2])
            sage: z = CDF(1,1); x = X0NPoint(z, 11, 30); y = X0NPoint(g.acton(z), 11, 30)
            sage: x == y
            True
            sage: z = CDF(1,1); x = X0NPoint(z, 11, 30); y = X0NPoint(z/2, 11, 30)
            sage: x == y
            False

        Points with different levels::

            sage: z = CDF(1,1); x=X0NPoint(z,5,30); y=X0NPoint(z,10,30)
            sage: x < y
            True
            sage: y < x
            False
        """
        c = cmp(self._N, right._N)
        if c: return c
        if self._N == 1:
            return cmp(self._sl2z_rep(), right.sl2z_rep())
        if is_gamma0N_equivalent(self._z, right._z, self._N, min(self._prec, right._prec)):
            return 0
        return cmp(self._z, right._z)

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import X0NPoint
            sage: X0NPoint(CDF(1,1), 11,30).__repr__()
            '[1.0000000 + 1.0000000*I]'
        """
        return "[%s]"%self._C(self._z)

    @cached_method
    def sl2z_rep(self):
        """
        Canonical element of fundamental domain for SL2Z that is
        equivalent to self.

        EXAMPLES::

            sage: import sage.schemes.elliptic_curves.chow_heegner as h
            sage: h.X0NPoint(CDF(2,.1), 11,30).sl2z_rep()
            10.000000*I
            sage: h.sl2z_rep_in_fundom(CDF(2,.1))[0]
            10.0*I
        """
        a = self._C(sl2z_rep_in_fundom(self._z)[0])
        b = canonicalize_sl2z(a)[0]
        return b

    def atkin_lehner(self, q=None):
        """
        Return image of this point under the Atkin-Lehner involution W_q.

        INPUT:

        - `q` -- integer that exactly divides the level, or
          None, in which case q defaults to the level.

        OUTPUT:

        - X0NPoint

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import X0NPoint
            sage: x = X0NPoint(CDF(1,1), 12, 30); x.atkin_lehner()
            [-0.041666667 + 0.041666667*I]
            sage: x.atkin_lehner().atkin_lehner()
            [1.0000000 + 1.0000000*I]
            sage: x.atkin_lehner(4)
            [0.32000000 + 0.010000000*I]
            sage: x.atkin_lehner(3)
            [-0.26016260 + 0.0081300813*I]
            sage: x.atkin_lehner(3).atkin_lehner(3)
            [1.0000000 + 1.0000000*I]
            sage: x.atkin_lehner(3).atkin_lehner(3) == x
            True
            sage: x.atkin_lehner(2)
            Traceback (most recent call last):
            ...
            ValueError: q must exactly divide N
            sage: x.atkin_lehner(5)
            Traceback (most recent call last):
            ...
            ValueError: q must divide N

        The parameter q need not be a prime power::

            sage: z = CDF(1,1)
            sage: P = X0NPoint(z,90,30) # N = 90 = 10*9
            sage: P.atkin_lehner(10)
            [0.11049724 + 0.00055248619*I]
            sage: P.atkin_lehner(10).atkin_lehner(10)
            [0.10554020 + 0.000013888696*I]

        This is because applying Atkin-Lehner results in substantial precision loss::

            sage: P.atkin_lehner(10).atkin_lehner(10) == P
            False

        Using higher precision for the underlying point fixes the problem::

            sage: P = X0NPoint(ComplexField(200)(1,1),90,30)
            sage: P.atkin_lehner(10).atkin_lehner(10) == P
            True
        """
        if q is None:
            # Main involution
            z, N, prec = self._z, self._N, self._prec
            return X0NPoint(-1/(N*z), N, prec)
        else:
            z, N, prec = self._z, self._N, self._prec
            q = ZZ(q)
            if N%q != 0:
                raise ValueError, "q must divide N"
            g, x, y = xgcd(q, -N//q)
            if g != 1:
                raise ValueError, "q must exactly divide N"
            # Now q*x - (N//q)*y = 1
            # W_q = [q*x,  y;  N,  q]
            return X0NPoint((q*x*z + y)/(N*z + q), N, prec)

    def __hash__(self):
        """
        We return 0 for the hash, which is useful so that we can use
        the data structures such as sets in Python with these points,
        though with a severe performance penalty.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import X0NPoint
            sage: X0NPoint(CDF(1,1), 12, 30).__hash__()
            0
        """
        return 0

def disk_to_h(q):
    """
    Given a floating point complex number q in the open unit disk,
    return a point tau in the upper half plane such that q =
    exp(2*pi*I*tau).

    NOTE: This function is just a scaled version of log, so it works
    for any nonzero complex number.

    INPUT:

    - q -- complex number in unit disk

    OUTPUT:

    - complex number in upper half plane

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.chow_heegner import disk_to_h, h_to_disk
        sage: disk_to_h(CDF(.5,.3))
        0.0860104348113 + 0.0858489451313*I
        sage: disk_to_h(ComplexField(100)(.5,.3))
        0.086010434811315337269743555113 + 0.085848945131318196784916564927*I
        sage: h_to_disk(disk_to_h(CDF(.5,.3)))
        0.5 + 0.3*I
    """
    K = q.parent()
    return q.log() / (2*K.pi()*K.gen())

def h_to_disk(z):
    """
    Given a floating point complex number z in the upper half plane,
    return the image of z under the map exp(2*pi*I*z), which is in the
    unit disk.

    NOTE: This function is just a scaled version of exp, so it works
    for any complex number.

    INPUT:

    - z -- complex number in upper half plane

    OUTPUT:

    - complex number in unit disk

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.chow_heegner import disk_to_h, h_to_disk
        sage: h_to_disk(CDF(0,2))
        3.48734235621e-06
        sage: disk_to_h(h_to_disk(CDF(0,2)))
        2.0*I
        sage: h_to_disk(CDF(2,1))
        0.00186744273171 - ...e-19*I

    Notice that the two maps are not mutual inverses, since the result
    of disk_to_h is always (approximately) in the vertical strip from
    0 to 1::

        sage: disk_to_h(h_to_disk(CDF(2,1)))
        -...e-17 + 1.0*I
        sage: disk_to_h(h_to_disk(CDF(0,1)))
        1.0*I
    """
    K = z.parent()
    return (2*K.pi()*K.gen()*z).exp()

class CloseEqual:
    """
    Object used in the throw_away_close function.

    TESTS::

        sage: from sage.schemes.elliptic_curves.chow_heegner import CloseEqual
        sage: c = CloseEqual(CDF(1,1), ComplexField(30))
        sage: TestSuite(c).run()
    """
    def __init__(self, x, C):
        """
        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import CloseEqual
            sage: c = CloseEqual(CDF(1,1), ComplexField(30)); c
            <sage.schemes.elliptic_curves.chow_heegner.CloseEqual instance at 0x...>
            sage: c.x, c.y
            (1.0 + 1.0*I, 1.0000000 + 1.0000000*I)
        """
        self.x = x
        self.y = C(x)

    def __hash__(self):
        """
        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import CloseEqual
            sage: c = CloseEqual(CDF(1,1), ComplexField(30)); c.__hash__()
            1000004
            sage: hash(c.y)
            1000004
        """
        return hash(self.y)

    def __cmp__(self, right):
        """
        Compare the images of x in C.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import CloseEqual
            sage: c = CloseEqual(CDF(1,1), ComplexField(30))
            sage: c2 = CloseEqual(CDF(1,1)+2^(-32), ComplexField(30))
            sage: c == c2
            True
            sage: c.x == c2.x
            False
        """
        return cmp(self.y, right.y)

def throw_away_close(v, prec):
    """
    Return list of 'distinct' elements of v, where we consider two
    elements close if they agree when coerced to complex numbers to
    prec bits of precision.

    INPUT:

    - `v` -- list of floating point numbers
    - ``prec`` -- integer; bits of precision to do comparisons

    OUTPUT:

    - sorted list of unique elements of v, where unique is defined
      as equality in ComplexField(prec)

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.chow_heegner import throw_away_close
        sage: v = [CDF(1,1), CDF(1,1)+2^(-32), CDF(1,1)+2^(-20)]

    All are distinct::

        sage: throw_away_close(v, 53)
        [1.0 + 1.0*I, 1.00000000023 + 1.0*I, 1.00000095367 + 1.0*I]

    Now v[1] and v[0] considered same::

        sage: throw_away_close(v, 30)
        [1.0 + 1.0*I, 1.00000095367 + 1.0*I]

    Now all the same::

        sage: throw_away_close(v, 20)
        [1.0 + 1.0*I]

    But not to slightly higher precision::

        sage: throw_away_close(v, 21)
        [1.0 + 1.0*I, 1.00000095367 + 1.0*I]
    """
    C = ComplexField(prec)
    w = [a.x for a in set([CloseEqual(x, C) for x in v])]
    w.sort()
    return w

def newton(f, x, max_iter=1000, max_err=1e-14):
    """
    Use the given number of steps of the Newton's algorithm to refine
    the given approximate root(s) x of f, where x is either a single
    root or a list of roots.  Returns list of refined roots, number of
    iterations, and an error bound.

    If the precision of base field of f is at most 53 and the
    coefficients are in fact real, then this function is particularly
    fast.

    INPUT:

    - `f` -- polynomial over a floating point complex field.
    - `x` -- an approximate root or *list* of approximate roots.
    - ``max_iter`` -- positive integer (default: 1000)
    - ``max_err`` -- small real (default: 1e-14)

    OUTPUT:

    - list of triples (r, n, err), where r is the root found via n
      iterations, where the iteration stopped because the
      difference between the n and n-1 value was err < max_err.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.chow_heegner import newton
        sage: R.<x> = CDF[]; f = x^3 + x - 3

    Starting value of 1 (not a list)::

        sage: t = newton(f, 1); t
        [(1.21341166276, 6, 0.0)]

    Resulting value is an approximate root as required::

        sage: f(t[0][0])
        -4.4408920985e-16

    In the following, we input three different starting values.  The
    first two converge to the same root, and the third to a different
    imaginary root::

        sage: t = newton(f, [1, 5, -.6-1.4*I]); t
        [(1.21341166276, 6, 0.0), (1.21341166276, 9, 0.0), (-0.606705831381 - 1.45061224919*I, 5, 0.0)]

    An example in which f has precision less than 53 bits::

        sage: R.<x> = ComplexField(30)[]; f = x^3 + x - 3
        sage: t = newton(f, [1]); t
        [(1.2134117, 6, 0.0)]
        sage: f(t[0][0])
        0.00000000

    An example in which f has complex coefficients but precision 53 bits::

        sage: from sage.schemes.elliptic_curves.chow_heegner import newton
        sage: R.<x> = CDF[]; f = x^3 + CDF.0*x - 3; f
        x^3 + I*x - 3.0
        sage: t = newton(f, 1); t
        [(1.44257744162 - 0.233096549477*I, 7, 2.77555756156e-17)]
        sage: f(t[0][0])
        4.4408920985e-16 - 1.11022302463e-16*I

    An example in which f has higher precision::

        sage: R.<x> = ComplexField(100)[]; f = x^3 + x - 3
        sage: t = newton(f, [1]); t
        [(1.2134116627622296341321313774, 6, 2.1019987680519080567000919351e-26)]

    Example with higher degree and precision::

        sage: R.<x> = ComplexField(500)[]; f = R([1..100])
        sage: t = newton(f, [1]); t
        [(-0.9485943966031...82, 34, 5.042484221...964289e-15)]

    However, the root we get isn't much better, due to not shrinking max_err::

        sage: f(t[0][0])
        -3.64698728...1825369e-26

    By shrinking max_err from the default, we get a much better answer::

        sage: t = newton(f, [1], max_err=2.0^(-400))
        sage: f(t[0][0])
        -2.443949090...136e-150

    Note that max_err is not a bound on f(approx_root).  Instead it is
    a bound on the difference of successive Newton approximations to
    the root.  For example, with max_err=0, we get::

        sage: t = newton(f, [1], max_err=0)
        sage: t[0][1]  # number of iterations
        39
        sage: f(t[0][0])      # error greater than 0
        -2.443949090...49136e-150
        sage: t[0][2]         # difference in successive approximations is 0 to our prec 500
        0.0000000...000000000
    """
    C = f.base_ring()

    if C.prec() <= 53:
        try:
            g = f.change_ring(RDF)
        except TypeError:
            # coefficients not real
            if C == CDF:
                f = f.change_ring(ComplexField(53))
        else:
            t = Polynomial_RDF_gsl(g).newton(x, max_iter=max_iter, max_err=max_err)
            if C != CDF:
                t = [(C(v[0]), v[1], v[2]) for v in t]
            return t

    if not isinstance(x, list):
        x = [x]

    f_prime = f.derivative()

    # ComplexPolynomial(...) is roughly 10 times faster to evaluate
    # than a usual Sage polynomial over a higher precision complex
    # field, so we use it instead.
    tm = verbose("Running Newton refinement on degree %s polynomial to precision %s on %s roots"%(
        f.degree(), C.prec(), len(x)))
    f = ComplexPolynomial(f); f_prime = ComplexPolynomial(f_prime)
    ans = []
    for root in x:
        root = C(root)
        last_root = root
        g_root = C(root)
        for i in range(max_iter):
            root = root - f(root)/f_prime(root)
            err = abs(last_root - root)
            if  err <= max_err:
                break
            last_root = root
        ans.append((root, i+1, err))
        verbose("Newton iterations: %s"%i)
    verbose("Completed %s Newton refinements in %s seconds."%(len(x), cputime(tm)))
    return ans

class NumericalPoint(object):
    """
    A wrapper around a numerical approximation to a point on an
    elliptic curve over QQ with several extra convenience functions.

    EXAMPLES::

        sage: Q = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b')).numerical_approx(min_imag=1e-3, deg1=130)
        sage: Q
        (6.00000000000... : 14.0000000000... : 1.00000000000000)
        sage: Q._eps
        0.0001
        sage: Q._P
        (6.00000000000... : 14.0000000000... : 1.00000000000000)

    TESTS::

        sage: TestSuite(Q).run()
    """
    def __init__(self, P, eps):
        """
        INPUT:

        - `P` -- actual point on some elliptic curve over a floating point base field
        - ``eps`` -- small positive real

        EXAMPLES::

            sage: Q = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b')).numerical_approx(min_imag=1e-3, deg1=100)
            sage: type(Q)
            <class 'sage.schemes.elliptic_curves.chow_heegner.NumericalPoint'>
        """
        self._P = P
        self._eps = eps

    def curve(self):
        """
        Return the curve over a floating point base field that this
        point is on.

        EXAMPLES::

            sage: Q = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b')).numerical_approx(min_imag=1e-3, deg1=100)
            sage: Q.curve()
            Elliptic Curve defined by y^2 + 1.00000000000000*y = x^3 + (-1.00000000000000)*x over Complex Field with 53 bits of precision
        """
        return self._P.curve()

    def rational_curve(self):
        """
        Return the elliptic curve that this point is on, but over the
        rational numbers.

        EXAMPLES::

            sage: Q = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b')).numerical_approx(min_imag=1e-3, deg1=100)
            sage: Q.rational_curve()
            Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        """
        return EllipticCurve(QQ,[int(a.real()) for a in self._P.curve().a_invariants()])

    def close_points(self, search_bound, eps=1e-3):
        """
        Return points close to self, sorted by distance (with closed
        point first).

        The elliptic curve must have rank 0 or 1, or a
        NotImplementedError is raised.

        INPUT:

        - ``search_bound`` -- we search through points of the form n*P + t,
          where P is a generator for the Mordell-Weil, t is any torsion point,
          and -search_bound<=n<=search_bound.

        - ``eps`` -- (default: 1e-3)

        OUTPUT:

        - list

        EXAMPLES::

            sage: Q = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b')).numerical_approx(min_imag=1e-3, deg1=100)
            sage: Q.close_points(100)
            [(6 : 14 : 1)]

        We try some absurdly large values of eps just to illustrate
        that we may see multiple points within eps of Q::

            sage: Q.close_points(10, eps=11)
            [(6 : 14 : 1)]
            sage: Q.close_points(10, eps=12.5)
            [(6 : 14 : 1), (2 : 2 : 1)]
            sage: Q.close_points(10, eps=13.5)
            [(6 : 14 : 1), (2 : 2 : 1), (0 : 1 : 0)]

        A current implementation limitation is that the curve must
        have rank at most 1::

            sage: E = EllipticCurve('389a'); P = E.change_ring(RDF)([0,0])
            sage: P
            (0.0 : 0.0 : 1.0)
            sage: from sage.schemes.elliptic_curves.chow_heegner import NumericalPoint
            sage: z = NumericalPoint(P, 1e-5)
            sage: z.close_points(10)
            Traceback (most recent call last):
            ...
            NotImplementedError: curve must have rank 0 or 1
        """
        E = self.rational_curve()
        g = E.gens()
        if len(g) == 0:
            P = E(0)
            search_bound=0
        elif len(g) == 1:
            P = g[0]
        else:
            raise NotImplementedError, "curve must have rank 0 or 1"
        T = E.torsion_points()
        v = []
        from sage.groups.generic import multiples
        for Q in multiples(P, 2*search_bound+1, -search_bound*P):
            for t in T:
                R = Q + t
                if abs(R[0] - self[0]) < eps and abs(R[1] - self[1]) < eps:
                    # record pair (distance, point)
                    v.append(  ((R[0] - self[0])**2 + (R[1] - self[1])**2, R)  )
        v.sort()
        return [R for _, R in v]

    def __getitem__(self, i):
        """
        INPUT:

        - i -- integer (0, 1 or 2)

        OUTPUT:

        - floating point number (x,y, or z coordinate)

        EXAMPLES::

            sage: Q = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b')).numerical_approx(min_imag=1e-3, deg1=130)
            sage: Q[0]
            6.00000000000...
            sage: Q[1]
            14.0000000000...
            sage: Q[2]
            1.00000000000000
        """
        return self._P.__getitem__(i)

    def __hash__(self):
        """
        We also return 0 for the hash.  It's convenient having a hash
        so that we can use the data structures such as sets in Python
        with these points, but this is never the bottlekneck in
        algorithms.  We must have all comparison go through the
        __cmp__ method.  By making the hash of everything 0, we force
        hash collisions to always happen and __cmp__ to always be used.

        EXAMPLES::

            sage: Q = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b')).numerical_approx(min_imag=1e-3, deg1=100)
            sage: hash(Q)
            0
        """
        return 0

    def __add__(self, right):
        """
        EXAMPLES::

            sage: Q = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b')).numerical_approx(min_imag=1e-3, deg1=130)
            sage: Q+Q
            (1.61355529131... : 1.184468407... : 1.00000000000000)
        """
        if isinstance(right, int) and right == 0:
            return self
        return NumericalPoint(self._P + self._P.curve().point(
            (right._P[0],right._P[1],right._P[2]),check=False), self._eps)

    def __radd__(self, left):
        """
        EXAMPLES::

            sage: Q = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b')).numerical_approx(min_imag=1e-3, deg1=130)
            sage: 0r+Q
            (6.00000000000... : 14.0000000000... : 1.00000000000000)
        """
        if isinstance(left, int) and left == 0:
            return self
        raise NotImplementedError

    def __rmul__(self, left):
        """
        EXAMPLES::

            sage: Q = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b')).numerical_approx(min_imag=1e-3, deg1=130)
            sage: 2*Q
            (1.61355529131... : 1.18446840788... : 1.00000000000000)
        """
        return NumericalPoint(left*self._P, self._eps)

    def __mul__(self, right):
        """
        EXAMPLES::

            sage: Q = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b')).numerical_approx(min_imag=1e-3, deg1=130)
            sage: Q*2
            (1.61355529131... : 1.18446840788... : 1.00000000000000)
        """
        return NumericalPoint(self._P*right, self._eps)

    def __cmp__(self, right):
        """
        EXAMPLES::

            sage: Q = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b')).numerical_approx(min_imag=1e-3, deg1=130)
            sage: S = 2*Q
            sage: Q
            (6.00000000000... : 14.0000000000... : 1.00000000000000)
            sage: S
            (1.61355529131... : 1.18446840788... : 1.00000000000000)
            sage: Q == Q
            True
            sage: S == S
            True
            sage: Q == S
            False

        Now make a point that is within eps, but clearly not equal, to check
        cmp in a more subtle case::

            sage: Q._eps
            0.0001
            sage: Z = Q.__class__(Q.curve().point([Q[0]+10*Q._eps, Q[1]+10*Q._eps, 1], check=False), Q._eps)
            sage: Q == Z
            False
            sage: Z = Q.__class__(Q.curve().point([Q[0]+.5*Q._eps, Q[1]+.5*Q._eps, 1], check=False), Q._eps)
            sage: Q == Z
            True
        """
        if max([abs(self._P[i]-right._P[i]) for i in range(3)]) < min(self._eps, right._eps):
            return 0
        return cmp(self._P, right._P)

    def __repr__(self):
        """
        EXAMPLES::

            sage: Q = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b')).numerical_approx(min_imag=1e-3, deg1=130)
            sage: Q.__repr__()
            '(6.00000000000... : 14.0000000000... : 1.00000000000000)'
        """
        return repr(self._P)


    def identify(self, search=200, eps=1e-5, infinity=1e12):
        """
        INPUT:

        - ``search`` -- (default: 200)
        - ``eps`` - (default: 1e-5)
        - ``infinity`` -- (default: 1e12)

        OUTPUT:

        - rational point on elliptic curve

        EXAMPLES::

            sage: Q = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b')).numerical_approx(min_imag=1e-3, deg1=130); Q
            (6.00000000000... : 14.0000000000... : 1.00000000000000)
            sage: Q.identify()
            (6 : 14 : 1)

        The index is 6, so searching through multiples only up to 5 of
        the generator of the Mordell-Weil group fails::

            sage: Q.identify(search=5)
            Traceback (most recent call last):
            ...
            ValueError: unable to identify rational point

        An example in which the point is 0::

            sage: c = EllipticCurve('57b').chow_heegner_point(EllipticCurve('57a')); c
            Chow-Heegner point on 57b1 associated to 57a1
            sage: c.numerical_approx().identify()
            (0 : 1 : 0)
        """
        if max(abs(self[0]),abs(self[1])) >= infinity or self._P == 0:
            return self.rational_curve()(0)

        m = self.close_points(search, eps)
        if len(m) != 1:
            raise ValueError, "unable to identify rational point"
        else:
            return m[0]



########################################################################
# Some helper functions needed for evaluation of the
# modular parametrization
########################################################################

def B_bound(ymin, prec):
    """
    Return an integer B so that using B terms of the series of a
    modular parameterization the tail end of the sum is bounded in
    absolute value by 2**(-prec) for any point with imaginary part at
    least ymin.

    INPUT:

    - ``ymin`` -- positive real number; minimum y coordinate
    - ``prec`` -- positive integer (bits of precision)

    OUTPUT:

    - integer

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.chow_heegner import B_bound
        sage: B_bound(1e-3,53)
        6765
        sage: B_bound(1e-5,53)
        749700
        sage: B_bound(1e-5,200)
        2371372

    We illustrate that the first bound above appears to be enough in
    one example, but using a smaller bound is not::

        sage: from sage.schemes.elliptic_curves.chow_heegner import phi_poly
        sage: phi = phi_poly(EllipticCurve([1..5]), 6765)
        sage: phi(CDF(e^(2*pi*i*(1+1e-3*I))))
        0.635422328062 + ...e-14*I
        sage: phi = phi_poly(EllipticCurve([1..5]), 20000)
        sage: phi(CDF(e^(2*pi*i*(1+1e-3*I))))
        0.635422328062 + ...e-14*I
        sage: phi = phi_poly(EllipticCurve([1..5]), 1000)
        sage: phi(CDF(e^(2*pi*i*(1+1e-3*I))))
        0.634984898068 + ...e-14*I
    """
    # It is important to use RR instead of RDF in this function
    # in order to get a larger range of exponents.
    y = RR(ymin)
    epsilon = RR(2)**(-(prec+1))
    pi = RR.pi()
    return int((epsilon*(1 - (-2*pi*y).exp())).log() / (-2*pi*y)) + 1

def phi_poly(E, B, base_field=QQ):
    """
    Return a polynomial over ``base_field`` that approximates the
    modular parametrization map associated to E.  This is the degree
    `B` polynomial `\sum a_n/n T^n`.

    INPUT:

    - `E` -- elliptic curve over the rational field
    - `B` -- a positive integer
    - ``base_field`` -- (default: QQ); a field

    OUTPUT:

    - a polynomial over the base field

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.chow_heegner import phi_poly
        sage: E = EllipticCurve([1..5])
        sage: phi_poly(E, 10)
        -3/10*q^10 - 1/3*q^9 - 3/8*q^8 - 1/7*q^7 - 3/5*q^5 - 1/4*q^4 + 1/2*q^2 + q
        sage: R.<q> = QQ[]
        sage: sum(E.an(n)/n*q^n for n in [1..10])
        -3/10*q^10 - 1/3*q^9 - 3/8*q^8 - 1/7*q^7 - 3/5*q^5 - 1/4*q^4 + 1/2*q^2 + q
        sage: phi_poly(E, 6, RDF)
        -0.6*q^5 - 0.25*q^4 + 0.5*q^2 + q
        sage: phi_poly(E, 4, RealField(100))
        -0.25000000000000000000000000000*q^4 + 0.50000000000000000000000000000*q^2 + q
    """
    R = base_field['q']
    v = E.anlist(B+1)
    return R([0] + [v[n]/n for n in range(1,B+1)])

def label(E):
    """
    Return the Cremona label of E if it is in the installed database,
    and otherwise just return the string representation of E. This
    function is used by some of the other printing functions in this
    module.

    INPUT:

    - `E` -- elliptic curve over the rational numbers

    OUTPUT:

    - string

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.chow_heegner import label
        sage: label(EllipticCurve('11a3'))
        '11a3'
        sage: label(EllipticCurve([0,2012]))
        'Elliptic Curve defined by y^2 = x^3 + 2012 over Rational Field'
    """
    try:
        return E.cremona_label()
    except RuntimeError:
        return str(E)

def check_optimal(E):
    """
    If E is in the installed database, verify that E is actually
    an optimal curve, and do nothing if E is not in the databse.
    If E is not optimal, raise a ValueError.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.chow_heegner import check_optimal
        sage: check_optimal(EllipticCurve('11a1'))
        sage: check_optimal(EllipticCurve('11a2'))
        Traceback (most recent call last):
        ...
        ValueError: curve must be optimal

    At least one of these curves is not optimal, but they have big
    conductor so they are not in the database::

        sage: E = EllipticCurve([0,-10001,0,10000,0]); E.conductor()
        266640
        sage: C = E.isogeny_class()
        sage: len(C)
        8
        sage: [check_optimal(X) for X in C]
        [None, None, None, None, None, None, None, None]

        sage: check_optimal(EllipticCurve('990h3'))
        True
        sage: check_optimal(EllipticCurve('990h1'))
        False
    """
    try:
        lbl = E.cremona_label()
        if lbl.startswith('990h'):
            return lbl.endswith('3')
        if not lbl.endswith('1'):
            raise ValueError, "curve must be optimal"
    except RuntimeError:
        pass

class ModularParametrization(object):
    """
    Modular parametrization map designed for the needs of the
    Chow-Heegner point algorithm.  There is another general modular
    parametrization class in Sage (in modular_parametrization.py), but
    it is not suited for this computation (wrong methods and
    semantics).

    The input curve is assumed optimal.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.chow_heegner import ModularParametrization
        sage: phi = ModularParametrization(EllipticCurve('37b'))
        sage: type(phi)
        <class 'sage.schemes.elliptic_curves.chow_heegner.ModularParametrization'>
        sage: phi
        Modular parametrization of 37b1 having degree 2
        sage: phi.elliptic_curve()
        Elliptic Curve defined by y^2 + y = x^3 + x^2 - 23*x - 50 over Rational Field
        sage: phi.degree()
        2
        sage: phi.cusps_over_torsion()
        {1: [Infinity], 3: [0]}
        sage: phi.images_of_cusps()
        {(0, 1): [0], (0, 0): [Infinity]}
        sage: phi(CDF(1,1))
        0.00186744489643 - ...e-19*I
        sage: phi(CDF(2,1))
        0.00186744489643 - ...e-19*I

    Find points in the fiber over the point 1 mod the period lattice::

        sage: v = phi.points_in_h(1, 1e-3); v
        [[-0.3512 + 0.001870*I], [-0.04588 + 0.01936*I]]
        sage: phi(v[0])
        1.0 ...
        sage: phi(v[1])
        1.0
    """
    def __init__(self, E):
        """
        INPUT:

        - `E` - elliptic curve over QQ assumed to be optimal

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import ModularParametrization
            sage: ModularParametrization(EllipticCurve('37b1'))
            Modular parametrization of 37b1 having degree 2

        The input curve must be optimal::

            sage: ModularParametrization(EllipticCurve('11a2'))
            Traceback (most recent call last):
            ...
            ValueError: curve must be optimal

        The input must be an elliptic curve over the rational numbers, so this works::

            sage: ModularParametrization(EllipticCurve([1/2,19/4]))
            Modular parametrization of Elliptic Curve defined by y^2 = x^3 + 1/2*x + 19/4 over Rational Field having degree 127872

        But these don't work::

            sage: ModularParametrization('37b1')
            Traceback (most recent call last):
            ...
            TypeError: E must be an elliptic curve
            sage: ModularParametrization(EllipticCurve(GF(7),[3,4]))
            Traceback (most recent call last):
            ...
            TypeError: E must be over QQ
        """
        if not is_EllipticCurve(E):
            raise TypeError, "E must be an elliptic curve"
        if not is_RationalField(E.base_field()):
            raise TypeError, "E must be over QQ"
        check_optimal(E)
        self._E = E
        self._label = label(E)

    def elliptic_curve(self):
        """
        Return the elliptic curve of which this is the modular
        parametrization.

        OUTPUT:

        - an elliptic curve over the rational numbers

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import ModularParametrization
            sage: phi = ModularParametrization(EllipticCurve('389a'))
            sage: phi.elliptic_curve()
            Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field
        """
        return self._E

    def degree(self):
        """
        Return the degree of this modular parametrization.

        OUTPUT:

        - positive integer

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import ModularParametrization
            sage: phi = ModularParametrization(EllipticCurve('389a'))
            sage: phi.degree()
            40
            sage: type(phi.degree())
            <type 'sage.rings.integer.Integer'>
            sage: phi = ModularParametrization(EllipticCurve([19,13]))
            sage: phi.elliptic_curve().conductor()
            127996
            sage: phi.degree()
            11400
        """
        return self._E.modular_degree()

    def cusps_over_torsion(self):
        """
        Return dictionary with keys integers n and values the cusps
        that map to a torsion point on the elliptic curve of exact
        order n.  This is computed using modular symbols, so is
        definitely correct, but will be slow if the conductor is at
        all large.

        OUTPUT:

        - dictionary

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import ModularParametrization

        A well known example -- level 11, where the cusp 0 maps to a 5-torsion point::

            sage: phi = ModularParametrization(EllipticCurve('11a'))
            sage: phi.cusps_over_torsion()
            {1: [Infinity], 5: [0]}

        A higher level example, where there is no nontrivial rational torsion::

            sage: phi = ModularParametrization(EllipticCurve('57a'))
            sage: phi.cusps_over_torsion()
            {1: [0, 1/19, 1/3, Infinity]}
            sage: phi.elliptic_curve().torsion_order()
            1

        An example with 2-torsion::

            sage: phi = ModularParametrization(EllipticCurve('57b'))
            sage: phi.cusps_over_torsion()
            {1: [Infinity], 2: [0, 1/19, 1/3]}
            sage: phi.elliptic_curve().torsion_subgroup().invariants()
            (2, 2)
        """
        phi = self._E.modular_symbol_space(0).integral_period_mapping()
        v = {}
        for c in Gamma0(self._E.conductor()).cusps():
            i = phi([c,infinity])
            d = i.denominator()
            if v.has_key(d):
                v[d].append(c)
            else:
                v[d] = [c]
        return v

    def images_of_cusps(self):
        """
        Return dictionary with keys elements of a finite module
        isomorphic to E[d] = (Z/dZ)x(Z/dZ).  The corresponding values
        for the keys are the cusps that map to that element of E[d].
        Here d is the exponent of the finite group generated by
        images of cusps.

        OUTPUT:

        - dictionary

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import ModularParametrization
            sage: phi = ModularParametrization(EllipticCurve('57b'))
            sage: phi.images_of_cusps()
            {(0, 1): [0], (1, 0): [1/19], (0, 0): [Infinity], (1, 1): [1/3]}
            sage: Z = phi.images_of_cusps(); Z
            {(0, 1): [0], (1, 0): [1/19], (0, 0): [Infinity], (1, 1): [1/3]}
            sage: Z.keys()[0].parent()
            Finitely generated module V/W over Integer Ring with invariants (2, 2)
            sage: phi = ModularParametrization(EllipticCurve('11a'))
            sage: phi.images_of_cusps()
            {(0, 0): [Infinity], (4, 0): [0]}
            sage: phi = ModularParametrization(EllipticCurve('57a'))
            sage: phi.images_of_cusps()
            {(): [0, 1/19, 1/3, Infinity]}
        """
        phi = self._E.modular_symbol_space(0).integral_period_mapping()
        d = ZZ(1)
        ims = []
        for c in Gamma0(self._E.conductor()).cusps():
            i = phi([c,infinity])
            d = d.lcm(i.denominator())
            ims.append((c,i))
        V = (ZZ**2)
        Q = V / V.scale(d)
        v = {}
        for c, i in ims:
            P = Q(d*i)
            if v.has_key(P):
                v[P].append(c)
            else:
                v[P] = [c]
        return v

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import ModularParametrization
            sage: ModularParametrization(EllipticCurve('37a')).__repr__()
            'Modular parametrization of 37a1 having degree 2'
        """
        return "Modular parametrization of %s having degree %s"%(self._label, self.degree())

    def __call__(self, z):
        """
        Compute the image of z under this modular parametrization map.
        We allow z to be either a point in the upper half plane or a
        list of such points.

        INPUT:

        - a complex number or list of complex numbers in the upper half plane

        OUTPUT:

        - complex number or list of complex numbers

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import ModularParametrization
            sage: phi = ModularParametrization(EllipticCurve('11a'))
            sage: phi(CDF(1,1))
            0.00186395322463 - ...e-19*I
            sage: phi(ComplexField(100)(1,1))
            0.0018639532246330691303545022211 - 6.3095154969118369542516994420e-34*I
            sage: phi([CDF(1,1), CDF(2,1), CDF(2,1/2)])
            [0.00186395322463..., 0.00186395322463..., 0.0413213515947...]

        A consistency check with the other modular parametrization in Sage::

            sage: phi0 = EllipticCurve('11a').modular_parametrization()
            sage: P = phi0(CDF(1,1)); P
            (287826.305812822...: -1.54416950808282e8... : 1.00000000000000)
            sage: z = phi(CDF(1,1)); z
            0.00186395322463...
            sage: E = phi.elliptic_curve()
            sage: E.elliptic_exponential(z)
            (287826.305812822...: -1.54416950808282e8... : 1.00000000000000)
        """
        if isinstance(z, list):
            if len(z) == 0:
                return []
            z = [(x.z() if isinstance(x,X0NPoint) else x) for x in z]
            d = max([B_bound(x.imag(), x.prec()) for x in z])
            is_list = True
        else:
            if isinstance(z, X0NPoint):
                z = z.z()
            z = [z]
            d = B_bound(z[0].imag(), z[0].prec())
            is_list = False
        f = phi_poly(self._E, d, base_field=z[0].parent())
        if z[0].prec() > 53:
            f = ComplexPolynomial(f)
        else:
            f = Polynomial_RDF_gsl(f)
        w = []
        for x in z:
            q = h_to_disk(x)
            w.append(f(q))
        if not is_list:
            return w[0]
        return w

    def _points_in_h_double(self, z, min_imag, max_iter=100, deg1=500):
        """
        Find points in upper half plane over double precision z.

        INPUT:

        - `z` -- (default: 0.1); floating point number that
          defines element of period lattice. WARNING: function wil
          be dramatically slower if z has nonzero imaginary part.
        - ``min_imag`` -- (default: 1e-4) positive real number; we
          only find points in the upper half plane with imaginary
          part at least this
        - ``max_iter`` -- maximum number of iterations used in
          newton calls
        - ``deg1`` -- (default: 500) degree used for first double
          precision root finding; if too small then some points
          will be missed.  Reasonable values are between 100 and 3000.

        OUTPUT:

        - list of CDF elements

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import ModularParametrization
            sage: phi = ModularParametrization(EllipticCurve('57b'))
            sage: v = phi._points_in_h_double(z=.1, min_imag=1e-3); v
            [-0.448258381107 + 0.00603424481535*I, -0.381277139742 + 0.00147982187514*I, -0.303618734588 + 0.0143125900531*I, 0.374289971606*I, 0.303618734588 + 0.0143125900531*I, 0.381277139742 + 0.00147982187514*I, 0.448258381107 + 0.00603424481535*I]
            sage: [phi(t) for t in v]
            [0.1..., 0.1..., 0.1..., 0.1, 0.1..., 0.1..., 0.1...]
        """
        z = CDF(z)
        f = phi_poly(self._E, deg1, base_field=CDF) - z
        try:
            v = [x for x in cdf_roots_of_rdf_poly(f) if abs(x) < 1]
        except TypeError:
            v = [x for x,_ in f.roots() if abs(x) < 1]
        verbose('Number of double precision roots in upper half plane: %s'%len(v))
        if len(v) == 0:
            return []

        # use actual minimum imaginary part found
        w = [disk_to_h(x).imag() for x in v]
        w = [x for x in w if x >= min_imag]
        if len(w) == 0:
            return []
        min_imag = min(w)

        f = phi_poly(self._E, B_bound(min_imag, 53), base_field=CDF) - z
        t = verbose("deg of refinement poly = %s"%f.degree())
        w = []
        if z.imag() == 0:
            roots = Polynomial_RDF_gsl(f).newton(v, max_iter)
        else:
            roots = newton(f, v, max_iter=max_iter)

        for b,i,err in roots:
            if abs(b) < 1 and i < max_iter:
                w.append(b)

        verbose("found %s double prec roots that refined in %s seconds"%(len(w), cputime(t)))
        w.sort()
        w = throw_away_close(w, prec=45)

        # Put points in upper half plane.
        w = [disk_to_h(z) for z in w]
        w = [z for z in w if z.imag() >= min_imag]  # 1000 is effectively oo
        w.sort()
        return w

    def points_in_h(self, z=CDF(0.1), min_imag=1e-4,
                    max_iter=25, deg1=500,
                    max_iter1=100, equiv_prec=None):
        r"""
        Return `\Gamma_0(N)`-inequivalent points tau in `X_0(N)` to
        precision z.prec(), represented as points in the open upper
        half plane, with tau.imag() >= min_imag; these points map to z
        via the modular parametrization map.  There is no guarantee at
        that representatives for all such points have been returned.

        INPUT:

        - `z` -- (default: 0.1); floating point number that
          defines element of period lattice. WARNING: Function wil
          be dramatically slower if z has nonzero imaginary part.
        - ``min_imag`` -- (default: 1e-4) positive real number; we
          only find points in the upper half plane with imaginary
          part at least this
        - ``max_iter`` -- maximum number of iterations used in
          newton iteration to refine roots to higher than 53 bits
          precision, if z has higher than 53 bits precision.
        - ``deg1`` -- (default: 500) degree used for first double
          precision root finding; if too small then some points
          will be missed.  Reasonable values are between 100 and 3000.
        - ``max_iter1`` -- (default: 100) max number of iterations
          used in refining initial roots to double precision
        - ``equiv_prec`` -- (default: prec//3, where prec is
          precision of z) used in determining Gamma_0(N)
          equivalence of points

        OUTPUT:

        - list of inequivalent numerical points in the upper half
          plane that define elements of the modular curve `X_0(N)`

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner import ModularParametrization
            sage: phi = ModularParametrization(EllipticCurve('37b'))

        Defaults to points over 0.1::

            sage: v = phi.points_in_h(deg1=100); v
            [[-0.3441 + 0.008841*I], [-0.2982 + 0.001291*I]]
            sage: phi(v[0]), phi(v[1])
            (0.1 ..., 0.1 ...)

        Compute points over 1::

            sage: v = phi.points_in_h(RDF(1), deg1=100); v
            [[-0.3512 + 0.001870*I], [-0.1683 + 0.001018*I]]
            sage: phi(v[0]), phi(v[1])
            (1.0 ..., 1.0 ...)

        Points over I::

            sage: v = phi.points_in_h(CDF(0,1), deg1=100); v    # long time
            [[-0.3443 + 0.0008692*I], [0.1500 + 0.007604*I]]
            sage: phi(v[0]), phi(v[1])                          # long time
            (...e-1... + 1.0*I, ...e-1... + 1.0*I)

        To higher precision::

            sage: v = phi.points_in_h(RealField(100)(1), min_imag=1e-3); v
            [[-0.35115510977014 + 0.0018701266234958*I], [-0.045882897517034 + 0.019362756466571*I]]
            sage: phi(v[0])
            1.0000000000000000000000000000 + 3.3526588471893001729998464025e-29*I
            sage: phi(v[1])
            1.0000000000000000000000000000 + 2.7610131682735413189410499785e-30*I

        Consider points with smaller imaginary parts::

            sage: phi.points_in_h(RDF(1), min_imag=1e-4)
            [[-0.3574 + 0.0004481*I], [-0.3512 + 0.001870*I]]

        Change the degree of the initial polynomial approximation to
        the modular parametrization map::

            sage: phi.points_in_h(RDF(1), deg1=100, min_imag=1e-3)
            [[-0.3512 + 0.001870*I], [-0.04588 + 0.01936*I]]

        Use a different bound on Newton iterations::

            sage: phi.points_in_h(RDF(1), max_iter1=10, min_imag=1e-3)
            [[-0.3512 + 0.001870*I], [-0.04588 + 0.01936*I]]

        Consider points as equivalent if they are equivalent to very
        low precision (silly example)::

            sage: phi.points_in_h(RDF(1), equiv_prec=2, min_imag=1e-3)
            [[-0.38 + 0.0020*I], [-0.094 + 0.0015*I], [-0.047 + 0.016*I]]

        Consider points equivalent only if equivalent to high
        precision, which will result in falsely considering points as
        inequivalent::

            sage: len(phi.points_in_h(RDF(1), equiv_prec=50, min_imag=1e-3))
            7
            sage: len(phi.points_in_h(RDF(1), equiv_prec=40, min_imag=1e-3))  # correct, since moddeg=2
            2
        """
        C = parent(z)
        if not (C == CDF or is_ComplexField(C)):
            try:
                prec = z.prec()
            except AttributeError:
                prec = 53
            C = CDF if prec==53 else ComplexField(prec)
            z = C(z)

        N = self._E.conductor()

        v = self._points_in_h_double(z, min_imag=min_imag, max_iter=max_iter1, deg1=deg1)
        if len(v) == 0:
            return []

        m_E = self.degree()

        prec = z.prec()
        if equiv_prec is None:
            equiv_prec = prec//3
        if prec <= 53:
            # double precision already enough
            ans = list(set([X0NPoint(C(x),N,prec=equiv_prec) for x in v]))
            ans.sort()
            if len(ans) < m_E:
                raise RuntimeError, "did not find enough points in the preimage. "
            return ans


        # refine to higher precision, and keep only good roots
        verbose("Number of double precision roots to refine via Newton iteration: %s"%len(v))
        t = cputime()
        f = phi_poly(self._E, B_bound(min_imag, prec), base_field=C) - z

        # We try m_E points at a time, since the newton step below is expensive.
        w0 = []; w = []
        #while len(w) < m_E:
        while len(v) > 0:
        # changed so that it always terminates
            for b,i,err in newton(f, [h_to_disk(C(a)) for a in v[:m_E]],
                                            max_iter, max_err=C(2)**(2-prec)):
                if abs(b) < 1 and i < max_iter:
                    w0.append(b)
            v = v[m_E:]
            w0 = throw_away_close(w0, prec=prec-10)

            # transform and take distinct X0 points in upper half plane.
            w = [disk_to_h(z) for z in w0]
            w = [z for z in w if z.imag() >= min_imag and z.imag() <= 1]  # above 1 --> is point at oo
            w = [X0NPoint(z, N, prec=prec//2) for z in w]
            w = list(set(w))

        verbose("found %s roots to prec %s that refined in %s seconds"%(len(w), prec, cputime(t)))
        w.sort()
        if len(w) < m_E:
            raise RuntimeError, "did not find enough points in the preimage, try changing parameters"
        return w


class ChowHeegnerPoint(object):
    """
    A Chow-Heegner point associated to an ordered pair of optimal
    elliptic curves of the same conductor.

    EXAMPLES::

        sage: P = EllipticCurve('57a').chow_heegner_point(EllipticCurve('57b')); P
        Chow-Heegner point on 57a1 associated to 57b1
        sage: P.numerical_approx()
        (1.44444444444440...: -1.03703703703... : 1.00000000000000)
        sage: P.point_exact()
        (13/9 : -28/27 : 1)

    You can recover the two curves used to define the point::

        sage: P.curves()
        (Elliptic Curve defined by y^2 + y = x^3 - x^2 - 2*x + 2 over Rational Field, Elliptic Curve defined by y^2 + x*y + y = x^3 - 7*x + 5 over Rational Field)
    """
    def __init__(self, E, F):
        """
        Initialize Chow-Heegner point associated to a pair of curves
        of the same conductor.

        INPUT:

        - `E`, `F` -- optimal elliptic curve over the rational
          numbers, or Cremona labels for elliptic curves

        EXAMPLES::

            sage: E = EllipticCurve('37a'); F = EllipticCurve('37b')
            sage: E.chow_heegner_point(F)
            Chow-Heegner point on 37a1 associated to 37b1

            sage: from sage.schemes.elliptic_curves.chow_heegner import ChowHeegnerPoint
            sage: ChowHeegnerPoint(E, F)
            Chow-Heegner point on 37a1 associated to 37b1
            sage: ChowHeegnerPoint(F, E)
            Chow-Heegner point on 37b1 associated to 37a1

        We allow E and F to be Cremona labels::

            sage: ChowHeegnerPoint('37a', '37b')
            Chow-Heegner point on 37a1 associated to 37b1

        The two curves must have the same conductor.  This is an
        implementation issue, not a theoretical issue, which is why a
        NotImplementedError is raised::

            sage: ChowHeegnerPoint(E, EllipticCurve('57b'))
            Traceback (most recent call last):
            ...
            NotImplementedError: E and F must currently have the same conductor

        The curves must also be optimal::

            sage: ChowHeegnerPoint(EllipticCurve('37a1'), EllipticCurve('37b2'))
            Traceback (most recent call last):
            ...
            ValueError: curve must be optimal
            sage: ChowHeegnerPoint(EllipticCurve('37b2'), EllipticCurve('37a1'))
            Traceback (most recent call last):
            ...
            ValueError: curve must be optimal

        They must be distinct and non-isogenous::

            sage: ChowHeegnerPoint(EllipticCurve('11a1'), EllipticCurve('11a1'))
            Traceback (most recent call last):
            ...
            ValueError: E and F must not be isomorphic

        The curves must be over the rational numbers::

            sage: E = EllipticCurve(QuadraticField(-1),[1..5])
            sage: ChowHeegnerPoint(E, E)
            Traceback (most recent call last):
            ...
            ValueError: E and F must be elliptic curves over QQ
        """
        if isinstance(E, str): E = EllipticCurve(E)
        if isinstance(F, str): F = EllipticCurve(F)

        if not is_RationalField(E.base_ring()) or not is_RationalField(F.base_ring()):
            raise ValueError, "E and F must be elliptic curves over QQ"
        if E.conductor() != F.conductor():
            raise NotImplementedError, "E and F must currently have the same conductor"
        if E.is_isomorphic(F):
            raise ValueError, "E and F must not be isomorphic"
        check_optimal(E)
        check_optimal(F)
        self._E = E
        self._F = F
        self._fE = ModularParametrization(E)
        self._fF = ModularParametrization(F)

    def __cmp__(self, right):
        """
        We compare the ordered pair of underlying elliptic curves.

        EXAMPLES::

            sage: P = EllipticCurve('57a').chow_heegner_point(EllipticCurve('57b'))
            sage: P == P
            True
            sage: Q = EllipticCurve('57b').chow_heegner_point(EllipticCurve('57a'))
            sage: P == Q
            False
            sage: P < Q
            True
            sage: Q < P
            False
        """
        assert isinstance(right, ChowHeegnerPoint)
        return cmp((self._E, self._F), (right._E, right._F))

    def __repr__(self):
        """
        EXAMPLES::

            sage: P = EllipticCurve('446a').chow_heegner_point(EllipticCurve('446b'))
            sage: P.__repr__()
            'Chow-Heegner point on 446a1 associated to 446b1'
        """
        return "Chow-Heegner point on %s associated to %s"%(label(self._E), label(self._F))

    def curves(self):
        """
        Return the two curves that this Chow-Heegner point is associated to.

        OUTPUT:

        - 2-tuple of elliptic curves over QQ

        EXAMPLES::

            sage: P = EllipticCurve('446a').chow_heegner_point(EllipticCurve('446b'))
            sage: E, F = P.curves()
            sage: E.cremona_label(), F.cremona_label()
            ('446a1', '446b1')
        """
        return self._E, self._F

    def parametrizations(self):
        """
        Return the modular parametrization maps associated to E and F.

        OUTPUT:

        - 2-tuple of modular parametrization maps

        EXAMPLES::

            sage: P = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b'))
            sage: P.parametrizations()
            (Modular parametrization of 37a1 having degree 2, Modular parametrization of 37b1 having degree 2)
        """
        return self._fE, self._fF

    def points_over_F(self, *args, **kwds):
        """
        Return some points on `X_0(N)` over a given real point on F
        represented as a real number modulo the period lattice.

        The input and output options are the same as those for the
        method ``points_in_h`` of the modular parametrization maps.
        Use self.parametrizations()[1] to get one of the maps, and
        look at the docstring.

        EXAMPLES::

            sage: P = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b'))
            sage: v = P.points_over_F(1/4, deg1=100); v
            [[-0.3465 + 0.008016*I], [-0.2989 + 0.001326*I]]
            sage: phi = P.parametrizations()[1]
            sage: phi(v[0]), phi(v[1])
            (0.25..., 0.25...)

        Here's an example of getting the points, but to higher precision::

            sage: v = P.points_over_F(RealField(100)(1/4), deg1=100); v  # long time
            [[-0.34645302134160 + 0.0080165180402406*I], [-0.29890732134688 + 0.0013263630385348*I]]
            sage: phi(v[0]), phi(v[1])      # long time
            (0.25000000000000000000000000000 + 3.2540512340366736973233803318e-30*I, 0.25000000000000000000000000004 + 1.2621774483536188886587657045e-29*I)
        """
        return self._fF.points_in_h(*args, **kwds)

    def numerical_approx(self, *args, **kwds):
        """
        Return a numerical approximation to this Chow-Heegner point.

        This function takes exactly the same inputs as the compute
        method, so see the docstring for that method.

        EXAMPLES::

            sage: P = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b'))
            sage: P.numerical_approx(deg1=130)
            (6.00000000000... : 14.0000000000... : 1.00000000000000)
            sage: P.numerical_approx(60, deg1=130)
            (6.00000000000... : 14.000000000000... : 1.0000000000000000)
            sage: P = EllipticCurve('37b').chow_heegner_point(EllipticCurve('37a'))
            sage: P.numerical_approx(deg1=100)
            (8.00000000000... : 18.0000000000... : 1.00000000000000)
        """
        P, _, _ = self.compute(*args, **kwds)
        if P is not None:
            return P
        else:
            raise RuntimeError, "failed to compute"

    @cached_method
    def point_exact(self, search=100, eps=1e-5, infinity=1e12, **kwds):
        """
        Return one of the best exact approximations to this point that
        is a bounded multiple of the generators, and has x and y
        coordinate that is at least within eps of the numerical
        approximation of this point, or the point infinity (see
        below).

        INPUT:

        - ``search`` -- (default: 100) search through all points
          of the form t + sum n*g, where t is torsion, g is a
          generator, and abs(n) <= search.
        - ``eps`` -- (default: 1e-5) a small float; x and y
          coordinates must be at least this close
        - ``infinity`` -- (default: 1e12) a large float; if either
          x or y coordinate of numerical approximation is at least
          this big, then we just declare this the point at
          infinity, and do not do a search.
        - ``**kwds`` -- all other keyword arguments are passed to
          the ``numerical_approx`` method.

        EXAMPLES::

            sage: P = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b'))
            sage: P.point_exact(deg1=100)
            (6 : 14 : 1)

        Another example that illustrates precision issues when testing
        equivalence of points in the upper half plane::

            sage: P = EllipticCurve('99a').chow_heegner_point(EllipticCurve('99b'))
            sage: P.point_exact()
            Traceback (most recent call last):
            ...
            RuntimeError: Found too many points (13 > 12) -- try (good) increasing precision of z (now=53) or (bad) *decreasing* equiv_prec (now=17)

        As suggested, either increasing the precision (slower, but
        safer) or somewhat decreasing the number equiv_prec of bits
        used for checking equivalence (faster, but less safe) both
        work in this case::

            sage: P.point_exact(equiv_prec=12)
            (105/64 : -897/512 : 1)
            sage: P.point_exact(prec=60)
            (105/64 : -897/512 : 1)
            sage: P.index(equiv_prec=12)
            4
        """
        P = self.numerical_approx(**kwds)
        return P.identify(search=search, eps=eps, infinity=infinity)

    def index(self, *args, **kwds):
        """
        Index of this Chow-Heegner point in the Mordell-Weil group
        modulo torsion.

        INPUT:

        - exactly the same as for the point_exact method

        OUTPUT:

        - Integer

        EXAMPLES::

            sage: P = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b'))
            sage: P.index(min_imag=1e-3, deg1=100)
            6
        """
        P = self.point_exact(*args, **kwds)
        h = P.height()
        E = self._E
        return ZZ(int((h/E.regulator()).sqrt().round()))

    @cached_method
    def compute(self, prec=53, z='0.1',
                equiv_prec=None, deg1=500, min_imag=1e-4,
                fiber=None, number_to_stabilize=6):
        """
        INPUT:

        - ``prec`` -- (default: 53) precision of numerical
          computations (should be at least 53)
        - `z` -- (default: 0.1) base point representative
        - ``equiv_prec`` -- (default: prec//3) precision for
          checking Gamma0(N) equivalence
        - ``deg1`` -- (default: 500) degree used for first double
          precision root finding; if too small then some points
          will be missed.  Reasonable values are between 100 and 3000.
        - ``min_imag`` -- (default: 1e-4) positive real number; we
          only find points in the upper half plane with imaginary
          part at least this
        - ``fiber`` -- (default: None); if given as a list, then
          we assume that the points in this list were the output
          from a previous run of the compute command, and start
          with this list, instead of starting from scratch
        - ``number_to_stabilize`` -- (default: 6); if after this
          many attempts (using varying equivalent base points) we
          find no further points in fiber, then this function
          terminates; if not enough points were found, the output
          point P is None.

        OUTPUT:

        - 3-tuple with entries:

            - ``P`` -- None, or the numerically computed Chow-Heegner
              point if we found enough points in the fiber
            - ``fiber`` -- the inequivalent points we found in the fiber
            - ``base_points`` -- the real numbers that define points
              on F that we computed the fiber over; these points are
              all equivalent to z modulo the period lattice

        EXAMPLES::

            sage: P = EllipticCurve('37a').chow_heegner_point(EllipticCurve('37b'))
            sage: m, fiber, base_points = P.compute(deg1=130, min_imag=1e-3)
            sage: m
            (6.00000000000... : 14.0000000000... : 1.00000000000000)
            sage: fiber
            [[-0.3441 + 0.008841*I], [-0.2982 + 0.001291*I]]
            sage: base_points
            [0.1]

        We compute (and get the same answer) using a complex base
        point, which gives us increased confidence in our result::

            sage: P.compute(z=CDF(1,1), min_imag=1e-3, deg1=130)
            ((6.00000000000... : 14.0000000000... : 1.00000000000000), [[-0.04423 + 0.001498*I], [0.08721 + 0.001363*I]], [1.0 + 1.0*I])
        """
        if prec == 53:
            z = CDF(z)
        else:
            z = ComplexField(prec)(z)

        if equiv_prec is None:
            equiv_prec = z.prec()//3

        if fiber is None:
            fiber = []
        else:
            fiber = list(fiber)

        base_points = [z]
        mF = self._fF.degree()
        verbose("Modular degree of F (target fiber size) = %s"%mF)
        target = mF
        Omega = self._F.period_lattice().basis(z.prec())[0]
        n = 0
        if z == 0:
            # reduce target by number of cusps that map to 0
            C = self._fF.cusps_over_torsion()[1]
            target -= len(C)
            verbose("Reducing target by %s cusps to %s (TODO *** result may be off by torsion! *** )"%(len(C), target))

        same_in_a_row = 0
        while len(fiber) < target:
            b = len(fiber)
            tm = verbose("** deg1 = %s, min_imag = %s, equiv_prec = %s, len(fiber) = %s, target = %s"%(deg1, min_imag, equiv_prec, b, target))
            for t in self.points_over_F(z, min_imag=min_imag, deg1=deg1, equiv_prec=equiv_prec):
                fiber.append(t)
            fiber = list(set(fiber))
            verbose("** Found %s new points in %s seconds"%(len(fiber)-b, cputime(tm)))
            if len(fiber) == b:
                same_in_a_row += 1
                if same_in_a_row >= number_to_stabilize:
                    verbose("*"*80)
                    verbose("** ERROR: same result %s times in a row, so 'stabilized' (GIVING UP!)"%same_in_a_row)
                    verbose("*"*80)
                    stabilize = True
                    break
            else:
                same_in_a_row = 0
            if len(fiber) < target:
                verbose('*'*80)
                # change base point
                n += 1
                z += ((-1)**(n+1)) * n*Omega
                base_points.append(z)
                verbose("z (=%s) |--> %s"%(base_points[-1], z))
            elif len(fiber) == target:
                verbose('+'*80)
                verbose("*** Found all %s points in the fiber ***"%target)
            else:
                verbose('-'*80)
                raise RuntimeError, "Found too many points (%s > %s) -- try (good) increasing precision of z (now=%s) or (bad) *decreasing* equiv_prec (now=%s)"%(
                    len(fiber),target,z.prec(),equiv_prec)

        if len(fiber) == target:
            verbose("Mapping points on modular curve to E...")
            t = cputime()
            m = self._fE(fiber)
            verbose("Mapped to E in %s seconds"%cputime(t))
            s = sum(m)
            if isinstance(s, int) and s == 0:
                P = NumericalPoint(self._E(0), 1e-4)
            else:
                P = NumericalPoint(self._E.elliptic_exponential(s), 1e-4)
        else:
            P = None

        return P, fiber, base_points


