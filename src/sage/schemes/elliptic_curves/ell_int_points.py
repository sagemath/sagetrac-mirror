r"""
Integral points of elliptic curves over number fields

The following examples are from :trac:`10152`. Previously Magma was
finding more integral points on these curves than Sage. We resolved
this issue.

EXAMPLES::

    sage: E = EllipticCurve('2082a1')
    sage: E.integral_points(algorithm = "old")
    [(-11 : 29 : 1), (-2 : 29 : 1), (4 : 11 : 1), (13 : 29 : 1)]
    sage: E.integral_points()
    [(-11 : 29 : 1), (-2 : -28 : 1), (4 : 11 : 1),
    (13 : 29 : 1), (507525709 : -11433961056931 : 1)]

::

    sage: E = EllipticCurve('6104b1')
    sage: E.integral_points(algorithm = "old")
    [(-41 : 266 : 1), (-27 : 294 : 1), (-13 : 266 : 1), (19 : 74 : 1),
    (22 : 49 : 1), (27 : 6 : 1), (29 : 14 : 1), (43 : 154 : 1),
    (53 : 266 : 1), (71 : 490 : 1), (414 : 8379 : 1), (1639 : 66346 : 1),
    (1667 : 68054 : 1)]
    sage: E.integral_points()
    [(-41 : 266 : 1), (-27 : 294 : 1), (-13 : -266 : 1), (19 : 74 : 1),
    (22 : -49 : 1), (27 : 6 : 1), (29 : -14 : 1), (43 : 154 : 1),
    (53 : 266 : 1), (71 : -490 : 1), (414 : 8379 : 1), (1639 : 66346 : 1),
    (1667 : 68054 : 1), (23036693 : 110568192854 : 1)]

::

    sage: E = EllipticCurve('8470g1')
    sage: E.integral_points(algorithm = "old")
    [(-1 : 0 : 1), (0 : 6 : 1), (3 : 12 : 1), (10 : 33 : 1), (15 : 56 : 1),
    (24 : 110 : 1), (43 : 264 : 1), (98 : 924 : 1), (395 : 7656 : 1),
    (1681690 : 2179974753 : 1)]
    sage: E.integral_points()
    [(-1 : 0 : 1), (0 : 6 : 1), (3 : -16 : 1), (10 : 33 : 1), (15 : -72 : 1),
    (24 : 110 : 1), (43 : -308 : 1), (98 : -1023 : 1), (395 : 7656 : 1),
    (1681690 : -2181656444 : 1)]

::

    sage: E = EllipticCurve('8470g2')
    sage: E.integral_points(algorithm = "old")
    [(-14 : 17 : 1), (-12 : 33 : 1), (-7 : 38 : 1), (-1 : 22 : 1),
    (0 : 17 : 1), (13 : 8 : 1), (14 : 17 : 1), (21 : 66 : 1),
    (33 : 158 : 1), (44 : 257 : 1), (63 : 458 : 1), (98 : 913 : 1),
    (175 : 2222 : 1), (1533 : 59258 : 1), (1561 : 60896 : 1)]
    sage: E.integral_points()
    [(-14 : 17 : 1), (-12 : 33 : 1), (-7 : -32 : 1), (-1 : 22 : 1),
    (0 : -18 : 1), (13 : -22 : 1), (14 : 17 : 1), (21 : -88 : 1),
    (33 : -192 : 1), (44 : -302 : 1), (63 : 458 : 1), (98 : -1012 : 1),
    (175 : 2222 : 1), (1533 : -60792 : 1), (1561 : 60896 : 1),
    (18888968 : -82103620687 : 1)]

::

    sage: from sage.schemes.elliptic_curves.ell_int_points import integral_points
    sage: len(EllipticCurve('15a1').integral_points(algorithm = "old")) == len(EllipticCurve('15a1').integral_points())
    True

AUTHORS:

- this file is based on work of Thotsaphon "Nook" Thongjunthug that
  was initially ported to Sage by Radoslav Kirov and Jackie Anderson.

- Radoslav Kirov, Jackie Anderson (2010): first version

- Justin Walker, Gagan Sekhon (2011): clean up

- Alyson Deines, Jen Balakrishnan (2011): Documentation, clean up,
  fixing bugs so this works for all examples in the paper [SMART]_

REFERENCE:

.. [SMART] Smart, N. P. and Stephens, N. M., Integral Points on Elliptic
   Curves over Number Fields. Math. Proc. Camb. Phil. Soc. (1997), 122, 9
"""


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

from copy import copy
from sage.functions.all import sqrt
from sage.matrix.all import zero_matrix
from sage.misc.all import verbose, prod
from sage.rings.all import polygen, ZZ, RealField, ComplexField
from sage.rings.all import QQ, integer_ceil, integer_floor
from sage.rings.real_mpfr import is_RealField


def abs_log_height(X, gcd_one=True, precision=None):
    r"""
    Computes the height of a point in a projective space over field `K`.

    It assumes the coordinates have gcd equal to 1. If not use
    gcd_one flag.

    INPUT:

    - ``X`` -- Point in projective space over a number field `K`

    - ``gcd_one`` -- (default: True) if false this throws in the
      places at the numerators in addition to the places at the denominator.

    - ``precision`` -- (default: None) bits of precision of the real and complex fields

    OUTPUT:

    The absolute logarithmic height of the point.

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: E = EllipticCurve('5077a')
        sage: [ellpts.abs_log_height(list(u)) for u in E.gens()]
        [1.09861228866811, 3.21887582486820, 0.000000000000000]

    ::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2-x-1)
        sage: F = EllipticCurve(K,[0,-a,1,-a-1,2*a+1])
        sage: P1, P2 = F.gens()
        sage: P2
        (-3/4*a + 1/4 : -5/4*a - 5/8 : 1)
        sage: ellpts.abs_log_height(list(P2))
        2.56625746464637
    """
    assert isinstance(X, list)
    K = X[0].parent()
    if precision is None:
        RR = RealField()
    else:
        RR = RealField(precision)
    places = set([])
    if K == QQ:
        embs = K.embeddings(RR)
        Xideal = X
    else:
        embs = K.places(precision)
        # skipping zero as it currently over K breaks Sage
        Xideal = [K.ideal(x) for x in X if not x.is_zero()]
    for x in Xideal:
        places = places.union(x.denominator().prime_factors())
        if not gcd_one:
            places = places.union(x.numerator().prime_factors())
    if K == QQ:
        non_arch = sum([max([RR(P)**(-x.valuation(P)) for x in X]).log()
                        for P in places])
    else:
        non_arch = sum([P.residue_class_degree() *
                        P.absolute_ramification_index() *
                        max([x.abs_non_arch(P, precision) for x in X]).log()
                        for P in places])
    arch = 0
    r, s = K.signature()
    for i, f in enumerate(embs):
        if i < r:
            arch += max([f(x).abs() for x in X]).log()
        else:
            arch += max([f(x).abs()**2 for x in X]).log()
    return (arch+non_arch)/K.degree()


def _compute_c6(E, emb):
    r"""
    Computes the `c_6` invariant from [SMART]_

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2+5)
        sage: E = EllipticCurve(K,[-112,400])
        sage: embs = K.embeddings(CC)
        sage: ellpts._compute_c6(E,embs[0])
        24.0994429807464
    """
    F = emb.codomain()
    if is_RealField(F):
        x = polygen(ComplexField(F.prec()))
    else:
        x = polygen(F)
    #f = x**3-27*emb(E.c4())*x-54*emb(E.c6())
    f = x**3-emb(E.c4()/48)*x-emb(E.c6()/864)
    R = f.roots(multiplicities=False)
    m = max([r.abs() for r in R])
    return 2*m


def _compute_c8(L):
    r"""
    Computes the `c_8` invariant from [SMART]_

    Original code transformed the lattice generators.
    Here we assume they are of standard form.

    INPUT:

    - `L` -- a basis for the period lattice of an elliptic curve over
      a number field via a given embedding

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2+5)
        sage: E = EllipticCurve(K,[-112,400])
        sage: embs = K.embeddings(CC)
        sage: Periods = E.period_lattice(embs[0]).normalised_basis();
        sage: ellpts._compute_c8(Periods)
        1.27486778253862
    """
    CC = ComplexField()
    w1, w2 = L
    return max(CC(w1).abs(), CC(w2).abs(), CC(w1+w2).abs())


def _c3(E, L):
    r"""
    Computes the `c_3` invariant of [SMART]_

    This is the least eigenvalue of the height pairing matrix of the
    Mordell-Weil basis.

    INPUT:

    - `E` -- an elliptic curve over a number field

    - `L` -- a list of points of `E` which generate the Mordell-Weil group

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2+5)
        sage: E = EllipticCurve(K,[-112,400])
        sage: ellpts._c3(E,[E((4,-4)),E((0,-20)),E((-4,-28)),
        ....: E((-89/5,-637*a/25))])
        0.366347735724...
    """
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return min(map(abs, E.height_pairing_matrix(L).eigenvalues()))


def _h_E(E):
    r"""
    Computes the height of the ellpitc curve `E`.

    INPUT:

    -`E` -- an elliptic curve

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2+5)
        sage: E = EllipticCurve(K,[-112,400])
        sage: ellpts._h_E(E)
        17.4513334798896
    """
    K = E.base_field()
    j = E.j_invariant()
    c4, c6 = E.c_invariants()
    g2, g3 = c4/12, c6/216
    return max(abs_log_height([K(1), g2, g3]), abs_log_height([K(1), j]))


def _h_m(E, P, ElogEmbedP, D7):
    r"""
    Computes the list of modified height of a point `P` on `E(K)` given
    a list of embeddings

    INPUT:

    - `E` -- an elliptic curve

    - `P` -- a point on `E`

    - `ElogEmbedP` -- elliptic logarithm of `P` in some embedding

    - `D7` -- constant `d_7`

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2+5)
        sage: embs = K.embeddings(CC)
        sage: E = EllipticCurve(K,[-112,400])
        sage: L = [E((4,-4)),E((0,-20)),E((-4,-28)),E((-89/5,-637*a/25))]
        sage: Elog = [ P.elliptic_logarithm(embs[0], precision = 300)
        ....: for P in L]
        sage: ellpts._h_m(E, L[0], Elog[0], 3.28)
        17.4513334798896
    """
    K = E.base_field()
    # M1 =  max([P.height(), _h_E(E), D7/K.degree()*abs((ElogEmbedP).log())**2])
    M2 = max([P.height(), _h_E(E), D7/K.degree()*abs(ElogEmbedP)**2])
    return M2


def _Extra_h_m(E, Periods, D7):
    r"""
    Computes two extra `h_m`'s given periods.

    INPUT:

    - `E` -- elliptic curve

    - `Periods` -- list of two periods of the fundamental parallelogram

    - `D7` -- constant `d_7`

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2+5)
        sage: embs = K.embeddings(CC)
        sage: E = EllipticCurve(K,[-112,400])
        sage: Periods = E.period_lattice(embs[0]).normalised_basis()
        sage: ellpts._Extra_h_m(E, Periods, 90.3)
        [48.6392854292010, 24.7424615832146]
    """
    D = E.base_field().degree()
    h = _h_E(E)
    return map(lambda x: max([0, h, D7/D*abs(x)**2]), Periods)


def _d8(E, L, Elog, Periods, D7):
    r"""
    The `c_8` constant (which we have renamed to `d_8` to avoid
    confusion with the other set of `c` constants) arising from
    David's lower bound for linear forms in elliptic logarithms.

    Notation as in Appendix A of Smart's book, "The Algorithmic
    Resolution of Diophantine Equations".

    INPUT:

    - `E` -- elliptic curve

    - `L` -- a list of points in `E(K)` (e.g. a basis for the
      Mordell-Weil group)

    - `Elog` -- a sequence of precomputed elliptic logs of each point in `L`

    - `D7` -- constant `d_7`

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2+5)
        sage: embs = K.embeddings(CC)
        sage: E = EllipticCurve(K,[-112,400])
        sage: Periods = E.period_lattice(embs[0]).normalised_basis()
        sage: L = [E((4,-4)),E((0,-20)),E((-4,-28)),E((-89/5,-637*a/25))]
        sage: Elog = [ P.elliptic_logarithm(embs[0], precision = 300)
        ....: for P in L]
        sage: ellpts._d8(E, L, Elog, Periods, 90.3)
        47.4376426807629
    """
    RR = RealField()
    C = [RR(1).exp() * _h_E(E)]
    D = E.base_field().degree()
    for i in range(len(L)):
        C.append(_h_m(E, L[i], Elog[i], D7) / D)
    C += [t / D for t in _Extra_h_m(E, Periods, D7)]
    return max(C)


def _d9(E, L, Elog, Periods, D7):
    r"""
    The `c_9` constant (which we have renamed to `d_9` to avoid
    confusion with the other set of `c` constants) arrising from
    David's lower bound for linear forms in elliptic logarithms.

    Notation as in Appendix A of Smart's book, "The Algorithmic
    Resolution of Diophantine Equations".

    INPUT:

    - `E` -- elliptic curve

    - `L` -- a list of points in `E(K)` (e.g. a basis for the
      Mordell-Weil group)

    - `Elog` -- a sequence of precomputed elliptic logs of each point in `L`

    - `D7` -- constant `d_7`

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2+5)
        sage: embs = K.embeddings(CC)
        sage: E = EllipticCurve(K,[-112,400])
        sage: Periods = E.period_lattice(embs[0]).normalised_basis()
        sage: L = [E((4,-4)),E((0,-20)),E((-4,-28)),E((-89/5,-637*a/25))]
        sage: Elog = [ P.elliptic_logarithm(embs[0], precision = 300)
        ....: for P in L]
        sage: ellpts._d9(E, L, Elog, Periods, 90.3)
        2.71828182845904
    """
    RR = RealField()
    D = E.base_field().degree()
    C = []
    for i in range(len(L)):
        tmp = RR(1.0).exp() * sqrt(D * _h_m(E, L[i], Elog[i], D7) / D7)
        C.append(tmp / abs(Elog[i]))
    Ehm = _Extra_h_m(E, Periods, D7)
    C += [RR(1).exp() * sqrt(D*Ehm[i]/D7) / abs(Periods[i]) for i in [0, 1]]
    verbose("d9: C = %s" % (C, ))
    return min(C)


def _d10(E, L, Elog, Periods, D7):
    r"""
    The `c_10` constant (which we have renamed to `d_10` to avoid
    confusion with the other set of `c` constants) arising from
    David's lower bound for linear forms in elliptic logarithms.

    Notation as in Appendix A of Smart's book, "The Algorithmic
    Resolution of Diophantine Equations".

    INPUT:

    - `E` -- elliptic curve

    - `L` -- a list of points in `E(K)` (e.g. a basis for the
      Mordell-Weil group)

    - `Elog` -- a sequence of precomputed elliptic logs of each point in `L`

    - `D7` -- constant `d_7`

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2+5)
        sage: embs = K.embeddings(CC)
        sage: E = EllipticCurve(K,[-112,400])
        sage: Periods = E.period_lattice(embs[0]).normalised_basis()
        sage: L = [E((4,-4)),E((0,-20)),E((-4,-28)),E((-89/5,-637*a/25))]
        sage: Elog = [ P.elliptic_logarithm(embs[0], precision = 300)
        ....: for P in L]
        sage: ellpts._d10(E, L, Elog, Periods, 90.3)
        1.19716657228519e226
    """
    RR = RealField()
    D = E.base_field().degree()
    n = len(L)+2
    scalar = 2 * 10**(8 + 7*n) * (2/RR(1).exp())**(2*n**2)
    scalar *= (n+1)**(4*n**2 + 10*n) * D**(2*n + 2)
    scalar *= (RR(_d9(E, L, Elog, Periods, D7)).log())**(-2*n-1)
    for i in range(len(L)):
        scalar *= _h_m(E, L[i], Elog[i], D7)
    scalar *= prod(_Extra_h_m(E, Periods, D7))
    return scalar


def _RHS(D, r, C9, C10, D9, D10, h, Q, expTors):
    r"""
    Computes the right hand side of the principal inequality.

    INPUT:

    - `D` -- degree of field extension

    - `r` -- Mordell-Weil rank

    - `C9` -- `c_9` invariant from [SMART]_

    - `C19` -- `c_{10}` invariant from [SMART]_

    - `D9` -- `D_9` invariant

    - `D10` -- `d_{10}` invariant

    - `h` -- `h_E(E)` where `E` is an elliptic curve

    - `Q` -- initial bound for the coefficients of `P_i`'s, where
      `P_i`'s are in the Mordell-Weil basis

    - `expTors` -- the exponent of the torsion subgroup of `E`

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2+5)
        sage: E = EllipticCurve(K,[-112,400])
        sage: L = [E((4,-4)),E((0,-20)),E((-4,-28)),E((-89/5,-637*a/25))]
        sage: C9 =  3425.58431073226
        sage: C10 = 0.183173867862383
        sage: Q0 =  5.12358476153997
        sage: D8 =  47.4376426807629
        sage: D9 =  6.19404412950310
        sage: D10 =  2.39513899842427e221
        sage: h =  17.4513334798896
        sage: expTors = 1
        sage: ellpts._RHS(K.degree(), len(L), C9, C10, D9, D10, h,
        ....: 1000000000000000000000, expTors)
        3.02137037238302e233
    """
    RR = RealField()
    bound = (RR(RR(Q*r*expTors).log()).log() + h + RR(D*D9).log())**(r+3)
    bound *= D10*(RR(Q*r*expTors).log() + RR(D*D9).log())
    bound += RR(C9*expTors).log()
    bound /= C10
    return bound


def _InitialQ(D, r, Q0, C9, C10, D8, D9, D10, h, expTors):
    r"""
    Computes the initial bound on `Q = max_{1 \leq i \leq r}(|q_i|)`.

    INPUT:

    - `D` -- degree of field extension

    - `r` -- Mordell-Weil rank

    - `Q0` -- the `K_0` constant in [SMART]_

    - `C9` -- `c_9` invariant from [SMART]_

    - `C10` -- `c_{10}` invariant from [SMART]_

    - `D9` -- `D_9` invariant

    - `D10` -- `d_{10}` invariant

    - `h` -- `h_E(E)` where `E` is an elliptic curve

    - `expTors` -- the exponent of the torsion subgroup of `E`

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2+5)
        sage: E = EllipticCurve(K,[-112,400])
        sage: L = [E((4,-4)),E((0,-20)),E((-4,-28)),E((-89/5,-637*a/25))]
        sage: C9 =  3425.58431073226
        sage: C10 = 0.183173867862383
        sage: Q0 =  5.12358476153997
        sage: D8 =  47.4376426807629
        sage: D9 =  6.19404412950310
        sage: D10 =  2.39513899842427e221
        sage: h =  17.4513334798896
        sage: expTors = 1
        sage: ellpts._InitialQ(K.degree(), len(L), Q0, C9, C10, D8, D9, D10,
        ....: h, expTors)
        118
    """
    RR = RealField()
    minQ = max(Q0, RR(D8).exp())
    Q = minQ + 1
    x = integer_ceil(Q.log(10))  # x = log_10(Q)
    verbose("x: %s, %s" % (type(Q), type(x)))
    exp_OK = 0    # the exponent that satisfies P.I.
    exp_fail = 0  # the exponent that fails P.I.
    verbose("Looping on %s < RHS(%s, %s, %s, %s %s, %s, %s, %s, %s"
            % (10**(2*x), D, r, C9, C10, D9, D10, h, 10**x, expTors))
    while 10**(2*x) < _RHS(D, r, C9, C10, D9, D10, h, 10**x, expTors):
        exp_OK = x  # Principal Inequality satisfied
        x *= 2      # double x, and retry
    exp_fail = x    # x that fails the Principal Inequality

    # So now x = log_10(Q) must lie between exp_OK and exp_fail
    # Refine x further using "binary search"
    while True:
        x = integer_floor((exp_OK + exp_fail)/2)
        if 10**(2*x) >= _RHS(D, r, C9, C10, D9, D10, h, 10**x, expTors):
            exp_fail = x
        else:
            exp_OK = x
        if (exp_fail - exp_OK) <= 1:
            break
    return exp_fail  # over-estimate


def _ReducedQ(E, f, LGen, Elog, C9, C10, Periods, expTors, initQ):
    r"""
    Internal reduction function which reduces the bound on `Q`.

    INPUT:

    - `E` -- elliptic curve over `K`

    - `f` -- embedding of `K` into `\mathbb{C}`

    - `LGen` -- a list of points in `E(K)` (e.g. a basis for the
      Mordell-Weil group)

    - `Elog` -- a sequence of precomputed elliptic logs of each point in `L`

    - `C9` -- `c_9` invariant from [SMART]_

    - `C10` -- `c_{10}` invariant from [SMART]_

    - `Periods` -- list of two periods of the fundamental parallelogram

    - `expTors` -- the exponent of the torsion subgroup of `E`

    - `initQ` -- inital guess for `Q` which is then reduced by LLL

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2+5)
        sage: E = EllipticCurve(K,[-112,400])
        sage: L = [E((4,-4)),E((0,-20)),E((-4,-28)),E((-89/5,-637*a/25))]
        sage: embs = K.embeddings(CC)
        sage: Elog = [ P.elliptic_logarithm(embs[0], precision = 300)
        ....: for P in L]
        sage: Periods = Periods = E.period_lattice(embs[0]).normalised_basis()
        sage: C9 =  3425.58431073226
        sage: C10 = 0.183173867862383
        sage: expTors = 1
        sage: initQ = 10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        sage: ellpts._ReducedQ(E, embs[0], L, Elog, C9, C10, Periods,
        ....: expTors, initQ)
        13
    """
    RR = RealField()
    w1, w2 = Periods
    r = len(LGen)
    newQ = initQ
    # Repeat LLL reduction until no further reduction is possible
    while True:
        Q = newQ
        S = r*(Q**2)*(expTors**2)
        T = 3*r*Q*expTors/sqrt(2.0)
        # Create the basis matrix
        C = 1
        while True:
            C *= Q**integer_ceil((r+2)/2)
            precbits = int(RR(C).log(2)+10)*4
            L = copy(zero_matrix(ZZ, r+2))
            # Elliptic logarithm should have precision "suitable to" C
            # e.g. If C = 10**100, then Re(Elog[i]) should be computed
            # correctly to at least 100 decimal places
            if precbits > Elog[0].prec():
                verbose("Need higher precision, recompute elliptic"
                        "logarithm ...")
                # Re-compute elliptic logarithm to the right precision
                verbose("precision in bits %i" % precbits)
                Elog = [P.elliptic_logarithm(f, precision=precbits)
                        for P in LGen]
                verbose("Elliptic logarithm recomputed")
                w1, w2 = E.period_lattice(f).normalised_basis(prec=precbits)
            # Assign all non-zero entries
            for i in range(r):
                L[i, i] = 1
                L[r, i] = (C*Elog[i].real_part()).trunc()
                L[r+1, i] = (C*Elog[i].imag_part()).trunc()
            L[r, r] = (C*w1.real_part()).trunc()
            L[r, r+1] = (C*w2.real_part()).trunc()
            L[r+1, r] = (C*w1.imag_part()).trunc()
            L[r+1, r+1] = (C*w2.imag_part()).trunc()
            # LLL reduction and constants
            L = L.transpose()
            L = L.LLL()
            b1 = L[0]  # 1st row of reduced basis
            #Norm(b1) = square of Euclidean norm
            normb1 = sum([i**2 for i in b1])
            lhs = RR(2)**(-r-1)*normb1 - S
            if (lhs >= 0) and (sqrt(lhs) > T):
                break
        newQ = (RR(C*C9*expTors).log() - RR(sqrt(lhs) - T).log()) / C10
        newQ = integer_floor(sqrt(newQ))
        verbose("After reduction, Q <= %f" % newQ)
        if ((Q - newQ) <= 1) or (newQ <= 1):
            break
    return newQ


def _cyc_iter(id, gens, mult, both_signs=False):
    r"""
    Cyclic group iterator

    INPUT:

    - `id` -- identity element

    - `gens` -- generators of the group

    - `mult` -- orders of the generators

    - ``both_signs`` --(Default: False) whether to return all elements
      or one per each {element, inverse} pair

    OUTPUT:

    Returns all elements of the cyclic group, remembering all intermediate steps

    EXAMPLE::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: E = EllipticCurve('57b')
        sage: ord = [P.order() for P in E.torsion_points()]
        sage: P_list = [i for i in ellpts._cyc_iter(E(0),
        ....: E.torsion_points(),ord)]
        sage: P_list
        [((0 : 1 : 0), [0, 0, 0, 0]), ((7/4 : -11/8 : 1), [0, 0, 0, 1]), ((1 :
        -1 : 1), [0, 0, 1, 0]), ((-3 : 1 : 1), [0, 0, 1, 1]), ((-3 : 1 : 1), [1,
        0, 0, 0]), ((1 : -1 : 1), [1, 0, 0, 1]), ((7/4 : -11/8 : 1), [1, 0, 1,
        0]), ((0 : 1 : 0), [1, 0, 1, 1])]
    """
    if len(gens) == 0:
        yield id, []
    else:
        P = gens[0]
        cur = id
        if both_signs:
            ran = xrange(mult[0])
        else:
            ran = xrange(mult[0]/2 + 1)
        for k in ran:
            for rest, coefs in _cyc_iter(id, gens[1:], mult[1:],
                                         both_signs or k != 0):
                yield cur + rest, [k] + coefs
            cur += P


#generates elements of form a_1G_1 + ... + a_nG_n
#where |a_i| <= bound and the first non-zero coefficient is positive
def _L_points_iter(id, gens, bound, all_zero=True):
    r"""
    Generates elements of form `a_1G_1 + ... + a_nG_n`
    where `|a_i| \leq` bound and the first non-zero coefficient is positive

    INPUT:

    - ``id`` -- identity element

    - ``gens`` -- a list of generators, `G_1,...,G_n`

    - ``bound`` -- a bound on the absolute value of the `a_i`'s, the
      coefficients of the elements in terms of the generators.

    - ``all_zero`` --(default: True)

    OUTPUT:

    An iterator that generates elements of the form `a_1G_1 + ... + a_nG_n`

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2+2)
        sage: E = EllipticCurve(K,[-16,16])
        sage: P_list = [[P,Pcoeff] for P,Pcoeff in
        ....: ellpts._L_points_iter(E(0),[E((0,4)),E((2,-2*a))],24)]
        sage: P_list[0:5]
        [[(0 : 1 : 0), [0, 0]], [(2 : -2*a : 1), [0, 1]], [(-9/2 : -5/4*a : 1),
        [0, 2]], [(418/169 : 4514/2197*a : 1), [0, 3]], [(-30241/200 :
        5257039/4000*a : 1), [0, 4]]]
    """
    if len(gens) == 0:
        yield id, []
    else:
        P = gens[0]
        cur = id
        for k in xrange(bound+1):
            for rest, coefs in _L_points_iter(id, gens[1:], bound,
                                              all_zero=(all_zero and k == 0)):
                yield cur + rest, [k] + coefs
                if k != 0 and not all_zero:
                    yield -cur + rest, [-k] + coefs
            cur += P


def _integral_points_with_Q(E, L, Q):
    r"""
    Given an elliptic curve `E` over a number field, its Mordell-Weil
    basis L,and a positive integer Q, return the sequence of all
    integral points modulo `[-1]` of the form `P = q_1*L[1] + ... +
    q_r*L[r] + T` with some torsion point `T` and `|q_i| \leq Q`.

    INPUT:

    - ``E`` -- an elliptic curve

    - ``L`` -- a basis for the Mordell-Weil group of `E`

    - ``Q`` -- positive integer, maximum for the absolute bound on all
      coefficients in the linear combination of points in `L`

    OUTPUT:

    Returns the sequence of all integral points modulo `[-1]`

    EXAMPLE::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2-x-1)
        sage: E = EllipticCurve(K,[2-a,-a-1,1-a,a,1])
        sage: L = E.gens()
        sage: ellpts._integral_points_with_Q(E,L,25)
        [(a : 1 : 1)]
    """
    if L == []:
        try:
            L = E.gens()
        except ValueError:
            raise ValueError("Generators of the Mordell-Weil group "
                             "must be supplied.")
    # Tors = E.torsion_subgroup()
    # expTors = Tors.exponent()
    Tpoints = E.torsion_points()
    id = E([0, 1, 0])
    verbose("Generators = %s" % L)
    PtsList = []

    verbose("Exact arithmetic")
    for P, Pcoeff in _L_points_iter(id, L, Q):
            for T in Tpoints:
                R = P+T
                if R[0].is_integral() and R[1].is_integral() and R != id:
                    if R not in PtsList and -R not in PtsList:
                        PtsList.append(R)
    verbose("*"*45)
    return PtsList


def _calculate_Q(E, L):
    r"""
    This calculates `Q`, one of the bounds in [SMART]_.

    `Q` is a positive integer, maximum for the absolute bound on all
    coefficients in the linear combination of points in `L`

    INPUT:

    - ``E`` -- Elliptic Curve

    - ``L`` -- a basis for the Mordell-Weil group of `E`

    OUTPUT:

    A positive integer, maximum for the absolute bound on all
    coefficients in the linear combination of points in `L`

    EXAMPLE::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: K.<a> = NumberField(x^2+5)
        sage: E = EllipticCurve(K,[-112,400])
        sage: L = [E((4,-4)),E((0,-20)),E((-4,-28)),E((-89/5,-637*a/25))]
        sage: ellpts._calculate_Q(E,L)
        13
    """
    from math import pi as Rpi
    CC = ComplexField()

    if len(L) == 0:
        return 0
    K = E.base_ring()
    expTors = E.torsion_subgroup().exponent()
    r, s = K.signature()
    RR = RealField()
    pi = RR(Rpi)
    b2 = E.b_invariants()[0]
    # Global constants (independent to the embedding of E)
    gl_con = {}
    #C2 = E.silverman_height_bounds[1] matches with [SMART]
    #C2 = -E.silverman_height_bound[0] matches with Nook's implementation in Magma
    #both return correct sets of integral points . . .
    gl_con['c2'] = C2 = -E.silverman_height_bounds()[0]
    gl_con['c3'] = C3 = _c3(E, L)
    gl_con['h_E'] = h = _h_E(E)
    verbose("Global constants")
    for name, val in gl_con.iteritems():
        verbose('%s = %s' % (name, val))
    verbose("-"*45)
    Q = []
    # Find the most reduced bound on each embedding of E
    # Sage bug, QQ.places() is not implemented
    # and K.embeddings gives each complex places twice
    if K is QQ:
        infplaces = K.embeddings(CC)
    else:
        infplaces = K.places()
    for i, f in enumerate(infplaces):
        # Bug in P.complex_logarithm(QQ.embeddings(CC)[0])
        if K is QQ:
            emb = None
        else:
            emb = f
        if i < r:
            nv = 1
            verbose("Real embedding #%i" % i)
        else:
            nv = 2
            verbose("Complex embedding #%i" % (i-r))
        # Create complex embedding of E
        # Embedding of all points in Mordell-Weil basis
        # Find complex elliptic logarithm of each embedded point
        precbits = integer_floor(100*RR(10).log(2))
        Elog = [P.elliptic_logarithm(emb, precision=precbits) for P in L]
        Periods = E.period_lattice(emb).normalised_basis()
        # Local constants (depending on embedding)
        # C9, C10 are used for the upper bound on the linear form in logarithm
        loc_con = {}
        loc_con['c4'] = C4 = C3 * K.degree() / RR(nv*(r+s))
        loc_con['c5'] = C5 = C2 * K.degree() / RR(nv*(r+s))
        loc_con['c6'] = C6 = _compute_c6(E, f)
        loc_con['delta'] = delta = 1 + (nv-1)*pi
        loc_con['c8'] = C8 = _compute_c8(Periods)
        loc_con['c7'] = C7 = RR(8)*(delta**2) + (C8**2)*abs(f(b2))/RR(12)
        #C9 = sqrt(C7)*RR(C5/2).exp() is a tigher bound, but [SMART] only uses
        #C9 = C7*RR(C5/2).exp() . . .
        #additionally, we could never get our C9 values to match with [SMART]
        loc_con['c9'] = C9 = C7*RR(C5/2).exp()
        #C10 was often correct to 3 or so digits
        loc_con['c10'] = C10 = C4/2
        #Q0 (K0 in [SMART]) is often correct to 1 or 2 digits
        #These descrepancies don't seem to matter . . .
        loc_con['Q0'] = Q0 = sqrt((RR(C6+abs(f(b2))/12).log() + C5) / C4)
        # Constants for David's lower bound on the linear form in logarithm
        # Magma and Sage output periods in different order, need to swap w1 and w2
        w2, w1 = Periods  # N.B. periods are already in "standard" form
        loc_con['d7'] = D7 = 3*pi / ((abs(w2)**2) * (w2/w1).imag_part())
        loc_con['d8'] = D8 = _d8(E, L, Elog, Periods, D7)
        loc_con['d9'] = D9 = _d9(E, L, Elog, Periods, D7)
        loc_con['d10'] = D10 = _d10(E, L, Elog, Periods, D7)
        for name, val in loc_con.iteritems():
            verbose("{0} = {1}".format(name, val))
        # Find the reduced bound for the coefficients in the linear logarithmic form
        loginitQ = _InitialQ(K.degree(), len(L), Q0, C9, C10, D8, D9, D10,
                             h, expTors)
        verbose("Initial Q <= 10^%f" % loginitQ)
        initQ = 10**loginitQ
        Q.append(_ReducedQ(E, emb, L, Elog, C9, C10, Periods, expTors, initQ))
        verbose("-"*45)
    Q = max(Q)
    verbose("Maximum absolute bound on coefficients = %i" % Q)
    return Q


def integral_points(self, L=None, both_signs=False, algorithm="new"):
    r"""
    Given an elliptic curve over a number field and its Mordell-Weil
    basis, return the sequence of all integral points modulo [-1].

    both_signs = True gives the set of all integral points.

    INPUT:

    - ``L`` -- (default: `L` = self.gens()) a basis for the
      Mordell-Weil group of self

    - ``both_signs`` -- (default: False) This gives the entire set of
      points, `P` and `-P` for each integral point.  If set to True,
      this just gives the `P` and not `-P`.

    - ``algorithm`` -- (default: "new") When "new", this runs the
      generic code in ell_int_points.py, when "old" this runs the old
      (and sometimes incorrect) code in ell_rational_field.py.

    OUTPUT:

    Returns the sequence of all integral points modulo `[-1]`

    .. note::

        The complexity increases exponentially in the rank of curve
        E. The computation time (but not the output!) depends on
        the Mordell-Weil basis. If mw_base is given but is not a
        basis for the Mordell-Weil group (modulo torsion), integral
        points which are not in the subgroup generated by the given
        points will almost certainly not be listed.

        All examples from [SMART]_ have been tested and return the
        correct integral points.


    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_int_points as ellpts
        sage: E = EllipticCurve('5077a')
        sage: ellpts.integral_points(E,E.gens(),both_signs = True) # long time
        [(-3 : -1 : 1), (-3 : 0 : 1), (-2 : -4 : 1), (-2 : 3 : 1),
        (-1 : -4 : 1), (-1 : 3 : 1), (0 : -3 : 1),
        (0 : 2 : 1), (1 : -1 : 1), (1 : 0 : 1), (2 : -1 : 1),
        (2 : 0 : 1), (3 : -4 : 1), (3 : 3 : 1), (4 : -7 : 1),
        (4 : 6 : 1), (8 : -22 : 1), (8 : 21 : 1), (11 : -36 : 1),
        (11 : 35 : 1), (14 : -52 : 1), (14 : 51 : 1), (21 : -96 : 1),
        (21 : 95 : 1), (37 : -225 : 1), (37 : 224 : 1), (52 : -375 : 1),
        (52 : 374 : 1), (93 : -897 : 1), (93 : 896 : 1), (342 : -6325 : 1),
        (342 : 6324 : 1), (406 : -8181 : 1),
        (406 : 8180 : 1), (816 : -23310 : 1), (816 : 23309 : 1)]

    Example 1 from [SMART]_::

        sage: K.<a> = NumberField(x^2+5)
        sage: E = EllipticCurve(K,[-112,400])
        sage: E.integral_points(L = [E((4,-4)),E((0,-20)),E((-4,-28)),E((-89/5,-637*a/25))])  # long time
        [(-12 : 4 : 1), (-8 : 28 : 1), (-7 : -29 : 1), (-4 : -28 : 1),
        (0 : -20 : 1), (1 : 17 : 1), (4 : -4 : 1), (8 : 4 : 1), (9 : -11 : 1),
        (12 : -28 : 1), (16 : 52 : 1), (25 : -115 : 1), (32 : 172 : 1),
        (44 : 284 : 1), (56 : -412 : 1), (84 : -764 : 1), (148 : 1796 : 1),
        (208 : 2996 : 1), (372 : -7172 : 1), (1368 : -50596 : 1),
        (1624 : -65444 : 1), (3264 : 186476 : 1)]

    Example 2 from [SMART]_::

        sage: K.<a> = NumberField(x^2+2)
        sage: E = EllipticCurve(K,[-16,16])
        sage: E.integral_points([E((0,4)),E((2,-2*a))])
        [(-4*a - 10 : -18*a + 28 : 1), (4*a - 10 : 18*a + 28 : 1),
        (-4*a - 4 : 20 : 1), (-4 : -4 : 1), (4*a - 4 : 20 : 1),
        (-4*a : -8*a - 12 : 1), (0 : 4 : 1), (4*a : 8*a - 12 : 1),
        (1 : -1 : 1), (2 : -2*a : 1), (4 : 4 : 1), (8 : -20 : 1),
        (24 : 116 : 1), (-40*a + 60 : -480*a + 316 : 1),
        (40*a + 60 : 480*a + 316 : 1)]

    Example 3 from [SMART]_, also an example where we can not compute
    generators using the simon two descent::

        sage: K.<a> = NumberField(x^3+3*x-1)
        sage: E = EllipticCurve(K,[3,0])
        sage: E.integral_points(L = [E((1,2)),E((a,1))])
        [(0 : 0 : 1), (3*a^2 + 9 : -9*a^2 - 3*a - 27 : 1), (3 : -6 : 1),
        (6*a^2 + 2*a + 17 : -24*a^2 - 6*a - 74 : 1), (12 : 42 : 1),
        (204*a^2 + 65*a + 632 : -5286*a^2 - 1704*a - 16405 : 1),
        (a : 1 : 1), (1 : 2 : 1)]
        sage: E.integral_points()
        Traceback (most recent call last):
        ...
        ValueError: Cannot compute a provable set of generators,
        please supply generators of the Mordell-Weil group.
    """
    if L is None or L == [] or L == 'auto':
        try:
            L = self.gens()
            r = self.rank()
        except ValueError:
            raise ValueError("Cannot compute a provable set of generators, "
                             "please "
                             "supply generators of the Mordell-Weil group.")
        assert len(L) == r

    if algorithm == "old":
        try:
            return self._integral_points_old(L, both_signs)
        except:
            pass

    Q = _calculate_Q(self, L)
    int_points = _integral_points_with_Q(self, L, Q)
    if both_signs:
        int_points += [-P for P in int_points]

    #need this so doc tests will always have points in the same order
    int_points.sort()
    return int_points
