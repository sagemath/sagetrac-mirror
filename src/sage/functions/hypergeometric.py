r"""
Hypergeometric Functions

This module implements manipulation of infinite hypergeometric series
represented in standard parametric form (as `\,_pF_q` functions).

AUTHORS:

- Fredrik Johansson (2010): initial version

- Eviatar Bach (2013): major changes

EXAMPLES:

Examples from :trac:`9908`::

    sage: maxima('integrate(bessel_j(2, x), x)').sage()
    1/24*x^3*hypergeometric((3/2,), (5/2, 3), -1/4*x^2)
    sage: sum(((2*I)^x/(x^3 + 1)*(1/4)^x), x, 0, oo)
    hypergeometric((1, 1, -1/2*I*sqrt(3) - 1/2, 1/2*I*sqrt(3) - 1/2),...
    (2, -1/2*I*sqrt(3) + 1/2, 1/2*I*sqrt(3) + 1/2), 1/2*I)
    sage: sum((-1)^x/((2*x + 1)*factorial(2*x + 1)), x, 0, oo)
    hypergeometric((1/2,), (3/2, 3/2), -1/4)

Simplification (note that ``simplify_full`` does not yet call
``simplify_hypergeometric``)::

    sage: hypergeometric([-2], [], x).simplify_hypergeometric()
    x^2 - 2*x + 1
    sage: hypergeometric([], [], x).simplify_hypergeometric()
    e^x
    sage: a = hypergeometric((hypergeometric((), (), x),), (),
    ....:                    hypergeometric((), (), x))
    sage: a.simplify_hypergeometric()
    1/((-e^x + 1)^e^x)
    sage: a.simplify_hypergeometric(algorithm='sage')
    (-e^x + 1)^(-e^x)

Equality testing::

    sage: bool(hypergeometric([], [], x).derivative(x) ==
    ....:      hypergeometric([], [], x))  # diff(e^x, x) == e^x
    True
    sage: bool(hypergeometric([], [], x) == hypergeometric([], [1], x))
    False

Computing terms and series::

    sage: z = var('z')
    sage: hypergeometric([], [], z).series(z, 0)
    Order(1)
    sage: hypergeometric([], [], z).series(z, 1)
    1 + Order(z)
    sage: hypergeometric([], [], z).series(z, 2)
    1 + 1*z + Order(z^2)
    sage: hypergeometric([], [], z).series(z, 3)
    1 + 1*z + 1/2*z^2 + Order(z^3)

    sage: hypergeometric([-2], [], z).series(z, 3)
    1 + (-2)*z + 1*z^2
    sage: hypergeometric([-2], [], z).series(z, 6)
    1 + (-2)*z + 1*z^2
    sage: hypergeometric([-2], [], z).series(z, 6).is_terminating_series()
    True
    sage: hypergeometric([-2], [], z).series(z, 2)
    1 + (-2)*z + Order(z^2)
    sage: hypergeometric([-2], [], z).series(z, 2).is_terminating_series()
    False

    sage: hypergeometric([1], [], z).series(z, 6)
    1 + 1*z + 1*z^2 + 1*z^3 + 1*z^4 + 1*z^5 + Order(z^6)
    sage: hypergeometric([], [1/2], -z^2/4).series(z, 11)
    1 + (-1/2)*z^2 + 1/24*z^4 + (-1/720)*z^6 + 1/40320*z^8 +...
    (-1/3628800)*z^10 + Order(z^11)

    sage: hypergeometric([1], [5], x).series(x, 5)
    1 + 1/5*x + 1/30*x^2 + 1/210*x^3 + 1/1680*x^4 + Order(x^5)

    sage: sum(hypergeometric([1, 2], [3], 1/3).terms(6)).n()
    1.29788359788360
    sage: hypergeometric([1, 2], [3], 1/3).n()
    1.29837194594696
    sage: hypergeometric([], [], x).series(x, 20)(x=1).n() == e.n()
    True

Plotting::

    sage: plot(hypergeometric([1, 1], [3, 3, 3], x), x, -30, 30)
    Graphics object consisting of 1 graphics primitive
    sage: complex_plot(hypergeometric([x], [], 2), (-1, 1), (-1, 1))
    Graphics object consisting of 1 graphics primitive

Numeric evaluation::

    sage: hypergeometric([1], [], 1/10).n()  # geometric series
    1.11111111111111
    sage: hypergeometric([], [], 1).n()  # e
    2.71828182845905
    sage: hypergeometric([], [], 3., hold=True)
    hypergeometric((), (), 3.00000000000000)
    sage: hypergeometric([1, 2, 3], [4, 5, 6], 1/2).n()
    1.02573619590134
    sage: hypergeometric([1, 2, 3], [4, 5, 6], 1/2).n(digits=30)
    1.02573619590133865036584139535
    sage: hypergeometric([5 - 3*I], [3/2, 2 + I, sqrt(2)], 4 + I).n()
    5.52605111678805 - 7.86331357527544*I
    sage: hypergeometric((10, 10), (50,), 2.)
    -1705.75733163554 - 356.749986056024*I

Conversions::

    sage: maxima(hypergeometric([1, 1, 1], [3, 3, 3], x))
    hypergeometric([1,1,1],[3,3,3],_SAGE_VAR_x)
    sage: hypergeometric((5, 4), (4, 4), 3)._sympy_()
    hyper((5, 4), (4, 4), 3)
    sage: hypergeometric((5, 4), (4, 4), 3)._mathematica_init_()
    'HypergeometricPFQ[{5,4},{4,4},3]'

Arbitrary level of nesting for conversions::

    sage: maxima(nest(lambda y: hypergeometric([y], [], x), 3, 1))
    1/(1-_SAGE_VAR_x)^(1/(1-_SAGE_VAR_x)^(1/(1-_SAGE_VAR_x)))
    sage: maxima(nest(lambda y: hypergeometric([y], [3], x), 3, 1))._sage_()
    hypergeometric((hypergeometric((hypergeometric((1,), (3,), x),), (3,),...
    x),), (3,), x)
    sage: nest(lambda y: hypergeometric([y], [], x), 3, 1)._mathematica_init_()
    'HypergeometricPFQ[{HypergeometricPFQ[{HypergeometricPFQ[{1},{},x]},...

The confluent hypergeometric functions can arise as solutions to second-order
differential equations (example from `here <http://ask.sagemath.org/question/
1168/how-can-one-use-maxima-kummer-confluent-functions>`_)::

    sage: m = var('m')
    sage: y = function('y')(x)
    sage: desolve(diff(y, x, 2) + 2*x*diff(y, x) - 4*m*y, y,
    ....:         contrib_ode=true, ivar=x)
    [y(x) == _K1*hypergeometric_M(-m, 1/2, -x^2) +...
     _K2*hypergeometric_U(-m, 1/2, -x^2)]

Series expansions of confluent hypergeometric functions::

    sage: hypergeometric_M(2, 2, x).series(x, 3)
    1 + 1*x + 1/2*x^2 + Order(x^3)
    sage: hypergeometric_U(2, 2, x).series(x == 3, 100).subs(x=1).n()
    0.403652637676806
    sage: hypergeometric_U(2, 2, 1).n()
    0.403652637676806
"""

#*****************************************************************************
#       Copyright (C) 2010 Fredrik Johansson <fredrik.johansson@gmail.com>
#       Copyright (C) 2013 Eviatar Bach <eviatarbach@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from six.moves import range

from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import Infinity
from sage.arith.all import binomial, rising_factorial, factorial
from sage.functions.other import sqrt, gamma, real_part
from sage.functions.log import exp, log
from sage.functions.trig import sin
from sage.functions.hyperbolic import cosh, sinh
from sage.functions.other import erf
from sage.symbolic.constants import pi
from sage.symbolic.all import I
from sage.symbolic.function import BuiltinFunction
from sage.symbolic.ring import SR
from sage.structure.element import get_coercion_model
from sage.misc.latex import latex
from sage.misc.misc_c import prod
from sage.libs.mpmath import utils as mpmath_utils
from sage.symbolic.expression import Expression
from sage.calculus.functional import derivative
from functools import reduce


def rational_param_as_tuple(x):
    r"""
    Utility function for converting rational `\,_pF_q` parameters to
    tuples (which mpmath handles more efficiently).

    EXAMPLES::

        sage: from sage.functions.hypergeometric import rational_param_as_tuple
        sage: rational_param_as_tuple(1/2)
        (1, 2)
        sage: rational_param_as_tuple(3)
        3
        sage: rational_param_as_tuple(pi)
        pi
    """
    try:
        x = x.pyobject()
    except AttributeError:
        pass
    try:
        if x.parent() is QQ:
            p = int(x.numer())
            q = int(x.denom())
            return p, q
    except AttributeError:
        pass
    return x


class Hypergeometric(BuiltinFunction):
    r"""
    Represents a (formal) generalized infinite hypergeometric series. It is
    defined as

    .. MATH::

        \,_pF_q(a_1, \ldots, a_p; b_1, \ldots, b_q; z)
        = \sum_{n=0}^{\infty} \frac{(a_1)_n \cdots (a_p)_n}{(b_1)_n
        \cdots(b_q)_n} \, \frac{z^n}{n!},

    where `(x)_n` is the rising factorial.
    """
    def __init__(self):
        """
        Initialize class.
        
        EXAMPLES::

            sage: maxima(hypergeometric)
            hypergeometric
        """
        BuiltinFunction.__init__(self, 'hypergeometric', nargs=3,
                                 conversions={'mathematica':
                                              'HypergeometricPFQ',
                                              'maxima': 'hypergeometric',
                                              'sympy': 'hyper'})

    def __call__(self, a, b, z, **kwargs):
        """
        Return symbolic hypergeometric function expression.
         
        INPUT:
    
        - ``a`` -- a list or tuple of parameters
        - ``b`` -- a list or tuple of parameters
        - ``z`` -- a number or symbolic expression
    
        EXAMPLES::

            sage: hypergeometric([], [], 1)
            hypergeometric((), (), 1)
            sage: hypergeometric([], [1], 1)
            hypergeometric((), (1,), 1)
            sage: hypergeometric([2, 3], [1], 1)
            hypergeometric((2, 3), (1,), 1)
            sage: hypergeometric([], [], x)
            hypergeometric((), (), x)
            sage: hypergeometric([x], [], x^2)
            hypergeometric((x,), (), x^2)
    
        The only simplification that is done automatically is returning 1
        if ``z`` is 0. For other simplifications use the
        ``simplify_hypergeometric`` method.
        """
        return BuiltinFunction.__call__(self,
                                        SR._force_pyobject(a),
                                        SR._force_pyobject(b),
                                        z, **kwargs)

    def _print_latex_(self, a, b, z):
        r"""
        TESTS::

            sage: latex(hypergeometric([1, 1], [2], -1))
            \,_2F_1\left(\begin{matrix} 1,1 \\ 2 \end{matrix} ; -1 \right)

        """
        aa = ",".join(latex(c) for c in a)
        bb = ",".join(latex(c) for c in b)
        z = latex(z)
        return (r"\,_{}F_{}\left(\begin{{matrix}} {} \\ {} \end{{matrix}} ; "
                r"{} \right)").format(len(a), len(b), aa, bb, z)

    def _eval_(self, a, b, z, **kwargs):
        """
        EXAMPLES::

            sage: hypergeometric([], [], 0)
            1
        """
        if not isinstance(a,tuple) or not isinstance(b,tuple):
            raise TypeError("The first two parameters must be of type list")
        if not isinstance(z, Expression) and z == 0:  # Expression is excluded
            return Integer(1)                         # to avoid call to Maxima

    def _evalf_try_(self, a, b, z):
        """
        Call :meth:`_evalf_` if one of the arguments is numerical and none
        of the arguments are symbolic.

        OUTPUT:

        - ``None`` if we didn't succeed to call :meth:`_evalf_` or if
          the input wasn't suitable for it.

        - otherwise, a numerical value for the function.

        EXAMPLES::

            sage: hypergeometric._evalf_try_((1.0,), (2.0,), 3.0)
            6.36184564106256
            sage: hypergeometric._evalf_try_((1.0, 1), (), 3.0)
            -0.0377593153441588 + 0.750349833788561*I
            sage: hypergeometric._evalf_try_((1, 1), (), 3)    # exact input
            sage: hypergeometric._evalf_try_((x,), (), 1.0)    # symbolic
            sage: hypergeometric._evalf_try_(1.0, 2.0, 3.0)    # not tuples
        """
        # We need to override this for hypergeometric functions since
        # the first 2 arguments are tuples and the generic _evalf_try_
        # cannot handle that.
        if not isinstance(a,tuple) or not isinstance(b,tuple):
            return None
        args = list(a) + list(b) + [z]
        if any(self._is_numerical(x) for x in args):
            if not any(isinstance(x, Expression) for x in args):
                p = get_coercion_model().common_parent(*args)
                return self._evalf_(a, b, z, parent=p)

    def _evalf_(self, a, b, z, parent, algorithm=None):
        """
        TESTS::

            sage: hypergeometric([1, 1], [2], -1).n()
            0.693147180559945
            sage: hypergeometric([], [], RealField(100)(1))
            2.7182818284590452353602874714

        """
        if not isinstance(a,tuple) or not isinstance(b,tuple):
            raise TypeError("The first two parameters must be of type list")
        from mpmath import hyper
        aa = [rational_param_as_tuple(c) for c in a]
        bb = [rational_param_as_tuple(c) for c in b]
        return mpmath_utils.call(hyper, aa, bb, z, parent=parent)

    def _tderivative_(self, a, b, z, *args, **kwargs):
        """
        EXAMPLES::

            sage: hypergeometric([1/3, 2/3], [5], x^2).diff(x)
            4/45*x*hypergeometric((4/3, 5/3), (6,), x^2)
            sage: hypergeometric([1, 2], [x], 2).diff(x)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative of hypergeometric function with...
            respect to parameters. Try calling .simplify_hypergeometric()...
            first.
            sage: hypergeometric([1/3, 2/3], [5], 2).diff(x)
            0
        """
        diff_param = kwargs['diff_param']
        if diff_param in hypergeometric(a, b, 1).variables():  # ignore z
            raise NotImplementedError("derivative of hypergeometric function "
                                      "with respect to parameters. Try calling"
                                      " .simplify_hypergeometric() first.")
        t = (reduce(lambda x, y: x * y, a, 1) *
             reduce(lambda x, y: x / y, b, Integer(1)))
        return (t * derivative(z, diff_param) *
                hypergeometric([c + 1 for c in a], [c + 1 for c in b], z))

    class EvaluationMethods(object):
        def _fast_float_(cls, self, *args):
            """
            Do not support the old ``fast_float``.

            OUTPUT:

            This method raises ``NotImplementedError``; use the newer
            ``fast_callable`` implementation.

            EXAMPLES::

                sage: f = hypergeometric([], [], x)
                sage: f._fast_float_()
                Traceback (most recent call last):
                ...
                NotImplementedError
            """
            raise NotImplementedError

        def _fast_callable_(cls, self, a, b, z, etb):
            """
            Override the ``fast_callable`` method.

            OUTPUT:

            A :class:`~sage.ext.fast_callable.ExpressionCall` representing the
            hypergeometric function in the expression tree.

            EXAMPLES::

                sage: h = hypergeometric([], [], x)
                sage: from sage.ext.fast_callable import ExpressionTreeBuilder
                sage: etb = ExpressionTreeBuilder(vars=['x'])
                sage: h._fast_callable_(etb)
                {hypergeometric((), (), x)}(v_0)

                sage: y = var('y')
                sage: fast_callable(hypergeometric([y], [], x),
                ....:               vars=[x, y])(3, 4)
                hypergeometric((4,), (), 3)
            """
            return etb.call(self, *map(etb.var, etb._vars))

        def sorted_parameters(cls, self, a, b, z):
            """
            Return with parameters sorted in a canonical order.

            EXAMPLES::

                sage: hypergeometric([2, 1, 3], [5, 4],
                ....:                1/2).sorted_parameters()
                hypergeometric((1, 2, 3), (4, 5), 1/2)
            """
            return hypergeometric(sorted(a), sorted(b), z)

        def eliminate_parameters(cls, self, a, b, z):
            """
            Eliminate repeated parameters by pairwise cancellation of identical
            terms in ``a`` and ``b``.

            EXAMPLES::

                sage: hypergeometric([1, 1, 2, 5], [5, 1, 4],
                ....:                1/2).eliminate_parameters()
                hypergeometric((1, 2), (4,), 1/2)
                sage: hypergeometric([x], [x], x).eliminate_parameters()
                hypergeometric((), (), x)
                sage: hypergeometric((5, 4), (4, 4), 3).eliminate_parameters()
                hypergeometric((5,), (4,), 3)
            """
            aa = list(a)  # tuples are immutable
            bb = list(b)
            p = pp = len(aa)
            q = qq = len(bb)
            i = 0
            while i < qq and aa:
                bbb = bb[i]
                if bbb in aa:
                    aa.remove(bbb)
                    bb.remove(bbb)
                    pp -= 1
                    qq -= 1
                else:
                    i += 1
            if (pp, qq) != (p, q):
                return hypergeometric(aa, bb, z)
            return self

        def is_termwise_finite(cls, self, a, b, z):
            """
            Determine whether all terms of ``self`` are finite. Any infinite
            terms or ambiguous terms beyond the first zero, if one exists,
            are ignored.

            Ambiguous cases (where a term is the product of both zero
            and an infinity) are not considered finite.

            EXAMPLES::

                sage: hypergeometric([2], [3, 4], 5).is_termwise_finite()
                True
                sage: hypergeometric([2], [-3, 4], 5).is_termwise_finite()
                False
                sage: hypergeometric([-2], [-3, 4], 5).is_termwise_finite()
                True
                sage: hypergeometric([-3], [-3, 4],
                ....:                5).is_termwise_finite()  # ambiguous
                False

                sage: hypergeometric([0], [-1], 5).is_termwise_finite()
                True
                sage: hypergeometric([0], [0],
                ....:                5).is_termwise_finite()  # ambiguous
                False
                sage: hypergeometric([1], [2], Infinity).is_termwise_finite()
                False
                sage: (hypergeometric([0], [0], Infinity)
                ....:  .is_termwise_finite())  # ambiguous
                False
                sage: (hypergeometric([0], [], Infinity)
                ....:  .is_termwise_finite())  # ambiguous
                False
            """
            if z == 0:
                return 0 not in b
            if abs(z) == Infinity:
                return False
            if abs(z) == Infinity:
                return False
            for bb in b:
                if bb in ZZ and bb <= 0:
                    if any((aa in ZZ) and (bb < aa <= 0) for aa in a):
                        continue
                    return False
            return True

        def is_terminating(cls, self, a, b, z):
            r"""
            Determine whether the series represented by self terminates
            after a finite number of terms, i.e. whether any of the
            numerator parameters are nonnegative integers (with no
            preceding nonnegative denominator parameters), or `z = 0`.

            If terminating, the series represents a polynomial of `z`.

            EXAMPLES::

                sage: hypergeometric([1, 2], [3, 4], x).is_terminating()
                False
                sage: hypergeometric([1, -2], [3, 4], x).is_terminating()
                True
                sage: hypergeometric([1, -2], [], x).is_terminating()
                True
            """
            if z == 0:
                return True
            for aa in a:
                if (aa in ZZ) and (aa <= 0):
                    return self.is_termwise_finite()
            return False

        def is_absolutely_convergent(cls, self, a, b, z):
            r"""
            Determine whether ``self`` converges absolutely as an infinite
            series. ``False`` is returned if not all terms are finite.

            EXAMPLES:

            Degree giving infinite radius of convergence::

                sage: hypergeometric([2, 3], [4, 5],
                ....:                6).is_absolutely_convergent()
                True
                sage: hypergeometric([2, 3], [-4, 5],
                ....:                6).is_absolutely_convergent()  # undefined
                False
                sage: (hypergeometric([2, 3], [-4, 5], Infinity)
                ....:  .is_absolutely_convergent())  # undefined
                False

            Ordinary geometric series (unit radius of convergence)::

                sage: hypergeometric([1], [], 1/2).is_absolutely_convergent()
                True
                sage: hypergeometric([1], [], 2).is_absolutely_convergent()
                False
                sage: hypergeometric([1], [], 1).is_absolutely_convergent()
                False
                sage: hypergeometric([1], [], -1).is_absolutely_convergent()
                False
                sage: hypergeometric([1], [], -1).n()  # Sum still exists
                0.500000000000000

            Degree `p = q+1` (unit radius of convergence)::

                sage: hypergeometric([2, 3], [4], 6).is_absolutely_convergent()
                False
                sage: hypergeometric([2, 3], [4], 1).is_absolutely_convergent()
                False
                sage: hypergeometric([2, 3], [5], 1).is_absolutely_convergent()
                False
                sage: hypergeometric([2, 3], [6], 1).is_absolutely_convergent()
                True
                sage: hypergeometric([-2, 3], [4],
                ....:                5).is_absolutely_convergent()
                True
                sage: hypergeometric([2, -3], [4],
                ....:                5).is_absolutely_convergent()
                True
                sage: hypergeometric([2, -3], [-4],
                ....:                5).is_absolutely_convergent()
                True
                sage: hypergeometric([2, -3], [-1],
                ....:                5).is_absolutely_convergent()
                False

            Degree giving zero radius of convergence::

                sage: hypergeometric([1, 2, 3], [4],
                ....:                2).is_absolutely_convergent()
                False
                sage: hypergeometric([1, 2, 3], [4],
                ....:                1/2).is_absolutely_convergent()
                False
                sage: (hypergeometric([1, 2, -3], [4], 1/2)
                ....:  .is_absolutely_convergent())  # polynomial
                True
            """
            p, q = len(a), len(b)
            if not self.is_termwise_finite():
                return False
            if p <= q:
                return True
            if self.is_terminating():
                return True
            if p == q + 1:
                if abs(z) < 1:
                    return True
                if abs(z) == 1:
                    if real_part(sum(b) - sum(a)) > 0:
                        return True
            return False

        def terms(cls, self, a, b, z, n=None):
            """
            Generate the terms of ``self`` (optionally only ``n`` terms).

            EXAMPLES::

                sage: list(hypergeometric([-2, 1], [3, 4], x).terms())
                [1, -1/6*x, 1/120*x^2]
                sage: list(hypergeometric([-2, 1], [3, 4], x).terms(2))
                [1, -1/6*x]
                sage: list(hypergeometric([-2, 1], [3, 4], x).terms(0))
                []
            """
            if n is None:
                n = Infinity
            t = Integer(1)
            k = 1
            while k <= n:
                yield t
                for aa in a:
                    t *= (aa + k - 1)
                for bb in b:
                    t /= (bb + k - 1)
                t *= z
                if t == 0:
                    break
                t /= k
                k += 1

        def deflated(cls, self, a, b, z):
            r"""
            Rewrite as a linear combination of functions of strictly lower
            degree by eliminating all parameters ``a[i]`` and ``b[j]`` such
            that ``a[i]`` = ``b[i]`` + ``m`` for nonnegative integer ``m``.

            EXAMPLES::

                sage: x = hypergeometric([6, 1], [3, 4, 5], 10)
                sage: y = x.deflated()
                sage: y
                1/252*hypergeometric((4,), (7, 8), 10)
                 + 1/12*hypergeometric((3,), (6, 7), 10)
                 + 1/2*hypergeometric((2,), (5, 6), 10)
                 + hypergeometric((1,), (4, 5), 10)
                sage: x.n(); y.n()
                2.87893612686782
                2.87893612686782

                sage: x = hypergeometric([6, 7], [3, 4, 5], 10)
                sage: y = x.deflated()
                sage: y
                25/27216*hypergeometric((), (11,), 10)
                 + 25/648*hypergeometric((), (10,), 10)
                 + 265/504*hypergeometric((), (9,), 10)
                 + 181/63*hypergeometric((), (8,), 10)
                 + 19/3*hypergeometric((), (7,), 10)
                 + 5*hypergeometric((), (6,), 10)
                 + hypergeometric((), (5,), 10)
                sage: x.n(); y.n()
                63.0734110716969
                63.0734110716969
            """
            return sum(map(prod, self._deflated()))

        def _deflated(cls, self, a, b, z):
            """
            Private helper to return list of deflated terms.
            
            EXAMPLES::

                sage: x = hypergeometric([5], [4], 3)
                sage: y = x.deflated()
                sage: y
                7/4*hypergeometric((), (), 3)
                sage: x.n(); y.n()
                35.1496896155784
                35.1496896155784
            """
            new = self.eliminate_parameters()
            aa = new.operands()[0].operands()
            bb = new.operands()[1].operands()
            for i, aaa in enumerate(aa):
                for j, bbb in enumerate(bb):
                    m = aaa - bbb
                    if m in ZZ and m > 0:
                        aaaa = aa[:i] + aa[i + 1:]
                        bbbb = bb[:j] + bb[j + 1:]
                        terms = []
                        for k in range(m + 1):
                            # TODO: could rewrite prefactors as recurrence
                            term = binomial(m, k)
                            for c in aaaa:
                                term *= rising_factorial(c, k)
                            for c in bbbb:
                                term /= rising_factorial(c, k)
                            term *= z ** k
                            term /= rising_factorial(aaa - m, k)
                            F = hypergeometric([c + k for c in aaaa],
                                               [c + k for c in bbbb], z)
                            unique = []
                            counts = []
                            for c, f in F._deflated():
                                if f in unique:
                                    counts[unique.index(f)] += c
                                else:
                                    unique.append(f)
                                    counts.append(c)
                            Fterms = zip(counts, unique)
                            terms += [(term * termG, G) for (termG, G) in
                                      Fterms]
                        return terms
            return ((1, new),)

hypergeometric = Hypergeometric()


def closed_form(hyp):
    """
    Try to evaluate ``hyp`` in closed form using elementary
    (and other simple) functions.

    It may be necessary to call :meth:`Hypergeometric.deflated` first to
    find some closed forms.

    EXAMPLES::

        sage: from sage.functions.hypergeometric import closed_form
        sage: var('a b c z')
        (a, b, c, z)
        sage: closed_form(hypergeometric([1], [], 1 + z))
        -1/z
        sage: closed_form(hypergeometric([], [], 1 + z))
        e^(z + 1)
        sage: closed_form(hypergeometric([], [1/2], 4))
        cosh(4)
        sage: closed_form(hypergeometric([], [3/2], 4))
        1/4*sinh(4)
        sage: closed_form(hypergeometric([], [5/2], 4))
        3/16*cosh(4) - 3/64*sinh(4)
        sage: closed_form(hypergeometric([], [-3/2], 4))
        19/3*cosh(4) - 4*sinh(4)
        sage: closed_form(hypergeometric([-3, 1], [var('a')], z))
        -3*z/a + 6*z^2/((a + 1)*a) - 6*z^3/((a + 2)*(a + 1)*a) + 1
        sage: closed_form(hypergeometric([-3, 1/3], [-4], z))
        7/162*z^3 + 1/9*z^2 + 1/4*z + 1
        sage: closed_form(hypergeometric([], [], z))
        e^z
        sage: closed_form(hypergeometric([a], [], z))
        (-z + 1)^(-a)
        sage: closed_form(hypergeometric([1, 1, 2], [1, 1], z))
        (z - 1)^(-2)
        sage: closed_form(hypergeometric([2, 3], [1], x))
        -1/(x - 1)^3 + 3*x/(x - 1)^4
        sage: closed_form(hypergeometric([1/2], [3/2], -5))
        1/10*sqrt(5)*sqrt(pi)*erf(sqrt(5))
        sage: closed_form(hypergeometric([2], [5], 3))
        4
        sage: closed_form(hypergeometric([2], [5], 5))
        48/625*e^5 + 612/625
        sage: closed_form(hypergeometric([1/2, 7/2], [3/2], z))
        1/5*z^2/(-z + 1)^(5/2) + 2/3*z/(-z + 1)^(3/2) + 1/sqrt(-z + 1)
        sage: closed_form(hypergeometric([1/2, 1], [2], z))
        -2*(sqrt(-z + 1) - 1)/z
        sage: closed_form(hypergeometric([1, 1], [2], z))
        -log(-z + 1)/z
        sage: closed_form(hypergeometric([1, 1], [3], z))
        -2*((z - 1)*log(-z + 1)/z - 1)/z
        sage: closed_form(hypergeometric([1, 1, 1], [2, 2], x))
        hypergeometric((1, 1, 1), (2, 2), x)
    """
    if hyp.is_terminating():
        return sum(hyp.terms())

    new = hyp.eliminate_parameters()

    def _closed_form(hyp):
        a, b, z = hyp.operands()
        a, b = a.operands(), b.operands()
        p, q = len(a), len(b)

        if z == 0:
            return Integer(1)
        if p == q == 0:
            return exp(z)
        if p == 1 and q == 0:
            return (1 - z) ** (-a[0])

        if p == 0 and q == 1:
            # TODO: make this require only linear time
            def _0f1(b, z):
                F12 = cosh(2 * sqrt(z))
                F32 = sinh(2 * sqrt(z)) / (2 * sqrt(z))
                if 2 * b == 1:
                    return F12
                if 2 * b == 3:
                    return F32
                if 2 * b > 3:
                    return ((b - 2) * (b - 1) / z * (_0f1(b - 2, z) -
                            _0f1(b - 1, z)))
                if 2 * b < 1:
                    return (_0f1(b + 1, z) + z / (b * (b + 1)) *
                            _0f1(b + 2, z))
                raise ValueError
            # Can evaluate 0F1 in terms of elementary functions when
            # the parameter is a half-integer
            if 2 * b[0] in ZZ and b[0] not in ZZ:
                return _0f1(b[0], z)

        # Confluent hypergeometric function
        if p == 1 and q == 1:
            aa, bb = a[0], b[0]
            if aa * 2 == 1 and bb * 2 == 3:
                t = sqrt(-z)
                return sqrt(pi) / 2 * erf(t) / t
            if a == 1 and b == 2:
                return (exp(z) - 1) / z
            n, m = aa, bb
            if n in ZZ and m in ZZ and m > 0 and n > 0:
                rf = rising_factorial
                if m <= n:
                    return (exp(z) * sum(rf(m - n, k) * (-z) ** k /
                            factorial(k) / rf(m, k) for k in
                            range(n - m + 1)))
                else:
                    T = sum(rf(n - m + 1, k) * z ** k /
                            (factorial(k) * rf(2 - m, k)) for k in
                            range(m - n))
                    U = sum(rf(1 - n, k) * (-z) ** k /
                            (factorial(k) * rf(2 - m, k)) for k in
                            range(n))
                    return (factorial(m - 2) * rf(1 - m, n) *
                            z ** (1 - m) / factorial(n - 1) *
                            (T - exp(z) * U))

        if p == 2 and q == 1:
            R12 = QQ('1/2')
            R32 = QQ('3/2')

            def _2f1(a, b, c, z):
                """
                Evaluation of 2F1(a, b; c; z), assuming a, b, c positive
                integers or half-integers
                """
                if b == c:
                    return (1 - z) ** (-a)
                if a == c:
                    return (1 - z) ** (-b)
                if a == 0 or b == 0:
                    return Integer(1)
                if a > b:
                    a, b = b, a
                if b >= 2:
                    F1 = _2f1(a, b - 1, c, z)
                    F2 = _2f1(a, b - 2, c, z)
                    q = (b - 1) * (z - 1)
                    return (((c - 2 * b + 2 + (b - a - 1) * z) * F1 +
                            (b - c - 1) * F2) / q)
                if c > 2:
                    # how to handle this case?
                    if a - c + 1 == 0 or b - c + 1 == 0:
                        raise NotImplementedError
                    F1 = _2f1(a, b, c - 1, z)
                    F2 = _2f1(a, b, c - 2, z)
                    r1 = (c - 1) * (2 - c - (a + b - 2 * c + 3) * z)
                    r2 = (c - 1) * (c - 2) * (1 - z)
                    q = (a - c + 1) * (b - c + 1) * z
                    return (r1 * F1 + r2 * F2) / q

                if (a, b, c) == (R12, 1, 2):
                    return (2 - 2 * sqrt(1 - z)) / z
                if (a, b, c) == (1, 1, 2):
                    return -log(1 - z) / z
                if (a, b, c) == (1, R32, R12):
                    return (1 + z) / (1 - z) ** 2
                if (a, b, c) == (1, R32, 2):
                    return 2 * (1 / sqrt(1 - z) - 1) / z
                if (a, b, c) == (R32, 2, R12):
                    return (1 + 3 * z) / (1 - z) ** 3
                if (a, b, c) == (R32, 2, 1):
                    return (2 + z) / (2 * (sqrt(1 - z) * (1 - z) ** 2))
                if (a, b, c) == (2, 2, 1):
                    return (1 + z) / (1 - z) ** 3
                raise NotImplementedError
            aa, bb = a
            cc, = b
            if z == 1:
                return (gamma(cc) * gamma(cc - aa - bb) / gamma(cc - aa) /
                        gamma(cc - bb))
            if ((aa * 2) in ZZ and (bb * 2) in ZZ and (cc * 2) in ZZ and
                aa > 0 and bb > 0 and cc > 0):
                try:
                    return _2f1(aa, bb, cc, z)
                except NotImplementedError:
                    pass
        return hyp
    return sum([coeff * _closed_form(pfq) for coeff, pfq in new._deflated()])

class Hypergeometric_M(BuiltinFunction):
    r"""
    The confluent hypergeometric function of the first kind,
    `y = M(a,b,z)`, is defined to be the solution to Kummer's differential
    equation

    .. math::

             zy'' + (b-z)y' - ay = 0.

    This is not the same as Kummer's `U`-hypergeometric function, though it
    satisfies the same DE that `M` does.

    .. warning::

       In the literature, both are called "Kummer confluent
       hypergeometric" functions.

    EXAMPLES::

        sage: hypergeometric_M(1, 1, 1)
        hypergeometric_M(1, 1, 1)
        sage: hypergeometric_M(1, 1, 1.)
        2.71828182845905
        sage: hypergeometric_M(1, 1, 1).n(70)
        2.7182818284590452354
        sage: hypergeometric_M(1, 1, 1).simplify_hypergeometric()
        e
        sage: hypergeometric_M(1, 1/2, x).simplify_hypergeometric()
        (-I*sqrt(pi)*x*erf(I*sqrt(-x))*e^x + sqrt(-x))/sqrt(-x)
        sage: hypergeometric_M(1, 3/2, 1).simplify_hypergeometric()
        1/2*sqrt(pi)*erf(1)*e
    """
    def __init__(self):
        """
        TESTS::
        
            sage: maxima(hypergeometric_M(1,1,x))
            kummer_m(1,1,_SAGE_VAR_x)
            sage: latex(hypergeometric_M(1,1,x))
            M\left(1, 1, x\right)
        """
        BuiltinFunction.__init__(self, 'hypergeometric_M', nargs=3,
                                 conversions={'mathematica':
                                              'Hypergeometric1F1',
                                              'maxima': 'kummer_m'},
                                 latex_name='M')

    def _eval_(self, a, b, z, **kwargs):
        """
        TESTS::
        
            sage: (a,b)=var('a,b')
            sage: hypergeometric_M(a,b,0)
            1
        """
        if not isinstance(z, Expression) and z == 0:
            return Integer(1)
        return

    def _evalf_(self, a, b, z, parent, algorithm=None):
        """
        TESTS::
        
            sage: hypergeometric_M(1,1,1).n()
            2.71828182845905
        """
        from mpmath import hyp1f1
        return mpmath_utils.call(hyp1f1, a, b, z, parent=parent)

    def _derivative_(self, a, b, z, diff_param):
        """
        TESTS::
        
            sage: diff(hypergeometric_M(1,1,x),x,3)
            hypergeometric_M(4, 4, x)
            sage: diff(hypergeometric_M(x,1,1),x,3)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative of hypergeometric function with respect to parameters
        """
        if diff_param == 2:
            return (a / b) * hypergeometric_M(a + 1, b + 1, z)
        raise NotImplementedError('derivative of hypergeometric function '
                                  'with respect to parameters')

    class EvaluationMethods(object):
        def generalized(cls, self, a, b, z):
            """
            Return as a generalized hypergeometric function

            EXAMPLES::

                sage: a, b, z = var('a b z')
                sage: hypergeometric_M(a, b, z).generalized()
                hypergeometric((a,), (b,), z)

            """
            return hypergeometric([a], [b], z)

hypergeometric_M = Hypergeometric_M()

class Hypergeometric_U(BuiltinFunction):
    r"""
    The confluent hypergeometric function of the second kind,
    `y = U(a,b,z)`, is defined to be the solution to Kummer's differential
    equation

    .. math::

             zy'' + (b-z)y' - ay = 0.

    This satisfies `U(a,b,z) \sim z^{-a}`, as
    `z\rightarrow \infty`, and is sometimes denoted
    `z^{-a}{}_2F_0(a,1+a-b;;-1/z)`. This is not the same as Kummer's
    `M`-hypergeometric function, denoted sometimes as
    `_1F_1(\alpha,\beta,z)`, though it satisfies the same DE that
    `U` does.

    .. warning::

       In the literature, both are called "Kummer confluent
       hypergeometric" functions.

    EXAMPLES::

        sage: hypergeometric_U(1, 1, 1)
        hypergeometric_U(1, 1, 1)
        sage: hypergeometric_U(1, 1, 1.)
        0.596347362323194
        sage: hypergeometric_U(1, 1, 1).n(70)
        0.59634736232319407434
        sage: hypergeometric_U(10^4, 1/3, 1).n()
        6.60377008885811e-35745
        sage: hypergeometric_U(2 + I, 2, 1).n()
        0.183481989942099 - 0.458685959185190*I
        sage: hypergeometric_U(1, 3, x).simplify_hypergeometric()
        (x + 1)/x^2
        sage: hypergeometric_U(1, 2, 2).simplify_hypergeometric()
        1/2

    """
    def __init__(self):
        """
        TESTS::
        
            sage: maxima(hypergeometric_U(1,1,x))
            kummer_u(1,1,_SAGE_VAR_x)
            sage: latex(hypergeometric_U(1,1,x))
            U\left(1, 1, x\right)
        """
        BuiltinFunction.__init__(self, 'hypergeometric_U', nargs=3,
                                 conversions={'mathematica':
                                              'HypergeometricU',
                                              'maxima': 'kummer_u'},
                                 latex_name='U')

    def _eval_(self, a, b, z, **kwargs):
        return

    def _evalf_(self, a, b, z, parent, algorithm=None):
        """
        TESTS::
        
            sage: hypergeometric_U(1,1,1).n()
            0.596347362323194
        """
        from mpmath import hyperu
        return mpmath_utils.call(hyperu, a, b, z, parent=parent)

    def _derivative_(self, a, b, z, diff_param):
        """
        TESTS::
        
            sage: diff(hypergeometric_U(1,1,x),x,3)
            -6*hypergeometric_U(4, 4, x)
            sage: diff(hypergeometric_U(x,1,1),x,3)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative of hypergeometric function with respect to parameters
        """
        if diff_param == 2:
            return -a * hypergeometric_U(a + 1, b + 1, z)
        raise NotImplementedError('derivative of hypergeometric function '
                                  'with respect to parameters')

    class EvaluationMethods(object):
        def generalized(cls, self, a, b, z):
            """
            Return in terms of the generalized hypergeometric function

            EXAMPLES::

                sage: a, b, z = var('a b z')
                sage: hypergeometric_U(a, b, z).generalized()
                z^(-a)*hypergeometric((a, a - b + 1), (), -1/z)
                sage: hypergeometric_U(1, 3, 1/2).generalized()
                2*hypergeometric((1, -1), (), -2)
                sage: hypergeometric_U(3, I, 2).generalized()
                1/8*hypergeometric((3, -I + 4), (), -1/2)

            """
            return z ** (-a) * hypergeometric([a, a - b + 1], [], -z ** (-1))

hypergeometric_U = Hypergeometric_U()
