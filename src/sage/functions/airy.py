r"""
Airy Functions

This module implements Airy functions and their generalized derivatives. It
supports symbolic functionality through Maxima and numeric evaluation through
mpmath and scipy.

Airy functions are solutions to the differential equation `f''(x) - x f(x) = 0`.

Four global function symbols are immediately available, please see

- :func:`airy_ai`: for the Airy Ai function

- :func:`airy_ai_prime()<FunctionAiryAiPrime>`: for the first differential
  of the Airy Ai function

- :func:`airy_bi`: for the Airy Bi function

- :func:`airy_bi_prime()<FunctionAiryBiPrime>`: for the first differential
   of the Airy Bi function

AUTHORS:

- Oscar Gerardo Lazo Arjona (2010): initial version

- Douglas McNeil (2012): rewrite

EXAMPLES:

Verify that the Airy functions are solutions to the differential equation::

    sage: diff(airy_ai(x), x, 2) - x * airy_ai(x)
    0
    sage: diff(airy_bi(x), x, 2) - x * airy_bi(x)
    0
"""

#*****************************************************************************
#      Copyright (C) 2010 Oscar Gerardo Lazo Arjona <algebraicamente@gmail.com>
#      Copyright (C) 2012 Douglas McNeil <dsm054@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.symbolic.function import BuiltinFunction
from sage.symbolic.expression import Expression
from sage.symbolic.ring import SR
from sage.structure.coerce import parent as sage_structure_coerce_parent
from sage.functions.other import gamma
from sage.rings.integer_ring import ZZ
from sage.rings.real_double import RDF
from sage.rings.rational import Rational as R
from sage.calculus.functional import derivative


class FunctionAiryAiGeneral(BuiltinFunction):
    def __init__(self):
        r"""
        The generalized derivative of the Airy Ai function

        INPUT:

        - ``alpha`` -- Return the `\alpha`-th order fractional derivative with
          respect to `z`.
          For `\alpha = n = 1,2,3,\ldots` this gives the derivative
          `\operatorname{Ai}^{(n)}(z)`, and for `\alpha = -n = -1,-2,-3,\ldots`
          this gives the `n`-fold iterated integral.

        .. math ::

            f_0(z) = \operatorname{Ai}(z)

            f_n(z) = \int_0^z f_{n-1}(t) dt

        - ``x`` -- The argument of the function

        EXAMPLES::

            sage: from sage.functions.airy import airy_ai_general
            sage: x, n = var('x n')
            sage: airy_ai_general(-2, x)
            airy_ai(-2, x)
            sage: derivative(airy_ai_general(-2, x), x)
            airy_ai(-1, x)
            sage: airy_ai_general(n, x)
            airy_ai(n, x)
            sage: derivative(airy_ai_general(n, x), x)
            airy_ai(n + 1, x)
        """
        BuiltinFunction.__init__(self, "airy_ai", nargs=2,
                                 latex_name=r"\operatorname{Ai}")

    def _derivative_(self, alpha, x, diff_param=None):
        """
        EXAMPLES::

            sage: from sage.functions.airy import airy_ai_general
            sage: x, n = var('x n')
            sage: derivative(airy_ai_general(n, x), x)
            airy_ai(n + 1, x)
            sage: derivative(airy_ai_general(n, x), n)
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot differentiate airy_ai in the first parameter
        """
        if diff_param == 0:
            raise NotImplementedError("cannot differentiate airy_ai in the"
                                      " first parameter")
        return airy_ai_general(alpha + 1, x)

    def _eval_(self, alpha, x):
        """
        EXAMPLES::

            sage: from sage.functions.airy import airy_ai_general
            sage: x, n = var('x n')
            sage: airy_ai_general(-2, 1.0)
            0.136645379421096
            sage: airy_ai_general(n, 1.0)
            airy_ai(n, 1.00000000000000)
        """
        if not isinstance(x, Expression) and \
                not isinstance(alpha, Expression):
            if self._is_numerical(x):
                return self._evalf_(alpha, x)
            if alpha == 0:
                return airy_ai_simple(x)
            if alpha == 1:
                return airy_ai_prime(x)
            if alpha == 2:
                return x*airy_ai_simple(x)
        else:
            return None

    def _evalf_(self, alpha, x, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: from sage.functions.airy import airy_ai_general
            sage: airy_ai_general(-2, 1.0)
            0.136645379421096
        """
        import mpmath
        from sage.libs.mpmath import utils as mpmath_utils
        return mpmath_utils.call(mpmath.airyai, x, derivative=alpha,
                                 parent=parent)


class FunctionAiryAiSimple(BuiltinFunction):
    def __init__(self):
        """
        The class for the Airy Ai function.

        EXAMPLES::

            sage: from sage.functions.airy import airy_ai_simple
            sage: f = airy_ai_simple(x); f
            airy_ai(x)
        """
        BuiltinFunction.__init__(self, "airy_ai",
                                 latex_name=r'\operatorname{Ai}',
                                 conversions=dict(mathematica='AiryAi',
                                                  maxima='airy_ai'))

    def _derivative_(self, x, diff_param=None):
        """
        EXAMPLES::

            sage: from sage.functions.airy import airy_ai_simple
            sage: derivative(airy_ai_simple(x), x)
            airy_ai_prime(x)
        """
        return airy_ai_prime(x)

    def _eval_(self, x):
        """
        EXAMPLES::

            sage: from sage.functions.airy import airy_ai_simple
            sage: airy_ai_simple(0)
            1/3*3^(1/3)/gamma(2/3)
            sage: airy_ai_simple(0.0)
            0.355028053887817
            sage: airy_ai_simple(I)
            airy_ai(I)
            sage: airy_ai_simple(1.0 * I)
            0.331493305432141 - 0.317449858968444*I
        """
        if x == 0:
            r = ZZ(2) / 3
            return 1 / (3 ** (r) * gamma(r))

    def _evalf_(self, x, **kwargs):
        """
        EXAMPLES::

            sage: from sage.functions.airy import airy_ai_simple
            sage: airy_ai_simple(0.0)
            0.355028053887817
            sage: airy_ai_simple(1.0 * I)
            0.331493305432141 - 0.317449858968444*I

        We can use several methods for numerical evaluation::

            sage: airy_ai_simple(3).n(algorithm='mpmath')
            0.00659113935746072
            sage: airy_ai_simple(3).n(algorithm='mpmath', prec=100)
            0.0065911393574607191442574484080
            sage: airy_ai_simple(3).n(algorithm='scipy')  # rel tol 1e-10
            0.006591139357460719
            sage: airy_ai_simple(I).n(algorithm='scipy')  # rel tol 1e-10
            0.33149330543214117 - 0.3174498589684438*I
            
        TESTS::

            sage: parent(airy_ai_simple(3).n(algorithm='scipy'))                                          
            Real Field with 53 bits of precision
            sage: airy_ai_simple(3).n(algorithm='scipy', prec=200)
            Traceback (most recent call last):
            ...
            NotImplementedError: airy_ai not implemented for precision > 53
        """
        algorithm = kwargs.get('algorithm', 'mpmath') or 'mpmath'
        parent = kwargs.get('parent')
        if algorithm == 'scipy':
            if hasattr(parent, 'prec') and parent.prec() > 53:
                raise NotImplementedError("%s not implemented for precision > 53"%self.name())
            from sage.rings.all import RR, CC
            from sage.functions.other import real,imag
            from scipy.special import airy as airy
            if x in RR:
                y = airy(real(x))[0]
                if parent is None:
                    return RR(y)
            else:
                y = airy(complex(real(x),imag(x)))[0]
                if parent is None:
                    return CC(y)
            return parent(y)
        elif algorithm == 'mpmath':
            import mpmath
            from sage.libs.mpmath import utils as mpmath_utils
            return mpmath_utils.call(mpmath.airyai, x, parent=parent)
        else:
            raise ValueError("unknown algorithm '%s'" % algorithm)


class FunctionAiryAiPrime(BuiltinFunction):
    def __init__(self):
        """
        The derivative of the Airy Ai function; see :func:`airy_ai`
        for the full documentation.

        EXAMPLES::

            sage: x, n = var('x n')
            sage: airy_ai_prime(x)
            airy_ai_prime(x)
            sage: airy_ai_prime(0)
            -1/3*3^(2/3)/gamma(1/3)
        """
        BuiltinFunction.__init__(self, "airy_ai_prime",
                                 latex_name=r"\operatorname{Ai}'",
                                 conversions=dict(mathematica='AiryAiPrime',
                                                  maxima='airy_dai'))

    def _derivative_(self, x, diff_param=None):
        """
        EXAMPLES::

           sage: derivative(airy_ai_prime(x), x)
            x*airy_ai(x)
        """
        return x * airy_ai_simple(x)

    def _eval_(self, x):
        """
        EXAMPLES::

            sage: airy_ai_prime(0)
            -1/3*3^(2/3)/gamma(1/3)
            sage: airy_ai_prime(0.0)
            -0.258819403792807
        """
        if x == 0:
            r = ZZ(1) / 3
            return -1 / (3 ** (r) * gamma(r))

    def _evalf_(self, x, **kwargs):
        """
        EXAMPLES::

            sage: airy_ai_prime(0.0)
            -0.258819403792807

        We can use several methods for numerical evaluation::

            sage: airy_ai_prime(4).n(algorithm='mpmath')
            -0.00195864095020418
            sage: airy_ai_prime(4).n(algorithm='mpmath', prec=100)
            -0.0019586409502041789001381409184
            sage: airy_ai_prime(4).n(algorithm='scipy')    # rel tol 1e-10
            -0.00195864095020418
            sage: airy_ai_prime(I).n(algorithm='scipy')    # rel tol 1e-10
            -0.43249265984180707 + 0.09804785622924324*I
            
        TESTS::

            sage: parent(airy_ai_prime(3).n(algorithm='scipy'))                                          
            Real Field with 53 bits of precision
            sage: airy_ai_prime(3).n(algorithm='scipy', prec=200)
            Traceback (most recent call last):
            ...
            NotImplementedError: airy_ai_prime not implemented for precision > 53
        """
        algorithm = kwargs.get('algorithm', 'mpmath') or 'mpmath'
        parent = kwargs.get('parent', None)
        if algorithm == 'scipy':
            if hasattr(parent, 'prec') and parent.prec() > 53:
                raise NotImplementedError("%s not implemented for precision > 53"%self.name())
            from sage.rings.all import RR, CC
            from sage.functions.other import real,imag
            from scipy.special import airy as airy
            if x in RR:
                y = airy(real(x))[1]
                if parent is None:
                    return RR(y)
            else:
                y = airy(complex(real(x),imag(x)))[1]
                if parent is None:
                    return CC(y)
            return parent(y)
        elif algorithm == 'mpmath':
            import mpmath
            from sage.libs.mpmath import utils as mpmath_utils
            return mpmath_utils.call(mpmath.airyai, x, derivative=1,
                                     parent=parent)
        else:
            raise ValueError("unknown algorithm '%s'" % algorithm)

airy_ai_general = FunctionAiryAiGeneral()
airy_ai_simple = FunctionAiryAiSimple()
airy_ai_prime = FunctionAiryAiPrime()


def airy_ai(alpha, x=None, hold_derivative=True, **kwds):
    r"""
    The Airy Ai function

    The Airy Ai function `\operatorname{Ai}(x)` is (along with
    `\operatorname{Bi}(x)`) one of the two linearly independent standard
    solutions to the Airy differential equation `f''(x) - x f(x) = 0`. It is
    defined by the initial conditions:

    .. math::

        \operatorname{Ai}(0)=\frac{1}{2^{2/3} \Gamma\left(\frac{2}{3}\right)},

        \operatorname{Ai}'(0)=-\frac{1}{2^{1/3}\Gamma\left(\frac{1}{3}\right)}.

    Another way to define the Airy Ai function is:

    .. math::

        \operatorname{Ai}(x)=\frac{1}{\pi}\int_0^\infty
        \cos\left(\frac{1}{3}t^3+xt\right) dt.

    INPUT:

    - ``alpha`` -- Return the `\alpha`-th order fractional derivative with
      respect to `z`.
      For `\alpha = n = 1,2,3,\ldots` this gives the derivative
      `\operatorname{Ai}^{(n)}(z)`, and for `\alpha = -n = -1,-2,-3,\ldots`
      this gives the `n`-fold iterated integral.

    .. math::

        f_0(z) = \operatorname{Ai}(z)

        f_n(z) = \int_0^z f_{n-1}(t) dt

    - ``x`` -- The argument of the function

    - ``hold_derivative`` -- Whether or not to stop from returning higher
      derivatives in terms of `\operatorname{Ai}(x)` and
      `\operatorname{Ai}'(x)`

    .. SEEALSO:: :func:`airy_bi`

    EXAMPLES::

        sage: n, x = var('n x')
        sage: airy_ai(x)
        airy_ai(x)

    It can return derivatives or integrals::

        sage: airy_ai(2, x)
        airy_ai(2, x)
        sage: airy_ai(1, x, hold_derivative=False)
        airy_ai_prime(x)
        sage: airy_ai(2, x, hold_derivative=False)
        x*airy_ai(x)
        sage: airy_ai(-2, x, hold_derivative=False)
        airy_ai(-2, x)
        sage: airy_ai(n, x)
        airy_ai(n, x)

    It can be evaluated symbolically or numerically for real or complex
    values::

        sage: airy_ai(0)
        1/3*3^(1/3)/gamma(2/3)
        sage: airy_ai(0.0)
        0.355028053887817
        sage: airy_ai(I)
        airy_ai(I)
        sage: airy_ai(1.0*I)
        0.331493305432141 - 0.317449858968444*I

    The functions can be evaluated numerically either using mpmath. which
    can compute the values to arbitrary precision, and scipy::

        sage: airy_ai(2).n(prec=100)
        0.034924130423274379135322080792
        sage: airy_ai(2).n(algorithm='mpmath', prec=100)
        0.034924130423274379135322080792
        sage: airy_ai(2).n(algorithm='scipy')  # rel tol 1e-10
        0.03492413042327323

    And the derivatives can be evaluated::

        sage: airy_ai(1, 0)
        -1/3*3^(2/3)/gamma(1/3)
        sage: airy_ai(1, 0.0)
        -0.258819403792807

    Plots::

        sage: plot(airy_ai(x), (x, -10, 5)) + plot(airy_ai_prime(x),
        ....:  (x, -10, 5), color='red')
        Graphics object consisting of 2 graphics primitives

    **References**

    - Abramowitz, Milton; Stegun, Irene A., eds. (1965), "Chapter 10"

    - :wikipedia:`Airy_function`
    """
    # We catch the case with no alpha
    if x is None:
        x = alpha
        return airy_ai_simple(x, **kwds)

    # We take care of all other cases.
    if not alpha in ZZ and not isinstance(alpha, Expression):
        return airy_ai_general(alpha, x, **kwds)
    if hold_derivative:
        return airy_ai_general(alpha, x, **kwds)
    elif alpha == 0:
        return airy_ai_simple(x, **kwds)
    elif alpha == 1:
        return airy_ai_prime(x, **kwds)
    elif alpha > 1:
        # We use a different variable here because if x is a
        # particular value, we would be differentiating a constant
        # which would return 0. What we want is the value of
        # the derivative at the value and not the derivative of
        # a particular value of the function.
        v = SR.symbol()
        return derivative(airy_ai_simple(v, **kwds), v, alpha).subs({v: x})
    else:
        return airy_ai_general(alpha, x, **kwds)

########################################################################
########################################################################


class FunctionAiryBiGeneral(BuiltinFunction):
    def __init__(self):
        r"""
        The generalized derivative of the Airy Bi function.

        INPUT:

        - ``alpha`` -- Return the `\alpha`-th order fractional derivative with
          respect to `z`.
          For `\alpha = n = 1,2,3,\ldots` this gives the derivative
          `\operatorname{Bi}^{(n)}(z)`, and for `\alpha = -n = -1,-2,-3,\ldots`
          this gives the `n`-fold iterated integral.

        .. math::

            f_0(z) = \operatorname{Bi}(z)

            f_n(z) = \int_0^z f_{n-1}(t) dt

        - ``x`` -- The argument of the function

        EXAMPLES::

            sage: from sage.functions.airy import airy_bi_general
            sage: x, n = var('x n')
            sage: airy_bi_general(-2, x)
            airy_bi(-2, x)
            sage: derivative(airy_bi_general(-2, x), x)
            airy_bi(-1, x)
            sage: airy_bi_general(n, x)
            airy_bi(n, x)
            sage: derivative(airy_bi_general(n, x), x)
            airy_bi(n + 1, x)
        """
        BuiltinFunction.__init__(self, "airy_bi", nargs=2,
                                 latex_name=r"\operatorname{Bi}")

    def _derivative_(self, alpha, x, diff_param=None):
        """
        EXAMPLES::

            sage: from sage.functions.airy import airy_bi_general
            sage: x, n = var('x n')
            sage: derivative(airy_bi_general(n, x), x)
            airy_bi(n + 1, x)
            sage: derivative(airy_bi_general(n, x), n)
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot differentiate airy_bi in the first parameter
        """
        if diff_param == 0:
            raise NotImplementedError("cannot differentiate airy_bi in the"
                                      " first parameter")
        return airy_bi_general(alpha + 1, x)

    def _eval_(self, alpha, x):
        """
        EXAMPLES::

            sage: from sage.functions.airy import airy_bi_general
            sage: x, n = var('x n')
            sage: airy_bi_general(-2, 1.0)
            0.388621540699059
            sage: airy_bi_general(n, 1.0)
            airy_bi(n, 1.00000000000000)
        """
        if not isinstance(x, Expression) and \
                not isinstance(alpha, Expression):
            if alpha == 0:
                return airy_bi_simple(x)
            if alpha == 1:
                return airy_bi_prime(x)
            if alpha == 2:
                return x*airy_bi_simple(x)

    def _evalf_(self, alpha, x, **kwargs):
        """
        EXAMPLES::

            sage: from sage.functions.airy import airy_bi_general
            sage: airy_bi_general(-2, 1.0)
            0.388621540699059

        """
        parent = kwargs.get('parent')
        import mpmath
        from sage.libs.mpmath import utils as mpmath_utils
        return mpmath_utils.call(mpmath.airybi, x, derivative=alpha,
                                 parent=parent)


class FunctionAiryBiSimple(BuiltinFunction):
    def __init__(self):
        """
        The class for the Airy Bi function.

        EXAMPLES::

            sage: from sage.functions.airy import airy_bi_simple
            sage: f = airy_bi_simple(x); f
            airy_bi(x)
        """
        BuiltinFunction.__init__(self, "airy_bi",
                                 latex_name=r'\operatorname{Bi}',
                                 conversions=dict(mathematica='AiryBi',
                                                  maxima='airy_bi'))

    def _derivative_(self, x, diff_param=None):
        """
        EXAMPLES::

            sage: from sage.functions.airy import airy_bi_simple
            sage: derivative(airy_bi_simple(x), x)
            airy_bi_prime(x)
        """
        return airy_bi_prime(x)

    def _eval_(self, x):
        """
        EXAMPLES::

            sage: from sage.functions.airy import airy_bi_simple
            sage: airy_bi_simple(0)
            1/3*3^(5/6)/gamma(2/3)
            sage: airy_bi_simple(0.0)
            0.614926627446001
            sage: airy_bi_simple(0).n() == airy_bi(0.0)
            True
            sage: airy_bi_simple(I)
            airy_bi(I)
            sage: airy_bi_simple(1.0 * I)
            0.648858208330395 + 0.344958634768048*I
        """
        if x == 0:
            one_sixth = ZZ(1) / 6
            return 1 / (3 ** (one_sixth) * gamma(4 * one_sixth))

    def _evalf_(self, x, **kwargs):
        """
        EXAMPLES::

            sage: from sage.functions.airy import airy_bi_simple
            sage: airy_bi_simple(0.0)
            0.614926627446001
            sage: airy_bi_simple(1.0 * I)
            0.648858208330395 + 0.344958634768048*I

        We can use several methods for numerical evaluation::

            sage: airy_bi_simple(3).n(algorithm='mpmath')
            14.0373289637302
            sage: airy_bi_simple(3).n(algorithm='mpmath', prec=100)
            14.037328963730232031740267314
            sage: airy_bi_simple(3).n(algorithm='scipy')  # rel tol 1e-10
            14.037328963730136
            sage: airy_bi_simple(I).n(algorithm='scipy')  # rel tol 1e-10
            0.648858208330395 + 0.34495863476804844*I
            
        TESTS::

            sage: parent(airy_bi_simple(3).n(algorithm='scipy'))                                          
            Real Field with 53 bits of precision
            sage: airy_bi_simple(3).n(algorithm='scipy', prec=200)
            Traceback (most recent call last):
            ...
            NotImplementedError: airy_bi not implemented for precision > 53
        """
        algorithm = kwargs.get('algorithm', 'mpmath') or 'mpmath'
        parent = kwargs.get('parent', None)
        if algorithm == 'scipy':
            if hasattr(parent, 'prec') and parent.prec() > 53:
                raise NotImplementedError("%s not implemented for precision > 53"%self.name())
            from sage.rings.all import RR, CC
            from sage.functions.other import real,imag
            from scipy.special import airy as airy
            if x in RR:
                y = airy(real(x))[2]
                if parent is None:
                    return RR(y)
            else:
                y = airy(complex(real(x),imag(x)))[2]
                if parent is None:
                    return CC(y)
            return parent(y)
        elif algorithm == 'mpmath':
            import mpmath
            from sage.libs.mpmath import utils as mpmath_utils
            return mpmath_utils.call(mpmath.airybi, x, parent=parent)
        else:
            raise ValueError("unknown algorithm '%s'" % algorithm)


class FunctionAiryBiPrime(BuiltinFunction):
    def __init__(self):
        """
        The derivative of the Airy Bi function; see :func:`airy_bi`
        for the full documentation.

        EXAMPLES::

            sage: x, n = var('x n')
            sage: airy_bi_prime(x)
            airy_bi_prime(x)
            sage: airy_bi_prime(0)
            3^(1/6)/gamma(1/3)
        """
        BuiltinFunction.__init__(self, "airy_bi_prime",
                                 latex_name=r"\operatorname{Bi}'",
                                 conversions=dict(mathematica='AiryBiPrime',
                                                  maxima='airy_dbi'))

    def _derivative_(self, x, diff_param=None):
        """
        EXAMPLES::

            sage: derivative(airy_bi_prime(x), x)
            x*airy_bi(x)
        """
        return x * airy_bi_simple(x)

    def _eval_(self, x):
        """
        EXAMPLES::

            sage: airy_bi_prime(0)
            3^(1/6)/gamma(1/3)
            sage: airy_bi_prime(0.0)
            0.448288357353826
        """
        if x == 0:
            one_sixth = ZZ(1) / 6
            return 3 ** (one_sixth) / gamma(2 * one_sixth)

    def _evalf_(self, x, **kwargs):
        """
        EXAMPLES::

            sage: airy_bi_prime(0.0)
            0.448288357353826

        We can use several methods for numerical evaluation::

            sage: airy_bi_prime(4).n(algorithm='mpmath')
            161.926683504613
            sage: airy_bi_prime(4).n(algorithm='mpmath', prec=100)
            161.92668350461340184309492429
            sage: airy_bi_prime(4).n(algorithm='scipy')  # rel tol 1e-10
            161.92668350461398
            sage: airy_bi_prime(I).n(algorithm='scipy')  # rel tol 1e-10
            0.135026646710819 - 0.1288373867812549*I
            
        TESTS::

            sage: parent(airy_bi_prime(3).n(algorithm='scipy'))                                          
            Real Field with 53 bits of precision
            sage: airy_bi_prime(3).n(algorithm='scipy', prec=200)
            Traceback (most recent call last):
            ...
            NotImplementedError: airy_bi_prime not implemented for precision > 53
        """
        algorithm = kwargs.get('algorithm', 'mpmath') or 'mpmath'
        parent = kwargs.get('parent', None)
        if algorithm == 'scipy':
            if hasattr(parent, 'prec') and parent.prec() > 53:
                raise NotImplementedError("%s not implemented for precision > 53"%self.name())
            from sage.rings.all import RR, CC
            from sage.functions.other import real,imag
            from scipy.special import airy as airy
            if x in RR:
                y = airy(real(x))[3]
                if parent is None:
                    return RR(y)
            else:
                y = airy(complex(real(x),imag(x)))[3]
                if parent is None:
                    return CC(y)
            return parent(y)
        elif algorithm == 'mpmath':
            import mpmath
            from sage.libs.mpmath import utils as mpmath_utils
            return mpmath_utils.call(mpmath.airybi, x, derivative=1,
                                     parent=parent)
        else:
            raise ValueError("unknown algorithm '%s'" % algorithm)

airy_bi_general = FunctionAiryBiGeneral()
airy_bi_simple = FunctionAiryBiSimple()
airy_bi_prime = FunctionAiryBiPrime()


def airy_bi(alpha, x=None, hold_derivative=True, **kwds):
    r"""
    The Airy Bi function

    The Airy Bi function `\operatorname{Bi}(x)` is (along with
    `\operatorname{Ai}(x)`) one of the two linearly independent standard
    solutions to the Airy differential equation `f''(x) - x f(x) = 0`. It is
    defined by the initial conditions:

    .. math::

        \operatorname{Bi}(0)=\frac{1}{3^{1/6} \Gamma\left(\frac{2}{3}\right)},

        \operatorname{Bi}'(0)=\frac{3^{1/6}}{ \Gamma\left(\frac{1}{3}\right)}.

    Another way to define the Airy Bi function is:

    .. math::

        \operatorname{Bi}(x)=\frac{1}{\pi}\int_0^\infty
        \left[ \exp\left( xt -\frac{t^3}{3} \right)
        +\sin\left(xt + \frac{1}{3}t^3\right) \right ] dt.

    INPUT:

    - ``alpha`` -- Return the `\alpha`-th order fractional derivative with
      respect to `z`.
      For `\alpha = n = 1,2,3,\ldots` this gives the derivative
      `\operatorname{Bi}^{(n)}(z)`, and for `\alpha = -n = -1,-2,-3,\ldots`
      this gives the `n`-fold iterated integral.

    .. math ::

        f_0(z) = \operatorname{Bi}(z)

        f_n(z) = \int_0^z f_{n-1}(t) dt

    - ``x`` -- The argument of the function

    - ``hold_derivative`` -- Whether or not to stop from returning higher
      derivatives in terms of `\operatorname{Bi}(x)` and
      `\operatorname{Bi}'(x)`

    .. SEEALSO:: :func:`airy_ai`

    EXAMPLES::

        sage: n, x = var('n x')
        sage: airy_bi(x)
        airy_bi(x)

    It can return derivatives or integrals::

        sage: airy_bi(2, x)
        airy_bi(2, x)
        sage: airy_bi(1, x, hold_derivative=False)
        airy_bi_prime(x)
        sage: airy_bi(2, x, hold_derivative=False)
        x*airy_bi(x)
        sage: airy_bi(-2, x, hold_derivative=False)
        airy_bi(-2, x)
        sage: airy_bi(n, x)
        airy_bi(n, x)

    It can be evaluated symbolically or numerically for real or complex
    values::

        sage: airy_bi(0)
        1/3*3^(5/6)/gamma(2/3)
        sage: airy_bi(0.0)
        0.614926627446001
        sage: airy_bi(I)
        airy_bi(I)
        sage: airy_bi(1.0*I)
        0.648858208330395 + 0.344958634768048*I

    The functions can be evaluated numerically using mpmath,
    which can compute the values to arbitrary precision, and scipy::

        sage: airy_bi(2).n(prec=100)
        3.2980949999782147102806044252
        sage: airy_bi(2).n(algorithm='mpmath', prec=100)
        3.2980949999782147102806044252
        sage: airy_bi(2).n(algorithm='scipy')  # rel tol 1e-10
        3.2980949999782134

    And the derivatives can be evaluated::

        sage: airy_bi(1, 0)
        3^(1/6)/gamma(1/3)
        sage: airy_bi(1, 0.0)
        0.448288357353826

    Plots::

        sage: plot(airy_bi(x), (x, -10, 5)) + plot(airy_bi_prime(x),
        ....:  (x, -10, 5), color='red')
        Graphics object consisting of 2 graphics primitives

    **References**

    - Abramowitz, Milton; Stegun, Irene A., eds. (1965), "Chapter 10"

    - :wikipedia:`Airy_function`
    """
    # We catch the case with no alpha
    if x is None:
        x = alpha
        return airy_bi_simple(x, **kwds)

    # We take care of all other cases.
    if not alpha in ZZ and not isinstance(alpha, Expression):
        return airy_bi_general(alpha, x, **kwds)
    if hold_derivative:
        return airy_bi_general(alpha, x, **kwds)
    elif alpha == 0:
        return airy_bi_simple(x, **kwds)
    elif alpha == 1:
        return airy_bi_prime(x, **kwds)
    elif alpha > 1:
        # We use a different variable here because if x is a
        # particular value, we would be differentiating a constant
        # which would return 0. What we want is the value of
        # the derivative at the value and not the derivative of
        # a particular value of the function.
        v = SR.symbol()
        return derivative(airy_bi_simple(v, **kwds), v, alpha).subs({v: x})
    else:
        return airy_bi_general(alpha, x, **kwds)
