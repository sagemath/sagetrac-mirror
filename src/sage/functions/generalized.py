r"""
Generalized Functions

Sage implements several generalized functions (also known as
distributions) such as Dirac delta, Heaviside step functions. These
generalized functions can be manipulated within Sage like any other
symbolic functions.


AUTHORS:

- Golam Mortuza Hossain (2009-06-26): initial version


EXAMPLES:

Dirac delta function::

    sage: dirac_delta(x)
    dirac_delta(x)

Heaviside step function::

    sage: heaviside(x)
    heaviside(x)

Unit step function::

    sage: unit_step(x)
    unit_step(x)

Signum (sgn) function::

    sage: sgn(x)
    sgn(x)

Kronecker delta function::

    sage: m,n=var('m,n')
    sage: kronecker_delta(m,n)
    kronecker_delta(m, n)

"""

##############################################################################
#
#       Copyright (C) 2009 Golam Mortuza Hossain <gmhossain@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL v2+)
#                  http://www.gnu.org/licenses/
#
##############################################################################

from sage.symbolic.function import BuiltinFunction
from sage.rings.all import ComplexIntervalField, ZZ

class FunctionDiracDelta(BuiltinFunction):
    r"""
    The Dirac delta (generalized) function, `\delta(x)` (``dirac_delta(x)``).

    INPUT:

    -  ``x`` - a real number or a symbolic expression

    DEFINITION:

    Dirac delta function `\delta(x)`, is defined in Sage as:

        `\delta(x) = 0` for real `x \ne 0` and
        `\int_{-\infty}^{\infty} \delta(x) dx = 1`

    Its alternate definition with respect to an arbitrary test
    function `f(x)` is

        `\int_{-\infty}^{\infty} f(x) \delta(x-a) dx = f(a)`

    EXAMPLES::

        sage: dirac_delta(1)
        0
        sage: dirac_delta(0)
        dirac_delta(0)
        sage: dirac_delta(x)
        dirac_delta(x)
        sage: integrate(dirac_delta(x), x, -1, 1, algorithm='sympy')
        1

    REFERENCES:

    -  http://en.wikipedia.org/wiki/Dirac_delta_function

    """
    def __init__(self):
        r"""
        The Dirac delta (generalized) function, ``dirac_delta(x)``.

        INPUT:

        -  ``x`` - a real number or a symbolic expression

        EXAMPLES::

            sage: dirac_delta(1)
            0
            sage: dirac_delta(0)
            dirac_delta(0)
            sage: dirac_delta(x)
            dirac_delta(x)
            sage: latex(dirac_delta(x))
            \delta\left(x\right)

            sage: loads(dumps(dirac_delta(x)))
            dirac_delta(x)
            sage: dirac_delta(x)._sympy_()
            DiracDelta(x)
        """
        BuiltinFunction.__init__(self, "dirac_delta", latex_name=r"\delta",
                                   conversions=dict(maxima='delta',
                                                  mathematica='DiracDelta',
                                                  sympy='DiracDelta'))

    def _eval_(self, x):
        """
        INPUT:

        -  ``x`` - a real number or a symbolic expression

        EXAMPLES::

            sage: dirac_delta(1)
            0
            sage: dirac_delta(0)
            dirac_delta(0)
            sage: dirac_delta(x)
            dirac_delta(x)
            sage: dirac_delta(exp(-10000000000000000000))
            0

        Evaluation test::

            sage: dirac_delta(x).subs(x=1)
            0
        """
        try:
            return self._evalf_(x)
        except (TypeError,ValueError):      # x is symbolic
            pass
        return None

    def _evalf_(self, x, **kwds):
        """
        TESTS::

            sage: h(x) = dirac_delta(x)
            sage: h(pi).numerical_approx()
            0.000000000000000
        """
        approx_x = ComplexIntervalField()(x)
        if bool(approx_x.imag() == 0):      # x is real
            if bool(approx_x.real() == 0):  # x is zero
                return None
            else:
                return 0
        raise ValueError("Numeric evaluation of symbolic expression")

dirac_delta = FunctionDiracDelta()

class FunctionHeaviside(BuiltinFunction):
    r"""
    The Heaviside step function, `H(x)` (``heaviside(x)``).

    INPUT:

    -  ``x`` - a real number or a symbolic expression

    DEFINITION:

    The Heaviside step function, `H(x)` is defined in Sage as:

        `H(x) = 0` for `x < 0` and `H(x) = 1` for `x > 0`

    EXAMPLES::

        sage: heaviside(-1)
        0
        sage: heaviside(1)
        1
        sage: heaviside(0)
        heaviside(0)
        sage: heaviside(x)
        heaviside(x)

    TESTS::

        sage: heaviside(x)._sympy_()
        Heaviside(x)

    REFERENCES:

    -  :wikipedia:`Heaviside_function`

    """
    def __init__(self):
        r"""
        The Heaviside step function, ``heaviside(x)``.

        INPUT:

        -  ``x`` - a real number or a symbolic expression

        EXAMPLES::

            sage: heaviside(-1)
            0
            sage: heaviside(1)
            1
            sage: heaviside(0)
            heaviside(0)
            sage: heaviside(x)
            heaviside(x)
            sage: latex(heaviside(x))
            H\left(x\right)
            sage: heaviside(x)._sympy_()
            Heaviside(x)
        """
        BuiltinFunction.__init__(self, "heaviside", latex_name="H",
                                 conversions=dict(maxima='hstep',
                                                  mathematica='HeavisideTheta',
                                                  sympy='Heaviside'))

    def _eval_(self, x):
        """
        INPUT:

        -  ``x`` - a real number or a symbolic expression

        EXAMPLES::

            sage: heaviside(-1/2)
            0
            sage: heaviside(1)
            1
            sage: heaviside(0)
            heaviside(0)
            sage: heaviside(x)
            heaviside(x)
            sage: heaviside(exp(-1000000000000000000000))
            1

        Evaluation test::

            sage: heaviside(x).subs(x=1)
            1
            sage: heaviside(x).subs(x=-1)
            0

        ::

            sage: ex = heaviside(x)+1
            sage: t = loads(dumps(ex)); t
            heaviside(x) + 1
            sage: bool(t == ex)
            True
            sage: t.subs(x=1)
            2
        """
        try:
            return self._evalf_(x)
        except (TypeError,ValueError):      # x is symbolic
            pass
        return None

    def _evalf_(self, x, **kwds):
        """
        TESTS::

            sage: h(x) = heaviside(x)
            sage: h(pi).numerical_approx()
            1.00000000000000
        """
        approx_x = ComplexIntervalField()(x)
        if bool(approx_x.imag() == 0):      # x is real
            if bool(approx_x.real() == 0):  # x is zero
                return None
            # Now we have a non-zero real
            if bool((approx_x**(0.5)).imag() == 0): # Check: x > 0
                return 1
            else:
                return 0
        raise ValueError("Numeric evaluation of symbolic expression")

    def _derivative_(self, x, diff_param=None):
        """
        Derivative of Heaviside step function

        EXAMPLES::

            sage: heaviside(x).diff(x)
            dirac_delta(x)
        """
        return dirac_delta(x)

heaviside = FunctionHeaviside()

class FunctionUnitStep(BuiltinFunction):
    r"""
    The unit step function, `\mathrm{u}(x)` (``unit_step(x)``).

    INPUT:

    -  ``x`` - a real number or a symbolic expression

    DEFINITION:

    The unit step function, `\mathrm{u}(x)` is defined in Sage as:

        `\mathrm{u}(x) = 0` for `x < 0` and `\mathrm{u}(x) = 1` for `x \geq 0`

    EXAMPLES::

        sage: unit_step(-1)
        0
        sage: unit_step(1)
        1
        sage: unit_step(0)
        1
        sage: unit_step(x)
        unit_step(x)
    """
    def __init__(self):
        r"""
        The unit step function, ``unit_step(x)``.

        INPUT:

        -  ``x`` - a real number or a symbolic expression

        EXAMPLES:

            sage: unit_step(-1)
            0
            sage: unit_step(1)
            1
            sage: unit_step(0)
            1
            sage: unit_step(x)
            unit_step(x)
            sage: latex(unit_step(x))
            \mathrm{u}\left(x\right)

        TESTS::

            sage: t = loads(dumps(unit_step(x)+1)); t
            unit_step(x) + 1
            sage: t.subs(x=0)
            2
        """
        BuiltinFunction.__init__(self, "unit_step", latex_name=r"\mathrm{u}",
                                   conversions=dict(mathematica='UnitStep'))

    def _eval_(self, x):
        """
        INPUT:

        -  ``x`` - a real number or a symbolic expression

        EXAMPLES::

            sage: unit_step(-1)
            0
            sage: unit_step(1)
            1
            sage: unit_step(0)
            1
            sage: unit_step(x)
            unit_step(x)
            sage: unit_step(-exp(-10000000000000000000))
            0

        Evaluation test::

            sage: unit_step(x).subs(x=1)
            1
            sage: unit_step(x).subs(x=0)
            1
        """
        try:
            return self._evalf_(x)
        except (TypeError,ValueError):      # x is symbolic
            pass
        return None

    def _evalf_(self, x, **kwds):
        """
        TESTS::

            sage: h(x) = unit_step(x)
            sage: h(pi).numerical_approx()
            1.00000000000000
        """
        approx_x = ComplexIntervalField()(x)
        if bool(approx_x.imag() == 0):      # x is real
            if bool(approx_x.real() == 0):  # x is zero
                return 1
            # Now we have a non-zero real
            if bool((approx_x**(0.5)).imag() == 0): # Check: x > 0
                return 1
            else:
                return 0
        raise ValueError("Numeric evaluation of symbolic expression")

    def _derivative_(self, x, diff_param=None):
        """
        Derivative of unit step function

        EXAMPLES::

            sage: unit_step(x).diff(x)
            dirac_delta(x)
        """
        return dirac_delta(x)

unit_step = FunctionUnitStep()

class FunctionSignum(BuiltinFunction):
    r"""
    The signum or sgn function `\mathrm{sgn}(x)` (``sgn(x)``).

    INPUT:

    -  ``x`` - a real number or a symbolic expression

    DEFINITION:

    The sgn function, `\mathrm{sgn}(x)` is defined as:

        `\mathrm{sgn}(x) =  1` for `x > 0`,
        `\mathrm{sgn}(x) =  0` for `x = 0` and
        `\mathrm{sgn}(x) = -1` for `x < 0`

    EXAMPLES::

        sage: sgn(-1)
        -1
        sage: sgn(1)
        1
        sage: sgn(0)
        0
        sage: sgn(x)
        sgn(x)

    We can also use ``sign``::

        sage: sign(1)
        1
        sage: sign(0)
        0
        sage: a = AA(-5).nth_root(7)
        sage: sign(a)
        -1

    TESTS:

    Check if conversion to sympy works :trac:`11921`::

        sage: sgn(x)._sympy_()
        sign(x)

    REFERENCES:

    -  http://en.wikipedia.org/wiki/Sign_function

    """
    def __init__(self):
        r"""
        The sgn function, ``sgn(x)``.

        EXAMPLES:

            sage: sgn(-1)
            -1
            sage: sgn(1)
            1
            sage: sgn(0)
            0
            sage: sgn(x)
            sgn(x)
            sage: sgn(x)._sympy_()
            sign(x)
        """
        BuiltinFunction.__init__(self, "sgn", latex_name=r"\mathrm{sgn}",
                conversions=dict(maxima='signum',mathematica='Sign',sympy='sign'),
                alt_name="sign")

    def _eval_(self, x):
        """

        EXAMPLES::

            sage: sgn(-1)
            -1
            sage: sgn(1)
            1
            sage: sgn(0)
            0
            sage: sgn(x)
            sgn(x)
            sage: sgn(-exp(-10000000000000000000))
            -1

        Evaluation test::

            sage: sgn(x).subs(x=1)
            1
            sage: sgn(x).subs(x=0)
            0
            sage: sgn(x).subs(x=-1)
            -1

        More tests::

            sage: sign(RR(2))
            1
            sage: sign(RDF(2))
            1
            sage: sign(AA(-2))
            -1
            sage: sign(AA(0))
            0
        """
        try:
            return self._evalf_(x)
        except (TypeError,ValueError):      # x is symbolic
            pass
        return None

    def _evalf_(self, x, **kwds):
        """
        TESTS:

        Check that :trac:`16587` is fixed::

            sage: M = sgn(3/2, hold=True); M
            sgn(3/2)
            sage: M.n()
            1
            sage: h(x) = sgn(x)
            sage: h(pi).numerical_approx()
            1.00000000000000
        """
        if hasattr(x,'sign'): # First check if x has a sign method
            return x.sign()
        if hasattr(x,'sgn'): # or a sgn method
            return x.sgn()
        approx_x = ComplexIntervalField()(x)
        if bool(approx_x.imag() == 0):      # x is real
            if bool(approx_x.real() == 0):  # x is zero
                return ZZ(0)
            # Now we have a non-zero real
            if bool((approx_x**(0.5)).imag() == 0): # Check: x > 0
                return ZZ(1)
            else:
                return ZZ(-1)
        raise ValueError("Numeric evaluation of symbolic expression")

    def _derivative_(self, x, diff_param=None):
        """
        Derivative of sgn function

        EXAMPLES::

            sage: sgn(x).diff(x)
            2*dirac_delta(x)
        """
        assert diff_param == 0
        return 2*dirac_delta(x)

sgn = FunctionSignum()
sign = sgn

class FunctionKroneckerDelta(BuiltinFunction):
    r"""
    The Kronecker delta function `\delta_{m,n}` (``kronecker_delta(m, n)``).

    INPUT:

    -  ``m`` - a number or a symbolic expression
    -  ``n`` - a number or a symbolic expression

    DEFINITION:

    Kronecker delta function `\delta_{m,n}` is defined as:

        `\delta_{m,n} = 0` for `m \ne n` and
        `\delta_{m,n} = 1` for `m = n`

    EXAMPLES::

        sage: kronecker_delta(1,2)
        0
        sage: kronecker_delta(1,1)
        1
        sage: m,n=var('m,n')
        sage: kronecker_delta(m,n)
        kronecker_delta(m, n)

    REFERENCES:

    - http://en.wikipedia.org/wiki/Kronecker_delta

    """
    def __init__(self):
        r"""
        The Kronecker delta function.

        EXAMPLES::

            sage: kronecker_delta(1,2)
            0
            sage: kronecker_delta(1,1)
            1
            sage: y = var('y')
            sage: kronecker_delta(x, y)._sympy_()
            KroneckerDelta(x, y)
        """
        BuiltinFunction.__init__(self, "kronecker_delta", nargs=2,
                                        conversions=dict(maxima='kron_delta',
                                        mathematica='KroneckerDelta',
                                        sympy='KroneckerDelta'))

    def _eval_(self, m, n):
        """
        The Kronecker delta function.

        EXAMPLES::

            sage: kronecker_delta(1,2)
            0
            sage: kronecker_delta(1,1)
            1

        Kronecker delta is a symmetric function. We keep arguments sorted to
        ensure that k_d(m, n) - k_d(n, m) cancels automatically::

            sage: x,y=var('x,y')
            sage: kronecker_delta(x, y)
            kronecker_delta(x, y)
            sage: kronecker_delta(y, x)
            kronecker_delta(x, y)
            sage: kronecker_delta(x,2*x)
            kronecker_delta(2*x, x)

        Evaluation test::

            sage: kronecker_delta(1,x).subs(x=1)
            1
        """
        try:
            return self._evalf_(m,n)
        except (TypeError,ValueError):      # x is symbolic
            pass
        return None

    def _evalf_(self, m, n, **kwds):
        """
        TESTS::

            sage: h(x) = kronecker_delta(3,x)
            sage: h(pi).numerical_approx()
            0.000000000000000
        """
        if bool(repr(m) > repr(n)):
            return kronecker_delta(n, m)

        x = m - n
        approx_x = ComplexIntervalField()(x)
        if bool(approx_x.imag() == 0):      # x is real
            if bool(approx_x.real() == 0):  # x is zero
                return 1
            else:
                return 0
        else:
            return 0            # x is complex
        raise ValueError("Numeric evaluation of symbolic expression")

    def _derivative_(self, *args, **kwds):
        """
        Derivative of Kronecker delta

        EXAMPLES::

            sage: kronecker_delta(x,1).diff(x)
            0
        """
        # Kronecker delta is non-zero (but finite) only in the set of
        # zero-measure unlike Dirac delta. Consequently, it is null
        # for the purpose of integration/differentiation. For *discrete sum*
        # Kronecker delta is however non-trivial.
        return 0

    def _print_latex_(self, m, n, **kwds):
        """
        Return latex expression

        EXAMPLES::

            sage: from sage.misc.latex import latex
            sage: m,n=var('m,n')
            sage: latex(kronecker_delta(m,n))
            \delta_{m,n}

        """
        from sage.misc.latex import latex
        return "\\delta_{%s,%s}"%(latex(m), latex(n))

kronecker_delta = FunctionKroneckerDelta()
