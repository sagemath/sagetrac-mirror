#   sympy.py - SymPy --> Sage interface
#   Author: (2017) Ralf Stephan <ralf@ark.in-berin.de>
#   Distributed under GNU GPL3, see www.gnu.org
#
#   The file consists of sage() methods that are added to SymPy objects.
#   Only Function objects where the names differ need their own sage() method.
#   Please be very specific with imports to minimalize startup time.
################################################################
from __future__ import absolute_import

#################         numbers and constants      ##############

def _sympysage_float(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import RealNumber as RN
        sage: assert SR(-1.34)._sympy_() == RN('-1.34')
        sage: assert SR(-1.34) == RN('-1.34')._sage_()
    """
    from sage.rings.real_mpfr import create_RealNumber
    return create_RealNumber(str(self))

def _sympysage_rational(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import Rational
        sage: assert SR(-5/7)._sympy_() == Rational(int(-5),int(7))
        sage: assert SR(-5/7) == Rational(int(-5),int(7))._sage_()
    """
    from sage.rings.integer import Integer
    from sage.rings.rational import Rational
    return Rational((Integer(self.p), Integer(self.q)))

def _sympysage_pinfty(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import oo as sinf
        sage: assert SR(oo)._sympy_() == sinf
        sage: assert SR(oo) == sinf._sage_()
    """
    from sage.rings.infinity import PlusInfinity
    return PlusInfinity()

def _sympysage_ninfty(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import oo as sinf
        sage: assert SR(-oo)._sympy_() == -sinf
        sage: assert SR(-oo) == (-sinf)._sage_()
    """
    from sage.rings.infinity import MinusInfinity
    return MinusInfinity()

def _sympysage_uinfty(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import zoo
        sage: assert unsigned_infinity._sympy_() == zoo
        sage: assert unsigned_infinity == zoo._sage_()
    """
    from sage.rings.infinity import unsigned_infinity
    return unsigned_infinity

def _sympysage_nan(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import nan as snan
        sage: assert NaN._sympy_() == snan
        sage: assert NaN == snan._sage_()
    """
    from sage.symbolic.constants import NaN
    return NaN

def _sympysage_e(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import E
        sage: assert e._sympy_() == E
        sage: assert e == E._sage_()
    """
    from sage.symbolic.constants import e
    return e

def _sympysage_pi(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import pi as spi
        sage: assert pi._sympy_() == spi
        sage: assert pi == spi._sage_()
    """
    from sage.symbolic.constants import pi
    return pi

def _sympysage_golden_ratio(self):
    """
    EXAMPLES::

        sage: from sympy.core.singleton import S
        sage: assert golden_ratio._sympy_() == S.GoldenRatio
        sage: assert golden_ratio == S.GoldenRatio._sage_()
    """
    from sage.symbolic.constants import golden_ratio
    return golden_ratio

def _sympysage_eulerg(self):
    """
    EXAMPLES::

        sage: from sympy.core.singleton import S
        sage: assert euler_gamma._sympy_() == S.EulerGamma
        sage: assert euler_gamma == S.EulerGamma._sage_()
    """
    from sage.symbolic.constants import euler_gamma
    return euler_gamma

def _sympysage_catalan(self):
    """
    EXAMPLES::

        sage: from sympy.core.singleton import S
        sage: assert catalan._sympy_() == S.Catalan
        sage: assert catalan == S.Catalan._sage_()
    """
    from sage.symbolic.constants import catalan
    return catalan

def _sympysage_i(self):
    """
    EXAMPLES::

        sage: from sympy.core.singleton import S
        sage: assert I._sympy_() == S.ImaginaryUnit
        sage: assert I == S.ImaginaryUnit._sage_()
    """
    from sage.symbolic.constants import I
    return I

##################       basic operators         ##############

def _sympysage_add(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol
        sage: from sympy.core.singleton import S
        sage: assert (x-pi+1)._sympy_() == Symbol('x')-S.Pi+1
        sage: assert x-pi+1 == (Symbol('x')-S.Pi+1)._sage_()
    """
    s = 0
    for x in self.args:
        s += x._sage_()
    return s

def _sympysage_mul(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol
        sage: from sympy.core.singleton import S
        sage: assert (-x*pi*5)._sympy_() == -Symbol('x')*S.Pi*5
        sage: assert -x*pi*5 == (-Symbol('x')*S.Pi*5)._sage_()
    """
    s = 1
    for x in self.args:
        s *= x._sage_()
    return s

def _sympysage_pow(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol
        sage: from sympy.core.singleton import S
        sage: assert (x^pi^5)._sympy_() == Symbol('x')**S.Pi**5
        sage: assert x^pi^5 == (Symbol('x')**S.Pi**5)._sage_()
    """
    return self.args[0]._sage_()**self.args[1]._sage_()

def _sympysage_symbol(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol
        sage: assert x._sympy_() == Symbol('x')
        sage: assert x == Symbol('x')._sage_()
    """
    from sage.symbolic.ring import SR
    return SR.var(self.name)

##############       functions       ###############

def _sympysage_function(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, Function, sin as Sin
        sage: assert sin(x)._sympy_() == Sin(Symbol('x'))
        sage: assert sin(x) == Sin(Symbol('x'))._sage_()

        sage: f = function('f')
        sage: F = Function('f')
        sage: assert f(x)._sympy_() == F(x)
        sage: assert f(x) == F(x)._sage_()

    Test that functions unknown to Sage raise an exception::

        sage: from sympy.functions.combinatorial.numbers import lucas
        sage: lucas(Symbol('x'))._sage_()
        Traceback (most recent call last):
        ...
        AttributeError...
        """
    from sage.functions import all as sagefuncs
    fname = self.func.__name__
    func = getattr(sagefuncs, fname, None)
    args = [arg._sage_() for arg in self.args]

    # In the case the function is not known in sage:
    if func is None:
        import sympy
        if getattr(sympy, fname, None) is None:
            # abstract function
            from sage.calculus.var import function
            return function(fname)(*args)

        else:
            # the function defined in sympy is not known in sage
            raise AttributeError

    return func(*args)

def _sympysage_integral(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, Integral
        sage: sx = Symbol('x')
        sage: assert integral(x, x, hold=True)._sympy_() == Integral(sx, sx)
        sage: assert integral(x, x, hold=True) == Integral(sx, sx)._sage_()
        sage: assert integral(x, x, 0, 1, hold=True)._sympy_() == Integral(sx, (sx,int(0),int(1)))
        sage: assert integral(x, x, 0, 1, hold=True) == Integral(sx, (sx,0,1))._sage_()
    """
    from sage.misc.functional import integral
    f, limits = self.function._sage_(), list(self.limits)
    for limit in limits:
        if len(limit) == 1:
            x = limit[0]
            f = integral(f, x._sage_(), hold=True)
        elif len(limit) == 2:
            x, b = limit
            f = integral(f, x._sage_(), b._sage_(), hold=True)
        else:
            x, a, b = limit
            f = integral(f, (x._sage_(), a._sage_(), b._sage_()), hold=True)
    return f

def _sympysage_derivative(self):
    from sage.calculus.functional import derivative
    args = [arg._sage_() for arg in self.args]
    return derivative(*args)

def _sympysage_order(self):
    from sage.rings.integer import Integer
    return Integer(0)

def _sympysage_rf(self):
    from sage.arith.all import rising_factorial
    return rising_factorial(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_ff(self):
    from sage.arith.all import falling_factorial
    return falling_factorial(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_lgamma(self):
    from sage.functions.other import log_gamma
    return log_gamma(self.args[0]._sage_())

def _sympysage_dirac_delta(self):
    from sage.functions.generalized import dirac_delta
    return dirac_delta(self.args[0]._sage_())

def _sympysage_heaviside(self):
    from sage.functions.generalized import heaviside
    return heaviside(self.args[0]._sage_())

def _sympysage_expint(self):
    from sage.functions.exp_integral import exp_integral_e
    return exp_integral_e(self.args[0]._sage_(), self.args[1]._sage_())

# are si,ci,shi,chi really necessary?
def _sympysage_si(self):
    from sage.functions.exp_integral import sin_integral
    return sin_integral(self.args[0]._sage_())

def _sympysage_ci(self):
    from sage.functions.exp_integral import cos_integral
    return cos_integral(self.args[0]._sage_())

def _sympysage_shi(self):
    from sage.functions.exp_integral import sinh_integral
    return sinh_integral(self.args[0]._sage_())

def _sympysage_chi(self):
    from sage.functions.exp_integral import cosh_integral
    return cosh_integral(self.args[0]._sage_())

def _sympysage_hyp(self):
    from sage.functions.hypergeometric import hypergeometric
    ap = [arg._sage_() for arg in self.args[0]]
    bq = [arg._sage_() for arg in self.args[1]]
    return hypergeometric(ap, bq, self.argument._sage_())

def _sympysage_elliptic_k(self):
    from sage.functions.special import elliptic_kc
    return elliptic_kc(self.args[0]._sage_())

def _sympysage_kronecker_delta(self):
    from sage.functions.generalized import kronecker_delta
    return kronecker_delta(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_besselj(self):
    from sage.functions.bessel import bessel_J
    return bessel_J(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_bessely(self):
    from sage.functions.bessel import bessel_Y
    return bessel_Y(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_besseli(self):
    from sage.functions.bessel import bessel_I
    return bessel_I(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_besselk(self):
    from sage.functions.bessel import bessel_K
    return bessel_K(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_ynm(self):
    from sage.functions.special import spherical_harmonic
    return spherical_harmonic(self.args[0]._sage_(),
                              self.args[1]._sage_(),
                              self.args[2]._sage_(),
                              self.args[3]._sage_())

#really necessary?
def _sympysage_csch(self):
    from sage.functions.hyperbolic import csch
    return csch(self.args[0]._sage_())

#really necessary?
def _sympysage_sech(self):
    from sage.functions.hyperbolic import sech
    return sech(self.args[0]._sage_())

def _sympysage_re(self):
    from sage.functions.other import real_part
    return real_part(self.args[0]._sage_())

def _sympysage_im(self):
    from sage.functions.other import imag_part
    return imag_part(self.args[0]._sage_())

#really necessary?
def _sympysage_sign(self):
    from sage.functions.generalized import sign
    return sign(self.args[0]._sage_())

def _sympysage_abs(self):
    from sage.functions.generalized import abs_symbolic
    return abs_symbolic(self.args[0]._sage_())


#------------------------------------------------------------------
from sage.repl.ipython_extension import run_once

@run_once
def sympy_init():
    from sympy import Add, Mul, Pow, Symbol
    from sympy.core.function import (Function, Derivative)
    from sympy.core.numbers import (Float, Integer, Rational, Infinity,
            NegativeInfinity, ComplexInfinity, Exp1, Pi, GoldenRatio,
            EulerGamma, Catalan, ImaginaryUnit)
    from sympy.core.numbers import NaN as sympy_nan
    from sympy.functions.combinatorial.factorials import (RisingFactorial,
            FallingFactorial)
    from sympy.functions.elementary.complexes import (re, im, sign, Abs)
    from sympy.functions.elementary.hyperbolic import (csch, sech)
    from sympy.functions.special.bessel import (besselj, bessely, besseli, besselk)
    from sympy.functions.special.delta_functions import (DiracDelta, Heaviside)
    from sympy.functions.special.error_functions import (expint, Si, Ci, Shi, Chi)
    from sympy.functions.special.elliptic_integrals import elliptic_k
    from sympy.functions.special.gamma_functions import loggamma
    from sympy.functions.special.hyper import hyper
    from sympy.functions.special.spherical_harmonics import Ynm
    from sympy.functions.special.tensor_functions import KroneckerDelta
    from sympy.integrals.integrals import Integral
    from sympy.series.order import Order

    Float._sage_ = _sympysage_float
    Rational._sage_ = _sympysage_rational
    Infinity._sage_ = _sympysage_pinfty
    NegativeInfinity._sage_ = _sympysage_ninfty
    ComplexInfinity._sage_ = _sympysage_uinfty
    sympy_nan._sage_ = _sympysage_nan
    Exp1._sage_ = _sympysage_e
    Pi._sage_ = _sympysage_pi
    GoldenRatio._sage_ = _sympysage_golden_ratio
    EulerGamma._sage_ = _sympysage_eulerg
    Catalan._sage_ = _sympysage_catalan
    ImaginaryUnit._sage_ = _sympysage_i
    Add._sage_ = _sympysage_add
    Mul._sage_ = _sympysage_mul
    Pow._sage_ = _sympysage_pow
    Symbol._sage_ = _sympysage_symbol
    Function._sage_ = _sympysage_function
    Integral._sage_ = _sympysage_integral
    Derivative._sage_ = _sympysage_derivative
    Order._sage_ = _sympysage_order
    RisingFactorial._sage_ = _sympysage_rf
    FallingFactorial._sage_ = _sympysage_ff
    loggamma._sage_ = _sympysage_lgamma
    DiracDelta._sage_ = _sympysage_dirac_delta
    Heaviside._sage_ = _sympysage_heaviside
    expint._sage_ = _sympysage_expint
    Si._sage_ = _sympysage_si
    Ci._sage_ = _sympysage_ci
    Shi._sage_ = _sympysage_shi
    Chi._sage_ = _sympysage_chi
    hyper._sage_ = _sympysage_hyp
    elliptic_k._sage_ = _sympysage_elliptic_k
    KroneckerDelta._sage_ = _sympysage_kronecker_delta
    besselj._sage_ = _sympysage_besselj
    bessely._sage_ = _sympysage_bessely
    besseli._sage_ = _sympysage_besseli
    besselk._sage_ = _sympysage_besselk
    Ynm._sage_ = _sympysage_ynm
    csch._sage_ = _sympysage_csch
    sech._sage_ = _sympysage_sech
    re._sage_ = _sympysage_re
    im._sage_ = _sympysage_im
    sign._sage_ = _sympysage_sign
    Abs._sage_ = _sympysage_sign

def check_expression(expr, var_symbols, only_from_sympy=False):
    """
    Does eval(expr) both in Sage and SymPy and does other checks.
    """

    # evaluate the expression in the context of Sage:
    if var_symbols:
        sage.var(var_symbols)
    a = globals().copy()
    # safety checks...
    a.update(sage.__dict__)
    assert "sin" in a
    is_different = False
    try:
        e_sage = eval(expr, a)
        assert not isinstance(e_sage, sympy.Basic)
    except (NameError, TypeError):
        is_different = True
        pass

    # evaluate the expression in the context of SymPy:
    if var_symbols:
        sympy_vars = sympy.var(var_symbols)
    b = globals().copy()
    b.update(sympy.__dict__)
    assert "sin" in b
    b.update(sympy.__dict__)
    e_sympy = eval(expr, b)
    assert isinstance(e_sympy, sympy.Basic)

    # Sympy func may have specific _sage_ method
    if is_different:
        _sage_method = getattr(e_sympy.func, "_sage_")
        e_sage = _sage_method(sympy.S(e_sympy))

    # Do the actual checks:
    if not only_from_sympy:
        assert sympy.S(e_sage) == e_sympy
    assert e_sage == sage.SR(e_sympy)


def test_basics():
    check_expression("x", "x")
    check_expression("x**2", "x")
    check_expression("x**2+y**3", "x y")
    check_expression("1/(x+y)**2-x**3/4", "x y")


def test_complex():
    check_expression("I", "")
    check_expression("23+I*4", "x")


def test_complex_fail():
    # Sage doesn't properly implement _sympy_ on I
    check_expression("I*y", "y")
    check_expression("x+I*y", "x y")


def test_integer():
    check_expression("4*x", "x")
    check_expression("-4*x", "x")


def test_real():
    check_expression("1.123*x", "x")
    check_expression("-18.22*x", "x")



def test_functions():
    # Test at least one Function without own _sage_ method
    assert not "_sage_" in sympy.factorial.__dict__
    check_expression("factorial(x)", "x")
    check_expression("sin(x)", "x")
    check_expression("cos(x)", "x")
    check_expression("tan(x)", "x")
    check_expression("cot(x)", "x")
    check_expression("asin(x)", "x")
    check_expression("acos(x)", "x")
    check_expression("atan(x)", "x")
    check_expression("atan2(y, x)", "x, y")
    check_expression("acot(x)", "x")
    check_expression("sinh(x)", "x")
    check_expression("cosh(x)", "x")
    check_expression("tanh(x)", "x")
    check_expression("coth(x)", "x")
    check_expression("asinh(x)", "x")
    check_expression("acosh(x)", "x")
    check_expression("atanh(x)", "x")
    check_expression("acoth(x)", "x")
    check_expression("exp(x)", "x")
    check_expression("log(x)", "x")
    check_expression("re(x)", "x")
    check_expression("im(x)", "x")
    check_expression("sign(x)", "x")
    check_expression("abs(x)", "x")
    check_expression("arg(x)", "x")
    check_expression("conjugate(x)", "x")

    # The following tests differently named functions
    check_expression("besselj(y, x)", "x, y")
    check_expression("bessely(y, x)", "x, y")
    check_expression("besseli(y, x)", "x, y")
    check_expression("besselk(y, x)", "x, y")
    check_expression("DiracDelta(x)", "x")
    check_expression("KroneckerDelta(x, y)", "x, y")
    check_expression("expint(y, x)", "x, y")
    check_expression("Si(x)", "x")
    check_expression("Ci(x)", "x")
    check_expression("Shi(x)", "x")
    check_expression("Chi(x)", "x")
    check_expression("loggamma(x)", "x")
    check_expression("Ynm(n,m,x,y)", "n, m, x, y")
    check_expression("hyper((n,m),(m,n),x)", "n, m, x")

def test_issue_4023():
    sage.var("a x")
    log = sage.log
    i = sympy.integrate(log(x)/a, (x, a, a + 1))
    i2 = sympy.simplify(i)
    s = sage.SR(i2)
    assert s == (a*log(1 + a) - a*log(a) + log(1 + a) - 1)/a

def test_integral():
    #test Sympy-->Sage
    check_expression("Integral(x, (x,))", "x", only_from_sympy=True)
    check_expression("Integral(x, (x, 0, 1))", "x", only_from_sympy=True)
    check_expression("Integral(x*y, (x,), (y, ))", "x,y", only_from_sympy=True)
    check_expression("Integral(x*y, (x,), (y, 0, 1))", "x,y", only_from_sympy=True)
    check_expression("Integral(x*y, (x, 0, 1), (y,))", "x,y", only_from_sympy=True)
    check_expression("Integral(x*y, (x, 0, 1), (y, 0, 1))", "x,y", only_from_sympy=True)
    check_expression("Integral(x*y*z, (x, 0, 1), (y, 0, 1), (z, 0, 1))", "x,y,z", only_from_sympy=True)

def test_integral_failing():
    # Note: sage may attempt to turn this into Integral(x, (x, x, 0))
    check_expression("Integral(x, (x, 0))", "x", only_from_sympy=True)
    check_expression("Integral(x*y, (x,), (y, 0))", "x,y", only_from_sympy=True)
    check_expression("Integral(x*y, (x, 0, 1), (y, 0))", "x,y", only_from_sympy=True)

def test_undefined_function():
    f = sympy.Function('f')
    sf = sage.function('f')
    x = sympy.symbols('x')
    sx = sage.var('x')
    assert bool(sf(sx) == f(x)._sage_())
    #assert bool(f == sympy.sympify(sf))

def test_abstract_function():
    from sage.symbolic.expression import Expression
    x,y = sympy.symbols('x y')
    f = sympy.Function('f')
    expr =  f(x,y)
    sexpr = expr._sage_()
    assert isinstance(sexpr,Expression), "converted expression %r is not sage expression" % sexpr
    # This test has to be uncommented in the future: it depends on the sage ticket #22802 (https://trac.sagemath.org/ticket/22802)
    # invexpr = sexpr._sympy_()
    # assert invexpr == expr, "inverse coversion %r is not correct " % invexpr

def test_relational():
    x = sympy.symbols('x')
    sx = sage.var('x')
    assert sympy.sympify(sx == 0) == sympy.Eq(x, 0)
    assert (sx == 0) == sage.SR(sympy.Eq(x, 0))
    assert sympy.sympify(sx != 0) == sympy.Ne(x, 0)
    assert (sx != 0) == sage.SR(sympy.Ne(x, 0))
    assert sympy.sympify(sx >= 0) == sympy.Ge(x, 0)
    assert (sx >= 0) == sage.SR(sympy.Ge(x, 0))
    assert sympy.sympify(sx <= 0) == sympy.Le(x, 0)
    assert (sx <= 0) == sage.SR(sympy.Le(x, 0))
    assert sympy.sympify(sx > 0) == sympy.Gt(x, 0)
    assert (sx > 0) == sage.SR(sympy.Gt(x, 0))
    assert sympy.sympify(sx < 0) == sympy.Lt(x, 0)
    assert (sx < 0) == sage.SR(sympy.Lt(x, 0))


