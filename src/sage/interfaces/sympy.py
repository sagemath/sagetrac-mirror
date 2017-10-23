"""
SymPy --> Sage conversion

The file consists of ``_sage_()`` methods that are added lazily to
the respective SymPy objects. Any call of the ``_sympy_()`` method
of a symbolic expression will trigger the addition. See
:class:`sage.symbolic.expression_conversion.SymPyConverter` for the
conversion to SymPy.

Only ``Function`` objects where the names differ need their own ``_sage()_``
method. There are several functions with differing name that have an alias
in Sage that is the same as the name in SymPy, so no explicit translation
is needed for them::

    sage: from sympy import Symbol, Si, Ci, Shi, Chi, sign
    sage: sx = Symbol('x')
    sage: assert sin_integral(x)._sympy_() == Si(sx)
    sage: assert sin_integral(x) == Si(sx)._sage_()
    sage: assert sinh_integral(x)._sympy_() == Shi(sx)
    sage: assert sinh_integral(x) == Shi(sx)._sage_()
    sage: assert cos_integral(x)._sympy_() == Ci(sx)
    sage: assert cos_integral(x) == Ci(sx)._sage_()
    sage: assert cosh_integral(x)._sympy_() == Chi(sx)
    sage: assert cosh_integral(x) == Chi(sx)._sage_()
    sage: assert sgn(x)._sympy_() == sign(sx)
    sage: assert sgn(x) == sign(sx)._sage_()

AUTHORS:

- Ralf Stephan (2017-10)
"""
################################################################
#   Distributed under GNU GPL3, see www.gnu.org
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

def _sympysage_Subs(self):
     """
     EXAMPLES::

         sage: from sympy import Symbol
         sage: from sympy.core.singleton import S
     """

     args = self.args
     substi = dict([(args[1][i]._sage_(),args[2][i]._sage_()) for i in range(len(args[1]))])

     return args[0]._sage_().subs(substi)


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
        sage: assert f(x+3)._sympy_() == F(x+3)
        sage: assert f(x+3) == F(x+3)._sage_()
        sage: assert (3*f(x))._sympy_() == 3*F(x)
        sage: assert 3*f(x) == (3*F(x))._sage_()

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
        sage: assert integral(x, x, 0, 1, hold=True)._sympy_() == Integral(sx, (sx,0,1)) # known bug
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
    """
    EXAMPLES::

        sage: from sympy import Derivative
        sage: f = function('f')
        sage: sympy_diff = Derivative(f(x)._sympy_(), x._sympy_())
        sage: assert diff(f(x),x)._sympy_() == sympy_diff
        sage: assert diff(f(x),x) == sympy_diff._sage_()
    """
    from sage.calculus.functional import derivative
    args = [arg._sage_() for arg in self.args]
    return derivative(*args)

def _sympysage_order(self):
    """
    EXAMPLES::

        sage: from sage.functions.other import Order
        sage: from sympy.series import Order as SOrder
        sage: assert Order(1)._sympy_() == SOrder(1)
        sage: assert Order(1) == SOrder(1)._sage_()
    """
    from sage.functions.other import Order
    return Order(self.args[0])._sage_()

def _sympysage_rf(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, rf
        sage: _ = var('x, y')
        sage: rfxy = rf(Symbol('x'), Symbol('y'))
        sage: assert rising_factorial(x,y)._sympy_() == rfxy.rewrite('gamma')
        sage: assert rising_factorial(x,y) == rfxy._sage_()
    """
    from sage.arith.all import rising_factorial
    return rising_factorial(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_ff(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, ff
        sage: _ = var('x, y')
        sage: ffxy = ff(Symbol('x'), Symbol('y'))
        sage: assert falling_factorial(x,y)._sympy_() == ffxy.rewrite('gamma') # known bug
        sage: assert falling_factorial(x,y) == ffxy._sage_()
    """
    from sage.arith.all import falling_factorial
    return falling_factorial(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_lgamma(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, loggamma
        sage: assert log_gamma(x)._sympy_() == loggamma(Symbol('x'))
        sage: assert log_gamma(x) == loggamma(Symbol('x'))._sage_()
    """
    from sage.functions.other import log_gamma
    return log_gamma(self.args[0]._sage_())

def _sympysage_dirac_delta(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, DiracDelta
        sage: assert dirac_delta(x)._sympy_() == DiracDelta(Symbol('x'))
        sage: assert dirac_delta(x) == DiracDelta(Symbol('x'))._sage_()
    """
    from sage.functions.generalized import dirac_delta
    return dirac_delta(self.args[0]._sage_())

def _sympysage_heaviside(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, Heaviside
        sage: assert heaviside(x)._sympy_() == Heaviside(Symbol('x'))
        sage: assert heaviside(x) == Heaviside(Symbol('x'))._sage_()
    """
    from sage.functions.generalized import heaviside
    return heaviside(self.args[0]._sage_())

def _sympysage_expint(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, expint
        sage: _ = var('x, y')
        sage: sy = expint(Symbol('x'), Symbol('y'))
        sage: assert exp_integral_e(x,y)._sympy_() == sy
        sage: assert exp_integral_e(x,y) == sy._sage_()
    """
    from sage.functions.exp_integral import exp_integral_e
    return exp_integral_e(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_hyp(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, hyper
        sage: _ = var('a,b,p,q,x')
        sage: sy = hyper((Symbol('a'), Symbol('b')), (Symbol('p'), Symbol('q')), Symbol('x'))
        sage: assert hypergeometric((a,b),(p,q),x)._sympy_() == sy
        sage: assert hypergeometric((a,b),(p,q),x) == sy._sage_()
    """
    from sage.functions.hypergeometric import hypergeometric
    ap = [arg._sage_() for arg in self.args[0]]
    bq = [arg._sage_() for arg in self.args[1]]
    return hypergeometric(ap, bq, self.argument._sage_())

def _sympysage_elliptic_k(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, elliptic_k
        sage: assert elliptic_kc(x)._sympy_() == elliptic_k(Symbol('x'))
        sage: assert elliptic_kc(x) == elliptic_k(Symbol('x'))._sage_()
    """
    from sage.functions.special import elliptic_kc
    return elliptic_kc(self.args[0]._sage_())

def _sympysage_kronecker_delta(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, KroneckerDelta
        sage: _ = var('x, y')
        sage: sy = KroneckerDelta(Symbol('x'), Symbol('y'))
        sage: assert kronecker_delta(x,y)._sympy_() == sy
        sage: assert kronecker_delta(x,y) == sy._sage_()
    """
    from sage.functions.generalized import kronecker_delta
    return kronecker_delta(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_ceiling(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, ceiling
        sage: assert ceil(x)._sympy_() == ceiling(Symbol('x'))
        sage: assert ceil(x) == ceiling(Symbol('x'))._sage_()
        sage: integrate(ceil(x), x, 0, infinity, algorithm='sympy')
        integrate(ceil(x), x, 0, +Infinity)
    """
    from sage.functions.other import ceil
    return ceil(self.args[0]._sage_())

def _sympysage_piecewise(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, pi as spi, Eq, Lt, Piecewise
        sage: sx = Symbol('x')
        sage: sp = Piecewise((spi, Lt(sx,0)), (1, Eq(sx,1)), (0, True))
        sage: ex = cases(((x<0, pi), (x==1, 1), (True, 0)))
        sage: assert ex._sympy_() == sp
        sage: assert ex == sp._sage_()

        sage: _ = var('y, z')
        sage: (x^y - z).integrate(y, algorithm="sympy")
        -y*z + cases(((log(x) == 0, y), (1, x^y/log(x))))
    """
    from sage.functions.other import cases
    return cases([(p.cond._sage_(),p.expr._sage_()) for p in self.args])

def _sympysage_besselj(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, besselj
        sage: _ = var('x, y')
        sage: sy = besselj(Symbol('x'), Symbol('y'))
        sage: assert bessel_J(x,y)._sympy_() == sy
        sage: assert bessel_J(x,y) == sy._sage_()
    """
    from sage.functions.bessel import bessel_J
    return bessel_J(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_bessely(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, bessely
        sage: _ = var('x, y')
        sage: sy = bessely(Symbol('x'), Symbol('y'))
        sage: assert bessel_Y(x,y)._sympy_() == sy
        sage: assert bessel_Y(x,y) == sy._sage_()
    """
    from sage.functions.bessel import bessel_Y
    return bessel_Y(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_besseli(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, besseli
        sage: _ = var('x, y')
        sage: sy = besseli(Symbol('x'), Symbol('y'))
        sage: assert bessel_I(x,y)._sympy_() == sy
        sage: assert bessel_I(x,y) == sy._sage_()
    """
    from sage.functions.bessel import bessel_I
    return bessel_I(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_besselk(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, besselk
        sage: _ = var('x, y')
        sage: sy = besselk(Symbol('x'), Symbol('y'))
        sage: assert bessel_K(x,y)._sympy_() == sy
        sage: assert bessel_K(x,y) == sy._sage_()
    """
    from sage.functions.bessel import bessel_K
    return bessel_K(self.args[0]._sage_(), self.args[1]._sage_())

def _sympysage_ynm(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, Ynm
        sage: _ = var('n,m,t,p')
        sage: sy = Ynm(Symbol('n'), Symbol('m'), Symbol('t'), Symbol('p'))
        sage: assert spherical_harmonic(n,m,t,p)._sympy_() == sy
        sage: assert spherical_harmonic(n,m,t,p) == sy._sage_()
    """
    from sage.functions.special import spherical_harmonic
    return spherical_harmonic(self.args[0]._sage_(),
                              self.args[1]._sage_(),
                              self.args[2]._sage_(),
                              self.args[3]._sage_())

def _sympysage_re(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, re
        sage: assert real_part(x)._sympy_() == re(Symbol('x'))
        sage: assert real_part(x) == re(Symbol('x'))._sage_()
    """
    from sage.functions.other import real_part
    return real_part(self.args[0]._sage_())

def _sympysage_im(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, im
        sage: assert imag_part(x)._sympy_() == im(Symbol('x'))
        sage: assert imag_part(x) == im(Symbol('x'))._sage_()
    """
    from sage.functions.other import imag_part
    return imag_part(self.args[0]._sage_())

def _sympysage_abs(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, Abs
        sage: assert abs(x)._sympy_() == Abs(Symbol('x'))
        sage: assert abs(x) == Abs(Symbol('x'))._sage_()
    """
    from sage.functions.other import abs_symbolic
    return abs_symbolic(self.args[0]._sage_())

def _sympysage_relational(self):
    """
    EXAMPLES::

        sage: from sympy import Eq, Ne, Gt, Ge, Lt, Le, Symbol
        sage: sx = Symbol('x')
        sage: assert (x == 0)._sympy_() == Eq(sx, 0)
        sage: assert (x == 0) == Eq(x, 0)._sage_()
        sage: assert (x != 0)._sympy_() == Ne(sx, 0)
        sage: assert (x != 0) == Ne(x, 0)._sage_()
        sage: assert (x > 0)._sympy_() == Gt(sx, 0)
        sage: assert (x > 0) == Gt(x, 0)._sage_()
        sage: assert (x >= 0)._sympy_() == Ge(sx, 0)
        sage: assert (x >= 0) == Ge(x, 0)._sage_()
        sage: assert (x < 0)._sympy_() == Lt(sx, 0)
        sage: assert (x < 0) == Lt(x, 0)._sage_()
        sage: assert (x <= 0)._sympy_() == Le(sx, 0)
        sage: assert (x <= 0) == Le(x, 0)._sage_()
     """
    from operator import eq, ne, gt, lt, ge, le
    from sympy import Eq, Ne, Gt, Ge, Lt, Le
    ops = {Eq : eq, Ne : ne, Gt : gt, Lt : lt, Ge : ge, Le : le}
    return ops.get(self.func)(self.lhs._sage_(), self.rhs._sage_())

def _sympysage_false(self):
    """
    EXAMPLES::

        sage: from sympy.logic.boolalg import BooleanFalse
        sage: assert SR(False)._sympy_() == BooleanFalse()      # known bug
        sage: assert SR(False) == BooleanFalse()._sage_()
    """
    from sage.symbolic.ring import SR
    return SR(False)

def _sympysage_true(self):
    """
    EXAMPLES::

        sage: from sympy.logic.boolalg import BooleanTrue
        sage: assert SR(True)._sympy_() == BooleanTrue()      # known bug
        sage: assert SR(True) == BooleanTrue()._sage_()
    """
    from sage.symbolic.ring import SR
    return SR(True)


#------------------------------------------------------------------
from sage.repl.ipython_extension import run_once

@run_once
def sympy_init():
    """
    Add ``_sage_()`` methods to SymPy objects where needed.

    This gets called with every call to ``Expression._sympy_()``
    so there is only need to call it if you bypass ``_sympy_()`` to
    create SymPy objects. Note that SymPy objects have ``_sage_()``
    methods hard installed but having them inside Sage as
    one file makes them easier to maintain for Sage developers.

    EXAMPLES::

        sage: from sage.interfaces.sympy import sympy_init
        sage: from sympy import Symbol, Abs
        sage: sympy_init()
        sage: assert abs(x) == Abs(Symbol('x'))._sage_()
    """
    from sympy import Add
    if Add._sage_ == _sympysage_add:
        return

    from sympy import Mul, Pow, Symbol, Subs
    from sympy.core.function import (Function, AppliedUndef, Derivative)
    from sympy.core.numbers import (Float, Integer, Rational, Infinity,
            NegativeInfinity, ComplexInfinity, Exp1, Pi, GoldenRatio,
            EulerGamma, Catalan, ImaginaryUnit)
    from sympy.core.numbers import NaN as sympy_nan
    from sympy.core.relational import Relational
    from sympy.functions.combinatorial.factorials import (RisingFactorial,
            FallingFactorial)
    from sympy.functions.elementary.complexes import (re, im, Abs)
    from sympy.functions.elementary.integers import ceiling
    from sympy.functions.elementary.piecewise import Piecewise
    from sympy.functions.special.bessel import (besselj, bessely, besseli, besselk)
    from sympy.functions.special.delta_functions import (DiracDelta, Heaviside)
    from sympy.functions.special.error_functions import expint
    from sympy.functions.special.elliptic_integrals import elliptic_k
    from sympy.functions.special.gamma_functions import loggamma
    from sympy.functions.special.hyper import hyper
    from sympy.functions.special.spherical_harmonics import Ynm
    from sympy.functions.special.tensor_functions import KroneckerDelta
    from sympy.logic.boolalg import BooleanTrue, BooleanFalse
    from sympy.integrals.integrals import Integral
    from sympy.series.order import Order

    Float._sage_ = _sympysage_float
    Rational._sage_ = _sympysage_rational
    Infinity._sage_ = _sympysage_pinfty
    NegativeInfinity._sage_ = _sympysage_ninfty
    ComplexInfinity._sage_ = _sympysage_uinfty
    sympy_nan._sage_ = _sympysage_nan
    Relational._sage_ = _sympysage_relational
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
    Subs._sage_ = _sympysage_Subs
    Function._sage_ = _sympysage_function
    AppliedUndef._sage_ = _sympysage_function
    Integral._sage_ = _sympysage_integral
    Derivative._sage_ = _sympysage_derivative
    Order._sage_ = _sympysage_order
    RisingFactorial._sage_ = _sympysage_rf
    FallingFactorial._sage_ = _sympysage_ff
    loggamma._sage_ = _sympysage_lgamma
    DiracDelta._sage_ = _sympysage_dirac_delta
    Heaviside._sage_ = _sympysage_heaviside
    expint._sage_ = _sympysage_expint
    hyper._sage_ = _sympysage_hyp
    elliptic_k._sage_ = _sympysage_elliptic_k
    KroneckerDelta._sage_ = _sympysage_kronecker_delta
    Piecewise._sage_ = _sympysage_piecewise
    besselj._sage_ = _sympysage_besselj
    bessely._sage_ = _sympysage_bessely
    besseli._sage_ = _sympysage_besseli
    besselk._sage_ = _sympysage_besselk
    Ynm._sage_ = _sympysage_ynm
    re._sage_ = _sympysage_re
    im._sage_ = _sympysage_im
    Abs._sage_ = _sympysage_abs
    BooleanFalse._sage_ = _sympysage_false
    BooleanTrue._sage_ = _sympysage_true
    ceiling._sage_ = _sympysage_ceiling

def check_expression(expr, var_symbols, only_from_sympy=False):
    """
    Does ``eval(expr)`` both in Sage and SymPy and does other checks.

    EXAMPLES::

        sage: from sage.interfaces.sympy import check_expression
        sage: check_expression("1.123*x", "x")
    """
    from sage import __dict__ as sagedict
    from sage.symbolic.ring import SR
    from sympy import (__dict__ as sympydict, Basic, S, var as svar)
    # evaluate the expression in the context of Sage:
    if var_symbols:
        SR.var(var_symbols)
    is_different = False
    try:
        e_sage = SR(expr)
        assert not isinstance(e_sage, Basic)
    except (NameError, TypeError):
        is_different = True
        pass

    # evaluate the expression in the context of SymPy:
    if var_symbols:
        sympy_vars = svar(var_symbols)
    b = globals().copy()
    b.update(sympydict)
    assert "sin" in b
    b.update(sympydict)
    e_sympy = eval(expr, b)
    assert isinstance(e_sympy, Basic)

    # Sympy func may have specific _sage_ method
    if is_different:
        _sage_method = getattr(e_sympy.func, "_sage_")
        e_sage = _sage_method(S(e_sympy))

    # Do the actual checks:
    if not only_from_sympy:
        assert S(e_sage) == e_sympy
    assert e_sage == SR(e_sympy)

def test_all():
    """
    Call some tests that were originally in SymPy.

    EXAMPLES::

        sage: from sage.interfaces.sympy import test_all
        sage: test_all()
    """
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
        from sympy import factorial
        assert not "_sage_" in factorial.__dict__
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
        check_expression("abs(x)", "x")
        check_expression("arg(x)", "x")
        check_expression("conjugate(x)", "x")

    def test_issue_4023():
        from sage.symbolic.ring import SR
        from sage.functions.all import log
        from sympy import integrate, simplify
        a,x = SR.var("a x")
        i = integrate(log(x)/a, (x, a, a + 1))
        i2 = simplify(i)
        s = SR(i2)
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
        from sage.symbolic.ring import SR
        from sage.calculus.var import function
        from sympy import Symbol, Function
        f = function('f')
        sf = Function('f')
        x,y = SR.var('x y')
        sx = Symbol('x')
        sy = Symbol('y')
        assert f(x)._sympy_() == sf(sx)
        assert f(x) == sf(sx)._sage_()
        assert f(x,y)._sympy_() == sf(sx, sy)
        assert f(x,y) == sf(sx, sy)._sage_()
        assert f._sympy_() == sf
        assert f == sf._sage_()

    test_basics()
    test_complex()
    test_complex_fail()
    test_integer()
    test_real()
    test_functions()
    test_issue_4023()
    test_integral()
    #test_integral_failing()
    test_undefined_function()
