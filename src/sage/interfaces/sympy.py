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
    from sage.rings.real_mpfr import RealNumber
    return RealNumber(str(self))

def _sympysage_rational(self):
    from sage.rings.integer import Integer
    from sage.rings.rational import Rational
    return Rational((Integer(self.p), Integer(self.q)))

def _sympysage_pinfty(self):
    from sage.rings.infinity import PlusInfinity
    return PlusInfinity()

def _sympysage_ninfty(self):
    from sage.rings.infinity import MinusInfinity
    return MinusInfinity()

def _sympysage_uinfty(self):
    from sage.rings.infinity import unsigned_infinity
    return unsigned_infinity

def _sympysage_nan(self):
    from sage.symbolic.constants import NaN
    return NaN

def _sympysage_e(self):
    from sage.symbolic.constants import e
    return e

def _sympysage_pi(self):
    from sage.symbolic.constants import pi
    return pi

def _sympysage_golden_ratio(self):
    from sage.symbolic.constants import golden_ratio
    return golden_ratio

def _sympysage_eulerg(self):
    from sage.symbolic.constants import euler_gamma
    return euler_gamma

def _sympysage_catalan(self):
    from sage.symbolic.constants import catalan
    return catalan

def _sympysage_i(self):
    from sage.symbolic.constants import I
    return I

##################       basic operators         ##############

def _sympysage_add(self):
    s = 0
    for x in self.args:
        s += x._sage_()
    return s

def _sympysage_mul(self):
    s = 1
    for x in self.args:
        s *= x._sage_()
    return s

def _sympysage_pow(self):
    return self.args[0]._sage_()**self.args[1]._sage_()

def _sympysage_symbol(self):
    from sage.symbolic.ring import SR
    return SR.var(self.name)

##############       functions       ###############

def _sympysage_function(self):
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
            # this exception is catched in sage
            raise AttributeError

    return func(*args)

def _sympysage_integral(self):
    from sage.misc.functional import integral
    f, limits = self.function._sage_(), list(self.limits)
    for limit in limits:
        if len(limit) == 1:
            x = limit[0]
            f = sintegral(f, x._sage_(), hold=True)
        elif len(limit) == 2:
            x, b = limit
            f = sintegral(f, x._sage_(), b._sage_(), hold=True)
        else:
            x, a, b = limit
            f = sintegral(f, (x._sage_(), a._sage_(), b._sage_()), hold=True)
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

    print('xy')

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
