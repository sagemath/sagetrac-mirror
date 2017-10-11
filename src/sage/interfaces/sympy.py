#   sympy.py - SymPy --> Sage interface
#   Author: (2017) Ralf Stephan <ralf@ark.in-berin.de>
#   Distributed under GNU GPL3, see www.gnu.org
#
#   The file consists of sage() methods that are added to SymPy objects.
#   Only Function objects where the names differ need their own sage() method.
#   Please be very specific with imports to minimalize startup time.
################################################################
from __future__ import absolute_import

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


#################         numbers and constants      ##############

def _sympysage_float(self):
    from sage.rings.real_mpfr import RealNumber
    return RealNumber(str(self))

Float.sage = _sympysage_float

def _sympysage_rational(self):
    from sage.rings.integer import Integer
    from sage.rings.rational import Rational
    return Rational((Integer(self.p), Integer(self.q)))

Rational.sage = _sympysage_rational

def _sympysage_pinfty(self):
    from sage.rings.infinity import oo
    return oo

Infinity.sage = _sympysage_pinfty

def _sympysage_ninfty(self):
    from sage.rings.infinity import oo
    return -oo

NegativeInfinity.sage = _sympysage_ninfty

def _sympysage_uinfty(self):
    from sage.rings.infinity import unsigned_infinity
    return unsigned_infinity

ComplexInfinity.sage = _sympysage_uinfty

def _sympysage_nan(self):
    from sage.symbolic.constants import NaN
    return NaN

sympy_nan.sage = _sympysage_nan

def _sympysage_e(self):
    from sage.symbolic.constants import e
    return e

Exp1.sage = _sympysage_e

def _sympysage_pi(self):
    from sage.symbolic.constants import pi
    return pi

Pi.sage = _sympysage_pi

def _sympysage_golden_ratio(self):
    from sage.symbolic.constants import golden_ratio
    return golden_ratio

GoldenRatio.sage = _sympysage_golden_ratio

def _sympysage_eulerg(self):
    from sage.symbolic.constants import euler_gamma
    return euler_gamma

EulerGamma.sage = _sympysage_eulerg

def _sympysage_catalan(self):
    from sage.symbolic.constants import catalan
    return catalan

Catalan.sage = _sympysage_catalan

def _sympysage_i(self):
    from sage.symbolic.constants import I
    return I

ImaginaryUnit.sage = _sympysage_i

##################       basic operators         ##############

def _sympysage_add(self):
    s = 0
    for x in self.args:
        s += x.sage()
    return s

Add.sage = _sympysage_add

def _sympysage_mul(self):
    s = 1
    for x in self.args:
        s *= x.sage()
    return s

Mul.sage = _sympysage_mul

def _sympysage_pow(self):
    return self.args[0].sage()**self.args[1].sage()

Pow.sage = _sympysage_pow

def _sympysage_symbol(self):
    from sage.symbolic.ring import SR
    return SR.var(self.name)

Symbol.sage = _sympysage_symbol

##############       functions       ###############

def _sympysage_function(self):
    from sage.functions import all as sagefuncs
    fname = self.func.__name__
    func = getattr(sagefuncs, fname, None)
    args = [arg.sage() for arg in self.args]

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

Function.sage = _sympysage_function

def _sympysage_integral(self):
    from sage.misc.functional import integral
    f, limits = self.function.sage(), list(self.limits)
    for limit in limits:
        if len(limit) == 1:
            x = limit[0]
            f = sintegral(f, x.sage(), hold=True)
        elif len(limit) == 2:
            x, b = limit
            f = sintegral(f, x.sage(), b.sage(), hold=True)
        else:
            x, a, b = limit
            f = sintegral(f, (x.sage(), a.sage(), b.sage()), hold=True)
    return f

Integral.sage = _sympysage_integral

def _sympysage_derivative(self):
    from sage.calculus.functional import derivative
    args = [arg.sage() for arg in self.args]
    return derivative(*args)

Derivative.sage = _sympysage_derivative

def _sympysage_order(self):
    from sage.rings.integer import Integer
    return Integer(0)

Order.sage = _sympysage_order

def _sympysage_rf(self):
    from sage.arith.all import rising_factorial
    return rising_factorial(self.args[0].sage(), self.args[1].sage())

RisingFactorial.sage = _sympysage_rf

def _sympysage_ff(self):
    from sage.arith.all import falling_factorial
    return falling_factorial(self.args[0].sage(), self.args[1].sage())

FallingFactorial.sage = _sympysage_ff

def _sympysage_lgamma(self):
    from sage.functions.other import log_gamma
    return log_gamma(self.args[0].sage())

loggamma.sage = _sympysage_lgamma

def _sympysage_dirac_delta(self):
    from sage.functions.generalized import dirac_delta
    return dirac_delta(self.args[0].sage())

DiracDelta.sage = _sympysage_dirac_delta

def _sympysage_heaviside(self):
    from sage.functions.generalized import heaviside
    return heaviside(self.args[0].sage())

Heaviside.sage = _sympysage_heaviside

def _sympysage_expint(self):
    from sage.functions.exp_integral import exp_integral_e
    return exp_integral_e(self.args[0].sage(), self.args[1].sage())

expint.sage = _sympysage_expint

# are si,ci,shi,chi really necessary?
def _sympysage_si(self):
    from sage.functions.exp_integral import sin_integral
    return sin_integral(self.args[0].sage())

Si.sage = _sympysage_si

def _sympysage_ci(self):
    from sage.functions.exp_integral import cos_integral
    return cos_integral(self.args[0].sage())

Ci.sage = _sympysage_ci

def _sympysage_shi(self):
    from sage.functions.exp_integral import sinh_integral
    return sinh_integral(self.args[0].sage())

Shi.sage = _sympysage_shi

def _sympysage_chi(self):
    from sage.functions.exp_integral import cosh_integral
    return cosh_integral(self.args[0].sage())

Chi.sage = _sympysage_chi

def _sympysage_hyp(self):
    from sage.functions.hypergeometric import hypergeometric
    ap = [arg.sage() for arg in self.args[0]]
    bq = [arg.sage() for arg in self.args[1]]
    return hypergeometric(ap, bq, self.argument.sage())

hyper.sage = _sympysage_hyp

def _sympysage_elliptic_k(self):
    from sage.functions.special import elliptic_kc
    return elliptic_kc(self.args[0].sage())

elliptic_k.sage = _sympysage_elliptic_k

def _sympysage_kronecker_delta(self):
    from sage.functions.generalized import kronecker_delta
    return kronecker_delta(self.args[0].sage(), self.args[1].sage())

KroneckerDelta.sage = _sympysage_kronecker_delta

def _sympysage_besselj(self):
    from sage.functions.bessel import bessel_J
    return bessel_J(self.args[0].sage(), self.args[1].sage())

besselj.sage = _sympysage_besselj

def _sympysage_bessely(self):
    from sage.functions.bessel import bessel_Y
    return bessel_Y(self.args[0].sage(), self.args[1].sage())

bessely.sage = _sympysage_bessely

def _sympysage_besseli(self):
    from sage.functions.bessel import bessel_I
    return bessel_I(self.args[0].sage(), self.args[1].sage())

besseli.sage = _sympysage_besseli

def _sympysage_besselk(self):
    from sage.functions.bessel import bessel_K
    return bessel_K(self.args[0].sage(), self.args[1].sage())

besselk.sage = _sympysage_besselk

def _sympysage_ynm(self):
    from sage.functions.special import spherical_harmonic
    return spherical_harmonic(self.args[0].sage(),
                              self.args[1].sage(),
                              self.args[2].sage(),
                              self.args[3].sage())

Ynm.sage = _sympysage_ynm

#really necessary?
def _sympysage_csch(self):
    from sage.functions.hyperbolic import csch
    return csch(self.args[0].sage())

csch.sage = _sympysage_csch

#really necessary?
def _sympysage_sech(self):
    from sage.functions.hyperbolic import sech
    return sech(self.args[0].sage())

sech.sage = _sympysage_sech

def _sympysage_re(self):
    from sage.functions.other import real_part
    return real_part(self.args[0].sage())

re.sage = _sympysage_re

def _sympysage_im(self):
    from sage.functions.other import imag_part
    return imag_part(self.args[0].sage())

im.sage = _sympysage_im

#really necessary?
def _sympysage_sign(self):
    from sage.functions.generalized import sign
    return sign(self.args[0].sage())

sign.sage = _sympysage_sign

def _sympysage_abs(self):
    from sage.functions.generalized import abs_symbolic
    return abs_symbolic(self.args[0].sage())

Abs.sage = _sympysage_sign







