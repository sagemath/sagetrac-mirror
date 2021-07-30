r"""
Double Exponential Integration for Vector-Valued Functions

A routine to use double exponential (otherwise known as tanh-sinh) quadrature
for numerical integration. 

EXAMPLES:

We verify that `\int_0^1 \frac{1}{2t^{1/2}} \, dt =1`::

    sage: from sage.numerical.double_exponential import integrate_vector_unbounded
    sage: P = 100
    sage: K = RealField(P)
    sage: integrand = lambda t: V([t^(-1/2)/2])
    sage: I = integrate_vector_unbounded(integrand, P)
    sage: error = (I-V([1])).norm()
    sage: bool(error<1e-31)
    True

AUTHORS:

 - Linden Disney-Hogg (2021-07-28): initial version
"""

# ****************************************************************************
#       Copyright (C) 2021 Linden Disney-Hogg <disneyhogg@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.real_mpfr import RealField
from sage.schemes.riemann_surfaces.riemann_surface import ConvergenceError

def integrate_vector_unbounded(f, prec, epsilon=None, n=1):
    r"""
    Integrate a one-argument vector valued function numerically using double
    exponential quadrature. 

    This functions returns an approximation to the integral `\int_0^1 f(t) \, dt`.
    By the nature of the double exponential scheme and the implementation, `f`
    may by singular at `t=0` provided the integral is well defined. A heuristic
    error bound, as described in [Bai2006]_, is used to determine when 
    convergence is achieved. Some optimisation of the method mentioned by Bailey
    are implemented, importantly the use of a higher precision field to 
    calculate the nodes and weights, see the input section. 

    INPUT:

     - ``f`` -- callable. Vector-valued integrand.

     - ``prec`` -- integer. Binary precision to be used.

     - ``epsilon`` -- multiprecision float (default: `2^{(-\text{prec}+3)}`). Target error bound.

     - ``n`` -- positive integer (default: 1). Multiplier of precision to be
       used when calculating nodes and weights. 

    OUTPUT:

    Vector approximating value of the integral.

    EXAMPLES::

    We verify that `\int_0^1 \frac{1}{10t^{9/10}} \, dt =1`::

        sage: from sage.numerical.double_exponential import integrate_vector_unbounded
        sage: P = 100
        sage: K = RealField(P)
        sage: epsilon = K(2)^(-P+3)
        sage: a = 9/10
        sage: integrand = lambda t: V([(1-a)*t^(-a)])
        sage: I = integrate_vector_unbounded(integrand, P)
        sage: error = (I-V([1])).norm()
        sage: bool(error<epsilon)
        True

    The error estimate is only a heuristic and this can be seen by considering
    poles of higher order::

        sage: from sage.numerical.double_exponential import integrate_vector_unbounded
        sage: P = 100
        sage: K = RealField(P)
        sage: epsilon = K(2)^(-P+3)
        sage: a = 99/100
        sage: integrand = lambda t: V([(1-a)*t^(-a)])
        sage: I = integrate_vector_unbounded(integrand, P)
        sage: error = (I-V([1])).norm()
        sage: bool(error<epsilon)
        False

    Numerical investigations suggest that the heuristic error bound is achieved
    for integrands of the form `t^{-a}` for `a \lessapprox 0.9625`.
    """
    R1 = RealField(prec)
    ONE = R1(1)
    if epsilon is None:
        epsilon = R1(2)**(-prec+3)

    Rn = RealField(n*prec)
    tau = Rn(epsilon**n)
    LAMBDA = Rn.pi()/2
    S0 = LAMBDA*f(1/2)/2

    h = ONE
    Nh = (-lambert_w(-1,-tau/2)/LAMBDA).log().ceil()
    hjs = [h*j for j in range(1,Nh+1)]
    u1 = [LAMBDA*hj.cosh() for hj in hjs]
    u2 = [LAMBDA*hj.sinh() for hj in hjs]
    extra_nodes = [(1/(2*U2.exp()*U2.cosh()),U1/(2*U2.cosh()**2)) for U1, U2 in zip(u1,u2)]
    vals = [w*f(y) for y, w in reversed(extra_nodes)] + [S0] + [w*f(1-y) for y, w in extra_nodes]
    results = [sum(vals)]
    D3 = epsilon*max([v.norm(oo) for v in vals])

    h /= 2
    Nh *= 2
    hjs = [h*j for j in range(1,Nh+1,2)]
    u1 = [LAMBDA*hj.cosh() for hj in hjs]
    u2 = [LAMBDA*hj.sinh() for hj in hjs]
    extra_nodes = [(1/(2*U2.exp()*U2.cosh()),U1/(2*U2.cosh()**2)) for U1, U2 in zip(u1,u2)]
    vals = [w*f(y) for y, w in reversed(extra_nodes)] + [w*f(1-y) for y, w in extra_nodes]
    results.append(h*sum(vals)+results[-1]/2)
    D3 = max(D3,epsilon*max([v.norm(oo) for v in vals]))

    for k in range(2,n*(prec-3)-1):
        h /= 2
        Nh *= 2
        hjs = [h*j for j in range(1,Nh+1,2)]
        # u notation from "A comparison of three high-precision quadrature schemes"
        u1 = [LAMBDA*hj.cosh() for hj in hjs]
        u2 = [LAMBDA*hj.sinh() for hj in hjs]
        extra_nodes = [(1/(2*U2.exp()*U2.cosh()),U1/(2*U2.cosh()**2)) for U1, U2 in zip(u1,u2)]
        vals = [w*f(y) for y, w in reversed(extra_nodes)] + [w*f(1-y) for y, w in extra_nodes]
        results.append(h*sum(vals)+results[-1]/2)
        D3 = max(D3,epsilon*max([v.norm(oo) for v in vals]))
        D4 = max([vals[j].norm(oo) for j in [-1,0]])

        if results[-1] == results[-2] or results[2] == results[-3]:
            D = epsilon
        else:
            D1 = (results[-1]-results[-2]).norm(oo)
            D2 = (results[-1]-results[-3]).norm(oo)
            D = min(ONE,max(D1**(D1.log()/D2.log()),D2**2,D3,D4,epsilon))

        if D <= epsilon:
            return results[-1]

    raise ConvergenceError("Integrator failed to converge in {} steps".format(k+1))
