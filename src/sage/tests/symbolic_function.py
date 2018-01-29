"""
Symbolic function tests

The tests :func:`func_test1` and :func:`func_test2` check most one
and two argument functions for matching return type.

EXAMPLES::

    sage: from sage.tests.symbolic_function import func_test1, func_test2
    sage: func_test1()
    sage: func_test2()
"""
from sage.functions.all import *
from sage.rings.all import RR, CC, RDF, CDF, RIF, CIF, RBF, CBF
from sage.misc.all import random
from sage.structure.all import parent

funcs1=[
    sin,
    cos,
    tan,
    cot,
    sec,
    csc,
    (arcsin, 0.7),
    (arccos, 0.7),
    arctan,
    arccot,
    (arcsec, 1.7),
    (arccsc, 1.7),
    zeta,
    sinh,
    cosh,
    tanh,
    (coth, 0.7),
    sech,
    (csch, 0.7),
    arcsinh,
    (arccosh, 1.7),
    (arctanh, 0.7),
    (arccoth, 0.7),
    (arcsech, 0.7),
    (arccsch, 0.7),
    exp,
    (log, 0.7),
    (dilog, 0.7),
    (log_gamma, 1.7),
    (gamma, 1.7),
    (psi, 1.7),
    (lambert_w, 1.7),
    exp_polar,
    (exp_integral_e1, 0.7),
    (log_integral, 0.7),
    (log_integral_offset, 0.7),
    sin_integral,
    (cos_integral, 0.7),
    sinh_integral,
    (cosh_integral, 0.7),
    Ei,
    dickman_rho,
    (elliptic_ec, -.7),
    (elliptic_kc, -.7),
    airy_ai,
    airy_ai,
    airy_ai_prime,
    airy_bi,
    airy_bi,
    airy_bi_prime,
    erf,
    erfi,
    erfc,
    (erfinv, 0.7),]

factories = {RR, CC, RDF, CDF, RIF, CIF, RBF, CBF}

def func_test1():
    for fu in funcs1:
        if isinstance(fu, tuple):
            f = fu[0]
            val = fu[1]
        else:
            f = fu
            val = 5-10*random()
        for g in factories:
            ret = f(g(val))
            p = parent(ret)
            if not p is g:
                print(f,val,g,p)
                #return False
        if not repr(parent(f(float(val)))) == "<type 'float'>":
            print(f,val,'float',p)
            #return False
        if not repr(parent(f(complex(val)))) == "<type 'complex'>":
            print(f,val,'complex',p)
            #return False

funcs2=[
    arctan2,
    zetaderiv,
    (polylog, 7, -5.),
    (hermite, 7, 1.7),
    harmonic_number,
    (exp_integral_e, 7, 1.7),
    gamma_inc_lower,
    hurwitz_zeta,
    bessel_J,
    (bessel_Y, 7, 0.7),
    bessel_I,
    (bessel_K, 7, 0.7),
    struve_H,
    struve_L,
    spherical_bessel_J,
    spherical_bessel_Y,
    (elliptic_e, 7, -0.7),
    (elliptic_eu, 7, 1.7),
    (elliptic_f, 7, -0.7),
    (jacobi_nd, 7, -0.7),
    (jacobi_ns, 7, -0.7),
    (jacobi_nc, 7, -0.7),
    (jacobi_dn, 7, -0.7),
    (jacobi_ds, 7, -0.7),
    (jacobi_dc, 7, -0.7),
    (jacobi_sn, 7, -0.7),
    (jacobi_sd, 7, -0.7),
    (jacobi_sc, 7, -0.7),
    (jacobi_cn, 7, -0.7),
    (jacobi_cd, 7, -0.7),
    (jacobi_cs, 7, -0.7),
    (jacobi_am, 7, 1.7),
    legendre_P,
    laguerre,]

def func_test2():
    for fu in funcs2:
        if isinstance(fu, tuple):
            f = fu[0]
            val1 = fu[1]
            val2 = fu[2]
        else:
            f = fu
            val1 = 7
            val2 = 5.0-10*random()
        for g in factories:
            ret = f(val1, g(val2))
            p = parent(ret)
            if not p is g:
                print(f, val1, val2, g, p)
                #return False
        if not repr(parent(f(val1, float(val2)))) == "<type 'float'>":
            print(f,val1,val2,'float',p)
            #return False
        if not repr(parent(f(val1, complex(val2)))) == "<type 'complex'>":
            print(f,val1,val2,'complex',p)
            #return False

# these are mostly complex
cfuncs2 = [
    hankel1,
    hankel2,
    spherical_hankel1,
    spherical_hankel2,
    inverse_jacobi_nd,
    inverse_jacobi_ns,
    inverse_jacobi_nc,
    inverse_jacobi_dn,
    inverse_jacobi_ds,
    inverse_jacobi_dc,
    inverse_jacobi_sn,
    inverse_jacobi_sd,
    inverse_jacobi_sc,
    inverse_jacobi_cn,
    inverse_jacobi_cd,
    inverse_jacobi_cs,
    legendre_Q,
    ]

"""
not_tested=[imag_part,
abs,
unit_step,
heaviside,
cases,
stieltjes,
    gegenbauer,
    elliptic_pi,
factorial,
binomial,
ceil,
    beta,
floor,
frac,
    chebyshev_T,
    chebyshev_U,
arg,
sum,
product,
    spherical_harmonic,
limit,
complex_root_of,
integrate,
laplace,
    gen_legendre_P,
    gen_legendre_Q,
    jacobi_P,
    gen_laguerre,
prime_pi,
dirac_delta,
sgn,
kronecker_delta,
max,
min,
hypergeometric,
hypergeometric_M,
hypergeometric_U,]
"""
