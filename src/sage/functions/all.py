from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import

lazy_import('sage.functions.piecewise_old', 'Piecewise')   # deprecated
lazy_import('sage.functions.piecewise', 'piecewise')
lazy_import('sage.functions.error', ['erf', 'erfc', 'erfi', 'erfinv'])

from .trig import ( sin, cos, sec, csc, cot, tan,
                   asin, acos, atan,
                   acot, acsc, asec,
                   arcsin, arccos, arctan,
                   arccot, arccsc, arcsec,
                   arctan2, atan2)

from .hyperbolic import ( tanh, sinh, cosh, coth, sech, csch,
                         asinh, acosh, atanh, acoth, asech, acsch,
                         arcsinh, arccosh, arctanh, arccoth, arcsech, arccsch )

reciprocal_trig_functions = {'sec': cos, 'csc': sin, 'cot': tan, 'sech': cosh, 'csch': sinh, 'coth': tanh}



from .other import ( ceil, floor, gamma, psi, factorial, beta, binomial,
                    abs_symbolic, sqrt, log_gamma,
                    gamma_inc, incomplete_gamma, gamma_inc_lower,
                    arg, real_part, real, frac,
                    imag_part, imag, imaginary, conjugate)

from .log import (exp, exp_polar, log, ln, polylog, dilog, lambert_w, harmonic_number)


from .transcendental import (zeta, zetaderiv, zeta_symmetric, hurwitz_zeta,
                            dickman_rho, stieltjes)

from .bessel import (bessel_I, bessel_J, bessel_K, bessel_Y,
                    Bessel, struve_H, struve_L, hankel1, hankel2,
                                   spherical_bessel_J, spherical_bessel_Y,
                                   spherical_hankel1, spherical_hankel2)

from .special import (spherical_harmonic, elliptic_e,
                     elliptic_f, elliptic_ec, elliptic_eu,
                     elliptic_kc, elliptic_pi, elliptic_j, error_fcn)

from .jacobi import (jacobi, inverse_jacobi, jacobi_nd, jacobi_ns, jacobi_nc,
                    jacobi_dn, jacobi_ds, jacobi_dc, jacobi_sn, jacobi_sd,
                    jacobi_sc, jacobi_cn, jacobi_cd, jacobi_cs, jacobi_am,
                    inverse_jacobi_nd, inverse_jacobi_ns, inverse_jacobi_nc,
                    inverse_jacobi_dn, inverse_jacobi_ds, inverse_jacobi_dc,
                    inverse_jacobi_sn, inverse_jacobi_sd, inverse_jacobi_sc,
                    inverse_jacobi_cn, inverse_jacobi_cd, inverse_jacobi_cs)

from .orthogonal_polys import (chebyshev_T,
                              chebyshev_U,
                              gen_laguerre,
                              gen_legendre_P,
                              gen_legendre_Q,
                              hermite,
                              jacobi_P,
                              laguerre,
                              legendre_P,
                              legendre_Q,
                              ultraspherical,
                              gegenbauer)

from .spike_function import spike_function

from .prime_pi import legendre_phi, partial_sieve_function, prime_pi

from .wigner import (wigner_3j, clebsch_gordan, racah, wigner_6j,
                    wigner_9j, gaunt)

from .generalized import (dirac_delta, heaviside, unit_step, sgn, sign,
                         kronecker_delta)

from .min_max import max_symbolic, min_symbolic

from .airy import airy_ai, airy_ai_prime, airy_bi, airy_bi_prime

from .exp_integral import (exp_integral_e, exp_integral_e1, log_integral, li, Li,
                          log_integral_offset,
                          sin_integral, cos_integral, Si, Ci,
                          sinh_integral, cosh_integral, Shi, Chi,
                          exponential_integral_1, Ei, exp_integral_ei)

from .hypergeometric import hypergeometric, hypergeometric_M, hypergeometric_U

