"""
Periods

Fast Cython code needed for computation of period integrals associated
to spaces of modular symbols.
"""

#############################################################################
#       Copyright (C) 2012 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

cdef extern:
    complex cexp(complex)

import math
cdef complex c0 = complex(0,2*math.pi)  # 2*pi*i

cpdef complex exp_z_integral(complex alpha, unsigned long n, unsigned int m):
    """
    Numerically compute the integral

       integral_{alpha}^{infty} exp(2*pi*i*n*z) z^m dz

    where alpha is complex, n>=1 an integer, and m>=0 an integer.
    This uses Lemma 10.4 from Stein's 2007 AMS modular forms book.

    Note that m is typically small but n gets large.
    """
    cdef complex t = n*c0                  # = 2*pi*i*n
    cdef complex one_over_t = 1/t
    cdef double j_prod = 1                 # = prod_{(m+1)-s}^m j
    cdef complex denom = 1/t               # = 1/(2*pi*i*n)^(s+1)
    cdef complex alpha_pow = alpha**m      # = alpha^(m-s)
    cdef complex alpha_inv = 1/alpha
    cdef int sgn = 1                       # = (-1)^s

    cdef unsigned int s

    cdef complex summation = sgn*alpha_pow*denom*j_prod

    for s in range(1,m+1):
        j_prod *= m+1-s
        denom *= one_over_t
        sgn *= -1
        alpha_pow *= alpha_inv
        summation += sgn*alpha_pow*denom*j_prod

    return cexp(t*alpha) * summation

cpdef complex extended_period_integral(unsigned int m, complex alpha, list v):
    """
    Entries of v = [a0,a1,a2,...] are assumed to be complex.
    """
    # There are many nicer ways that this code could be written, e.g., using
    # sum, range, enumerate, etc. -- don't bother, as they are all way slower,
    # at least with Cython 0.15.1.
    cdef complex summation = 0
    cdef unsigned long n = 1
    cdef complex z
    for z in v[1:]:
        summation += z * exp_z_integral(alpha, n, m)
        n += 1
    return summation * c0
