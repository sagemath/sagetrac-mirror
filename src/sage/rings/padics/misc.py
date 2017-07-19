r"""
Miscellaneous Functions

This file contains several miscellaneous functions used by `p`-adics.

- ``gauss_sum`` -- Computes gauss sums using the Gross-Koblitz formula.
- ``min`` -- a version of ``min`` that returns `\infty` on empty input.
- ``max`` -- a version of ``max`` that returns `-\infty` on empty input.

AUTHORS:

- David Roe
- Adriana Salerno
- Ander Steele
"""

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from six.moves.builtins import min as python_min
from six.moves.builtins import max as python_max
from sage.rings.infinity import infinity
from sage.rings.padics.factory import Zp
from sage.rings.all import PolynomialRing

def gauss_sum(a, p, f, prec = 20):
    """
    Return the gauss sum `g_q(a)` as a `p`-adic number.

    `g_q(a)` is defined by `g_q(a)= \sum_{u\in F_q^*} \omega(u)^{(-a)} \zeta_q^u`
    where `q = p^f`, `\omega` is the Teichmuller character and `\zeta_q` is some arbitrary 
    choice of primitive `q`-th root of unity


   INPUT:

    - ``a`` -- integer.

    - ``p`` -- prime.

    - ``f`` -- positive integer.

    - ``prec`` -- positive integer. set to 20 by default.

    OUTPUT:

    - a `p`-adic number in an Eisenstein extension of Q_p

    .. NOTE::

        This is based on GP code written by Adriana Salerno

    EXAMPLES:

    In this example, we verify that g_3(0) = -1

    ::

        sage: from sage.rings.padics.misc import gauss_sum
        sage: -gauss_sum(0,3,1)
        1 + O(pi^40)

    Next, we verify that g_5(a)*g_5(-a) = 5*(-1)^a 

    ::

        sage: from sage.rings.padics.misc import gauss_sum
        sage: gauss_sum(2,5,1)^2-5
        O(pi^84)
        sage: gauss_sum(1,5,1)*gauss_sum(3,5,1)+5
        O(pi^84)

    Finally, we compute a non-trivial value

    ::

        sage: from sage.rings.padics.misc import gauss_sum
        sage: gauss_sum(2,13,2)
        6*pi^2 + 7*pi^14 + 11*pi^26 + 3*pi^62 + 6*pi^74 + 3*pi^86 + 5*pi^98 + pi^110 + 7*pi^134 + 9*pi^146 + 4*pi^158 + 6*pi^170 + 4*pi^194 + pi^206 + 6*pi^218 + 9*pi^230 + O(pi^242)

    """
    a = a % (p**f)
    R = Zp(p, prec)
    R_poly = PolynomialRing(R,name='X')
    X = R_poly.gen()
    F = R.ext(X**(p-1)+p, names='pi')
    pi = F.gen()
    digits = Zp(p)(a).list(start_val = 0)
    n = len(digits)
    digits = digits + [0]*(f-n)
    s = sum(digits)
    out = -pi**(s)
    for i in range(0,f):
        a_i = R(sum([digits[k]*p**((i+k)%f) for k in range(f)]))        
        if a_i:
            out = out*R((a_i/(p**f-1)).gamma())                                
    return out

def min(*L):
    r"""
    Return the minimum of the inputs, where the minimum of the empty
    list is `\infty`.

    EXAMPLES::

        sage: from sage.rings.padics.misc import min
        sage: min()
        +Infinity
        sage: min(2,3)
        2
    """
    if len(L) == 1 and isinstance(L[0], (list, tuple)):
        L = L[0]
    try:
        return python_min(L)
    except ValueError:
        return infinity


def max(*L):
    r"""
    Return the maximum of the inputs, where the maximum of the empty
    list is `-\infty`.

    EXAMPLES::

        sage: from sage.rings.padics.misc import max
        sage: max()
        -Infinity
        sage: max(2,3)
        3
    """
    if len(L) == 1 and isinstance(L[0], (list, tuple)):
        L = L[0]
    try:
        return python_max(L)
    except ValueError:
        return -infinity
