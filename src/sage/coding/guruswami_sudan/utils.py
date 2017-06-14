r"""
Guruswami-Sudan utility methods

AUTHORS:

- Johan S. R. Nielsen, original implementation (see [Nie]_ for details)
- David Lucas, ported the original implementation in Sage
"""

#*****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#                     2015 Johan S. R. Nielsen <jsrn@jsrn.dk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.functions.other import binomial, floor, sqrt
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.arith.all import lcm
from sage.combinat.permutation import Permutation

def polynomial_to_list(p, len):
    r"""
    Returns ``p`` as a list of its coefficients of length ``len``.

    INPUT:

    - ``p`` -- a polynomial

    - ``len`` -- an integer. If ``len`` is smaller than the degree of ``p``, the
      returned list will be of size degree of ``p``, else it will be of size ``len``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import polynomial_to_list
        sage: F.<x> = GF(41)[]
        sage: p = 9*x^2 + 8*x + 37
        sage: polynomial_to_list(p, 4)
        [37, 8, 9, 0]
    """
    return list(p) + [0]*max(0, len-p.degree()-1)

def johnson_radius(n, d):
    r"""
    Returns the Johnson-radius for the code length `n` and the minimum distance `d`.

    The Johnson radius is defined as `n - \sqrt(n(n-d))`.

    INPUT:

    - ``n`` -- an integer, the length of the code
    - ``d`` -- an integer, the minimum distance of the code

    EXAMPLES::

        sage: sage.coding.guruswami_sudan.utils.johnson_radius(250, 181)
        -5*sqrt(690) + 250
    """
    return n - sqrt(n*(n-d))

def ligt(x):
    r"""
    Returns the least integer greater than ``x``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import ligt
        sage: ligt(41)
        42

    It works with any type of numbers (not only integers)::

        sage: ligt(41.041)
        42
    """
    return floor(x+1)

def gilt(x):
    r"""
    Returns the greatest integer smaller than ``x``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import gilt
        sage: gilt(43)
        42

    It works with any type of numbers (not only integers)::

        sage: gilt(43.041)
        43
    """
    if x in ZZ:
        return Integer(x-1)
    else:
        return floor(x)

def solve_degree2_to_integer_range(a,b,c):
    r"""
    Returns the greatest integer range `[i_1, i_2]` such that
    `i_1 > x_1` and `i_2 < x_2` where `x_1, x_2` are the two zeroes of the equation in `x`:
    `ax^2+bx+c=0`.

    If there is no real solution to the equation, it returns an empty range with negative coefficients.

    INPUT:

    - ``a``, ``b`` and ``c`` -- coefficients of a second degree equation, ``a`` being the coefficient of
      the higher degree term.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import solve_degree2_to_integer_range
        sage: solve_degree2_to_integer_range(1, -5, 1)
        (1, 4)

    If there is no real solution::

        sage: solve_degree2_to_integer_range(50, 5, 42)
        (-2, -1)
    """
    D = b**2 - 4*a*c
    if D < 0:
        return (-2,-1)
    sD = float(sqrt(D))
    minx, maxx = (-b-sD)/2.0/a , (-b+sD)/2.0/a
    mini, maxi = (ligt(minx), gilt(maxx))
    if mini > maxi:
        return (-2,-1)
    else:
        return (mini,maxi)

def _degree_of_vector(v, shifts=None):
    r"""
    Returns the greatest degree among the entries of the polynomial vector `v`.

    INPUT:

    - ``v`` -- a vector of polynomials.

    - ``shifts`` -- (default: ``None``) a list of integer shifts to consider
      ``v`` under, i.e. compute `\max(\deg v_i + s_i)`, where `s_1,\ldots, s_n`
      is the list of shifts.

      If ``None``, all shifts are considered as ``0``.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.utils import _degree_of_vector
        sage: F.<x> = GF(7)[]
        sage: v = vector(F, [0, 1, x, x^2])
        sage: _degree_of_vector(v)
        2
        sage: _degree_of_vector(v, shifts=[10, 1, 0, -3])
        1
    """
    if not shifts:
        return max( vi.degree() for vi in v )
    else:
        if v.is_zero():
            -1
        return max( degi + si for (degi,si) in zip([ vi.degree() for vi in v ],shifts) if degi > -1 )
