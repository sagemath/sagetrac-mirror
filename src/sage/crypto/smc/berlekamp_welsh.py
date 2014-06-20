# coding: UTF-8
r"""
Berlekamp-Welsh algorithm

The Berlekamp-Welsh algorithm serves for the purpose of reconstructing 
polynomials with erroneous points.
The algorithm is well explained at
`Wikipedia <https://en.wikipedia.org/wiki/Berlekamp–Welch_algorithm>`_.

AUTHORS:

- Thomas Loruenser (2013): initial version

"""
###############################################################################
# Copyright 2013, Thomas Loruenser <thomas.loruenser@ait.ac.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

from sage.structure.sage_object import SageObject


def berlekamp_welsh(deg, points):
    r"""
    Reconstruct polynomial with Berlekamp-Welsh algorithm.

    INPUT:

    - ``deg``    --  degree of polynomial to reconstruct.
    - ``points`` --  array of points (list of (x,y)-tuples).

    OUTPUT:

    Reconstructed polynomial.

    EXAMPLES::

        sage: from sage.crypto.smc.berlekamp_welsh import berlekamp_welsh
        sage: from sage.rings.finite_rings.constructor import FiniteField
        sage: from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

    Reconstruction with errors::

        sage: order = 2**8
        sage: F = FiniteField(order, 'a')
        sage: P = PolynomialRing(F, 'x')
        sage: n = 7
        sage: deg = 2
        sage: poly = F.fetch_int(42)
        sage: for i in range(1, deg+1): poly += F.random_element() * P.gen()**i

        sage: # evaluate polynomial at different points (shares)
        sage: points = [(F.fetch_int(i), poly(F.fetch_int(i))) for i in range(1, n+1)]
        sage: poly == berlekamp_welsh(deg, points)
        True

        sage: # introduce error
        sage: points[0] = (points[0][0], points[0][1] + F.fetch_int(9))
        sage: poly == berlekamp_welsh(deg, points)
        True

    """
    # check input vector
    F = points[0][0].parent()
    if not F.is_field():
        raise TypeError("points must be of field type.")
    for x, y in points:
        if x.parent() != F or y.parent() != F:
            raise TypeError("all points must be from same field.")
        
    # generate and solve system of linear equations
    from sage.functions.all import floor
    deg_E = floor((len(points) - (deg + 1)) / 2.)
    deg_Q = deg_E + deg
    from sage.matrix.all import Matrix
    from sage.all import vector
    sys_size = deg_Q + 1 + deg_E
    A = Matrix(F, sys_size)
    b = vector(F, sys_size)
    for n, (x, y) in enumerate(points):
        A[n] = ([x**i for i in range(deg_Q+1)]+[-y * x**i for i in range(deg_E)])
        b[n] = (y * x**deg_E)
        QE = A.solve_right(b)

    # reconstruct polynomial
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    P = PolynomialRing(F, 'x')
    Q = sum([coeff * P.gen()**i for i, coeff in enumerate(QE[:deg_Q+1])]);
    E = P.gen()**deg_E
    E += sum([coeff * P.gen()**i for i, coeff in enumerate(QE[deg_Q+1:])]);
    P = Q.quo_rem(E)[0]
    return P

# vim: set fileencoding=UTF-8 filetype=python :
