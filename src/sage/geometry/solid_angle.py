r"""
Solid angle of a polyhedral cone

This module implements the *normalized solid angle measure* of polyhedral
cones. The normalized solid angle measure of a cone is the amount of space
the cone takes up in reference to some given set.

EXAMPLES:

Compute the normalized solid angle measure of the first quadrant::

    sage: from sage.geometry.solid_angle import solid_angle_simplicial_2d
    sage: RDF(solid_angle_simplicial_2d(matrix([[0,1],[1,0]])))
    0.25

AUTHORS:

- Yuan Zhou (2022)

- Allison Fitisone (2022)
"""
# ****************************************************************************
#       Copyright (C) 2022 Yuan Zhou <yuan.zhou@uky.edu>
#       Copyright (C) 2022 Allison Fitisone <allison.fitisone@uky.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.symbolic.constants import pi
from sage.matrix.constructor import matrix
from sage.structure.element import is_Matrix
from sage.functions.trig import arccos
from sage.symbolic.subring import SymbolicSubring
from sage.rings.integer_ring import ZZ


def solid_angle_simplicial_2d(A):
    r"""
    Return the normalized solid angle measure of the cone generated by the
    row vector(s) of ``A``.

    INPUT:

    - ``A`` -- 2x2 matrix in the form of ``matrix([[a,b],[c,d]])`` or
      simply ``[[a,b],[c,d]]`` where the nonzero linearly independent
      rows of ``A`` generate the cone in `\RR^2`

    OUTPUT:

    - the normalized solid angle measure of the cone generated by the
      vectors

    EXAMPLES:

    This example shows the normalized measure of the solid angle spanned
    by the rows of the matrix::

        sage: from sage.geometry.solid_angle import solid_angle_simplicial_2d
        sage: solid_angle_simplicial_2d(matrix([[0, 1], [1, 0]]))
        1/4

    The input can be a list of vectors instead of a matrix as well::

        sage: solid_angle_simplicial_2d([[1, 0], [-1, sqrt(3)]])
        1/3

        sage: solid_angle_simplicial_2d(matrix([[0, 1], [4.5, 0]]))
        0.785398163397448/pi

        sage: RDF(solid_angle_simplicial_2d([[2, 13], [-1, 7]]))  # abs tol 1e-15
        0.04687851282419763

    This example illustrates how the solid angle measure will not be
    greater than 0.5 as the function always outputs the minimal angle
    between the two rays::

        sage: RDF(solid_angle_simplicial_2d([[1, 0], [-1, -1]]))  # abs tol 1e-15
        0.375

    .. NOTE::

        This function uses the dot product of two vectors to determine the
        angle between them.

    The following tests check for corner cases where the vectors are
    antiparallel, parallel and perpendicular, respectively::

        sage: solid_angle_simplicial_2d([[1, 1], [-1, -1]])
        0

        sage: solid_angle_simplicial_2d([[1, 2], [2, 4]])
        0

        sage: solid_angle_simplicial_2d([[2, 2], [-1, 1]])
        1/4

    The following tests check the assumptions of the input::

        sage: solid_angle_simplicial_2d([[-3, 2]])
        Traceback (most recent call last):
        ...
        ValueError: input matrix has incorrect dimension

        sage: solid_angle_simplicial_2d([[1, 4], [0, 0]])
        Traceback (most recent call last):
        ...
        ValueError: input matrix has a row that is zero

    TESTS:

    In the following examples, we check the parent of the output corresponding
    to an input whose parent is exact.::

        sage: A = matrix([[2/3, 0], [-1/3, 5]])
        sage: solid_angle_simplicial_2d(A)
        1/2*arccos(-1/226*sqrt(226))/pi
        sage: solid_angle_simplicial_2d(A).parent()
        Symbolic Ring

        sage: A = matrix([[2/3, 0], [-1/3, 0]])
        sage: solid_angle_simplicial_2d(A)
        0
        sage: solid_angle_simplicial_2d(A).parent()
        Symbolic Constants Subring

        sage: A = matrix([[1, 1], [-4, -4]])
        sage: solid_angle_simplicial_2d(A)
        0
        sage: solid_angle_simplicial_2d(A).parent()
        Symbolic Constants Subring

    In the following examples, we check the parent of the output corresponding
    to an input whose parent is not exact::

        sage: A = matrix([[12, 0], [RBF(1.5), 9]])
        sage: solid_angle_simplicial_2d(A)
        ([0.702823824690135 +/- 4.26e-16])/pi
        sage: solid_angle_simplicial_2d(A).parent()
        Symbolic Ring

        sage: A = matrix([[RBF(12.7), 0], [-1, 0]])
        sage: solid_angle_simplicial_2d(A)
        0
        sage: solid_angle_simplicial_2d(A).parent()
        Real ball field with 53 bits of precision

        sage: A = matrix([[12, 0], [sqrt(17), 9]])
        sage: solid_angle_simplicial_2d(A)
        1/2*arccos(1/14*sqrt(17)*sqrt(2))/pi
        sage: solid_angle_simplicial_2d(A).parent()
        Symbolic Ring

        sage: A = matrix([[0, sqrt(18)], [0, -sqrt(pi)]])
        sage: solid_angle_simplicial_2d(A)
        0
        sage: solid_angle_simplicial_2d(A).parent()
        Symbolic Ring

        sage: A = matrix([[0, 1], [RDF(pi), RDF(pi)]])
        sage: solid_angle_simplicial_2d(A)
        0.39269908169872414/pi
        sage: solid_angle_simplicial_2d(A).parent()
        Symbolic Ring

        sage: A = matrix([[-1, -1], [RDF(pi), RDF(pi)]])
        sage: solid_angle_simplicial_2d(A)
        0.0
        sage: solid_angle_simplicial_2d(A).parent()
        Real Double Field

        sage: A = matrix([[0, 1], [RR(-pi), RR(0)]])
        sage: solid_angle_simplicial_2d(A)
        0.785398163397448/pi
        sage: solid_angle_simplicial_2d(A).parent()
        Symbolic Ring

        sage: A = matrix([[1, 0], [RR(-pi), RR(0)]])
        sage: solid_angle_simplicial_2d(A)
        0.000000000000000
        sage: solid_angle_simplicial_2d(A).parent()
        Real Field with 53 bits of precision
    """
    if not is_Matrix(A):
        A = matrix(A)
    P = A.base_ring()
    if A.nrows() != 2 or A.ncols() != 2:
        raise ValueError("input matrix has incorrect dimension")
    if any(r == 0 for r in A.rows()):
        raise ValueError("input matrix has a row that is zero")
    if A.rank() < 2:
        import sage.rings.abc
        if P.is_exact():
            return SymbolicSubring(no_variables=True)(ZZ(0))
        else:
            return P.zero()
    u = A.row(0)
    v = A.row(1)
    p = u.dot_product(v)
    a = u.norm()
    b = v.norm()
    cs = p/(a*b)
    final_calc = arccos(cs) / (2*pi)
    return final_calc
