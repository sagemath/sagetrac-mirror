r"""
Solid angle of a polyhedral cone

This module implements the *normalized solid angle measure* of polyhedral
cones. The normalized solid angle measure of a cone is the amount of space
the cone takes up in reference to some given set.

EXAMPLES:

Compute the normalized solid angle measure of the first quadrant::

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

        sage: solid_angle_simplicial_2d(matrix([[0, 1], [1, 0]]))
        1/4

    The input can be a list of vectors instead of a matrix as well::

        sage: solid_angle_simplicial_2d([[1, 0], [-1, sqrt(3)]])
        1/3

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
        ValueError: input matrix has incorrect dimension.

        sage: solid_angle_simplicial_2d([[1, 4], [0, 0]])
        Traceback (most recent call last):
        ...
        ValueError: input matrix has a row that is zero.
    """
    if not is_Matrix(A):
        A = matrix(A)
    if A.nrows() != 2 or A.ncols() != 2:
        raise ValueError("input matrix has incorrect dimension.")
    if any(r == 0 for r in A.rows()):
        raise ValueError("input matrix has a row that is zero.")
    if A.rank() < 2:
        return 0
    u = A.row(0)
    v = A.row(1)
    p = u.dot_product(v)
    a = u.norm()
    b = v.norm()
    cs = p/(a*b)
    final_calc = arccos(cs) / (2*pi)
    return final_calc
